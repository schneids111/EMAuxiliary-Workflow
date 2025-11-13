## ===== helpers: within/between markers (kept for back-compat; show 1x note) =====
within  <- function(x) x
between <- function(x) x
w <- function(x) x
b <- function(x) x
.__blimpwrap_markers_warned <- FALSE

## ===== small utilities =====
.bw_strip_markers <- function(s) {
  s <- gsub("(within|w)\\(([^)]+)\\)", "\\2", s, perl = TRUE)
  s <- gsub("(between|b)\\(([^)]+)\\)", "\\2", s, perl = TRUE)
  s
}
.bw_terms_to_tokens <- function(vec) {
  out <- gsub("(within|w)\\(([^)]+)\\)", "\\2", vec, perl = TRUE)     # within -> x
  out <- gsub("(between|b)\\(([^)]+)\\)", "\\2.mean", out, perl = TRUE) # between -> x.mean
  out
}
.bw_has_markers <- function(ftxt) {
  grepl("(within|w)\\(|(between|b)\\(", ftxt, perl = TRUE)
}
.bw_is_level2_var <- function(v, id, data) {
  if (!v %in% names(data)) return(FALSE)
  sds <- tapply(data[[v]], data[[id]], function(z) stats::sd(z, na.rm = TRUE))
  all(is.na(sds) | sds < 1e-8)
}
.bw_icc <- function(v, id, data) {
  x <- data[[v]]
  grp <- data[[id]]
  m_g <- tapply(x, grp, function(z) mean(z, na.rm = TRUE))
  m_g <- m_g[match(grp, names(m_g))]
  x_c <- x - m_g
  wvar <- stats::var(x_c, na.rm = TRUE)
  bvar <- stats::var(m_g[!duplicated(grp)], na.rm = TRUE)
  if (!is.finite(wvar) || !is.finite(bvar) || (wvar + bvar) <= 0) return(NA_real_)
  bvar / (bvar + wvar)
}

## Parse last PSR from BLIMP output
.bw_last_psr <- function(rblimp_fit) {
  out <- utils::capture.output(rblimp::output(rblimp_fit))
  idx <- grep("^\\s*Comparing iterations across", out)
  if (!length(idx)) return(NA_real_)
  block <- out[idx:length(out)]
  lines <- grep("^\\s*\\d+ to \\d+\\s+\\S+\\s+\\d+\\s*$", block, value = TRUE)
  if (!length(lines)) return(NA_real_)
  last <- tail(lines, 1)
  psr <- sub(".*?\\s([0-9]+\\.[0-9]+)\\s+\\d+\\s*$", "\\1", last)
  suppressWarnings(as.numeric(psr))
}

## ===== printing helpers (optional) =====
blimp_print_psr <- function(fit_obj) {
  cat("---- PSR section ----\n")
  out <- utils::capture.output(rblimp::output(fit_obj$fit))
  from <- grep("^\\s*BURN-IN POTENTIAL SCALE REDUCTION", out)
  to   <- grep("^\\s*METROPOLIS-HASTINGS ACCEPTANCE RATES", out)
  toalt <- grep("^\\s*DATA INFORMATION", out)
  if (length(from) && length(to)) {
    cat(paste(out[from:to], collapse = "\n"), "\n")
  } else if (length(from) && length(toalt)) {
    cat(paste(out[from:toalt], collapse = "\n"), "\n")
  } else {
    cat("PSR section not found in output.\n")
  }
}

blimp_print_focal <- function(obj, y = NULL) {
  # obj can be rblimp fit or our blimp_wrap_fit (with $fit)
  fit <- if (inherits(obj, "blimp_wrap_fit")) obj$fit else obj
  
  out <- utils::capture.output(rblimp::output(fit))
  if (!length(out)) {
    message("No output captured from rblimp::output().")
    return(invisible(NULL))
  }
  
  re_escape <- function(s) gsub("([][(){}^$.|*+?\\-\\\\])", "\\\\\\1", s)
  
  # If y not provided and obj has meta$y, use it; otherwise fall back to first outcome
  if (is.null(y) && !is.null(obj$meta) && !is.null(obj$meta$y)) {
    y <- obj$meta$y
  }
  
  # Find start: "Outcome Variable: <y>" if known; else first "Outcome Variable:"
  start_idx <- NA_integer_
  if (!is.null(y) && nzchar(y)) {
    pat <- paste0("^\\s*Outcome\\s+Variable:\\s+", re_escape(y), "\\b")
    hit <- which(grepl(pat, out))
    if (length(hit)) start_idx <- hit[0 + 1L]
  }
  if (is.na(start_idx)) {
    hit <- which(grepl("^\\s*Outcome\\s+Variable:\\s+", out))
    if (length(hit)) start_idx <- hit[0 + 1L]
  }
  if (is.na(start_idx)) {
    message("Could not find an 'Outcome Variable:' section.")
    return(invisible(NULL))
  }
  
  # End just before "auxiliary.model block:" or "PREDICTOR MODEL ESTIMATES:"
  aux_idx  <- which(grepl("^\\s*auxiliary\\.model\\s+block:\\s*$", out))
  pred_idx <- which(grepl("^\\s*PREDICTOR\\s+MODEL\\s+ESTIMATES:\\s*$", out))
  candidates <- c(aux_idx[aux_idx > start_idx], pred_idx[pred_idx > start_idx])
  end_idx <- if (length(candidates)) min(candidates) - 1L else length(out)
  
  focal <- out[start_idx:end_idx]
  cat(paste(focal, collapse = "\n"), "\n")
  invisible(focal)
}


## ===== main wrapper =====
EMAuxiliary <- function(formula, data, id,
                       aux = NULL,
                       center_group = NULL,
                       center_grand = NULL,
                       ordinal = NULL,
                       nominal = NULL,
                       burn = 10000, iter = 10000, chains = 2, seed = 12345) {

  stopifnot(is.data.frame(data))
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Install lme4.")
  if (!requireNamespace("rblimp", quietly = TRUE)) stop("Install rblimp.")

  ftxt <- paste(deparse(formula), collapse = "")
  if (.bw_has_markers(ftxt) && !.__blimpwrap_markers_warned) {
    message("[EMAuxilary] Note: within()/between() (and w()/b()) are supported for convenience, ",
            "but we recommend BLIMP-native syntax instead:\n",
            "  • Put bases in center_group = c(\"x\", ...)\n",
            "  • Refer to latent means as x.mean in the formula.\n")
    .__blimpwrap_markers_warned <<- TRUE
  }

  fixed_form <- lme4::nobars(formula)
  fb <- lme4::findbars(formula)
  rand_slopes_raw <- character(0)
  if (length(fb)) {
    inside_terms <- unique(unlist(lapply(fb, function(b) {
      left <- deparse(b[[2]])
      attr(terms(as.formula(paste("~", left))), "term.labels")
    })))
    rand_slopes_raw <- setdiff(inside_terms, "1")
  }
  
  tl_raw    <- attr(terms(fixed_form), "term.labels")
  tl_tokens <- .bw_terms_to_tokens(tl_raw)
  base_raw  <- unique(.bw_strip_markers(unique(unlist(strsplit(tl_raw, ":", fixed = TRUE)))))
  base_raw  <- base_raw[nzchar(base_raw)]
  
  blimp_terms <- if (length(tl_tokens)) gsub(":", "*", tl_tokens, fixed = TRUE) else character(0)
  rhs_tokens  <- unique(blimp_terms)
  fixed_rhs   <- if (length(rhs_tokens)) paste(rhs_tokens, collapse = " ") else "1"
  
  y <- deparse(formula[[2L]])
  
  # Only auto-center y if we need y.mean for auxiliaries or if the focal RHS references y.mean
  need_y_mean <- (length(aux) > 0L) ||
    any(grepl(paste0("\\b", gsub("([][(){}^$.|*+?\\-\\\\])","\\\\\\1", y), "\\.mean\\b"), rhs_tokens))
  
  # Start from user-specified group mean list; add y only if needed
  cg <- unique(c(center_group, if (need_y_mean) y else NULL))
  cg <- cg[!is.na(cg) & nzchar(cg)]
  
  
  
  # tokens that appear on RHS (already with * instead of :)
  rhs_tokens  <- unique(blimp_terms)
  
  # --- robust extraction of .mean bases ---
  # 1) split RHS into atomic pieces (by space or '*')
  rhs_atoms <- unique(unlist(strsplit(rhs_tokens, "\\s+|\\*", perl = TRUE)))
  rhs_atoms <- rhs_atoms[nzchar(rhs_atoms)]
  
  # 2) keep only clean .mean atoms like name.mean (no operators)
  mean_atoms <- grep("^[A-Za-z]\\w*\\.mean$", rhs_atoms, value = TRUE)
  
  # 3) strip .mean to get base names
  means_in_rhs  <- unique(sub("\\.mean$", "", mean_atoms))
  
  # 4) enforce that bases are in center_group (plus y if need_y_mean)
  missing_bases <- setdiff(means_in_rhs, cg)
  if (length(missing_bases)) {
    stop(
      "You used ", paste0("'", paste0(paste0(missing_bases, ".mean"), collapse = "', '"), "'"),
      " in the formula but did not include the base name(s) in `center_group`.\n",
      "Add: center_group = c(", paste(sprintf('\"%s\"', missing_bases), collapse = ", "), ", ...).",
      call. = FALSE
    )
  }
  
  
  # Build CENTER: only if there is anything to center
  center_pieces <- character(0)
  if (length(cg))            center_pieces <- c(center_pieces, paste0("groupmean = ", paste(cg, collapse = " ")))
  if (length(center_grand))  center_pieces <- c(center_pieces, paste0("grandmean = ", paste(center_grand, collapse = " ")))
  center_block <- if (length(center_pieces)) paste0("CENTER: ", paste(center_pieces, collapse = " ; "), ";") else NULL
  
  # ---- robust '.mean' validation (handles interactions and [cluster] suffix) ----
  split_tokens <- function(v) {
    out <- unlist(strsplit(v, "\\*"))
    trimws(out[nzchar(out)])
  }
  atomic_tokens <- unique(unlist(lapply(rhs_tokens, split_tokens)))
  
  # tokens like x.mean or x.mean[id]
  is_mean_tok  <- grepl("\\.mean(\\[[^]]+\\])?$", atomic_tokens)
  mean_tokens  <- atomic_tokens[is_mean_tok]
  
  # strip ".mean" and optional "[cluster]" → base var
  means_in_rhs <- unique(gsub("\\.mean(\\[[^]]+\\])?$", "", mean_tokens))
  
  missing_bases <- setdiff(means_in_rhs, cg)
  if (length(missing_bases)) {
    stop(
      paste0(
        "You used ",
        paste(sprintf("'%s.mean'", missing_bases), collapse = ", "),
        " in the formula but did not include the base variable(s) in center_group.\n",
        "Add: center_group = c(",
        paste(sprintf('"%s"', missing_bases), collapse = ", "),
        if (length(center_group)) ", ...)" else ")"
      ),
      call. = FALSE
    )
  }
  

  if (length(rand_slopes_raw)) {
    if (any(grepl("\\.mean\\b", rand_slopes_raw))) {
      stop("Random slopes cannot include latent mean tokens (*.mean). Use only within-level bases.", call. = FALSE)
    }
  }
  rand_slopes <- .bw_strip_markers(rand_slopes_raw)

  center_lines <- character(0)
  if (length(cg)) {
    center_lines <- c(center_lines, paste0("groupmean = ", paste(cg, collapse = " "), ";"))
  }
  if (!is.null(center_grand) && length(center_grand)) {
    center_lines <- c(center_lines, paste0("grandmean = ", paste(center_grand, collapse = " "), ";"))
  }
  center_block <- if (length(center_lines)) paste0("CENTER: ", paste(center_lines, collapse = " ")) else NULL

  if (length(rand_slopes) == 0) {
    model_line      <- paste0("MODEL:\n  ", y, " ~ ", fixed_rhs, ";")
    blimp_model_arg <- paste0(y, " ~ ", fixed_rhs)
  } else {
    model_line      <- paste0("MODEL:\n  ", y, " ~ ", fixed_rhs, " | ", paste(rand_slopes, collapse = " "), ";")
    blimp_model_arg <- paste0(y, " ~ ", fixed_rhs, " | ", paste(rand_slopes, collapse = " "))
  }

  rhs_bet_means <- unique(sub("\\.mean$", "", grep("\\.mean\\b", rhs_tokens, value = TRUE)))
  pure_L2_obs <- base_raw[vapply(base_raw, function(v) .bw_is_level2_var(v, id, data), logical(1))]

  between_token_for <- function(x) {
    if (x %in% pure_L2_obs) return(x)
    if (x %in% cg) return(paste0(x, ".mean"))
    NA_character_
  }

  aux_block <- NULL
  aux_model_arg <- NULL
  aux_scaled <- character(0)

  if (!is.null(aux) && length(aux)) {

    # ---- Check for missing auxiliary variables ----
    aux_missing <- setdiff(aux, names(data))
        if (length(aux_missing)) {
      #warning(
      #  "[EMAuxiliary] The following auxiliary variable(s) were listed but not found in the data: ",
      #  paste(aux_missing, collapse = ", "),
      #  "\nThey will be ignored.",
      stop(
        "[EMAuxiliary] ERROR: Auxiliary variable(s) not found in the data: ",
            paste(aux_missing, collapse = ", "),
            "\nPlease correct these before running EMAuxiliary().",
        call. = FALSE
      )
    }
        # Keep only valid auxiliaries
    aux <- intersect(aux, names(data))
    
    # Make a safe local copy
    data2 <- data
    aux_local <- aux
    # Variables that must never be scaled
    no_scale <- unique(c(ordinal, nominal))
    aux_scaled <- character(0)
    # Loop over *all* auxiliary variables
    for (a in aux_local) {
    # Skip auxiliaries not in data
      if (!a %in% names(data2)) next
    # Skip if ordinal or nominal
      if (a %in% no_scale) next
        xa <- data2[[a]]
    # Skip non-numeric variables
      if (!is.numeric(xa)) next
    # Skip binaries or degenerate vars
      if (length(unique(xa[!is.na(xa)])) < 3) next
    # --- Apply standardization ---
      data2[[a]] <- as.numeric(scale(xa))
      aux_scaled <- c(aux_scaled, a)
    }
    # Now replace the working copy
    data <- data2
    
    
    main_bases <- base_raw
    full_rhs <- rhs_tokens
    aux_lines <- character(0)
    y_mean <- paste0(y, ".mean")

    for (a in aux) {
      if (!a %in% names(data2)) next
      a_is_L2 <- .bw_is_level2_var(a, id, data2)
      a_icc   <- if (!a_is_L2) .bw_icc(a, id, data2) else 1

      if (isTRUE(a_is_L2)) {
        bet_preds <- unique(Filter(Negate(is.na),
                                   vapply(main_bases, between_token_for, FUN.VALUE = character(1))))
        rhs <- unique(c(y_mean, bet_preds))
      } else if (is.finite(a_icc) && a_icc < 0.05) {
        raw_ok <- intersect(main_bases, names(data2))
        rhs <- unique(c(y, raw_ok))
      } else {
        rhs <- unique(c(y, y_mean, full_rhs))
      }

      rhs <- setdiff(rhs, a)
      if (length(rhs)) aux_lines <- c(aux_lines, paste0(a, " ~ ", paste(rhs, collapse = " "), ";"))
    }

    if (length(aux_lines)) {
      aux_block     <- paste0("auxiliary.model:\n  ", paste(aux_lines, collapse = "\n  "))
      aux_model_arg <- paste(aux_lines, collapse = "\n  ")
    }

    data <- data2
    if (length(aux_scaled)) {
      message("[EMAuxilary] Standardized continuous auxiliaries: ",
              paste(aux_scaled, collapse = ", "))
    }
  }

  obs_vars <- unique(c(id,
                       .bw_strip_markers(all.vars(fixed_form)),
                       aux, ordinal, nominal,
                       center_group, center_grand))
  obs_vars <- obs_vars[obs_vars %in% names(data)]
  obs_vars <- setdiff(obs_vars, c("within","between","w","b"))

  blimp_code <- paste(
    "DATA: <in-memory by rblimp>;",
    paste0("VARIABLES: ", paste(obs_vars, collapse = " "), ";"),
    paste0("CLUSTERID: ", id, ";"),
    if (!is.null(ordinal) && length(ordinal)) paste0("ORDINAL: ", paste(ordinal, collapse = " "), ";") else NULL,
    if (!is.null(nominal) && length(nominal)) paste0("NOMINAL: ", paste(nominal, collapse = " "), ";") else NULL,
    if (!is.null(center_block)) center_block else NULL,
    model_line,
    if (!is.null(aux_block)) aux_block else NULL,
    paste0("SEED: ", seed, ";"),
    paste0("BURN: ", burn, ";"),
    paste0("ITER: ", iter, ";"),
    paste0("CHAINS: ", chains, ";"),
    sep = "\n"
  )

  args <- list(
    data = data,
    clusterid = id,
    model = blimp_model_arg,
    seed = seed, burn = burn, iter = iter, chains = chains
  )
  if (length(center_lines))  args$center  <- paste(center_lines, collapse = " ")
  if (!is.null(ordinal) && length(ordinal)) args$ordinal <- paste(ordinal, collapse = " ")
  if (!is.null(nominal) && length(nominal)) args$nominal <- paste(nominal, collapse = " ")

  if (!is.null(aux_model_arg)) {
    fm <- names(formals(rblimp::rblimp))
    aux_names <- intersect(c("auxiliary.model", "auxmodel", "auxiliarymodel", "auxiliary"), fm)
    if (length(aux_names)) {
      args[[aux_names[1]]] <- aux_model_arg
      used_aux_route <- paste0("arg:", aux_names[1])
    } else {
      args$model <- paste0(blimp_model_arg, ";\nauxiliary.model:\n  ", aux_model_arg, ";")
      used_aux_route <- "inlined-into-model"
    }
  } else {
    used_aux_route <- "no-aux"
  }
  message("Auxiliaries route: ", used_aux_route)

  fit <- do.call(rblimp::rblimp, args)

  psr <- tryCatch(.bw_last_psr(fit), error = function(e) NA_real_)
  if (is.finite(psr) && psr > 1.05) {
    header <- sprintf(
      "Final PSR is %.3f. The final PSR should be < 1.05.\nThe results are NOT TRUSTWORTHY!\nConsider increasing burn/iter, more chains, simplifying the model, or reducing auxiliaries.\n",
      psr
    )
    ratio_thresh <- 6
    abs_thresh   <- 3
    base_preds <- base_raw
    focal_cols <- unique(c(y, base_preds))
    sds <- vapply(focal_cols, function(v) {
      if (!v %in% names(data) || !is.numeric(data[[v]])) return(NA_real_)
      stats::sd(data[[v]], na.rm = TRUE)
    }, numeric(1))
    sd_ok <- sds[is.finite(sds) & sds > 0]
    addendum <- ""
    if (length(sd_ok) >= 1) {
      max_sd <- max(sd_ok); min_sd <- min(sd_ok); ratio <- if (length(sd_ok) >= 2) max_sd/min_sd else 1
      trigger <- (length(sd_ok) >= 2 && ratio > ratio_thresh) || any(sd_ok > abs_thresh)
      if (trigger) {
        sd_tbl <- paste(sprintf("  %s: %.3f", names(sds), sds), collapse = "\n")
        addendum <- paste0(
          "Potential cause for poor PSR: large scale differences among focal variables.\n",
          sprintf("  max SD = %.3f, min SD = %.3f, max/min = %.2f (thresholds: ratio > %g OR any SD > %g)\n",
                  max_sd, min_sd, ratio, ratio_thresh, abs_thresh),
          "Focal SDs:\n", sd_tbl, "\n\n",
          "Tip: Consider z-scaling y and the largest-SD predictors BEFORE calling EMAuxilary().\n"
        )
      }
    }
    warning(paste0(header, addendum), call. = FALSE)
  }

  if (!is.null(aux_model_arg) && exists("used_aux_route") && used_aux_route == "inlined-into-model") {
    model_plus_aux <- paste0(gsub("^MODEL:\n\\s*", "", model_line), "\n",
                             "auxiliary.model:\n  ", aux_model_arg, ";")
    blimp_code <- sub("MODEL:[\\s\\S]*?;", paste0("MODEL:\n  ", model_plus_aux), blimp_code)
  }

  structure(list(
    blimp_code = blimp_code,
    fit = fit,
    meta = list(y = y,
                center_group = cg,
                center_grand = center_grand,
                aux_scaled = aux_scaled)
  ), class = "blimp_wrap_fit")
}

print.blimp_wrap_fit <- function(x, ...) {
  cat("BLIMP code (conventional ML dialect):\n",
      "--------------------------------------\n", x$blimp_code, "\n\n", sep = "")
  cat("BLIMP output (first lines; use rblimp::output() for full tables)\n",
      "----------------------------------------------------------------\n", sep = "")
  out <- utils::capture.output(rblimp::output(x$fit))
  cat(paste(utils::head(out, 25), collapse = "\n"), "\n")
  invisible(x)
}

