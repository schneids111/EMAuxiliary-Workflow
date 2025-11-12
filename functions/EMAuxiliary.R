# EMAuxiliary.R
# -------------------------------------------------------------------
# BLIMP wrapper for lme4-style formulas using explicit `.mean` tokens.
# - No w()/b() markers. You write x and/or x.mean explicitly.
# - Example focal model: y ~ x + x.mean + x:x.mean + (1 + x | id)
# - To *use* x.mean, you must also list the base in center_group (e.g., center_group="x").
# - Auxiliaries are standardized (z-scored) if numeric.
# - Aux RHS follows rules described in the comments below.
# - Adds: chains=, nominal=, center_group=, center_grand=
# - Provides: blimp_print_psr(), blimp_print_focal()
# -------------------------------------------------------------------

EMAuxiliary <- function(formula,
                       data,
                       id,
                       aux = NULL,
                       ordinal = NULL,
                       nominal = NULL,
                       center_group = NULL,   # e.g., c("x","time","y") to make x.mean, time.mean, y.mean available
                       center_grand = NULL,   # e.g., c("x.mean","time.mean") to GM-center those between parts
                       burn = 5000,
                       iter = 10000,
                       chains = 2,
                       seed = 12345) {
  stopifnot(is.data.frame(data))
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Install lme4.")
  if (!requireNamespace("rblimp", quietly = TRUE)) stop("Install rblimp.")
  
  # ----------------------- small helpers ----------------------------
  .is_num <- function(x) is.numeric(x) && !is.factor(x)
  .safe_sd <- function(v) if (is.null(v) || !is.numeric(v)) NA_real_ else stats::sd(v, na.rm = TRUE)
  .strip_space <- function(x) x[nzchar(trimws(x))]
  
  # escape a variable name for regex in blimp_print_focal
  .re_esc <- function(s) gsub("([][{}()+*.^$|?\\\\])", "\\\\\\1", s, perl = TRUE)
  
  # detect L2-only variables (constant within clusters)
  .is_L2 <- function(vname) {
    if (!vname %in% names(data)) return(FALSE)
    sds <- tapply(data[[vname]], data[[id]], function(z) stats::sd(z, na.rm = TRUE))
    all(is.na(sds) | sds < 1e-8)
  }
  
  # quick ICC for auxiliaries (returns NA if not computable)
  .icc <- function(vname) {
    if (!vname %in% names(data)) return(NA_real_)
    z <- data[[vname]]
    grp <- data[[id]]
    if (!is.numeric(z)) return(NA_real_)
    m_j <- tapply(z, grp, mean, na.rm = TRUE)
    v_b <- stats::var(m_j, na.rm = TRUE)
    v_w <- mean(tapply(z, grp, function(a) stats::var(a, na.rm = TRUE)), na.rm = TRUE)
    if (is.na(v_b) || is.na(v_w)) return(NA_real_)
    if ((v_b + v_w) <= 0) return(NA_real_)
    v_b / (v_b + v_w)
  }
  
  # ------------------ parse focal model (no random parts) -----------
  fixed_form <- lme4::nobars(formula)
  tl_raw <- attr(terms(fixed_form), "term.labels")  # tokens like c("x", "x.mean", "x:x.mean", ...)
  tl_raw <- if (length(tl_raw)) tl_raw else character(0)
  
  # Build BLIMP tokens by replacing ":" with "*"
  blimp_terms <- if (length(tl_raw)) gsub(":", "*", tl_raw, fixed = TRUE) else character(0)
  rhs_tokens <- unique(blimp_terms)
  fixed_rhs  <- if (length(rhs_tokens)) paste(rhs_tokens, collapse = " ") else "1"
  
  # outcome
  y <- deparse(formula[[2L]])
  
  # -------------------- centering logic (explicit) ------------------
  # Users must request cluster means for any base that appears as "base.mean" in the formula.
  # We do NOT silently add those bases (to avoid altering the user's raw x scaling).
  # However, if aux != NULL, we *do* ensure y.mean is available by adding y to group-mean list.
  
  base_names_in_means <- {
    m <- grep("\\.mean\\b", rhs_tokens, value = TRUE)
    if (length(m)) unique(sub("\\.mean$", "", m)) else character(0)
  }
  
  cg <- unique(.strip_space(center_group))
  gg <- unique(.strip_space(center_grand))
  
  # Require bases for any *.mean that appear in the focal RHS
  missing_bases <- setdiff(base_names_in_means, cg)
  if (length(missing_bases)) {
    stop(
      "You used '", paste0(missing_bases, collapse = ", "),
      ".mean' in the formula but did not include the base name(s) in `center_group`.\n",
      "Add: center_group = c(", paste(sprintf('"%s"', missing_bases), collapse = ", "), ", ...).",
      call. = FALSE
    )
  }
  
  # Ensure y.mean exists for auxiliary equations (only if aux are provided).
  if (!is.null(aux) && length(aux)) {
    if (!y %in% cg) cg <- c(cg, y)
  }
  
  # Grand-mean centering: user should pass tokens like "x.mean"
  # Validate that anything in center_grand actually ends with ".mean"
  wrong_gg <- gg[!grepl("\\.mean$", gg)]
  if (length(wrong_gg)) {
    stop("`center_grand` must list between tokens that end with '.mean'. Bad value(s): ",
         paste(wrong_gg, collapse = ", "), call. = FALSE)
  }
  
  # Render CENTER block
  center_lines <- character(0)
  if (length(cg)) center_lines <- c(center_lines, paste0("groupmean = ", paste(unique(cg), collapse = " "), ";"))
  if (length(gg)) center_lines <- c(center_lines, paste0("grandmean = ", paste(unique(gg), collapse = " "), ";"))
  center_block <- if (length(center_lines)) paste("CENTER:", paste(center_lines, collapse = " "), sep = "\n  ") else NULL
  
  # -------------------- random slopes (from original formula) -------
  fb <- lme4::findbars(formula)
  rand_slopes <- character(0)
  if (length(fb)) {
    inside_terms <- unique(unlist(lapply(fb, function(b) {
      left <- deparse(b[[2]])
      attr(terms(as.formula(paste("~", left))), "term.labels")
    })))
    rand_slopes <- setdiff(inside_terms, "1")
  }
  
  # Disallow ordinal/nominal as random-slope predictors
  bad_rs <- intersect(rand_slopes, c(ordinal, nominal))
  if (length(bad_rs)) {
    stop("Categorical (ordinal/nominal) variables cannot be random-slope predictors: ",
         paste(bad_rs, collapse = ", "), call. = FALSE)
  }
  
  # Assemble focal model line
  if (!length(rand_slopes)) {
    model_line <- paste0("MODEL:\n  ", y, " ~ ", fixed_rhs, ";")
    blimp_model_arg <- paste0(y, " ~ ", fixed_rhs)
  } else {
    model_line <- paste0("MODEL:\n  ", y, " ~ ", fixed_rhs, " | ", paste(rand_slopes, collapse = " "), ";")
    blimp_model_arg <- paste0(y, " ~ ", fixed_rhs, " | ", paste(rand_slopes, collapse = " "))
  }
  
  # --------------------- standardize numeric auxiliaries ------------
  data_work <- data
  if (!is.null(aux) && length(aux)) {
    for (a in aux) {
      if (a %in% names(data_work) && .is_num(data_work[[a]])) {
        data_work[[a]] <- as.numeric(scale(data_work[[a]], center = TRUE, scale = TRUE)[,1])
      }
    }
  }
  
  # ------------------ build auxiliary.model block -------------------
  aux_block <- NULL
  aux_model_arg <- NULL
  
  if (!is.null(aux) && length(aux)) {
    # We need:
    # - main "base" predictors (no interactions) as raw tokens
    # - full RHS tokens (include interactions, may contain *.mean)
    # - classification of auxiliaries via ICC
    main_bases <- unique(unlist(strsplit(tl_raw, ":", fixed = TRUE)))
    main_bases <- main_bases[nzchar(main_bases)]
    
    full_rhs <- rhs_tokens
    
    # Helper: choose mean token for a base var if needed in between-level RHS
    # If predictor is L2-only -> use raw name, else use base.mean
    mean_token <- function(xname) {
      if (.is_L2(xname)) xname else paste0(xname, ".mean")
    }
    
    # Rule buckets:
    # - pure-between aux: ICC ~ 1.0 (we'll use .is_L2() as a stronger condition)
    # - near-pure within: ICC < .05  => use y + main_bases (no *.mean)
    # - multi-level: otherwise => use y + y.mean + full RHS tokens
    
    aux_lines <- character(0)
    
    for (a in aux) {
      # Skip if aux not in data
      if (!a %in% names(data_work)) next
      
      a_is_L2 <- .is_L2(a)
      a_icc   <- .icc(a)
      
      if (isTRUE(a_is_L2)) {
        # PURE BETWEEN auxiliary:
        #   RHS: y.mean + main effects, with predictors mapped to *.mean unless they are pure L2 themselves
        rhs_means <- unique(c(paste0(y, ".mean"),
                              vapply(main_bases, mean_token, character(1))))
        rhs <- rhs_means
      } else if (is.finite(a_icc) && a_icc < 0.05) {
        # NEAR-PURE WITHIN auxiliary:
        #   RHS: y + main effects (raw), no *.mean even if present in focal
        rhs <- unique(c(y, main_bases))
      } else {
        # MULTI-LEVEL auxiliary:
        #   RHS: y + y.mean + FULL RHS tokens (including interactions and any *.mean the user requested)
        rhs <- unique(c(y, paste0(y, ".mean"), full_rhs))
      }
      
      # Drop self-prediction and any empties
      rhs <- rhs[!is.na(rhs) & nzchar(rhs) & rhs != a]
      
      if (length(rhs)) {
        aux_lines <- c(aux_lines, paste0(a, " ~ ", paste(rhs, collapse = " "), ";"))
      }
    }
    
    if (length(aux_lines)) {
      aux_block <- paste0("auxiliary.model:\n  ", paste(aux_lines, collapse = "\n  "))
      aux_model_arg <- paste(aux_lines, collapse = "\n  ")
    }
  }
  
  # ------------------------- VARIABLES block ------------------------
  # Include id, y, every raw base token (split on ':'), any *.mean tokens explicitly used,
  # plus any aux/ordinal/nominal variables.
  raw_bases <- unique(unlist(strsplit(tl_raw, ":", fixed = TRUE)))
  raw_bases <- raw_bases[nzchar(raw_bases)]
  # Collect explicit .mean tokens from formula (user responsibility)
  explicit_means <- unique(grep("\\.mean\\b", tl_raw, value = TRUE))
  explicit_means <- unique(sub("^(.+)$", "\\1", explicit_means)) # keep as written
  
  # VARIABLES must list observed names (no *.mean here).
  # So strip ".mean" for inclusion in VARIABLES (they are constructed by CENTER).
  vars_needed <- unique(c(
    id, y,
    raw_bases,
    sub("\\.mean$", "", explicit_means),
    aux, ordinal, nominal
  ))
  vars_needed <- vars_needed[vars_needed %in% names(data_work)]
  vars_needed <- .strip_space(vars_needed)
  
  blimp_code <- paste(
    "DATA: <in-memory by rblimp>;",
    paste0("VARIABLES: ", paste(vars_needed, collapse = " "), ";"),
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
  
  # ----------------------------- run BLIMP --------------------------
  args <- list(
    data = data_work,
    clusterid = id,
    model = blimp_model_arg,
    seed = seed,
    burn = burn,
    iter = iter,
    chains = chains
  )
  if (length(cg)) {
    args$center <- paste0(
      paste0("groupmean = ", paste(unique(cg), collapse = " ")),
      if (length(gg)) paste0("; grandmean = ", paste(unique(gg), collapse = " ")) else ""
    )
  } else if (length(gg)) {
    # grandmean without groupmean is allowed (for tokens already created by user elsewhere),
    # but in BLIMP practice, *.mean come from groupmean. We'll pass anyway if user insists.
    args$center <- paste0("grandmean = ", paste(unique(gg), collapse = " "))
  }
  fm <- names(formals(rblimp::rblimp))
  if (!is.null(ordinal) && length(ordinal) && "ordinal" %in% fm) args$ordinal <- paste(ordinal, collapse = " ")
  if (!is.null(nominal) && length(nominal) && "nominal" %in% fm) args$nominal <- paste(nominal, collapse = " ")
  if (!is.null(aux_model_arg)) {
    aux_names <- c("auxiliary.model", "auxmodel", "auxiliarymodel", "auxiliary")
    matched <- intersect(aux_names, fm)
    if (length(matched)) {
      args[[matched[1]]] <- aux_model_arg
    } else {
      # inline fallback
      args$model <- paste0(blimp_model_arg, ";\nauxiliary.model:\n  ", aux_model_arg, ";")
    }
  }
  
  fit <- do.call(rblimp::rblimp, args)
  
  # ---------------------- PSR + scale-audit warning -----------------
  psr <- tryCatch(blimp_last_psr(fit), error = function(e) NA_real_)
  if (is.finite(psr) && psr > 1.05) {
    # scale audit across focal *observed* columns (y and main bases)
    base_cols <- unique(.strip_space(unlist(strsplit(tl_raw, ":", fixed = TRUE))))
    focal_cols <- unique(c(y, base_cols))
    sds <- vapply(focal_cols, function(v) .safe_sd(data[[v]]), numeric(1))
    sd_ok <- sds[is.finite(sds) & sds > 0]
    ratio_thresh <- 6
    abs_thresh   <- 3
    
    msg <- sprintf("Final PSR is %.3f. The final PSR should be < 1.05.\n", psr)
    add <- ""
    if (length(sd_ok) >= 1) {
      max_sd <- max(sd_ok); min_sd <- min(sd_ok); ratio <- if (length(sd_ok) >= 2) max_sd / min_sd else 1
      trigger <- (length(sd_ok) >= 2 && ratio > ratio_thresh) || any(sd_ok > abs_thresh)
      if (trigger) {
        sd_tbl <- paste(sprintf("  %s: %.3f", names(sds), sds), collapse = "\n")
        add <- paste0(
          "Potential cause for poor PSR: large scale differences among focal variables.\n",
          sprintf("  max SD = %.3f, min SD = %.3f, max/min = %.2f (thresholds: ratio > %g OR any SD > %g)\n",
                  max_sd, min_sd, ratio, ratio_thresh, abs_thresh),
          "Focal SDs:\n", sd_tbl, "\n\n",
          "Tip: Consider z-scaling y and large-SD predictors BEFORE calling blimp_wrap().\n"
        )
      }
    }
    warning(paste0(msg, add, "Consider longer burn/iter, more chains, or simplifying the model."), call. = FALSE)
  }
  
  # ------------------------------ return ----------------------------
  structure(list(
    blimp_code = blimp_code,
    fit = fit,
    meta = list(y = y, id = id)
  ), class = "blimp_wrap_fit")
}

# ------------------------- Printers/helpers -------------------------

print.blimp_wrap_fit <- function(x, ...) {
  cat("BLIMP code (explicit .mean idiom):\n")
  cat("----------------------------------\n")
  cat(x$blimp_code, "\n\n")
  cat("BLIMP output (first lines; use rblimp::output() for full tables)\n")
  cat("----------------------------------------------------------------\n")
  out <- utils::capture.output(rblimp::output(x$fit))
  cat(paste(head(out, 25), collapse = "\n"), "\n")
  invisible(x)
}

# Extract the last reported PSR row's final value (highest PSR in the last block)
blimp_last_psr <- function(fit_obj) {
  out <- utils::capture.output(rblimp::output(fit_obj))
  idx <- grep("^\\s*5001 to 10000|^\\s*\\d+ to \\d+\\s+\\d+\\.\\d+\\s+\\d+\\s*$", out)
  if (!length(idx)) {
    # fallback: find all PSR rows and take the last numeric
    idx <- grep("\\bHighest PSR\\b", out)
    if (!length(idx)) return(NA_real_)
    line <- out[max(idx)]
  } else {
    line <- out[max(idx)]
  }
  # pull the first number that looks like a PSR (e.g., 1.029)
  num <- as.numeric(regmatches(line, regexpr("\\d+\\.\\d+", line)))
  num
}

# Print only the PSR block (burn-in PSR through just before "METROPOLIS-HASTINGS...")
blimp_print_psr <- function(obj) {
  # Accept wrapper or raw blimp_obj
  bl <- if (inherits(obj, "blimp_wrap_fit")) obj$fit else obj
  out <- utils::capture.output(rblimp::output(bl))
  
  s <- grep("^\\s*BURN-IN POTENTIAL SCALE REDUCTION", out, ignore.case = TRUE)
  e <- grep("^\\s*METROPOLIS-HASTINGS ACCEPTANCE RATES", out, ignore.case = TRUE)
  
  if (!length(s)) {
    cat("PSR section not found.\n")
    return(invisible())
  }
  if (!length(e)) e <- length(out) + 1L  # print to end if we don't find MH header
  
  cat(paste(out[s:(e[1] - 1L)], collapse = "\n"), "\n")
}

# Print focal model block only (from "Outcome Variable: <y>" up to before auxiliary.model block
# or before the next "Outcome Variable:" heading, whichever comes first)
blimp_print_focal <- function(obj) {
  # Accept wrapper or raw blimp_obj (needs y if wrapper)
  if (inherits(obj, "blimp_wrap_fit")) {
    bl <- obj$fit
    y  <- obj$meta$y
  } else {
    stop("Pass the wrapper object (class 'blimp_wrap_fit') so I know the focal outcome name.")
  }
  
  out <- utils::capture.output(rblimp::output(bl))
  y_esc <- gsub("([][{}()+*.^$|?\\\\])", "\\\\\\1", y, perl = TRUE)
  
  # Start at focal outcome header
  start <- grep(paste0("^\\s*Outcome Variable:\\s+", y_esc, "\\b"), out)
  if (!length(start)) {
    cat("Focal section not found for outcome '", y, "'.\n", sep = "")
    return(invisible())
  }
  start <- start[1]
  
  # Candidate end markers (first one after 'start' wins)
  aux_hdr       <- grep("^\\s*auxiliary\\.model block\\s*:\\s*$", out, ignore.case = TRUE)
  next_outcome  <- setdiff(grep("^\\s*Outcome Variable:\\s+", out, ignore.case = TRUE), start)
  next_outcome  <- next_outcome[next_outcome > start]
  predictor_hdr <- grep("^\\s*PREDICTOR MODEL ESTIMATES\\s*:\\s*$", out, ignore.case = TRUE)
  
  cand <- c(aux_hdr, next_outcome, predictor_hdr, length(out) + 1L)
  cand <- cand[cand > start]
  end  <- if (length(cand)) min(cand) - 1L else length(out)
  
  cat(paste(out[start:end], collapse = "\n"), "\n")
}
