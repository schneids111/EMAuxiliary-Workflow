# =========================================================================
# REPRODUCIBLE MONTE CARLO SIMULATION WITH BLIMP DIAGNOSTICS
# =========================================================================
# This script assumes that complete datasets (complete_data_*.csv)
# have already been generated (see 01_simulate_complete_datasets.R)

# Load necessary libraries
library(dplyr)
library(lme4)
library(rblimp)

# Set path to your local Blimp installation
blimp_path <- "C:/Program Files/Blimp/blimp.exe"
rblimp::set_blimp(blimp_path)
# Set your working directory here
# setwd("path/to/your/simulation/folder")

# -------------------------------------------------------------------------
# 1. PATHS AND SETTINGS
# -------------------------------------------------------------------------

# Newly generated files will be saved in this location
data_dir <- paste0(
  "C:/mypath"
)

output_dir <- file.path(data_dir, "monte_carlo")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results_file <- file.path(
  output_dir,
  "simulation_results_details.csv"
)

diagnostics_file <- file.path(
  output_dir,
  "blimp_diagnostics_details.csv"
)

# Dedicated temporary root.
# Use a drive with adequate free space if desired.
blimp_temp_root <- file.path(output_dir, "blimp_temp")
dir.create(blimp_temp_root, recursive = TRUE, showWarnings = FALSE)

n_sims <- 1000

blimp_burn <- 10000
blimp_iter <- 10000

# Base R seed controlling the induced missingness.
# Replication i will use base_seed + i.
base_seed <- 20260723

# BLIMP seed. This can remain fixed because each model is fitted to a
# different dataset/missingness realization. Alternatively, use a
# replication-specific BLIMP seed below.
base_blimp_seed <- 12345

# -------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# -------------------------------------------------------------------------

# Identify the final PSR row and summarize it.
#
# BLIMP documentation indicates:
#   - Rows are burn-in checkpoints.
#   - The final row represents the end of burn-in.
#   - Standardized coefficients and derived R2/SD quantities can sometimes
#     show elevated PSR values that are not the primary diagnostic focus.
summarize_psr <- function(model) {
  
  psr_table <- model@psr
  
  if (
    is.null(psr_table) ||
    nrow(psr_table) == 0L ||
    ncol(psr_table) == 0L
  ) {
    return(list(
      final_psr_max_all = NA_real_,
      final_psr_max_core = NA_real_,
      final_psr_parameter_all = NA_character_,
      final_psr_parameter_core = NA_character_,
      final_psr_focal_bw = NA_real_,
      final_psr_focal_bb = NA_real_
    ))
  }
  
  final_psr <- unlist(psr_table[nrow(psr_table), ], use.names = TRUE)
  final_psr <- as.numeric(final_psr) |>
    setNames(colnames(psr_table))
  
  # Exclude derived standardized coefficients, SDs, and R2 quantities.
  # Retain regression coefficients, intercepts, variances, and predictor
  # distribution parameters such as grand means.
  derived <- grepl(
    pattern = "standardized|r2:|\\bsd$",
    x = names(final_psr),
    ignore.case = TRUE
  )
  
  core_psr <- final_psr[!derived]
  
  max_all_index <- which.max(final_psr)
  max_core_index <- which.max(core_psr)
  
  # These two searches extract the focal outcome coefficients where present.
  bw_name <- grep(
    "^y(_mar)?\\.regressed on\\.x(_mar)?$",
    names(final_psr),
    value = TRUE
  )
  
  bb_name <- grep(
    "^y(_mar)?\\.regressed on\\.x(_mar)?\\.mean\\[id2\\]$",
    names(final_psr),
    value = TRUE
  )
  
  list(
    final_psr_max_all = unname(final_psr[max_all_index]),
    final_psr_max_core = unname(core_psr[max_core_index]),
    
    final_psr_parameter_all = names(final_psr)[max_all_index],
    final_psr_parameter_core = names(core_psr)[max_core_index],
    
    final_psr_focal_bw = if (length(bw_name) > 0L) {
      unname(final_psr[bw_name[1]])
    } else {
      NA_real_
    },
    
    final_psr_focal_bb = if (length(bb_name) > 0L) {
      unname(final_psr[bb_name[1]])
    } else {
      NA_real_
    }
  )
}


# Extract effective sample size from the BLIMP estimates table.
summarize_ess <- function(model) {
  
  est <- model@estimates
  
  if (is.null(est) || nrow(est) == 0L) {
    return(list(
      min_ess_all = NA_real_,
      min_ess_core = NA_real_,
      min_ess_parameter_all = NA_character_,
      min_ess_parameter_core = NA_character_
    ))
  }
  
  ess_column <- grep(
    pattern = "^N_?Eff$|effective",
    x = colnames(est),
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(ess_column) == 0L) {
    return(list(
      min_ess_all = NA_real_,
      min_ess_core = NA_real_,
      min_ess_parameter_all = NA_character_,
      min_ess_parameter_core = NA_character_
    ))
  }
  
  ess <- as.numeric(est[, ess_column[1]])
  names(ess) <- rownames(est)
  
  valid <- is.finite(ess)
  ess <- ess[valid]
  
  if (length(ess) == 0L) {
    return(list(
      min_ess_all = NA_real_,
      min_ess_core = NA_real_,
      min_ess_parameter_all = NA_character_,
      min_ess_parameter_core = NA_character_
    ))
  }
  
  derived <- grepl(
    pattern = "standardized|r2:|\\bsd$",
    x = names(ess),
    ignore.case = TRUE
  )
  
  core_ess <- ess[!derived]
  
  min_all_index <- which.min(ess)
  min_core_index <- which.min(core_ess)
  
  list(
    min_ess_all = unname(ess[min_all_index]),
    min_ess_core = unname(core_ess[min_core_index]),
    
    min_ess_parameter_all = names(ess)[min_all_index],
    min_ess_parameter_core = names(core_ess)[min_core_index]
  )
}


# Extract one estimate and its interval safely.
extract_blimp_parameter <- function(model, parameter) {
  
  if (!parameter %in% rownames(model@estimates)) {
    return(c(
      estimate = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      ess = NA_real_
    ))
  }
  
  tab <- model@estimates
  
  ess_column <- grep(
    "^N_?Eff$|effective",
    colnames(tab),
    ignore.case = TRUE,
    value = TRUE
  )
  
  c(
    estimate = tab[parameter, "Estimate"],
    lower = tab[parameter, "2.5%"],
    upper = tab[parameter, "97.5%"],
    ess = if (length(ess_column) > 0L) {
      tab[parameter, ess_column[1]]
    } else {
      NA_real_
    }
  )
}


# Fit one BLIMP model, extract everything needed, and then immediately
# remove both the object and its temporary directory.
fit_blimp_safely <- function(
    data,
    model_name,
    model_syntax,
    center_syntax = NULL,
    blimp_seed,
    focal_bw_name,
    focal_bb_name,
    variance_u0_name,
    variance_e_name
) {
  
  fit_dir <- tempfile(
    pattern = paste0("blimp_", model_name, "_"),
    tmpdir = blimp_temp_root
  )
  
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  
  old_wd <- getwd()
  
  start_time <- Sys.time()
  
  # Defaults returned if BLIMP produces an error.
  output <- list(
    model_name = model_name,
    status = "error",
    error_message = NA_character_,
    exitcode = NA_integer_,
    elapsed_seconds = NA_real_,
    
    b_0_est = NA_real_,
    b_w_est = NA_real_,
    b_b_est = NA_real_,
    b_w_lower = NA_real_,
    b_w_upper = NA_real_,
    b_b_lower = NA_real_,
    b_b_upper = NA_real_,
    var_u0_est = NA_real_,
    var_e_est = NA_real_,
    
    final_psr_max_all = NA_real_,
    final_psr_max_core = NA_real_,
    final_psr_parameter_all = NA_character_,
    final_psr_parameter_core = NA_character_,
    final_psr_focal_bw = NA_real_,
    final_psr_focal_bb = NA_real_,
    
    min_ess_all = NA_real_,
    min_ess_core = NA_real_,
    min_ess_parameter_all = NA_character_,
    min_ess_parameter_core = NA_character_,
    ess_focal_bw = NA_real_,
    ess_focal_bb = NA_real_
  )
  
  tryCatch({
    
    setwd(fit_dir)
    
    fit_args <- list(
      data = data,
      clusterid = "id2",
      model = model_syntax,
      seed = blimp_seed,
      burn = blimp_burn,
      iter = blimp_iter
    )
    
    if (!is.null(center_syntax)) {
      fit_args$center <- center_syntax
    }
    
    fit <- do.call(rblimp, fit_args)
    
    psr_info <- summarize_psr(fit)
    ess_info <- summarize_ess(fit)
    
    bw <- extract_blimp_parameter(fit, focal_bw_name)
    bb <- extract_blimp_parameter(fit, focal_bb_name)
    vu <- extract_blimp_parameter(fit, variance_u0_name)
    ve <- extract_blimp_parameter(fit, variance_e_name)
    
    intercept_name <- grep(
      pattern = "^(y|y_mar) ~ Intercept$",
      x = rownames(fit@estimates),
      value = TRUE
    )
    
    intercept <- if (length(intercept_name) > 0L) {
      fit@estimates[intercept_name[1], "Estimate"]
    } else {
      NA_real_
    }
    
    exitcode <- attr(fit@output, "exitcode")
    
    output <- modifyList(
      output,
      c(
        list(
          status = if (
            !is.null(exitcode) &&
            identical(as.integer(exitcode), 0L)
          ) {
            "completed"
          } else {
            "nonzero_exit"
          },
          
          exitcode = if (is.null(exitcode)) {
            NA_integer_
          } else {
            as.integer(exitcode)
          },
          
          b_0_est = unname(intercept),
          b_w_est = unname(bw["estimate"]),
          b_b_est = unname(bb["estimate"]),
          b_w_lower = unname(bw["lower"]),
          b_w_upper = unname(bw["upper"]),
          b_b_lower = unname(bb["lower"]),
          b_b_upper = unname(bb["upper"]),
          var_u0_est = unname(vu["estimate"]),
          var_e_est = unname(ve["estimate"]),
          ess_focal_bw = unname(bw["ess"]),
          ess_focal_bb = unname(bb["ess"])
        ),
        psr_info,
        ess_info
      )
    )
    
    # Delete the potentially very large object before the next fit.
    rm(fit)
    invisible(gc(verbose = FALSE))
    
  }, error = function(e) {
    
    output$error_message <<- conditionMessage(e)
    
  }, finally = {
    
    setwd(old_wd)
    
    # Repeated attempts help on Windows if a process briefly retains a handle.
    for (attempt in 1:5) {
      if (!dir.exists(fit_dir)) break
      
      unlink(
        fit_dir,
        recursive = TRUE,
        force = TRUE
      )
      
      if (dir.exists(fit_dir)) {
        Sys.sleep(1)
      }
    }
    
    invisible(gc(verbose = FALSE))
  })
  
  output$elapsed_seconds <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  )
  
  as.data.frame(output, stringsAsFactors = FALSE)
}


# Coverage helper.
covered <- function(lower, upper, truth) {
  as.integer(
    is.finite(lower) &&
      is.finite(upper) &&
      lower < truth &&
      upper > truth
  )
}

# -------------------------------------------------------------------------
# 3. RESUME SUPPORT
# -------------------------------------------------------------------------

if (file.exists(results_file)) {
  
  results <- read.csv(
    results_file,
    stringsAsFactors = FALSE
  )
  
  completed_sims <- unique(results$sim_id)
  
  message(
    "Resuming. Completed replications: ",
    length(completed_sims)
  )
  
} else {
  
  results <- data.frame()
  completed_sims <- integer()
}

if (file.exists(diagnostics_file)) {
  diagnostics <- read.csv(
    diagnostics_file,
    stringsAsFactors = FALSE
  )
} else {
  diagnostics <- data.frame()
}

# -------------------------------------------------------------------------
# 4. SIMULATION LOOP
# -------------------------------------------------------------------------

for (i in seq_len(n_sims)) {
  
  if (i %in% completed_sims) {
    next
  }
  
  message(
    "\n========================================",
    "\nStarting replication ", i, " of ", n_sims,
    "\n========================================"
  )
  
  file_name <- file.path(
    data_dir,
    paste0("complete_data_", i, ".csv")
  )
  
  if (!file.exists(file_name)) {
    warning("Missing input file: ", file_name)
    next
  }
  
  dat <- read.csv(file_name)
  
  # Deterministic replication-specific missingness.
  set.seed(base_seed + i)
  
  dat <- dat %>%
    mutate(
      prob_missing = plogis(-0.5 + 2.5 * aux),
      is_missing = rbinom(n(), 1, prob_missing),
      y_mar = ifelse(is_missing == 1, NA, y),
      x_mar = ifelse(is_missing == 1, NA, x)
    ) %>%
    group_by(id2) %>%
    mutate(
      x_mean_true = mean(x),
      x_cwc_true = x - x_mean_true,
      
      x_mean_biased = mean(x_mar, na.rm = TRUE),
      x_cwc_biased = x_mar - x_mean_biased
    ) %>%
    ungroup()
  
  # -----------------------------------------------------------------------
  # MODEL 1: lmer complete
  # -----------------------------------------------------------------------
  
  m1 <- lmer(
    y ~ x_cwc_true + x_mean_true + (1 | id2),
    data = dat
  )
  
  m1_sum <- summary(m1)$coefficients
  m1_vc <- as.data.frame(VarCorr(m1))
  
  m1_row <- data.frame(
    sim_id = i,
    model_type = "1_lmer_complete",
    
    b_0_est = m1_sum["(Intercept)", "Estimate"],
    b_w_est = m1_sum["x_cwc_true", "Estimate"],
    b_b_est = m1_sum["x_mean_true", "Estimate"],
    
    cov_b_w = covered(
      m1_sum["x_cwc_true", "Estimate"] -
        1.96 * m1_sum["x_cwc_true", "Std. Error"],
      m1_sum["x_cwc_true", "Estimate"] +
        1.96 * m1_sum["x_cwc_true", "Std. Error"],
      0.4
    ),
    
    cov_b_b = covered(
      m1_sum["x_mean_true", "Estimate"] -
        1.96 * m1_sum["x_mean_true", "Std. Error"],
      m1_sum["x_mean_true", "Estimate"] +
        1.96 * m1_sum["x_mean_true", "Std. Error"],
      0.6
    ),
    
    var_u0_est = m1_vc[
      m1_vc$grp == "id2",
      "vcov"
    ],
    
    var_e_est = m1_vc[
      m1_vc$grp == "Residual",
      "vcov"
    ]
  )
  
  rm(m1)
  invisible(gc(verbose = FALSE))
  
  # -----------------------------------------------------------------------
  # MODEL 2: BLIMP complete
  # -----------------------------------------------------------------------
  
  m2 <- fit_blimp_safely(
    data = dat,
    model_name = paste0("m2_sim_", i),
    
    center_syntax = "groupmean = x",
    
    model_syntax = "
      y ~ x x.mean;
    ",
    
    blimp_seed = base_blimp_seed + i * 10L + 2L,
    
    focal_bw_name = "y ~ x",
    focal_bb_name = "y ~ x.mean[id2]",
    variance_u0_name = "y level-2 intercept variance",
    variance_e_name = "y residual variance"
  )
  
  # -----------------------------------------------------------------------
  # MODEL 3: lmer naive
  # -----------------------------------------------------------------------
  
  m3 <- lmer(
    y_mar ~ x_cwc_biased + x_mean_biased + (1 | id2),
    data = dat
  )
  
  m3_sum <- summary(m3)$coefficients
  m3_vc <- as.data.frame(VarCorr(m3))
  
  m3_row <- data.frame(
    sim_id = i,
    model_type = "3_lmer_naive",
    
    b_0_est = m3_sum["(Intercept)", "Estimate"],
    b_w_est = m3_sum["x_cwc_biased", "Estimate"],
    b_b_est = m3_sum["x_mean_biased", "Estimate"],
    
    cov_b_w = covered(
      m3_sum["x_cwc_biased", "Estimate"] -
        1.96 * m3_sum["x_cwc_biased", "Std. Error"],
      m3_sum["x_cwc_biased", "Estimate"] +
        1.96 * m3_sum["x_cwc_biased", "Std. Error"],
      0.4
    ),
    
    cov_b_b = covered(
      m3_sum["x_mean_biased", "Estimate"] -
        1.96 * m3_sum["x_mean_biased", "Std. Error"],
      m3_sum["x_mean_biased", "Estimate"] +
        1.96 * m3_sum["x_mean_biased", "Std. Error"],
      0.6
    ),
    
    var_u0_est = m3_vc[
      m3_vc$grp == "id2",
      "vcov"
    ],
    
    var_e_est = m3_vc[
      m3_vc$grp == "Residual",
      "vcov"
    ]
  )
  
  rm(m3)
  invisible(gc(verbose = FALSE))
  
  # -----------------------------------------------------------------------
  # MODEL 4: BLIMP latent mean centering, no auxiliary
  # -----------------------------------------------------------------------
  
  m4 <- fit_blimp_safely(
    data = dat,
    model_name = paste0("m4_sim_", i),
    
    center_syntax = "groupmean = x_mar",
    
    model_syntax = "
      y_mar ~ x_mar x_mar.mean;
    ",
    
    blimp_seed = base_blimp_seed + i * 10L + 4L,
    
    focal_bw_name = "y_mar ~ x_mar",
    focal_bb_name = "y_mar ~ x_mar.mean[id2]",
    variance_u0_name = "y_mar level-2 intercept variance",
    variance_e_name = "y_mar residual variance"
  )
  
  # -----------------------------------------------------------------------
  # MODEL 5: BLIMP observed mean centering with auxiliary
  # -----------------------------------------------------------------------
  
  m5 <- fit_blimp_safely(
    data = dat,
    model_name = paste0("m5_sim_", i),
    
    model_syntax = "
      y_mar ~ x_cwc_biased x_mean_biased;
      aux ~ y_mar x_cwc_biased x_mean_biased;
    ",
    
    blimp_seed = base_blimp_seed + i * 10L + 5L,
    
    focal_bw_name = "y_mar ~ x_cwc_biased",
    focal_bb_name = "y_mar ~ x_mean_biased",
    variance_u0_name = "y_mar level-2 intercept variance",
    variance_e_name = "y_mar residual variance"
  )
  
  # -----------------------------------------------------------------------
  # MODEL 6: BLIMP latent mean centering with auxiliary
  # -----------------------------------------------------------------------
  
  m6 <- fit_blimp_safely(
    data = dat,
    model_name = paste0("m6_sim_", i),
    
    center_syntax = "groupmean = x_mar",
    
    model_syntax = "
      y_mar ~ x_mar x_mar.mean;
      aux ~ y_mar y_mar.mean x_mar x_mar.mean;
    ",
    
    blimp_seed = base_blimp_seed + i * 10L + 6L,
    
    focal_bw_name = "y_mar ~ x_mar",
    focal_bb_name = "y_mar ~ x_mar.mean[id2]",
    variance_u0_name = "y_mar level-2 intercept variance",
    variance_e_name = "y_mar residual variance"
  )
  
  # -----------------------------------------------------------------------
  # CONVERT BLIMP FITS TO ANALYSIS ROWS
  # -----------------------------------------------------------------------
  
  make_blimp_result_row <- function(x, model_type) {
    
    data.frame(
      sim_id = i,
      model_type = model_type,
      
      b_0_est = x$b_0_est,
      b_w_est = x$b_w_est,
      b_b_est = x$b_b_est,
      
      cov_b_w = covered(
        x$b_w_lower,
        x$b_w_upper,
        0.4
      ),
      
      cov_b_b = covered(
        x$b_b_lower,
        x$b_b_upper,
        0.6
      ),
      
      var_u0_est = x$var_u0_est,
      var_e_est = x$var_e_est
    )
  }
  
  simulation_rows <- bind_rows(
    m1_row,
    make_blimp_result_row(m2, "2_blimp_complete"),
    m3_row,
    make_blimp_result_row(m4, "4_blimp_lmc_noaux"),
    make_blimp_result_row(m5, "5_blimp_omc_aux"),
    make_blimp_result_row(m6, "6_blimp_rescue")
  )
  
  diagnostic_rows <- bind_rows(
    transform(
      m2,
      sim_id = i,
      model_type = "2_blimp_complete"
    ),
    transform(
      m4,
      sim_id = i,
      model_type = "4_blimp_lmc_noaux"
    ),
    transform(
      m5,
      sim_id = i,
      model_type = "5_blimp_omc_aux"
    ),
    transform(
      m6,
      sim_id = i,
      model_type = "6_blimp_rescue"
    )
  )
  
  results <- bind_rows(results, simulation_rows)
  diagnostics <- bind_rows(diagnostics, diagnostic_rows)
  
  # Save after every complete replication.
  write.csv(
    results,
    results_file,
    row.names = FALSE
  )
  
  write.csv(
    diagnostics,
    diagnostics_file,
    row.names = FALSE
  )
  
  rm(
    dat,
    m1_row,
    m2,
    m3_row,
    m4,
    m5,
    m6,
    simulation_rows,
    diagnostic_rows
  )
  
  invisible(gc(verbose = FALSE))
  
  message(
    "Completed replication ", i,
    "; available disk space should now be stable."
  )
}

# -------------------------------------------------------------------------
# 5. SUMMARIZE SIMULATION RESULTS
# -------------------------------------------------------------------------

summary_table <- results %>%
  group_by(model_type) %>%
  summarize(
    n_replications = sum(is.finite(b_w_est)),
    
    b_w_mean = mean(b_w_est, na.rm = TRUE),
    b_w_abs_bias = b_w_mean - 0.4,
    b_w_emp_se = sd(b_w_est, na.rm = TRUE),
    b_w_std_bias = b_w_abs_bias / b_w_emp_se,
    b_w_coverage = mean(cov_b_w, na.rm = TRUE),
    
    b_b_mean = mean(b_b_est, na.rm = TRUE),
    b_b_abs_bias = b_b_mean - 0.6,
    b_b_emp_se = sd(b_b_est, na.rm = TRUE),
    b_b_std_bias = b_b_abs_bias / b_b_emp_se,
    b_b_coverage = mean(cov_b_b, na.rm = TRUE),
    
    var_u0_mean = mean(var_u0_est, na.rm = TRUE),
    var_u0_abs_bias = var_u0_mean - 1.0,
    var_u0_std_bias =
      var_u0_abs_bias / sd(var_u0_est, na.rm = TRUE),
    
    var_e_mean = mean(var_e_est, na.rm = TRUE),
    var_e_abs_bias = var_e_mean - 1.0,
    var_e_std_bias =
      var_e_abs_bias / sd(var_e_est, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )

write.csv(
  summary_table,
  file.path(output_dir, "simulation_summary.csv"),
  row.names = FALSE
)

print(t(summary_table))

# -------------------------------------------------------------------------
# 6. SUMMARIZE BLIMP DIAGNOSTICS
# -------------------------------------------------------------------------

diagnostic_summary <- diagnostics %>%
  group_by(model_type) %>%
  summarize(
    n_fits = n(),
    
    n_completed = sum(status == "completed", na.rm = TRUE),
    n_errors = sum(status == "error", na.rm = TRUE),
    n_nonzero_exit = sum(
      status == "nonzero_exit",
      na.rm = TRUE
    ),
    
    median_final_psr_core = median(
      final_psr_max_core,
      na.rm = TRUE
    ),
    
    maximum_final_psr_core = max(
      final_psr_max_core,
      na.rm = TRUE
    ),
    
    proportion_final_psr_core_below_1_05 = mean(
      final_psr_max_core < 1.05,
      na.rm = TRUE
    ),
    
    proportion_final_psr_core_below_1_10 = mean(
      final_psr_max_core < 1.10,
      na.rm = TRUE
    ),
    
    median_min_ess_core = median(
      min_ess_core,
      na.rm = TRUE
    ),
    
    minimum_ess_core = min(
      min_ess_core,
      na.rm = TRUE
    ),
    
    proportion_min_ess_core_above_100 = mean(
      min_ess_core > 100,
      na.rm = TRUE
    ),
    
    median_elapsed_minutes = median(
      elapsed_seconds / 60,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

write.csv(
  diagnostic_summary,
  file.path(output_dir, "diagnostic_summary.csv"),
  row.names = FALSE
)

print(diagnostic_summary)