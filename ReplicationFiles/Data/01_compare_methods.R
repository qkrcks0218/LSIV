## Batch jobs: one cross-fitting repetition of the Job Corps analysis per job.
## BATCH is the seed. Merge the jobs with 02_density_ratio.R and 03_summary.R.

args <- (commandArgs(trailingOnly=TRUE))
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

library(SIV)
source("00_common.R")

## SIV: ATT, QTT(0.25, 0.5, 0.75) and the falsification test.
## Competing methods: ATT only.
## One repetition per job; SIV()'s repeated cross-fitting median adjustment is
## applied across jobs by 03_summary.R.

SL_LIST <- c(1:9)
NUM_CROSSFIT <- 1L
NUM_FOLD <- 2L
QTT_QUANTILE <- c(0.25, 0.5, 0.75)
BW_QUANTILE <- c(0.5, 0.75, 0.9, 0.95)
NUM_BOOT <- 1000L

SEED <- as.integer(BATCH)
d <- jc_load_data("JC.csv")

message(sprintf("BATCH=%d: N=%d, SL.list=%s, %d fold(s)",
                SEED, length(d$Y), paste(SL_LIST, collapse = ","), NUM_FOLD))

## ---------------------------------------------------------------------------
## SIV
## ---------------------------------------------------------------------------

set.seed(SEED)
siv_run <- jc_capture(SIV::SIV(
  d$Y, d$A, d$Z, d$X,
  estimand = c("ATT", "QTT"),
  level = 0.95,
  Overall.Estimation_SL.list = SL_LIST,
  Overall.Estimation_num.crossfit = NUM_CROSSFIT,
  Overall.Estimation_num.fold = NUM_FOLD,
  y_type = "continuous",
  y_grid = JC_Y_GRID,
  y_weights = JC_Y_WEIGHTS,
  y_ref = 0,
  outcome_density = "np_gss_cv",
  ratio_bandwidth = "cv_separate",
  ratio_bandwidth.cv.grid = BW_QUANTILE,
  QTT_quantile = QTT_QUANTILE,
  QTT_support = JC_Y_GRID,
  QTT_bounds = c(0, 52),
  falsification = TRUE,
  falsification_num.boot = NUM_BOOT,
  falsification_num.bin = c(20L, rep.int(5L, ncol(d$X))),
  falsification_categorical = colnames(d$X)
))
fit <- siv_run$value
fit

siv_estimands <- c("ATT", paste0("QTT_", vapply(QTT_QUANTILE, format,
                                                character(1L), trim = TRUE,
                                                scientific = FALSE)))
if (is.null(fit)) {
  siv_rows <- do.call(rbind, lapply(siv_estimands, function(nm) {
    jc_row("SIV", nm, seconds = siv_run$seconds,
           warning = siv_run$warning, error = siv_run$error)
  }))
  siv_test <- data.frame(
    seed = SEED, converged = NA, g_floor_frac = NA_real_,
    min_relevance = NA_real_, contraction_modulus = NA_real_,
    ratio_method = NA_character_,
    falsification_statistic = NA_real_,
    falsification_critical_value = NA_real_,
    falsification_p_value = NA_real_, falsification_reject = NA,
    falsification_dictionary = NA_integer_,
    seconds = siv_run$seconds, warning = siv_run$warning,
    error = siv_run$error, stringsAsFactors = FALSE
  )
} else {
  diagnostic_ok <- fit$converged &&
    fit$diagnostics[["g_floor_frac"]] <= 1e-3 &&
    fit$diagnostics[["min_relevance"]] >= 1e-3
  est <- stats::coef(fit)
  se <- SIV::siv_se(fit)
  ## ATT: Wald interval from the influence function (Theorem 5.1).
  ## QTT: Berger-Boos interval; there is no Wald SE.
  siv_rows <- do.call(rbind, lapply(siv_estimands, function(nm) {
    jc_row("SIV", nm, est[nm], se[nm],
           fit$ci_matrix[nm, "lower"], fit$ci_matrix[nm, "upper"],
           ok = diagnostic_ok, seconds = siv_run$seconds,
           warning = siv_run$warning)
  }))
  ft <- fit$falsification
  siv_test <- data.frame(
    seed = SEED, converged = fit$converged,
    g_floor_frac = fit$diagnostics[["g_floor_frac"]],
    min_relevance = fit$diagnostics[["min_relevance"]],
    contraction_modulus = fit$diagnostics[["contraction_modulus"]],
    ratio_method = fit$ratio_method,
    falsification_statistic = ft$statistic,
    falsification_critical_value = ft$critical_value,
    falsification_p_value = ft$p_value,
    falsification_reject = ft$reject,
    falsification_dictionary = ft$n_dictionary,
    seconds = siv_run$seconds, warning = siv_run$warning,
    error = "", stringsAsFactors = FALSE
  )

  ## Cross-fitted conditional densities for the overlap diagnostic.
  dens <- siv_conditional_densities(fit)
  dens$seed <- SEED
  saveRDS(dens, sprintf("JC_SIV_Density_B%0.5d.rds", SEED))
}

## ---------------------------------------------------------------------------
## Competing methods (ATT only). Settings match
## Simulation_Version2_based_on_SIV_Package/Main/01_compare_methods.R.
## ---------------------------------------------------------------------------

## A separate stream keeps these results independent of the SIV call.
set.seed(SEED + 100000L)
rows <- list()

rows[["SIV_Liu"]] <- run_method("SIV_Liu", "ATT", {
  fit_liu <- SIV::SIV_Liu(
    d$Y, d$A, d$Z, d$X,
    estimator = "dr", outcome_type = "continuous",
    basis = SIV::siv_basis_spec("linear"), n_nodes = 12L
  )
  list(estimate = stats::coef(fit_liu)["ETT"],
       se = SIV::siv_se(fit_liu)["ETT"],
       ok = fit_liu$converged && all(fit_liu$checks))
})

rows[["SIV_Sun"]] <- run_method("SIV_Sun", "ATT", {
  fit_sun <- SIV::SIV_Sun(
    d$Y, d$A, d$Z, d$X,
    outcome_type = "continuous",
    basis = SIV::siv_basis_spec("linear"), n_nodes = 12L
  )
  list(estimate = stats::coef(fit_sun)["ATT"],
       se = SIV::siv_se(fit_sun)["ATT"],
       ok = fit_sun$converged && all(fit_sun$checks))
})

## SIV_LATE estimates the LATE; it is reported in the ATT column as the
## conventional-IV proxy, as in the simulation.
rows[["SIV_LATE"]] <- run_method("SIV_LATE", "ATT", {
  fit_late <- SIV::SIV_LATE(
    d$Y, d$A, d$Z, d$X,
    learner = "SL", SL.list = SL_LIST,
    crossfit = 2L
  )
  list(estimate = stats::coef(fit_late)["LATE"],
       se = SIV::siv_se(fit_late)["LATE"],
       ok = fit_late$converged && abs(fit_late$first_stage) >= 1e-3)
})

rows[["2SLS"]] <- run_method("2SLS", "ATT", {
  ans <- siv_fit_2sls(d$Y, d$A, d$Z, d$X)
  list(estimate = ans["estimate"], se = ans["se"], ok = TRUE)
})

## Unadjusted 2SLS requested by Reviewer 1. With a binary instrument, its
## point estimate is the Wald estimand without conditioning on X.
rows[["2SLS_no_X"]] <- run_method("2SLS_no_X", "ATT", {
  X_empty <- matrix(numeric(0L), nrow = length(d$Y), ncol = 0L)
  ans <- siv_fit_2sls(d$Y, d$A, d$Z, X_empty)
  list(estimate = ans["estimate"], se = ans["se"], ok = TRUE)
})

rows[["Ign"]] <- run_method("Ign", "ATT", {
  ans <- siv_fit_ignorability(d$Y, d$A, d$Z, d$X, crossfit = 2L)
  list(estimate = ans["estimate"], se = ans["se"], ok = TRUE)
})

rows[["linear model"]] <- run_method("linear model", "ATT", {
  ans <- siv_fit_linear_model(d$Y, d$A, d$Z, d$X)
  list(estimate = ans["estimate"], se = ans["se"], ok = TRUE)
})

## ---------------------------------------------------------------------------
## Output: one file of each kind per job.
## ---------------------------------------------------------------------------

raw <- rbind(siv_rows, do.call(rbind, rows))
rownames(raw) <- NULL
raw <- cbind(seed = SEED, raw)

result_columns <- c(
  "seed", "method", "estimand", "estimate", "se",
  "ci95_lower", "ci95_upper"
)
test_columns <- c(
  "seed", "falsification_statistic", "falsification_critical_value",
  "falsification_p_value", "falsification_reject"
)

write.csv(raw[, result_columns, drop = FALSE],
          sprintf("JC_Result_B%0.5d.csv", SEED), row.names = FALSE)
write.csv(siv_test[, test_columns, drop = FALSE],
          sprintf("JC_SIV_Test_B%0.5d.csv", SEED), row.names = FALSE)
 
