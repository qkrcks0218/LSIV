## One estimated-nuisance falsification-test replication per batch job.

args <- (commandArgs(trailingOnly=TRUE))
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

library(SIV)
source("00_common.R")

NUM_BOOT <- 1000L
SL_LIST <- 1:9
NUM_CROSSFIT <- 3L
NUM_FOLDS <- 2L
Y_BINS <- 20L
X_BINS <- 5L
N_GRID <- 101L
BW_QUANTILE <- c(0.5,0.75,0.9,0.95)

DGP_VALUES <- c("SIV", "mild", "severe")
DGP_LABELS <- c(
  SIV = "SIV holds",
  mild = "Mild violation",
  severe = "Severe violation"
)

HYPERPARAMETER <- expand.grid(
  n = seq.int(500L, 4000L, by = 500L),
  dgp = DGP_VALUES,
  type = "continuous",
  seed = 1:510,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

BATCH <- as.integer(BATCH)

N <- HYPERPARAMETER$n[BATCH]
DGP <- HYPERPARAMETER$dgp[BATCH]
TYPE <- HYPERPARAMETER$type[BATCH]
SEED <- HYPERPARAMETER$seed[BATCH]

message(sprintf(
  "BATCH=%d/%d: N=%d, dgp=%s, type=%s, seed=%d",
  BATCH, nrow(HYPERPARAMETER), N, DGP, TYPE, SEED
))

## The test is run through the public API, `SIV(falsification = TRUE)`, so the
## simulation exercises exactly the code path a user gets: the Super Learner
## nuisance estimates, the ATT fixed point and the falsification test all come
## out of one call, on the cross-fitted nuisances the test requires.
## The package always uses the marginal X-cell dictionary.

one_rep <- function(seed, n, dgp, type) {
  violation <- tolower(dgp)
  data_seed <- as.integer(seed)
  bootstrap_seed <- data_seed + 500000000L
  warnings <- character(0)
  started <- proc.time()[[3L]]
  
  ans <- withCallingHandlers(
    tryCatch({
      set.seed(data_seed)
      d <- siv_dgp(n, type = type, violation = violation)
      
      ## Separating the data seed from the estimation seed keeps the fit
      ## invariant to the number of parallel workers.  The second seed governs
      ## the cross-fitting folds, the Super Learner cross-validation and the
      ## multiplier bootstrap alike, since SIV() draws them all in this call.
      set.seed(bootstrap_seed)
      fit <- SIV::SIV(
        d$Y, d$A, d$Z, d$X,
        y_type = type,
        Overall.Estimation_SL.list = SL_LIST,
        Overall.Estimation_num.crossfit = NUM_CROSSFIT,
        Overall.Estimation_num.fold = NUM_FOLDS,
        y_num.grid = N_GRID,
        outcome_include.observed = FALSE,
        outcome_density="np_gss_cv",
        ratio_bandwidth = "cv_common",
        ratio_bandwidth.cv.grid = BW_QUANTILE,
        falsification = TRUE,
        falsification_num.boot = NUM_BOOT,
        falsification_num.bin = c(Y_BINS, rep.int(X_BINS, ncol(d$X))),
        falsification_floor = NULL,
        falsification_categorical = NULL,
        level = 0.95
      )
      ## SIV() retains the nuisance fit from its last cross-fitting repetition.
      ## Label these auxiliary diagnostics accordingly. The reported test
      ## itself uses the median p-value across all repetitions.
      nuisance <- fit$nuisance
      
      list(
        ft = fit$falsification,
        att = unname(stats::coef(fit)["ATT"]),
        converged = fit$converged,
        g_floor_frac = unname(fit$diagnostics[["g_floor_frac"]]),
        min_relevance = unname(fit$diagnostics[["min_relevance"]]),
        marginal_first_stage = mean(d$pA1 - d$pA0),
        reversal_probability = mean(d$pA1 < d$pA0),
        reversal_magnitude = mean(pmax(d$pA0 - d$pA1, 0)),
        last_crossfit_min_conditional_first_stage =
          min(nuisance$pA1_z1 - nuisance$pA1_z0),
        last_crossfit_mean_conditional_first_stage =
          mean(nuisance$pA1_z1 - nuisance$pA1_z0)
      )
    }, error = function(e) {
      structure(conditionMessage(e), class = "siv_falsification_error")
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seconds <- proc.time()[[3L]] - started
  
  prefix <- data.frame(
    dgp = dgp,
    dgp_label = unname(DGP_LABELS[dgp]),
    n = n,
    type = type,
    seed = seed,
    data_seed = data_seed,
    bootstrap_seed = bootstrap_seed,
    num_crossfit = NUM_CROSSFIT,
    num_fold = NUM_FOLDS,
    stringsAsFactors = FALSE
  )
  if (inherits(ans, "siv_falsification_error")) {
    return(cbind(
      prefix,
      data.frame(
        reject = NA, p_value = NA_real_, statistic = NA_real_,
        critical_value = NA_real_, dictionary_size = NA_integer_,
        split_p_values = "", representative_crossfit = NA_integer_,
        att = NA_real_, converged = NA, g_floor_frac = NA_real_,
        min_relevance = NA_real_,
        marginal_first_stage = NA_real_, reversal_probability = NA_real_,
        reversal_magnitude = NA_real_,
        last_crossfit_min_conditional_first_stage = NA_real_,
        last_crossfit_mean_conditional_first_stage = NA_real_,
        seconds = unname(seconds),
        n_warnings = length(warnings),
        warning = paste(unique(warnings), collapse = " | "),
        error = as.character(ans), stringsAsFactors = FALSE
      )
    ))
  }
  
  cbind(
    prefix,
    data.frame(
      reject = ans$ft$reject,
      p_value = ans$ft$p_value,
      statistic = ans$ft$statistic,
      critical_value = ans$ft$critical_value,
      dictionary_size = ans$ft$n_dictionary,
      split_p_values = paste(
        signif(ans$ft$median_adjustment$p_values, 8L), collapse = ";"
      ),
      representative_crossfit = ans$ft$median_adjustment$representative,
      att = ans$att,
      converged = ans$converged,
      g_floor_frac = ans$g_floor_frac,
      min_relevance = ans$min_relevance,
      marginal_first_stage = ans$marginal_first_stage,
      reversal_probability = ans$reversal_probability,
      reversal_magnitude = ans$reversal_magnitude,
      last_crossfit_min_conditional_first_stage =
        ans$last_crossfit_min_conditional_first_stage,
      last_crossfit_mean_conditional_first_stage =
        ans$last_crossfit_mean_conditional_first_stage,
      seconds = unname(seconds),
      n_warnings = length(warnings),
      warning = paste(unique(warnings), collapse = " | "),
      error = "", stringsAsFactors = FALSE
    )
  )
}

## ---------------------------------------------------------------------------
## Run this batch's single cell and write it out, exactly as
## 01_compare_methods.R does: one CSV per job, merged afterwards.
## ---------------------------------------------------------------------------

row <- one_rep(
  seed = SEED,
  n = N,
  dgp = DGP,
  type = TYPE
)

output_columns <- c("dgp", "dgp_label", "n", "type", "seed", "reject")

write.csv(
  row[, output_columns, drop = FALSE],
  sprintf("Falsification_N%0.4d_%s_B%0.5d.csv", N, DGP, SEED),
  row.names = FALSE
)
