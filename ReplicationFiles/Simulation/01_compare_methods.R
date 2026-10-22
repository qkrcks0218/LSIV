## Batch jobs

args <- (commandArgs(trailingOnly=TRUE))
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

library(SIV)
source("00_common.R")

## Paired comparison of the seven requested methods.
## Defaults are intentionally a quick check: continuous DGPs only, one
## replication per batch, and all nine Super Learner groups. Set
## SIM_TYPE=binary or SIM_TYPE=both explicitly when those designs are wanted.

REPS <- 1
TRUTH_N <- 10^(5)
SL_LIST <- 1:9
NUM_CROSSFIT <- 3L
NUM_FOLDS <- 2L
BW_QUANTILE <- c(0.5,0.75,0.9,0.95)

HYPERPARAMETER <- expand.grid(
  n = c(1000,2000,4000),
  type = c("continuous","binary"),
  seed = 1:1010,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

if (length(BATCH) != 1L || !is.finite(BATCH) || BATCH != as.integer(BATCH) ||
    BATCH < 1L || BATCH > nrow(HYPERPARAMETER)) {
  stop(sprintf("BATCH must be an integer from 1 to %d for SIM_TYPE=%s.",
               nrow(HYPERPARAMETER), Sys.getenv("SIM_TYPE", "continuous")),
       call. = FALSE)
}
BATCH <- as.integer(BATCH)

N     <- HYPERPARAMETER$n[BATCH]
TYPES <- HYPERPARAMETER$type[BATCH]
SEED  <- HYPERPARAMETER$seed[BATCH]

message(sprintf("BATCH=%d: N=%d, type=%s, seed=%d",
                BATCH, N, TYPES, SEED))


run_method <- function(method, estimand, truth, code) {
  warnings <- character(0)
  started <- proc.time()[[3L]]
  value <- withCallingHandlers(
    tryCatch(force(code), error = function(e) {
      structure(list(message = conditionMessage(e)), class = "siv_run_error")
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seconds <- proc.time()[[3L]] - started
  if (inherits(value, "siv_run_error")) {
    return(data.frame(
      method = method, estimand = estimand, truth = truth,
      estimate = NA_real_, se = NA_real_,
      ci95_lower = NA_real_, ci95_upper = NA_real_,
      computed = FALSE, ok = FALSE, seconds = seconds,
      warning = paste(unique(warnings), collapse = " | "),
      error = value$message, stringsAsFactors = FALSE
    ))
  }
  estimate <- unname(value$estimate)
  se <- unname(value$se)
  computed <- is.finite(estimate) && is.finite(se)
  data.frame(
    method = method, estimand = estimand, truth = truth,
    estimate = estimate, se = se,
    ci95_lower = estimate - stats::qnorm(0.975) * se,
    ci95_upper = estimate + stats::qnorm(0.975) * se,
    computed = computed,
    ok = isTRUE(value$ok) && computed,
    seconds = seconds, warning = paste(unique(warnings), collapse = " | "),
    error = "", stringsAsFactors = FALSE
  )
}

fit_one_dataset <- function(d, type, truth) {
  holder <- new.env(parent = emptyenv())
  holder$siv <- NULL
  rows <- list()
  
  ## Package methods: the calls are kept explicit so the comparison settings
  ## can be audited at a glance. SL_LIST includes learner groups 1 through 9.
  rows[["SIV"]] <- run_method("SIV", "ATT", truth["ATT"], {
    fit <- SIV::SIV(
      d$Y, d$A, d$Z, d$X,
      y_type = type,
      Overall.Estimation_SL.list = SL_LIST,
      Overall.Estimation_num.crossfit = NUM_CROSSFIT,
      Overall.Estimation_num.fold = NUM_FOLDS,
      outcome_include.observed = FALSE,
      outcome_density="np_gss_cv",
      ratio_bandwidth = "cv_common",
      ratio_bandwidth.cv.grid = BW_QUANTILE
    )
    holder$siv <- fit
    diagnostic_ok <- fit$diagnostics[["g_floor_frac"]] <= 1e-3 &&
      fit$diagnostics[["min_relevance"]] >= 1e-3
    list(estimate = stats::coef(fit)["ATT"], se = SIV::siv_se(fit)["ATT"],
         ok = fit$converged && diagnostic_ok)
  })
  
  rows[["SIV_Liu"]] <- run_method("SIV_Liu", "ATT", truth["ATT"], {
    fit <- SIV::SIV_Liu(
      d$Y, d$A, d$Z, d$X,
      estimator = "dr", outcome_type = type,
      basis = SIV::siv_basis_spec("linear"), n_nodes = 12L
    )
    list(estimate = stats::coef(fit)["ETT"], se = SIV::siv_se(fit)["ETT"],
         ok = fit$converged && all(fit$checks))
  })
  
  rows[["SIV_Sun"]] <- run_method("SIV_Sun", "ATT", truth["ATT"], {
    fit <- SIV::SIV_Sun(
      d$Y, d$A, d$Z, d$X,
      outcome_type = type,
      basis = SIV::siv_basis_spec("linear"), n_nodes = 12L
    )
    list(estimate = stats::coef(fit)["ATT"], se = SIV::siv_se(fit)["ATT"],
         ok = fit$converged && all(fit$checks))
  })
  
  ## Deliberately score the LATE estimate against ATT: this quantifies the
  ## error from using conventional IV/LATE as a proxy for the ATT target.
  rows[["SIV_LATE"]] <- run_method("SIV_LATE", "ATT", truth["ATT"], {
    fit <- SIV::SIV_LATE(
      d$Y, d$A, d$Z, d$X,
      learner = "SL", SL.list = SL_LIST,
      crossfit = 2L
    )
    list(estimate = stats::coef(fit)["LATE"], se = SIV::siv_se(fit)["LATE"],
         ok = fit$converged && abs(fit$first_stage) >= 1e-3)
  })
  
  ## Reference methods. They are intentionally not added to SIV_Package.
  ## Likewise, score the 2SLS coefficient against ATT rather than its own
  ## population projection to measure its error as an ATT proxy.
  rows[["2SLS"]] <- run_method(
    "2SLS", "ATT", truth["ATT"], {
      ans <- siv_fit_2sls(d$Y, d$A, d$Z, d$X)
      list(estimate = ans["estimate"], se = ans["se"], ok = TRUE)
    })
  
  ## Ignorability benchmark: a self-contained cross-fitted EIF estimator with
  ## glm nuisances, so it no longer borrows SIV's fitted functions.
  rows[["Ign"]] <- run_method("Ign", "ATT", truth["ATT"], {
    ans <- siv_fit_ignorability(d$Y, d$A, d$Z, d$X, crossfit = 2L)
    list(estimate = ans["estimate"], se = ans["se"], ok = TRUE)
  })
  
  rows[["linear model"]] <- run_method("linear model", "ATT", truth["ATT"], {
    ans <- siv_fit_linear_model(d$Y, d$A, d$Z, d$X)
    list(estimate = ans["estimate"], se = ans["se"], ok = TRUE)
  })
  
  do.call(rbind, rows)
}

truths <- lapply(TYPES, function(type) {
  message(sprintf("Computing %s reference values with %s draws ...", type,
                  format(TRUTH_N, big.mark = ",", scientific = FALSE)))
  siv_compute_truths(type, n = TRUTH_N, seed = SEED)
})
names(truths) <- TYPES

raw <- list()
for (type in TYPES) {
  for (r in seq_len(REPS)) {
    message(sprintf("[%s] replication %d/%d", type, r, REPS))
    set.seed(SEED + 10000L * match(type, TYPES) + r)
    d <- siv_dgp(N, type)
    set.seed(SEED + 20000L * match(type, TYPES) + r)
    one <- fit_one_dataset(d, type, truths[[type]])
    one$type <- type
    one$n <- N
    one$rep <- r
    raw[[paste(type, r)]] <- one
  }
}
raw <- do.call(rbind, raw)
rownames(raw) <- NULL
summary <- siv_summarise_results(raw)

raw.seed <- cbind(SEED,raw)
colnames(raw.seed) <- c("seed",colnames(raw))

output_columns <- c(
  "seed", "method", "truth", "estimate", "se",
  "ci95_lower", "ci95_upper", "type", "n"
)

write.csv(raw.seed[, output_columns, drop = FALSE],
          sprintf("Result_N%0.4d_Type%s_B%0.5d.csv",N,substr(TYPES,1,1),SEED), 
          row.names = FALSE) 
