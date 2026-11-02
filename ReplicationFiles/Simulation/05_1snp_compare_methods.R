## Batch jobs

args <- (commandArgs(trailingOnly=TRUE))
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

library(SIV)
source("00_1snp_common.R")

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


## Script-local implementation of Supplement Section A.6.1, equation (S.29).
## Estimate e(X)=Pr(Z=1|X) and m0(X)=E(Y|Z=0,X) on training observations,
## then evaluate sum{(Z-e_hat)/(1-e_hat)*(Y-m0_hat)}/sum(A) out of fold.
siv_fit_supp_one_sided <- function(Y, A, Z, X,
                                   type = c("continuous", "binary"),
                                   learner = c("SL", "glm"), SL.list = 1:9,
                                   crossfit = 2L, num.crossfit = 1L,
                                   trim = 1e-2) {
  type <- match.arg(type)
  learner <- match.arg(learner)
  checked <- getFromNamespace("siv_check_data", "SIV")(Y, A, Z, X)
  Y <- checked$Y; A <- checked$A; Z <- checked$Z; X <- checked$X
  n <- length(Y)
  if (any(A > Z)) {
    stop("The one-sided ATT estimator requires A=0 whenever Z=0.",
         call. = FALSE)
  }
  if (any(!is.finite(Y)) || any(!is.finite(X))) {
    stop("Y and X must contain only finite values.", call. = FALSE)
  }
  if (type == "binary" && any(!Y %in% c(0, 1))) {
    stop("Binary Y must be coded as 0/1.", call. = FALSE)
  }
  if (length(crossfit) != 1L || !is.finite(crossfit) ||
      crossfit < 2L || crossfit > n || crossfit != as.integer(crossfit)) {
    stop("crossfit must be an integer between 2 and n.", call. = FALSE)
  }
  if (length(num.crossfit) != 1L || !is.finite(num.crossfit) ||
      num.crossfit < 1L || num.crossfit != as.integer(num.crossfit)) {
    stop("num.crossfit must be a positive integer.", call. = FALSE)
  }
  if (length(trim) != 1L || !is.finite(trim) || trim <= 0 || trim >= 0.5) {
    stop("trim must lie in (0, 0.5).", call. = FALSE)
  }

  make_folds <- getFromNamespace("siv_make_folds", "SIV")
  fit_nuisance <- getFromNamespace("siv_learn_fit", "SIV")
  predict_nuisance <- getFromNamespace("siv_learn_predict", "SIV")
  df <- data.frame(Y = Y, Z = Z, X, check.names = FALSE)
  xnames <- colnames(X)
  outcome_family <- if (type == "binary") stats::binomial() else
    stats::gaussian()
  rho_hat <- mean(A)

  splits <- lapply(seq_len(num.crossfit), function(s) {
    folds <- make_folds(n, A, Z, crossfit)
    e_hat <- m0_hat <- rep(NA_real_, n)
    for (te in folds) {
      tr <- setdiff(seq_len(n), te)
      tr0 <- tr[Z[tr] == 0]
      if (length(tr0) < 2L || !any(Z[tr] == 1)) {
        stop("Each training sample must contain both instrument arms.",
             call. = FALSE)
      }
      fit_e <- fit_nuisance(
        df[tr, , drop = FALSE], "Z", xnames, stats::binomial(),
        learner = learner, SL.list = SL.list, SL.N = n
      )
      fit_m0 <- fit_nuisance(
        df[tr0, , drop = FALSE], "Y", xnames, outcome_family,
        learner = learner, SL.list = SL.list, SL.N = n
      )
      e_hat[te] <- predict_nuisance(fit_e, df[te, , drop = FALSE])
      m0_hat[te] <- predict_nuisance(fit_m0, df[te, , drop = FALSE])
    }
    if (any(!is.finite(c(e_hat, m0_hat)))) {
      stop("Non-finite one-sided nuisance prediction.", call. = FALSE)
    }
    n_trimmed <- sum(e_hat < trim | e_hat > 1 - trim)
    e_hat <- pmin(pmax(e_hat, trim), 1 - trim)
    score <- (Z - e_hat) / (1 - e_hat) * (Y - m0_hat)
    estimate <- sum(score) / sum(A)
    eif <- (score - estimate * A) / rho_hat
    list(estimate = estimate, variance = mean(eif^2) / n,
         eif = eif, e_hat = e_hat, m0_hat = m0_hat,
         folds = folds, n_trimmed = n_trimmed)
  })

  ## Use the same repeated-split median adjustment as the SIV comparison.
  estimates <- vapply(splits, function(x) x$estimate, numeric(1))
  variances <- vapply(splits, function(x) x$variance, numeric(1))
  estimate <- stats::median(estimates)
  se <- sqrt(stats::median(variances + (estimates - estimate)^2))
  list(estimate = estimate, se = se,
       ok = all(is.finite(c(estimates, variances, se))),
       split_estimates = estimates, split_variances = variances,
       splits = splits)
}

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
  
  ## Constant assignment in this one-sided DGP makes LATE equal to ATT.
  rows[["SIV_LATE"]] <- run_method("SIV_LATE", "ATT", truth["ATT"], {
    fit <- SIV::SIV_LATE(
      d$Y, d$A, d$Z, d$X,
      learner = "SL", SL.list = SL_LIST,
      crossfit = 2L
    )
    list(estimate = stats::coef(fit)["LATE"], se = SIV::siv_se(fit)["LATE"],
         ok = fit$converged && abs(fit$first_stage) >= 1e-3)
  })
  
  rows[["SIV_one_sided"]] <- run_method(
    "SIV_one_sided", "ATT", truth["ATT"], {
      siv_fit_supp_one_sided(
        d$Y, d$A, d$Z, d$X, type = type,
        learner = "SL", SL.list = SL_LIST,
        crossfit = NUM_FOLDS, num.crossfit = NUM_CROSSFIT
      )
    })
  
  do.call(rbind, rows[c("SIV", "SIV_one_sided", "SIV_LATE")])
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
          sprintf("Result_1snp_N%0.4d_Type%s_B%0.5d.csv",N,substr(TYPES,1,1),SEED), 
          row.names = FALSE) 
