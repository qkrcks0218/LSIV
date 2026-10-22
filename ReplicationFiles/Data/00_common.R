## Shared utilities for the version-2 Job Corps analysis.

if (!requireNamespace("SIV", quietly = TRUE)) {
  stop(
    "Install SIV and all required dependencies first:\n  Rscript Code/Install.R",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(SIV))

## ---------------------------------------------------------------------------
## Data. The variable construction is identical to Data_Version1/1.SIV_JC.R.
## ---------------------------------------------------------------------------

jc_load_data <- function(path = "JC.csv") {
  JC.O <- read.csv(path)
  Z <- JC.O$assignment
  A <- round((JC.O$trainy1 + JC.O$trainy2) / 2)
  Y <- round(JC.O$pworky4 / 100 * 52)
  X <- cbind(
    X01 = JC.O$female,
    X02 = (JC.O$age - 16) / 8,
    X03 = JC.O$white,
    X04 = JC.O$black,
    X05 = JC.O$hispanic,
    X06 = JC.O$hsdegree,
    X07 = JC.O$english,
    X08 = JC.O$cohabmarried,
    X09 = JC.O$haschild,
    X10 = JC.O$everwkd
  )
  list(Y = Y, A = A, Z = Z, X = X)
}

## Weeks worked in year 4 take the values 0, ..., 52. As in Version 1, SIV
## uses the nonparametric continuous-outcome nuisances (np baseline density,
## KLIEP/uLSIF density ratio) evaluated on this integer grid with unit
## weights, so every integral is a sum over the 53 support points.
JC_Y_GRID <- 0:52
JC_Y_WEIGHTS <- rep(1, length(JC_Y_GRID))

## ---------------------------------------------------------------------------
## Running a method and recording the outcome as one row per estimand.
## ---------------------------------------------------------------------------

jc_capture <- function(code) {
  warnings <- character(0)
  started <- proc.time()[[3L]]
  value <- withCallingHandlers(
    tryCatch(force(code), error = function(e) {
      structure(list(message = conditionMessage(e)), class = "jc_run_error")
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  failed <- inherits(value, "jc_run_error")
  list(value = if (failed) NULL else value,
       error = if (failed) value$message else "",
       warning = paste(unique(warnings), collapse = " | "),
       seconds = proc.time()[[3L]] - started)
}

jc_row <- function(method, estimand, estimate = NA_real_, se = NA_real_,
                   lower = NA_real_, upper = NA_real_, ok = FALSE,
                   seconds = NA_real_, warning = "", error = "") {
  estimate <- unname(estimate)
  data.frame(
    method = method, estimand = estimand,
    estimate = estimate, se = unname(se),
    ci95_lower = unname(lower), ci95_upper = unname(upper),
    computed = is.finite(estimate),
    ok = isTRUE(ok) && is.finite(estimate),
    seconds = unname(seconds), warning = warning, error = error,
    stringsAsFactors = FALSE
  )
}

## `code` returns list(estimate, se, ok); the interval is the Wald interval.
run_method <- function(method, estimand, code) {
  res <- jc_capture(code)
  if (is.null(res$value)) {
    return(jc_row(method, estimand, seconds = res$seconds,
                  warning = res$warning, error = res$error))
  }
  estimate <- unname(res$value$estimate)
  se <- unname(res$value$se)
  zcrit <- stats::qnorm(0.975)
  jc_row(method, estimand, estimate, se,
         estimate - zcrit * se, estimate + zcrit * se,
         ok = isTRUE(res$value$ok) && is.finite(se),
         seconds = res$seconds, warning = res$warning)
}

## ---------------------------------------------------------------------------
## Conditional outcome densities implied by one SIV fit, for the overlap
## diagnostic of Supplement A.8.2 (02_density_ratio.R).
##
##   f(y | A=0, Z=z, X): the untreated densities entering the fixed point;
##   f(y | A=1, Z=z, X): the treatment-free outcome among the treated,
##       f(Y(0)=y | A=1, Z=z, X)
##         = delta(y,X) f(y|A=0,Z=z,X) / int delta(t,X) f(t|A=0,Z=z,X) dt,
##   because the odds Pr(A=1|Y(0)=y,Z=z,X) / Pr(A=0|Y(0)=y,Z=z,X) equal
##   delta(y,X) gamma(z,X) and Z is independent of Y(0) given X.  At the
##   fixed point this equals fYA1Z* of Data_Version1/1.SIV_JC.R ("Product").
##
## Every row is cross-fitted: unit i's densities come from the models fitted
## on the fold that excludes i. SIV() keeps only the nuisances of its last
## cross-fitting repetition, so the fit must use one repetition.
## ---------------------------------------------------------------------------

siv_conditional_densities <- function(fit) {
  if (!identical(as.integer(fit$num_crossfit), 1L)) {
    stop("`fit` must use Overall.Estimation_num.crossfit = 1.", call. = FALSE)
  }
  fp <- fit$fixed_point
  W <- matrix(fit$weights, nrow(fp$delta), length(fit$weights), byrow = TRUE)
  normalise <- function(f) f / rowSums(f * W)
  f_a0z0 <- normalise(fp$fY_z0)
  f_a0z1 <- normalise(fp$fY_z1)
  list(
    grid = fit$grid, weights = fit$weights,
    fY_a0z0 = f_a0z0,
    fY_a0z1 = f_a0z1,
    fY_a1z0 = normalise(fp$delta * f_a0z0),
    fY_a1z1 = normalise(fp$delta * f_a0z1),
    pZ1 = fit$nuisance$pZ1,
    pA1_z0 = fit$nuisance$pA1_z0,
    pA1_z1 = fit$nuisance$pA1_z1
  )
}

## ---------------------------------------------------------------------------
## Reference estimators. Copied unchanged from
## Simulation_Version2_based_on_SIV_Package/Main/00_common.R.
## ---------------------------------------------------------------------------

siv_fit_2sls <- function(Y, A, Z, X, compute_se = TRUE) {
  X <- as.matrix(X)
  D <- cbind(`(Intercept)` = 1, A = A, X)
  Q <- cbind(`(Intercept)` = 1, Z = Z, X)
  qd <- crossprod(Q, D)
  b <- as.vector(solve(qd, crossprod(Q, Y)))
  if (!compute_se) return(c(estimate = b[2L], se = NA_real_))
  resid <- Y - as.vector(D %*% b)
  s2 <- sum(resid^2) / (length(Y) - ncol(D))
  iq <- solve(qd)
  V <- s2 * iq %*% crossprod(Q) %*% t(iq)
  c(estimate = b[2L], se = sqrt(V[2L, 2L]))
}

siv_fit_linear_model <- function(Y, A, Z, X) {
  dat <- data.frame(Y = Y, A = A, Z = Z, as.data.frame(X))
  rhs <- paste(c("A", "Z", colnames(X)), collapse = " + ")
  fit <- lm(as.formula(paste("Y ~", rhs)), data = dat)
  c(estimate = unname(coef(fit)["A"]),
    se = unname(sqrt(diag(vcov(fit)))["A"]))
}

## The ignorability benchmark: the cross-fitted EIF estimator of tau_Ign with
## L = (Z, X) and glm nuisances fitted on the other folds. See the simulation
## copy of this function for the derivation.
siv_fit_ignorability <- function(Y, A, Z, X, crossfit = 2L, trim = 1e-2) {
  X <- as.matrix(X)
  n <- length(Y)
  crossfit <- max(1L, as.integer(crossfit))
  df <- data.frame(Y = Y, A = as.integer(A), Z = Z, as.data.frame(X))
  rhs <- paste(c("Z", colnames(X)), collapse = " + ")
  f_a <- stats::as.formula(paste("A ~", rhs))
  f_y <- stats::as.formula(paste("Y ~", rhs))

  if (crossfit <= 1L) {
    folds <- list(seq_len(n))
  } else {
    cell <- interaction(df$A, df$Z, drop = TRUE)
    lab <- integer(n)
    for (lv in levels(cell)) {
      idx <- which(cell == lv)
      lab[idx] <- sample(rep_len(seq_len(crossfit), length(idx)))
    }
    folds <- split(seq_len(n), lab)
  }

  pA <- numeric(n)
  m1 <- numeric(n)
  for (te in folds) {
    tr <- if (length(folds) == 1L) te else setdiff(seq_len(n), te)
    fa <- stats::glm(f_a, data = df[tr, , drop = FALSE],
                     family = stats::binomial())
    pA[te] <- stats::predict(fa, newdata = df[te, , drop = FALSE],
                             type = "response")
    tr1 <- tr[df$A[tr] == 1L]
    fy <- stats::glm(f_y, data = df[tr1, , drop = FALSE],
                     family = stats::gaussian())
    m1[te] <- stats::predict(fy, newdata = df[te, , drop = FALSE],
                             type = "response")
  }
  pA <- pmin(pmax(pA, trim), 1 - trim)

  num <- (A - (1 - A) * pA / (1 - pA)) * (Y - m1)
  estimate <- sum(num) / sum(A)
  eif <- (num - estimate * A) / mean(A)
  c(estimate = estimate, se = stats::sd(eif) / sqrt(n))
}
