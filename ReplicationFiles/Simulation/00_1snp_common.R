## Shared utilities for the R1-3 LATE comparison.
## X ~ N_p(0,I), W=sum(X), Z|X ~ Bernoulli(1/2), and
## A|Z,X ~ Bernoulli(Z*plogis(W)). The default p is 3.
## Treatment is generated independently of the potential outcomes given X,Z.
## Potential-outcome distributions retain the Final18_FF specifications.
## Consequently alpha=1, beta(0,X)=0, and beta(1,X)=exp(W).

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x
}

if (!requireNamespace("LSIV", quietly = TRUE)) {
  stop(
    "Install LSIV and its required R packages with:\n  devtools::install_github(\"qkrcks0218/LSIV\")",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(LSIV))

required_lsiv_functions <- c("lsiv_check_data", "lsiv_make_folds",
                             "lsiv_learn_fit", "lsiv_learn_predict",
                             "lsiv_late_eif")
missing_lsiv_functions <- required_lsiv_functions[!vapply(
  required_lsiv_functions, exists, logical(1), envir = asNamespace("LSIV"),
  inherits = FALSE)]
if (length(missing_lsiv_functions)) {
  stop("Update the installed LSIV package. Missing functions: ",
       paste(missing_lsiv_functions, collapse = ", "), call. = FALSE)
}

lsiv_env_int <- function(name, default, lower = 1L) {
  ans <- suppressWarnings(as.integer(Sys.getenv(name) %||% as.character(default)))
  if (length(ans) != 1L || is.na(ans) || ans < lower) {
    stop(sprintf("%s must be an integer >= %d.", name, lower), call. = FALSE)
  }
  ans
}

## Comma-separated list of positive integers, e.g. SIM_N="500,1000,2000".
lsiv_env_int_vector <- function(name, default, lower = 1L) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) value <- paste(default, collapse = ",")
  ans <- suppressWarnings(as.integer(trimws(
    strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(ans) || anyNA(ans) || any(ans < lower)) {
    stop(sprintf("%s must be a comma-separated list of integers >= %d.",
                 name, lower), call. = FALSE)
  }
  ans
}

lsiv_env_types <- function(name = "SIM_TYPE", default = "both") {
  value <- tolower(Sys.getenv(name) %||% default)
  if (identical(value, "both")) return(c("continuous", "binary"))
  value <- strsplit(value, ",", fixed = TRUE)[[1L]]
  value <- trimws(value)
  if (!length(value) || !all(value %in% c("continuous", "binary"))) {
    stop(sprintf("%s must be continuous, binary, or both.", name), call. = FALSE)
  }
  unique(value)
}

## ---------------------------------------------------------------------------
## Simulation DGPs. These live here, not in the LSIV package: they are
## properties of the simulation designs rather than of any estimator
## implementation.
## ---------------------------------------------------------------------------

## EDIT THE DGP HERE. All observed-data generation, potential-treatment draws,
## oracle nuisances, and Monte Carlo truths below use this specification.
lsiv_dgp_forall <- function(type = c("continuous", "binary"), p = 3L,
                            violation = "lsiv", one_sided = TRUE,
                            assignment = "randomized") {
  type <- match.arg(type)
  if (!identical(violation, "lsiv") || !isTRUE(one_sided) ||
      !identical(assignment, "randomized")) {
    stop("This comparison uses one_sided=TRUE, assignment='randomized', violation='lsiv'.",
         call. = FALSE)
  }
  if (length(p) != 1L || !is.finite(p) || p < 1 || p != as.integer(p)) {
    stop("p must be a positive integer.", call. = FALSE)
  }
  p <- as.integer(p)
  draw_x <- function(n) {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- sprintf("X%02d", seq_len(p))
    X
  }
  x_index <- function(X) rowSums(as.matrix(X))
  p_z1 <- function(W) rep(0.5, length(W))
  p_a1 <- function(W) plogis(0.5*W)
  log_gamma <- function(z, W) {
    ifelse(rep_len(z, length(W)) == 0, -Inf, 0.5*W)
  }
  alpha_shape <- function(y) numeric(length(y))
  log_delta <- function(y0, W, X = NULL) numeric(length(W))
  log_alpha <- function(y0, z, W, X = NULL) numeric(length(W))
  alpha <- function(y0, z, W, X = NULL) rep(1, length(W))
  if (type == "continuous") {
    mean_y0 <- function(W, X = NULL) 0.3 - 0.2 * W
    mean_y1 <- function(W, X = NULL) 1 + W
    draw_y0 <- function(W, X = NULL) mean_y0(W) + rnorm(length(W))
    draw_y1 <- function(W, X = NULL) mean_y1(W) + 0.25 * rnorm(length(W))
    density_y0 <- function(grid, W, X = NULL) {
      vapply(grid, function(y) dnorm(y, mean_y0(W), 1), numeric(length(W)))
    }
  } else {
    mean_y0 <- function(W, X = NULL) plogis(-W)
    mean_y1 <- function(W, X = NULL) plogis(2 + W)
    draw_y0 <- function(W, X = NULL) rbinom(length(W), 1, mean_y0(W))
    draw_y1 <- function(W, X = NULL) rbinom(length(W), 1, mean_y1(W))
    density_y0 <- function(grid, W, X = NULL) {
      prob <- mean_y0(W)
      vapply(grid, function(y) if (y == 1) prob else 1 - prob,
             numeric(length(W)))
    }
  }
  list(type = type, p = p, draw_x = draw_x, x_index = x_index,
       p_z1 = p_z1, p_a1 = p_a1, draw_y0 = draw_y0, draw_y1 = draw_y1,
       mean_y0 = mean_y0, mean_y1 = mean_y1, density_y0 = density_y0,
       alpha_shape = alpha_shape, log_delta = log_delta, log_gamma = log_gamma,
       violation = "lsiv", violation_parameters = list(kappa = 0),
       baseline_alpha_shape = alpha_shape, log_alpha = log_alpha, alpha = alpha,
       one_sided = TRUE, assignment = "randomized")
}

## Compatibility helper used by a few diagnostics.
lsiv_alpha_shape <- function(y, type = c("continuous", "binary")) {
  type <- match.arg(type)
  lsiv_dgp_forall(type)$alpha_shape(y)
}

## Generic sampler. Do not repeat DGP formulas here.
lsiv_dgp <- function(n, type = c("continuous", "binary"), p = 3L,
                     violation = "lsiv",
                     one_sided = TRUE,
                     assignment = "randomized") {
  type <- match.arg(type)
  violation <- match.arg(violation)
  assignment <- match.arg(assignment)
  spec <- lsiv_dgp_forall(type, p, violation, one_sided, assignment)
  X <- spec$draw_x(n)
  W <- spec$x_index(X)
  pZ1 <- spec$p_z1(W)
  Z <- rbinom(n, 1, pZ1)
  Y0 <- spec$draw_y0(W, X)
  Y1 <- spec$draw_y1(W, X)
  ## With `violation = "lsiv"` (the default) log_alpha == log_delta for both z,
  ## so this is exactly the separable design used by the estimation scripts.
  logalpha_z0 <- spec$log_alpha(Y0, 0, W, X)
  logalpha_z1 <- spec$log_alpha(Y0, 1, W, X)
  logdelta <- spec$log_delta(Y0, W, X)
  loggamma_z0 <- spec$log_gamma(0, W)
  loggamma_z1 <- spec$log_gamma(1, W)
  pA0 <- plogis(logalpha_z0 + loggamma_z0)
  pA1 <- plogis(logalpha_z1 + loggamma_z1)
  treatment_rank <- runif(n)
  A0 <- as.integer(treatment_rank <= pA0)
  A1 <- as.integer(treatment_rank <= pA1)
  A <- ifelse(Z == 1, A1, A0)
  pA <- ifelse(Z == 1, pA1, pA0)
  loggamma <- ifelse(Z == 1, loggamma_z1, loggamma_z0)
  Y <- (1 - A) * Y0 + A * Y1
  list(Y = Y, A = A, Z = Z, X = X, Y0 = Y0, Y1 = Y1, W = W,
       A0 = A0, A1 = A1, pA = pA, pA0 = pA0, pA1 = pA1, pZ1 = pZ1,
       logdelta = logdelta, loggamma = loggamma, type = type,
       violation = violation,
       alpha = exp(ifelse(Z == 1, logalpha_z1, logalpha_z0)),
       alpha_z0 = exp(logalpha_z0), alpha_z1 = exp(logalpha_z1),
       logalpha = ifelse(Z == 1, logalpha_z1, logalpha_z0),
       logalpha_z0 = logalpha_z0, logalpha_z1 = logalpha_z1,
       one_sided = one_sided, assignment = assignment,
       mean_y1 = spec$mean_y1(W, X))
}

## Generic oracle calculation from the same central DGP specification.
lsiv_oracle_nuisance <- function(X, grid, weights,
                                 type = c("continuous", "binary"),
                                 violation = "lsiv", one_sided = TRUE,
                                 assignment = "randomized") {
  type <- match.arg(type)
  X <- as.matrix(X)
  n <- nrow(X); M <- length(grid)
  spec <- lsiv_dgp_forall(type, ncol(X), violation, one_sided, assignment)
  W <- spec$x_index(X)
  fy <- matrix(spec$density_y0(grid, W, X), nrow = n, ncol = M)
  Wm <- matrix(weights, n, M, byrow = TRUE)
  fY0 <- fy / rowSums(fy * Wm)
  p1 <- spec$p_a1(W)
  zero <- matrix(0, n, M)
  one <- matrix(1, n, M)
  list(pZ1 = spec$p_z1(W), pA1_z0 = numeric(n), pA1_z1 = p1,
       fY_z0 = fY0, fY_z1 = fY0,
       fYA0_z0 = fY0, fYA0_z1 = (1 - p1) * fY0,
       delta = one, alpha_z0 = one, alpha_z1 = one,
       logalpha_z0 = zero, logalpha_z1 = zero,
       gamma_z0 = numeric(n), gamma_z1 = exp(W), fY0 = fY0,
       one_sided = TRUE, assignment = "randomized")
}

## ---------------------------------------------------------------------------
## Experimental one-sided boundary calculations, adapted from R1_3_OneSided.
## These support continuous X. No saturated four-cell regression is carried
## over from the earlier two-Bernoulli-covariate pilot.
## ---------------------------------------------------------------------------

## Analytic normalized update for compatible positive density inputs. The
## returned delta and gamma_z0=0 are boundary continuations of the interior
## parameters. Incompatible fitted densities are reported, not silently fixed.
lsiv_one_sided_fixed_point <- function(pA1_z1, fY_z0, fY_z1, grid, weights,
                                       y_ref = NULL) {
  f0 <- as.matrix(fY_z0)
  f1 <- as.matrix(fY_z1)
  n <- nrow(f0)
  M <- ncol(f0)
  if (!identical(dim(f0), dim(f1)) || length(pA1_z1) != n ||
      length(grid) != M || length(weights) != M ||
      any(!is.finite(c(f0, f1, pA1_z1, grid, weights))) ||
      any(f0 <= 0) || any(f1 <= 0) || any(weights <= 0) ||
      any(pA1_z1 <= 0 | pA1_z1 >= 1)) {
    stop("Positive normalized densities and a first stage in (0,1) are required.")
  }
  if (max(abs(as.vector(f0 %*% weights) - 1),
          abs(as.vector(f1 %*% weights) - 1)) > 1e-8) {
    stop("The supplied densities must integrate to one.")
  }
  RY <- f1 / f0
  slack <- f0 - (1 - pA1_z1) * f1
  if (any(slack <= 0)) {
    stop("The one-sided reconstruction requires f0-(1-p1)*f1>0 on the grid.")
  }
  h <- slack / (pA1_z1 * RY)
  g <- h / as.vector(h %*% weights)
  if (is.null(y_ref)) y_ref <- grid[which.min(abs(grid))]
  iref <- which.min(abs(grid - y_ref))
  if (abs(grid[iref] - y_ref) > 1e-8) stop("y_ref must be a grid point.")
  delta <- g / f0
  delta <- delta / delta[, iref]
  denom0 <- as.vector((delta * f0) %*% weights)
  denom1 <- as.vector((delta * f1) %*% weights)
  list(g = g, delta = delta, gamma_z0 = numeric(n),
       gamma_z1 = pA1_z1 / (1 - pA1_z1) / denom1,
       denom_z0 = denom0, denom_z1 = denom1,
       fY_z0 = f0, fY_z1 = f1, RY = RY,
       iter = 1L, converged = rep(TRUE, n), min_slack = min(slack),
       y_ref_index = iref)
}

## Fit arm-specific regressions on all covariates. The first stage is fitted
## only in Z=1, and its Z=0 counterpart is exactly zero. LSIV and LATE receive
## the same m0, instrument propensities, and cross-fitting splits.
lsiv_one_sided_nuisance <- function(Y, A, Z, X,
                                    type = c("continuous", "binary"),
                                    learner = c("SL", "glm"), SL.list = 1:9,
                                    crossfit = 2L, folds = NULL,
                                    known_pZ1 = NULL, trim = 1e-2) {
  type <- match.arg(type)
  learner <- match.arg(learner)
  checked <- getFromNamespace("lsiv_check_data", "LSIV")(Y, A, Z, X)
  Y <- checked$Y; A <- checked$A; Z <- checked$Z; X <- checked$X
  n <- length(Y)
  if (any(A > Z)) stop("One-sided estimation requires A=0 whenever Z=0.")
  if (type == "binary" && any(!Y %in% c(0, 1))) stop("Binary Y must be 0/1.")
  if (length(trim) != 1L || !is.finite(trim) || trim <= 0 || trim >= 0.5) {
    stop("trim must lie in (0,0.5).")
  }
  if (is.null(folds)) {
    folds <- getFromNamespace("lsiv_make_folds", "LSIV")(n, A, Z, crossfit)
  }
  if (!is.list(folds) || !length(folds) || any(lengths(folds) == 0) ||
      !identical(sort(as.integer(unlist(folds))), seq_len(n))) {
    stop("folds must partition the observations exactly once.")
  }
  if (!is.null(known_pZ1)) {
    if (!length(known_pZ1) %in% c(1L, n) ||
        any(!is.finite(known_pZ1)) || any(known_pZ1 <= 0 | known_pZ1 >= 1)) {
      stop("known_pZ1 must be one probability or one probability per observation.")
    }
    known_pZ1 <- rep_len(known_pZ1, n)
  }
  df <- data.frame(Y = Y, A = A, Z = Z, X, check.names = FALSE)
  xnames <- colnames(X)
  fit <- getFromNamespace("lsiv_learn_fit", "LSIV")
  predict_fit <- getFromNamespace("lsiv_learn_predict", "LSIV")
  guard <- function(x) pmin(pmax(x, trim), 1 - trim)
  family <- if (type == "binary") binomial() else gaussian()
  e <- m0 <- m1 <- p1 <- numeric(n)
  for (te in folds) {
    tr <- if (length(folds) == 1L) te else setdiff(seq_len(n), te)
    tr0 <- tr[Z[tr] == 0]
    tr1 <- tr[Z[tr] == 1]
    if (min(length(tr0), length(tr1)) < 2L) stop("A training arm is empty.")
    regress <- function(response, rows, fam) {
      model <- fit(df[rows, , drop = FALSE], response, xnames, fam,
                   learner = learner, SL.list = SL.list, SL.N = n)
      predict_fit(model, df[te, , drop = FALSE])
    }
    e[te] <- if (is.null(known_pZ1)) guard(regress("Z", tr, binomial())) else
      known_pZ1[te]
    m0[te] <- regress("Y", tr0, family)
    m1[te] <- regress("Y", tr1, family)
    p1[te] <- guard(regress("A", tr1, binomial()))
  }
  if (any(!is.finite(c(e, m0, m1, p1)))) stop("Non-finite nuisance prediction.")
  list(pZ1 = e, mY_z0 = m0, mY_z1 = m1,
       mA_z0 = numeric(n), mA_z1 = p1, folds = folds)
}

lsiv_one_sided_ratio <- function(numerator, denominator) {
  den <- mean(denominator)
  if (!is.finite(den) || abs(den) < 1e-6) stop("Degenerate score denominator.")
  estimate <- mean(numerator) / den
  eif <- (numerator - estimate * denominator) / den
  list(estimate = estimate, se = sd(eif) / sqrt(length(eif)), eif = eif,
       numerator = numerator, denominator = denominator,
       first_stage = den, ok = TRUE, converged = TRUE)
}

## Finite boundary limit of the unrestricted LSIV ATT score. It does not require
## fitting outcome densities or evaluating the singular interior B,C system.
lsiv_one_sided_score <- function(Y, A, Z, nuisance) {
  e <- nuisance$pZ1
  lsiv_one_sided_ratio((Z - e) / (1 - e) * (Y - nuisance$mY_z0), A)
}

## Projection for a known constant randomized instrument. This is exactly the
## LATE ratio score times e, including its denominator, when A<=Z.
lsiv_one_sided_projected_score <- function(Y, A, Z, nuisance) {
  e <- nuisance$pZ1
  if (any(A > Z) || max(abs(e - e[1L])) > 1e-12) {
    stop("The projected identity requires A<=Z and a constant instrument probability.")
  }
  delta <- nuisance$mY_z1 - nuisance$mY_z0
  lsiv_one_sided_ratio(
    (Z - e) / (1 - e) * (Y - nuisance$mY_z0) - (Z - e) * delta,
    A - (Z - e) * nuisance$mA_z1)
}

## Paired cross-fitted estimates. num.crossfit controls repeated splits; each
## method uses the median estimate and the usual median variance adjustment.
## known_pZ1=NULL estimates e(X); known_pZ1=0.5 uses the randomized design.
## score="boundary" retains the original unrestricted-model score limit.
## score="randomized" uses its projection in the known randomization model.
## With the same nuisance fits, the latter equals the LATE estimator exactly.
lsiv_compare_one_sided <- function(Y, A, Z, X,
                                   type = c("continuous", "binary"),
                                   learner = c("SL", "glm"), SL.list = 1:9,
                                   crossfit = 2L, num.crossfit = 1L,
                                   known_pZ1 = NULL, nuisance = NULL,
                                   score = c("boundary", "randomized")) {
  type <- match.arg(type)
  learner <- match.arg(learner)
  score <- match.arg(score)
  if (score == "randomized" && is.null(nuisance) && is.null(known_pZ1)) {
    stop("The randomized score requires a known constant instrument probability.")
  }
  if (length(num.crossfit) != 1L || !is.finite(num.crossfit) ||
      num.crossfit < 1 || num.crossfit != as.integer(num.crossfit)) {
    stop("num.crossfit must be a positive integer.")
  }
  if (!is.null(nuisance) && num.crossfit != 1L) {
    stop("A supplied nuisance object represents one cross-fitting split.")
  }
  if (any(A > Z)) stop("One-sided comparison requires A<=Z.")
  late_score <- getFromNamespace("lsiv_late_eif", "LSIV")
  splits <- lapply(seq_len(num.crossfit), function(s) {
    nu <- if (is.null(nuisance)) lsiv_one_sided_nuisance(
      Y, A, Z, X, type, learner, SL.list, crossfit,
      known_pZ1 = known_pZ1) else nuisance
    if (any(nu$mA_z0 != 0)) stop("The Z=0 treatment probability must be exactly zero.")
    lsiv <- if (score == "boundary") lsiv_one_sided_score(Y, A, Z, nu) else
      lsiv_one_sided_projected_score(Y, A, Z, nu)
    raw_late <- late_score(Y, A, Z, nu$pZ1, nu$mY_z0, nu$mY_z1,
                          nu$mA_z0, nu$mA_z1)
    late <- list(estimate = raw_late$late, se = sqrt(raw_late$var),
                 eif = raw_late$eif, first_stage = raw_late$first_stage,
                 ok = is.finite(raw_late$late), converged = TRUE)
    list(LSIV = lsiv, LSIV_LATE = late, nuisance = nu)
  })
  combine <- function(method) {
    estimates <- vapply(splits, function(x) x[[method]]$estimate, numeric(1))
    variances <- vapply(splits, function(x) x[[method]]$se^2, numeric(1))
    if (any(!is.finite(c(estimates, variances)))) stop("A split-specific fit failed.")
    estimate <- median(estimates)
    list(estimate = estimate,
         se = sqrt(median(variances + (estimates - estimate)^2)),
         first_stage = median(vapply(splits, function(x) x[[method]]$first_stage,
                                      numeric(1))),
         ok = TRUE, converged = TRUE, split_estimates = estimates,
         split_variances = variances)
  }
  list(LSIV = combine("LSIV"), LSIV_LATE = combine("LSIV_LATE"),
       splits = splits, score = score)
}

## Validation DGPs from the competing papers, used by 02_check_competitors.R.
lsiv_dgp_liu2021 <- function(n) {
  C1 <- rbinom(n, 1, 0.4)
  C2 <- rbinom(n, 1, 0.6)
  Z <- rbinom(n, 1, plogis(0.2 + 0.4 * C1 - 0.5 * C2))
  Y0 <- rbinom(n, 1, plogis(0.6 + 0.8 * C1 - 2 * C2))
  Y1 <- rbinom(n, 1, plogis(0.7 - 0.3 * C1))
  A <- rbinom(n, 1,
              plogis(0.4 + 2 * Z + 0.8 * C1 - 0.6 * Y0 - 1.6 * C1 * Z))
  Y <- Y0 * (1 - A) + Y1 * A
  list(Y = Y, A = A, Z = Z, X = cbind(C1 = C1, C2 = C2),
       Y0 = Y0, Y1 = Y1, psi_true = -0.6)
}

lsiv_dgp_liu2021_cont <- function(n) {
  C1 <- rbinom(n, 1, 0.4)
  C2 <- rbinom(n, 1, 0.6)
  Z <- rbinom(n, 1, plogis(0.7 + 0.8 * C1 - C2))
  Y0 <- rnorm(n, 0.5 + C1 + 3 * C2, 1)
  Y1 <- rnorm(n, 1.1 - 1.3 * C1, 1)
  A <- rbinom(n, 1,
              plogis(-0.2 - 3 * Z - 3 * C1 + 0.3 * Y0 + 4 * C1 * Z))
  Y <- Y0 * (1 - A) + Y1 * A
  list(Y = Y, A = A, Z = Z, X = cbind(C1 = C1, C2 = C2),
       Y0 = Y0, Y1 = Y1, psi_true = 0.3)
}

lsiv_dgp_sun2021 <- function(n) {
  X1 <- rbinom(n, 1, 0.4)
  X2 <- rbinom(n, 1, 0.6)
  Z <- rbinom(n, 1,
              plogis(0.4 + 0.9 * X1 - 0.7 * X2 - 0.8 * X1 * X2))
  Y <- rbinom(n, 1, plogis(1 - 1.2 * X1 + 1.5 * X2))
  R <- rbinom(n, 1,
              plogis(-1.5 + 2.5 * Z + 0.8 * X1 - 1.2 * X2 + 1.8 * Y))
  list(Y = Y, Yobs = ifelse(R == 1, Y, NA_real_), R = R, Z = Z,
       X = cbind(X1 = X1, X2 = X2), gamma_true = 1.8, mu_true = 0.769)
}

lsiv_dgp_sun2025 <- function(n, type = c("binary", "continuous"),
                             gamma_true = 2) {
  type <- match.arg(type)
  if (type == "binary") {
    L <- chol(matrix(c(1, 0.2, 0.2, 1), 2, 2))
    U <- matrix(rnorm(2 * n), n, 2) %*% L
  } else {
    rtn <- function(k) qnorm(runif(
      k, pnorm(-1), pnorm(1)))
    U <- cbind(rtn(n), rtn(n))
  }
  U1 <- U[, 1L]; U2 <- U[, 2L]
  Z <- rbinom(n, 1,
              plogis(1 + 2 * U1 - U2 - 0.8 * U1 * U2))
  if (type == "binary") {
    Y <- rbinom(n, 1, plogis(0.5 - 2 * U1 + U2))
    R <- rbinom(n, 1,
                plogis(2 - 3 * Z + 0.8 * U1 + U2 + gamma_true * Y))
  } else {
    Y <- rnorm(n, 0.5 - 2 * U1 + U2, sqrt(0.25))
    R <- rbinom(n, 1,
                plogis(2 - 4 * Z + 0.8 * U1 + U2 + gamma_true * Y))
  }
  list(Y = Y, Yobs = ifelse(R == 1, Y, NA_real_), R = R, Z = Z,
       X = cbind(U1 = U1, U2 = U2), gamma_true = gamma_true, type = type)
}

## Integrate out both assignment and treatment. In the one-sided design the
## ATT weight is e(X)*pA1, whereas the LATE weight is pA1. Use E(Y1|X) when
## available to remove independent Y1 simulation noise from the reference value.
lsiv_att_from_draw <- function(d) {
  effect <- (if (is.null(d$mean_y1)) d$Y1 else d$mean_y1) - d$Y0
  weight <- d$pZ1 * d$pA1 + (1 - d$pZ1) * d$pA0
  sum(weight * effect) / sum(weight)
}

lsiv_att_truth <- function(type = c("continuous", "binary"),
                           n = 1000000L, seed = 20260902L, p = 3L,
                           one_sided = TRUE, assignment = "randomized") {
  ## n and seed are retained for existing callers; this truth uses integration.
  unname(lsiv_compute_truths(match.arg(type), p = p, one_sided = one_sided,
                             assignment = assignment)["ATT"])
}

## Reference estimators. These deliberately remain outside the LSIV package.
lsiv_fit_2sls <- function(Y, A, Z, X, compute_se = TRUE) {
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

lsiv_fit_linear_model <- function(Y, A, Z, X) {
  dat <- data.frame(Y = Y, A = A, Z = Z, as.data.frame(X))
  rhs <- paste(c("A", "Z", colnames(X)), collapse = " + ")
  fit <- lm(as.formula(paste("Y ~", rhs)), data = dat)
  c(estimate = unname(coef(fit)["A"]),
    se = unname(sqrt(diag(vcov(fit)))["A"]))
}

## The ignorability benchmark: the cross-fitted EIF estimator of tau_Ign,
##
##   tau_Ign = sum_k sum_{i in I_k} {A_i - (1-A_i) Pr^(-k)(A_i=1|L_i)
##                                          / Pr^(-k)(A_i=0|L_i)}
##                                 {Y_i - E^(-k)(Y_i | A_i=1, L_i)}
##             / sum_i A_i ,
##
## with L = (Z, X) and both nuisances fitted by `glm()` on the other K-1 folds.
##
## It is the ATT one would report if A were assumed independent of Y(0) given L
## -- i.e. ignoring the selection on the treatment-free outcome that the
## logit-separable IV model allows.  Its bias against the true ATT is therefore
## the price of that assumption, which is exactly what the comparison is for.
##
## Note that the outcome regression cancels from the estimand: for any g(L),
## E[(1-A) {Pr(A=1|L)/Pr(A=0|L)} g(L)] = E[Pr(A=1|L) g(L)] = E[A g(L)],
## so the two g-terms cancel and only the treatment model has to be correct.
## The regression affects efficiency, not consistency; it is included as
## written above, fitted among the treated.
lsiv_fit_ignorability <- function(Y, A, Z, X, crossfit = 2L, trim = 1e-2) {
  X <- as.matrix(X)
  n <- length(Y)
  crossfit <- max(1L, as.integer(crossfit))
  df <- data.frame(Y = Y, A = as.integer(A), Z = Z, as.data.frame(X))
  rhs <- paste(c("Z", colnames(X)), collapse = " + ")
  f_a <- stats::as.formula(paste("A ~", rhs))
  f_y <- stats::as.formula(paste("Y ~", rhs))
  
  ## Folds stratified on (A, Z) so that every training set contains treated
  ## units (needed by the outcome regression) and both instrument arms.
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
  ## Pr(A=1|L) -> 1 makes the odds diverge; bound it as LSIV_LATE() bounds its
  ## own propensities.
  pA <- pmin(pmax(pA, trim), 1 - trim)
  
  num <- (A - (1 - A) * pA / (1 - pA)) * (Y - m1)
  estimate <- sum(num) / sum(A)
  eif <- (num - estimate * A) / mean(A)
  c(estimate = estimate, se = stats::sd(eif) / sqrt(n))
}

lsiv_compute_truths <- function(type, n = 300000L, seed = 20260902L, p = 3L,
                                one_sided = TRUE, assignment = "randomized") {
  ## W~N(0,p), E{plogis(W)}=1/2, and compliance is independent of Y0,Y1
  ## given X. Integrate the conditional treatment effect over complier weights.
  ## The same deterministic reference is used by every server batch.
  spec <- lsiv_dgp_forall(type, p, one_sided = one_sided, assignment = assignment)
  average <- function(fun) integrate(function(w) fun(w) * dnorm(w, sd = sqrt(p)),
                                     -Inf, Inf, rel.tol = 1e-10)$value
  numerator <- average(function(w) spec$p_a1(w) *
                         (spec$mean_y1(w) - spec$mean_y0(w)))
  target <- numerator / 0.5
  c(ATT = target, LATE = target, TSLS_projection = target)
}

lsiv_summarise_results <- function(raw) {
  groups <- split(raw, interaction(raw$type, raw$method, drop = TRUE))
  out <- lapply(groups, function(g) {
    ## Preserve every finite estimate in the numerical summary. Diagnostics
    ## are flags, not missing-data instructions: an estimate that failed A3/A5
    ## must remain visible alongside n_ok/n_fail. Only genuine computation
    ## errors (non-finite estimate) are omitted from numerical summaries.
    computed <- is.finite(g$estimate)
    with_se <- computed & is.finite(g$se)
    ok <- g$ok & computed
    est <- g$estimate[computed]
    se <- g$se[with_se]
    truth <- g$truth[1L]
    data.frame(
      type = g$type[1L], method = g$method[1L], estimand = g$estimand[1L],
      truth = truth, n_computed = sum(computed),
      n_ok = sum(ok), n_fail = sum(!g$ok),
      mean = if (length(est)) mean(est) else NA_real_,
      bias = if (length(est) && is.finite(truth)) mean(est) - truth else NA_real_,
      empirical_sd = if (length(est) > 1L) sd(est) else NA_real_,
      rmse = if (length(est) && is.finite(truth))
        sqrt(mean((est - truth)^2)) else NA_real_,
      average_se = if (length(se)) mean(se) else NA_real_,
      coverage95 = if (any(with_se) && is.finite(truth))
        mean(abs(g$estimate[with_se] - truth) <=
               qnorm(0.975) * g$se[with_se]) else NA_real_,
      mean_seconds = mean(g$seconds),
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans[order(ans$type, ans$estimand, ans$method), ]
}
