## ---------------------------------------------------------------------------
## The proposed estimator, verified separately for type = "binary" and
## type = "continuous" (Task 6 / Task 11).
## ---------------------------------------------------------------------------

test_that("SIV() exposes only the simplified public interface", {
  expect_identical(names(formals(SIV)), c(
    "Y", "A", "Z", "X", "estimand", "level",
    "Overall.Estimation_SL.list", "Overall.Estimation_num.crossfit",
    "Overall.Estimation_num.fold", "y_type", "y_grid", "y_weights",
    "y_num.grid", "y_ref", "outcome_density", "outcome_formula",
    "outcome_num.bins", "outcome_cv.folds", "outcome_gss.formula",
    "outcome_include.observed", "outcome_floor", "ratio_method",
    "ratio_bandwidth", "ratio_bandwidth.cv.grid",
    "ratio_kernels",
    "ratio_folds", "ratio_formula",
    "ratio_floor", "contraction_solver",
    "contraction_maxit", "contraction_tol", "QTT_quantile", "QTT_c1",
    "QTT_numBB", "falsification", "falsification_num.boot",
    "falsification_num.bin", "falsification_floor",
    "falsification_categorical", "diagnostic_warnings", "QTT_support",
    "QTT_bounds"
  ))
  expect_identical(formals(SIV)$outcome_floor, 1e-4)
  expect_identical(formals(SIV)$ratio_floor, 1e-4)
  expect_identical(formals(SIV)$Overall.Estimation_SL.list, 1L)
  expect_identical(formals(SIV)$Overall.Estimation_num.crossfit, 3L)
  expect_identical(formals(SIV)$Overall.Estimation_num.fold, 2L)
  expect_identical(formals(SIV)$y_num.grid, 101L)
  expect_null(formals(SIV)$y_ref)
  expect_identical(eval(formals(SIV)$ratio_bandwidth)[1L], "cv_separate")
  expect_identical(eval(formals(SIV)$ratio_bandwidth.cv.grid),
                   c(0.5, 0.75, 0.9))
  expect_false("ratio_ensemble.loss" %in% names(formals(SIV)))
  expect_identical(eval(formals(SIV)$outcome_density)[1L], "np_gss_cv")
  expect_identical(formals(siv_estimate)$density_floor, 1e-4)
  expect_identical(formals(siv_estimate)$dr_floor, 1e-4)
  expect_identical(eval(formals(siv_estimate)$ratio_bandwidth)[1L],
                   "cv_separate")
  expect_identical(eval(formals(siv_estimate)$ratio_bandwidth.cv.grid),
                   c(0.5, 0.75, 0.9))
  expect_false("ratio_ensemble.loss" %in% names(formals(siv_estimate)))
  expect_identical(eval(formals(siv_estimate)$outcome_density)[1L], "np_gss_cv")
  expect_identical(formals(siv_fixed_point)$density_floor, 1e-4)
  expect_identical(formals(siv_fixed_point)$dr_floor, 1e-4)

  set.seed(30)
  d <- test_dgp(100, "binary")
  expect_error(SIV(
    d$Y, d$A, d$Z, d$X,
    Overall.Estimation_num.crossfit = 0L
  ), "positive integer")
  expect_error(SIV(
    d$Y, d$A, d$Z, d$X,
    Overall.Estimation_num.fold = 1.5
  ), "positive integer")
})

test_that("median adjustment excludes non-finite split estimates", {
  estimates <- cbind(
    ATT = c(1, NA, 3),
    QTT = c(NA, NA, NA)
  )
  variances <- cbind(
    ATT = c(0.1, NA, 0.3),
    QTT = c(NA, NA, NA)
  )

  adjusted <- SIV:::siv_repeated_median_adjustment(estimates, variances)

  expect_equal(adjusted$estimate[["ATT"]], 2)
  expect_equal(adjusted$variance[["ATT"]], 1.2)
  expect_true(adjusted$adjustable[["ATT"]])
  expect_equal(colSums(adjusted$valid_estimate)[["ATT"]], 2)
  expect_true(is.na(adjusted$estimate[["QTT"]]))
  expect_true(is.na(adjusted$variance[["QTT"]]))
  expect_false(adjusted$adjustable[["QTT"]])
  expect_equal(SIV:::siv_finite_median(c(1, NA, Inf, 3)), 2)
})

test_that("SIV() defaults to binary and validates the outcome type", {
  set.seed(31)
  d <- test_dgp(200, "binary")
  ## match.arg rejects anything outside the documented set
  expect_identical(eval(formals(SIV)$y_type),
                   c("binary", "continuous", "discrete"))
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "ordinal"), "should be one of")
  ## the type is validated against the data
  dc <- test_dgp(200, "continuous")
  expect_error(SIV(dc$Y, dc$A, dc$Z, dc$X, y_type = "binary"), "not all 0/1")
  expect_warning(SIV(d$Y, d$A, d$Z, d$X, y_type = "continuous",
                     y_num.grid = 11,
                     Overall.Estimation_num.crossfit = 1L,
                     Overall.Estimation_num.fold = 1L,
                     outcome_density = "gaussian"),
                 "at most two distinct values")
})

test_that("binary outcome nuisance pools untreated observations and includes Z", {
  set.seed(311)
  d <- test_dgp(400, "binary")
  outcome_fit_calls <- 0L
  outcome_z <- integer(0)
  outcome_all_untreated <- outcome_uses_z <- logical(0)
  original_learn_fit <- siv_learn_fit
  testthat::local_mocked_bindings(
    siv_learn_fit = function(data, yname, xnames, ...) {
      if (identical(yname, "Y")) {
        outcome_fit_calls <<- outcome_fit_calls + 1L
        outcome_z <<- c(outcome_z, unique(data$Z))
        outcome_all_untreated <<- c(outcome_all_untreated, all(data$A == 0L))
        outcome_uses_z <<- c(outcome_uses_z, "Z" %in% xnames)
      }
      original_learn_fit(data, yname, xnames, ...)
    },
    .package = "SIV"
  )
  fit <- suppressWarnings(siv_estimate(
    d$Y, d$A, d$Z, d$X, type = "binary", learner = "glm", crossfit = 2L
  ))
  expect_identical(outcome_fit_calls, 2L)
  expect_identical(as.integer(table(outcome_z)), c(2L, 2L))
  expect_true(all(outcome_all_untreated))
  expect_true(all(outcome_uses_z))
  expect_equal(fit$nuisance$fY_z0[, 1L],
               1 - fit$nuisance$fY_z0[, 2L], tolerance = 1e-14)
  expect_equal(fit$nuisance$fY_z1[, 1L],
               1 - fit$nuisance$fY_z1[, 2L], tolerance = 1e-14)
})

test_that("SIV() derives integration weights from a supplied y_grid", {
  set.seed(310)
  n <- 220
  Z <- stats::rbinom(n, 1, 0.5)
  A <- stats::rbinom(n, 1, stats::plogis(-0.5 + Z))
  Y <- A + 0.4 * Z + stats::rnorm(n)
  grid <- sort(unique(c(0, seq(min(Y), max(Y), length.out = 21L))))
  fit <- suppressWarnings(SIV(
    Y, A, Z, X = NULL, y_type = "continuous", y_grid = grid,
    y_ref = 0,
    ratio_method = "logistic", Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 1L
  ))
  expect_identical(fit$grid, grid)
  expect_length(fit$weights, length(grid))
  expect_true(all(fit$weights > 0))
})

test_that("SIV(type = 'binary') is unbiased with oracle nuisance functions", {
  set.seed(32)
  est <- numeric(60)
  for (r in seq_along(est)) {
    d <- test_dgp(1500, "binary")
    g <- siv_grid(d$Y, d$A, "binary")
    nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "binary")
    f <- siv_estimate(d$Y, d$A, d$Z, d$X, type = "binary", nuisance = nu,
                      grid = g$grid, weights = g$weights)
    est[r] <- f$estimate
  }
  ## true ATT for the binary DGP is approximately 0.082 (manuscript, Section 7)
  expect_lt(abs(mean(est) - 0.082), 3 * sd(est) / sqrt(length(est)) + 0.005)
})

test_that("SIV(type = 'continuous') is unbiased with oracle nuisance functions", {
  set.seed(33)
  est <- se <- numeric(40)
  for (r in seq_along(est)) {
    d <- test_dgp(1500, "continuous")
    g <- siv_grid(d$Y, d$A, "continuous", n_grid = 121, include_observed = FALSE)
    nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "continuous")
    f <- siv_estimate(d$Y, d$A, d$Z, d$X, type = "continuous", nuisance = nu,
                      grid = g$grid, weights = g$weights)
    est[r] <- f$estimate; se[r] <- f$se
  }
  ## True ATT for the continuous DGP in Simulation_F_Both/0.DGP.R is
  ## approximately 1.043.  This locks the package DGP to the legacy simulation
  ## rather than to the short-lived smoothstep variant.
  expect_lt(abs(mean(est) - 1.043), 3 * sd(est) / sqrt(length(est)) + 0.02)
  ## the influence-function SE tracks the empirical SD
  expect_lt(abs(mean(se) / sd(est) - 1), 0.25)
})

test_that("SIV() returns a well-formed object with working methods", {
  set.seed(34)
  d <- test_dgp(400, "binary")
  f <- SIV(d$Y, d$A, d$Z, d$X, y_type = "binary")
  expect_identical(f$ratio_method, "plugin (automatic)")
  expect_s3_class(f, "SIV")
  expect_named(coef(f), "ATT")
  V <- vcov(f)
  expect_true(all(siv_check_vcov(V, p = 1)))
  expect_equal(sqrt(V[1, 1]), unname(siv_se(f)[1]))
  expect_equal(attr(siv_se(f), "type"), f$se_type)
  ci <- confint(f)
  expect_true(ci[1, 1] < coef(f) && coef(f) < ci[1, 2])
  s <- summary(f)
  expect_s3_class(s, "summary.siv_fit")
  expect_output(print(s), "influence function")
  expect_output(print(f), "Park & Tchetgen Tchetgen")
  ## the SE must NOT be advertised as a sandwich SE
  expect_false(grepl("sandwich", f$se_type))
})

test_that("SIV() applies the repeated-cross-fitting median adjustment", {
  set.seed(341)
  d <- test_dgp(400, "binary")
  f <- suppressWarnings(SIV(
    d$Y, d$A, d$Z, d$X, y_type = "binary",
    Overall.Estimation_num.crossfit = 3L,
    Overall.Estimation_num.fold = 2L
  ))
  ma <- f$median_adjustment
  tau <- ma$estimates[, "ATT"]
  sigma2 <- ma$variances[, "ATT"]
  tau_median <- stats::median(tau)
  expect_true(ma$applied)
  expect_identical(ma$num_crossfit, 3L)
  expect_identical(ma$num_fold, 2L)
  expect_length(ma$folds, 3L)
  expect_true(all(lengths(ma$folds) == 2L))
  expect_equal(unname(f$estimate["ATT"]), tau_median)
  expect_equal(unname(f$vcov["ATT", "ATT"]),
               stats::median(sigma2 + (tau - tau_median)^2))
})

test_that("SIV() resolves the automatic reference outcome", {
  set.seed(342)
  db <- test_dgp(250, "binary")
  fb <- suppressWarnings(SIV(
    db$Y, db$A, db$Z, db$X, y_type = "binary",
    Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 1L
  ))
  expect_identical(fb$y_ref, 0)

  dc <- test_dgp(250, "continuous")
  fc <- suppressWarnings(SIV(
    dc$Y, dc$A, dc$Z, dc$X, y_type = "continuous",
    Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 1L,
    y_num.grid = 31L, outcome_include.observed = FALSE,
    outcome_density = "gaussian", ratio_method = "logistic"
  ))
  expect_equal(fc$y_ref, mean(dc$Y[dc$A == 0]))
  expect_true(any(abs(fc$grid - fc$y_ref) < 1e-12))
})

test_that("SIV() supports discrete outcomes on their observed support", {
  set.seed(35)
  n <- 800
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("X1", "X2")
  W <- rowSums(X)
  Z <- rbinom(n, 1, expit(0.2 * W))
  Y0 <- sample(c(0, 1, 2), n, TRUE, c(0.4, 0.35, 0.25))
  A <- rbinom(n, 1, expit(0.5 * Y0 - 1 + 2 * Z + 0.1 * W))
  Y1 <- Y0 + 1
  Y <- (1 - A) * Y0 + A * Y1
  f <- SIV(Y, A, Z, X, y_type = "discrete")
  expect_equal(f$grid, c(0, 1, 2))
  expect_true(all(f$weights == 1))
  expect_true(is.finite(coef(f)))
  expect_true(all(siv_check_vcov(vcov(f), 1)))
})

test_that("SIV() cross-fitting runs and the solver choice is validated", {
  set.seed(36)
  d <- test_dgp(500, "binary")
  f <- SIV(d$Y, d$A, d$Z, d$X, y_type = "binary",
           Overall.Estimation_num.crossfit = 1L,
           Overall.Estimation_num.fold = 2L)
  expect_true(is.finite(coef(f)))
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "discrete",
                   contraction_solver = "closed-form"),
               "only available")
})

test_that("SIV() rejects malformed data", {
  set.seed(37)
  d <- test_dgp(200, "binary")
  expect_error(SIV(d$Y, d$A, rep(0, 200), d$X, y_type = "binary"), "both values")
  expect_error(SIV(d$Y, rep(1, 200), d$Z, d$X, y_type = "binary"), "both values")
  expect_error(SIV(d$Y, d$A, d$Z * 2, d$X, y_type = "binary"), "binary")
  expect_error(SIV(d$Y, d$A, d$Z, d$X[1:10, ], y_type = "binary"), "rows")
  nu <- list(pZ1 = rep(0.5, 200))
  expect_error(siv_estimate(d$Y, d$A, d$Z, d$X, type = "binary",
                            nuisance = nu), "missing element")
})

test_that("SIV() reports and warns about (A3) / (A5) diagnostics", {
  set.seed(38)
  d <- test_dgp(800, "binary")
  g <- siv_grid(d$Y, d$A, "binary")
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "binary")
  f <- siv_estimate(d$Y, d$A, d$Z, d$X, type = "binary", nuisance = nu,
                    grid = g$grid, weights = g$weights)
  expect_named(f$diagnostics,
               c("g_floor_frac", "min_relevance", "contraction_modulus"))
  ## with the true nuisance functions neither assumption is close to binding
  expect_equal(unname(f$diagnostics["g_floor_frac"]), 0)
  expect_gt(f$diagnostics[["min_relevance"]], 1e-2)
  ## the fixed point is far from the a -> 1 degeneracy
  expect_lt(f$diagnostics[["contraction_modulus"]], 0.99)

  ## A genuinely weak first stage trips the relevance diagnostic even after
  ## density normalization and the falsification-bound projection.
  bad <- nu
  bad$pA1_z1 <- bad$pA1_z0 + 1e-5
  ## ... but the warnings are OFF by default: only the diagnostics change
  quiet <- expect_silent(siv_estimate(
    d$Y, d$A, d$Z, d$X, type = "binary",
    nuisance = bad, grid = g$grid, weights = g$weights
  ))
  expect_lt(quiet$diagnostics[["min_relevance"]], 1e-3)

  ## ... and are recovered with `diagnostic_warnings = TRUE`
  loud <- NULL
  expect_warning(loud <- siv_estimate(
    d$Y, d$A, d$Z, d$X, type = "binary",
    nuisance = bad, grid = g$grid, weights = g$weights,
    diagnostic_warnings = TRUE
  ),
                 "Assumption|strong|modulus")
  ## opting into the warnings must not change the estimate
  expect_equal(stats::coef(loud), stats::coef(quiet))
})

test_that("the fixed point respects the (A3) lower bound", {
  set.seed(39)
  d <- test_dgp(300, "binary")
  g <- siv_grid(d$Y, d$A, "binary")
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "binary")
  ## a very aggressive floor must bind and must be respected
  fp <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1,
                        g$grid, g$weights, g_floor = 0.99)
  expect_gt(fp$g_floor_frac, 0)
  expect_true(all(fp$g >= 0.99 * apply(fp$g, 1, max) - 1e-10))
  expect_true(all(is.finite(log(fp$delta))))
  ## the default floor does not bind with well-behaved nuisance functions
  fp0 <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1,
                         g$grid, g$weights)
  expect_equal(fp0$g_floor_frac, 0)
})
