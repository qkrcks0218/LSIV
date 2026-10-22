test_that("np/gss and separate-bandwidth KLIEP/uLSIF are the continuous defaults", {
  skip_if_not_installed("np")
  skip_if_not_installed("gss")
  skip_if_not_installed("densityratio")
  skip_on_cran()

  set.seed(7302)
  n <- 240
  X <- cbind(X1 = stats::rnorm(n))
  Z <- stats::rbinom(n, 1, stats::plogis(0.2 * X[, 1]))
  Y0 <- 0.5 * X[, 1] + stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(-0.8 + 1.4 * Z + 0.2 * Y0))
  Y <- ifelse(A == 0, Y0, Y0 + 1)

  fit <- suppressWarnings(SIV(
    Y, A, Z, X, y_type = "continuous",
    Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 1L,
    y_num.grid = 21L, outcome_include.observed = FALSE,
    ratio_kernels = 20L
  ))

  expect_identical(fit$outcome_density, "np_gss_cv")
  expect_identical(fit$ratio_bandwidth, "cv_separate")
  expect_true(fit$ratio_method %in%
                c("kliep_ulsif",
                  "logistic (kliep_ulsif fallback)"))
  expect_identical(fit$nuisance$baseline_density$method,
                   "cross-fitted np/gss MSPE ensemble")
  expect_identical(fit$nuisance$baseline_density$candidates, c("np", "gss"))
  expect_identical(fit$nuisance$baseline_density$bandwidth,
                   "cv.ml")
  expect_gte(min(fit$nuisance$fY_z0), 1e-4 - 1e-12)
  if (identical(fit$ratio_method, "kliep_ulsif")) {
    expect_identical(fit$nuisance$density_ratio_ensemble$candidates,
                     c("kliep", "ulsif"))
    expect_equal(rowSums(fit$nuisance$density_ratio_ensemble$weights), 1,
                 tolerance = 1e-8)
    expect_true(all(is.finite(
      fit$nuisance$density_ratio_ensemble$mspe_weights
    )))
    expect_null(fit$nuisance$density_ratio_ensemble$kl_weights)
    expect_equal(fit$nuisance$density_ratio_ensemble$weights,
                 fit$nuisance$density_ratio_ensemble$mspe_weights)
    expect_identical(fit$nuisance$density_ratio_ensemble$ensemble_loss, "mspe")
  }
  expect_equal(fit$nuisance$baseline_density$weights,
               fit$nuisance$baseline_density$mspe_weights,
               tolerance = 1e-8)
  expect_equal(fit$nuisance$baseline_density$mse_weights,
               fit$nuisance$baseline_density$mspe_weights)
  expect_null(fit$nuisance$baseline_density$kl_weights)
  expect_identical(fit$nuisance$baseline_density$ensemble_loss, "mspe")
  expect_equal(rowSums(fit$nuisance$baseline_density$weights), 1,
               tolerance = 1e-8)
  W <- matrix(fit$weights, n, length(fit$weights), byrow = TRUE)
  expect_equal(rowSums(fit$nuisance$fY_z0 * W), rep(1, n),
               tolerance = 1e-8)
  expect_equal(rowSums(fit$nuisance$fY_z1 * W), rep(1, n),
               tolerance = 1e-8)
  expect_true(is.finite(unname(stats::coef(fit))))
})

test_that("the default gss formula contains main effects only", {
  f <- SIV:::siv_gss_default_formula(c("X1", "X2"))
  expect_identical(attr(stats::terms(f), "term.labels"), c("Y", "X1", "X2"))
  expect_identical(attr(stats::terms(f), "order"), rep.int(1L, 3L))
})

test_that("np_cv uses the likelihood-CV bandwidth without a gss candidate", {
  skip_if_not_installed("np")
  skip_on_cran()

  set.seed(7305)
  n <- 100
  X <- cbind(X1 = stats::rnorm(n))
  Z <- rep(rep(0:1, each = 25L), 2L)
  A <- rep(c(rep(0, 50L), rep(1, 50L)))
  Y <- 0.4 * X[, 1L] + stats::rnorm(n)
  grid <- seq(-2, 2, length.out = 11L)
  weights <- rep(1, length(grid))

  fit <- siv_baseline_density_np_cv(
    Y, A, Z, X, grid, weights, folds = list(seq_len(n))
  )
  expect_identical(fit$bandwidth, "cv.ml")
  expect_identical(fit$candidates, "np")
  expect_identical(fit$method, "cross-fitted np cv.ml estimator")
  expect_equal(rowSums(fit$density), rep(1, n), tolerance = 1e-8)
})

test_that("the slow Version1 np/gss and CV KLIEP/uLSIF path remains available", {
  skip_if_not_installed("np")
  skip_if_not_installed("gss")
  skip_if_not_installed("densityratio")
  skip_on_cran()

  set.seed(7304)
  n <- 180
  X <- cbind(X1 = stats::rnorm(n))
  Z <- stats::rbinom(n, 1, stats::plogis(0.2 * X[, 1]))
  Y0 <- 0.5 * X[, 1] + stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(-1 + Z + 0.2 * Y0))
  Y <- ifelse(A == 0, Y0, Y0 + 1)

  fit <- suppressWarnings(SIV(
    Y, A, Z, X, y_type = "continuous",
    Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 1L,
    outcome_density = "np_gss_cv", ratio_method = "kliep_ulsif",
    ratio_bandwidth = "cv_common", outcome_cv.folds = 2L,
    outcome_gss.formula = ~Y + X1,
    y_num.grid = 15L, outcome_include.observed = FALSE,
    ratio_kernels = 15L, ratio_folds = 2L
  ))

  expect_identical(fit$outcome_density, "np_gss_cv")
  expect_true(fit$ratio_method %in%
                c("kliep_ulsif", "logistic (kliep_ulsif fallback)"))
  expect_identical(fit$ratio_bandwidth, "cv_common")
  expect_identical(fit$nuisance$baseline_density$method,
                   "cross-fitted np/gss MSPE ensemble")
  expect_identical(fit$nuisance$baseline_density$candidates, c("np", "gss"))
  expect_identical(fit$nuisance$baseline_density$bandwidth, "cv.ml")
  expect_equal(fit$nuisance$baseline_density$weights,
               fit$nuisance$baseline_density$mspe_weights)
  expect_null(fit$nuisance$baseline_density$kl_weights)
  expect_identical(fit$nuisance$baseline_density$ensemble_loss, "mspe")
})

test_that("MSPE weights recover the interior minimum missed by the penalty optimizer", {
  H <- matrix(c(0.274649277529102, 0.264281863900480,
                0.264281863900480, 0.269356147720356), 2L)
  h <- c(0.278153223208423, 0.268002610678459)
  fit <- SIV:::siv_density_mspe_weights(H, h)
  expect_equal(unname(fit$weight), c(0.985960021598615, 0.014039978401385),
               tolerance = 1e-10)
  score <- function(w) as.numeric(.5 * t(w) %*% H %*% w - sum(w * h))
  expect_lt(score(fit$weight), score(c(1, 0)) - 1e-6)
  expect_equal(fit$weight, SIV:::siv_ratio_ensemble_weights(H, h)$weight)
})

test_that("MSPE weights include exact endpoints and reject nonfinite scores", {
  expect_identical(SIV:::siv_density_mspe_weights(diag(2), c(2, 0))$weight,
                   c(1, 0))
  expect_identical(SIV:::siv_density_mspe_weights(diag(2), c(0, 2))$weight,
                   c(0, 1))
  flat <- SIV:::siv_density_mspe_weights(matrix(1, 2, 2), c(1, 1))
  expect_equal(sum(flat$weight), 1)
  expect_true(all(is.finite(flat$weight)))
  expect_error(SIV:::siv_density_mspe_weights(diag(2), c(NA_real_, 0)),
               "finite scores")
})

test_that("the np/gss baseline path validates its CV settings", {
  set.seed(7303)
  d <- test_dgp(200, "continuous")
  expect_error(
    SIV(d$Y, d$A, d$Z, d$X, y_type = "continuous",
        outcome_density = "np_gss_cv", outcome_cv.folds = 1L),
    "density_cv_folds"
  )
  expect_identical(
    eval(formals(SIV)$outcome_density),
    c("np_gss_cv", "np_normal", "np_cv", "gaussian", "binned")
  )

  g <- siv_grid(d$Y, d$A, "continuous", n_grid = 21L,
                include_observed = FALSE)
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "continuous")
  alias_fit <- siv_estimate(
    d$Y, d$A, d$Z, d$X, type = "continuous", nuisance = nu,
    grid = g$grid, weights = g$weights, outcome_density = "np_gss"
  )
  expect_identical(alias_fit$outcome_density, "np_gss_cv")
  expect_false("lower_density" %in% names(formals(SIV)))
})
