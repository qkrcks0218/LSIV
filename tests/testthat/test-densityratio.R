## ---------------------------------------------------------------------------
## The KLIEP density-ratio engine and the R_Y estimator of eq. (15).
## ---------------------------------------------------------------------------

test_that("KLIEP and KLIEP/uLSIF use the requested numerical floor", {
  expect_identical(formals(lsiv_kliep_ratio_grid)$dr_floor, 1e-4)
  expect_identical(formals(lsiv_kliep_ulsif_median_ratio_grid)$dr_floor, 1e-4)
  expect_identical(formals(lsiv_kliep_ulsif_ratio_grid)$dr_floor, 1e-4)
  expect_identical(
    formals(lsiv_kliep_ulsif_median_ratio_grid)$density_floor, 1e-4
  )
  expect_identical(formals(lsiv_kliep_ulsif_ratio_grid)$density_floor, 1e-4)
})

test_that("lsiv_kliep() recovers a known Gaussian density ratio", {
  ## p1 = N(mu, 1), p2 = N(0, 1) => log r(y) = mu*y - mu^2/2
  set.seed(401)
  mu <- 0.5
  x1 <- matrix(stats::rnorm(4000, mu)); x2 <- matrix(stats::rnorm(4000, 0))
  fit <- lsiv_kliep(x1, x2, n_y_cols = 1, sigma = "median")
  yy <- seq(-1.5, 1.5, by = 0.5)
  est <- lsiv_kliep_predict(matrix(yy), fit)
  truth <- stats::dnorm(yy, mu) / stats::dnorm(yy, 0)
  expect_true(all(est > 0))
  expect_true(all(is.finite(est)))
  ## KLIEP is a weighted sum of Gaussian kernels and is not constrained to be
  ## monotone, so only the trend and the accuracy on the bulk of the support
  ## are asserted.
  expect_gt(stats::cor(log(est), log(truth)), 0.95)
  expect_gt(est[length(est)], est[1L])
  expect_lt(max(abs(log(est) - log(truth))), 0.5)
})

test_that("the 'auto' bandwidth search improves on the median heuristic", {
  skip_on_cran()
  set.seed(402)
  mu <- 1.0
  x1 <- matrix(stats::rnorm(4000, mu)); x2 <- matrix(stats::rnorm(4000, 0))
  yy <- seq(-1.5, 1.5, by = 0.5)
  truth <- log(stats::dnorm(yy, mu) / stats::dnorm(yy, 0))
  err <- function(sg) {
    set.seed(403)
    f <- lsiv_kliep(x1, x2, n_y_cols = 1, sigma = sg)
    max(abs(log(lsiv_kliep_predict(matrix(yy), f)) - truth))
  }
  ## the cross-validated search is the hyper-parameter tuning; it must not be
  ## worse than the median heuristic on a problem where the heuristic is poor
  expect_lt(err("auto"), err("median"))
})

test_that("the bandwidth search is reproducible under set.seed()", {
  ## The CV fold labels are drawn with R's RNG, not arma::randperm(), precisely
  ## so that this holds.
  set.seed(404); a <- matrix(stats::rnorm(300, 0.5)); b <- matrix(stats::rnorm(300))
  set.seed(405); f1 <- lsiv_kliep(a, b, sigma = "auto")
  set.seed(405); f2 <- lsiv_kliep(a, b, sigma = "auto")
  expect_identical(f1$kernel_info$sigma, f2$kernel_info$sigma)
  expect_equal(f1$kernel_weights, f2$kernel_weights)
  ## and the median path too (kernel centres are sampled with R's RNG)
  set.seed(406); g1 <- lsiv_kliep(a, b, sigma = "median")
  set.seed(406); g2 <- lsiv_kliep(a, b, sigma = "median")
  expect_equal(g1$kernel_weights, g2$kernel_weights)
})

test_that("lsiv_kliep() validates its inputs", {
  a <- matrix(stats::rnorm(50)); b <- matrix(stats::rnorm(50))
  expect_error(lsiv_kliep(a, cbind(b, b)), "same number of columns")
  expect_error(lsiv_kliep(a[0, , drop = FALSE], b), "non-empty")
  expect_error(lsiv_kliep(a, b, n_y_cols = 2), "n_y_cols")
})

test_that("the per-dimension kernel helpers behave as documented", {
  set.seed(407)
  x <- matrix(stats::rnorm(20), 10, 2)
  ctr <- matrix(stats::rnorm(8), 4, 2)
  sg <- c(0.7, 1.3)
  K <- LSIV:::lsiv_gaussian_kernel_matrix_bw(x, ctr, sg)
  expect_equal(dim(K), c(10L, 4L))
  ## reference computation
  ref <- outer(seq_len(10), seq_len(4), Vectorize(function(i, c)
    exp(-sum((x[i, ] - ctr[c, ])^2 / (2 * sg^2)))))
  expect_equal(K, ref, tolerance = 1e-12)
  ## the ratio is the kernel matrix times the weights
  w <- c(0.1, 0.4, 0.2, 0.3)
  expect_equal(as.numeric(LSIV:::lsiv_gaussian_kernel_ratio_bw(x, ctr, sg, w)),
               as.numeric(K %*% w), tolerance = 1e-12)
  ## scalar-bandwidth path agrees with the per-dimension one at equal sigmas
  expect_equal(as.numeric(LSIV:::lsiv_gaussian_kernel_ratio(x, ctr, 0.9, w)),
               as.numeric(LSIV:::lsiv_gaussian_kernel_ratio_bw(x, ctr, c(0.9, 0.9), w)),
               tolerance = 1e-12)
})

test_that("LSIV(ratio_method = 'kliep') runs and reports what it used", {
  skip_on_cran()
  set.seed(408)
  d <- test_dgp(1500, "continuous")
  f <- suppressWarnings(LSIV(
    d$Y, d$A, d$Z, d$X, y_type = "continuous",
    y_num.grid = 61L, outcome_include.observed = FALSE,
    Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 2L, outcome_density = "gaussian",
    ratio_method = "kliep"
  ))
  expect_true(is.finite(coef(f)))
  expect_true(is.finite(lsiv_se(f)))
  expect_true(f$ratio_method %in%
                c("kliep", "logistic (kliep fallback)"))
  if (identical(f$ratio_method, "kliep")) {
    expect_false(is.null(f$nuisance$RY))
    expect_true(all(is.finite(f$nuisance$RY)))
    expect_true(all(f$nuisance$RY > 0))
  }
})

test_that("the KLIEP/uLSIF ratio estimator ensembles both candidates", {
  skip_if_not_installed("densityratio")
  skip_on_cran()
  set.seed(413)
  d <- test_dgp(1000, "continuous")
  gw <- lsiv_grid(d$Y, d$A, "continuous", n_grid = 41L,
                  include_observed = FALSE)
  folds <- lsiv_make_folds(nrow(d$X), d$A, d$Z, 2L)
  oracle <- test_oracle_nuisance(d$X, gw$grid, gw$weights, "continuous")
  fit <- lsiv_kliep_ulsif_ratio_grid(
    d$Y, d$A, d$Z, d$X, gw$grid, gw$weights, folds,
    oracle$fY_z0, kernel_num = 40L, fold = 2L
  )
  expect_identical(fit$candidates, c("kliep", "ulsif"))
  expect_equal(rowSums(fit$weights),
               rep(1, 2), tolerance = 1e-8)
  expect_length(fit$selected_bandwidth, 2L)
  expect_true(all(vapply(fit$selected_bandwidth, function(bw) {
    identical(dim(bw), c(2L, 2L)) && all(is.finite(bw)) && all(bw > 0)
  }, logical(1L))))
  expect_true(all(vapply(fit$bandwidth_range, function(x) {
    length(x$h_y) == 3L && length(x$h_x) == 3L &&
      identical(as.numeric(x$probs), c(0.5, 0.75, 0.9))
  }, logical(1L))))
  expect_identical(fit$cv_grid, c(0.5, 0.75, 0.9))
  expect_true(all(is.finite(fit$ratio)))
  W <- matrix(gw$weights, nrow(d$X), length(gw$weights), byrow = TRUE)
  expect_equal(rowSums(fit$ratio * oracle$fY_z0 * W), rep(1, nrow(d$X)),
               tolerance = 1e-8)
})

test_that("separate-bandwidth preprocessing standardises every X column", {
  set.seed(415)
  x1 <- data.frame(Y = stats::rnorm(80), X1 = stats::rnorm(80, sd = 2),
                   X2 = stats::rnorm(80, sd = 10))
  x0 <- data.frame(Y = stats::rnorm(100), X1 = stats::rnorm(100, sd = 3),
                   X2 = stats::rnorm(100, sd = 20))
  tr <- lsiv_kliep_ulsif_transform_info(x1, x0)
  z <- lsiv_kliep_ulsif_transform(
    rbind(x1, x0), tr, h_y = 0.7, h_x = 2
  )
  expect_equal(unname(apply(z[c("X1", "X2")], 2L, stats::sd)), c(0.5, 0.5),
               tolerance = 1e-12)
  expect_equal(z$Y[1L], (x1$Y[1L] - tr$y_center) / 0.7)
})

test_that("common-bandwidth preprocessing standardises Y and every X column", {
  skip_if_not_installed("densityratio")
  set.seed(416)
  x1 <- data.frame(Y = stats::rnorm(80, sd = 7),
                   X1 = stats::rnorm(80, sd = 2),
                   X2 = stats::rnorm(80, sd = 10))
  x0 <- data.frame(Y = stats::rnorm(100, sd = 4),
                   X1 = stats::rnorm(100, sd = 3),
                   X2 = stats::rnorm(100, sd = 20))
  tr <- lsiv_kliep_ulsif_common_transform_info(x1, x0)
  z <- lsiv_kliep_ulsif_common_transform(rbind(x1, x0), tr)
  expect_equal(unname(colMeans(z)), rep(0, 3L), tolerance = 1e-12)
  expect_equal(unname(apply(z, 2L, stats::sd)), rep(1, 3L),
               tolerance = 1e-12)

  cv_fit <- lsiv_kliep_ulsif_YXcommonbw_fit(
    x1, x0, kernel_num = 20L, fold = 2L, methods = "ulsif"
  )
  median_fit <- lsiv_kliep_ulsif_median_fit(
    x1, x0, kernel_num = 20L, methods = "ulsif"
  )
  expect_s3_class(cv_fit$fits$ulsif, "lsiv_kliep_ulsif_common")
  expect_s3_class(median_fit$fits$ulsif, "lsiv_kliep_ulsif_common")
  expect_equal(cv_fit$transform$scale, tr$scale)
  expect_equal(median_fit$transform$scale, tr$scale)
  expect_length(cv_fit$sigma_candidates, 3L)
  expect_true(all(is.finite(cv_fit$sigma_candidates)))
  expect_true(all(cv_fit$sigma_candidates > 0))
  expect_true(all(is.finite(median_fit$sigma_candidates)))
  expect_true(all(median_fit$sigma_candidates > 0))
})

test_that("the fast ratio estimator CV-weights both median candidates", {
  skip_if_not_installed("densityratio")
  skip_on_cran()
  set.seed(414)
  d <- test_dgp(800, "continuous")
  gw <- lsiv_grid(d$Y, d$A, "continuous", n_grid = 31L,
                  include_observed = FALSE)
  folds <- lsiv_make_folds(nrow(d$X), d$A, d$Z, 2L)
  oracle <- test_oracle_nuisance(d$X, gw$grid, gw$weights, "continuous")
  fit <- lsiv_kliep_ulsif_median_ratio_grid(
    d$Y, d$A, d$Z, d$X, gw$grid, gw$weights, folds,
    oracle$fY_z0, kernel_num = 30L
  )
  expect_identical(fit$candidates, c("kliep", "ulsif"))
  expect_identical(fit$ratio_method, "kliep_ulsif")
  expect_identical(fit$ratio_bandwidth, "median_common")
  expect_identical(fit$ensemble_loss, "mspe")
  expect_equal(rowSums(fit$weights), rep(1, 2L), tolerance = 1e-8)
  expect_true(all(is.finite(fit$weights)))
  expect_true(all(is.finite(fit$mspe_weights)))
  expect_equal(fit$weights, fit$mspe_weights)
  expect_null(fit$kl_weights)
  expect_true(all(is.finite(fit$selected_sigma)))
  expect_true(all(fit$selected_sigma > 0))
  expect_true(all(is.finite(fit$ratio)))
})

test_that("the ratio interface validates its arguments", {
  set.seed(409)
  d <- test_dgp(300, "continuous")
  expect_error(LSIV(d$Y, d$A, d$Z, d$X, y_type = "continuous",
                    ratio_method = "nonesuch"), "should be one of")
  expect_error(LSIV(d$Y, d$A, d$Z, d$X, y_type = "continuous",
                    ratio_bandwidth = "nonesuch"), "should be one of")
  expect_error(LSIV(d$Y, d$A, d$Z, d$X, y_type = "continuous",
                    ratio_ensemble.loss = "ise"), "unused argument")
  expect_error(LSIV(
    d$Y, d$A, d$Z, d$X, y_type = "continuous",
    ratio_bandwidth.cv.grid = c(0, 0.5)
  ), "strictly between 0 and 1")
  expect_false("estimate_ratio" %in% names(formals(LSIV)))
  expect_false("plugin" %in% eval(formals(LSIV)$ratio_method))

  ## Binary outcomes always use the automatic plug-in ratio.
  db <- test_dgp(300, "binary")
  f <- LSIV(db$Y, db$A, db$Z, db$X, y_type = "binary",
            ratio_method = "ulsif")
  expect_identical(f$ratio_method, "plugin (automatic)")
})

test_that("ratio ensemble minimizes MSPE without a loss-selection option", {
  H <- diag(2)
  h <- c(0.8, 0.2)
  fit <- LSIV:::lsiv_ratio_ensemble_weights(H, h)

  expect_equal(fit$weight, c(0.8, 0.2), tolerance = 1e-12)
  expect_equal(fit$weight, fit$mspe)
  expect_equal(fit$mse, fit$mspe)
  expect_identical(fit$loss, "mspe")
  expect_error(LSIV:::lsiv_ratio_ensemble_weights(H, h, loss = "logscore"),
               "unused argument")
})

test_that("standalone uLSIF uses the selected bandwidth rule", {
  skip_if_not_installed("densityratio")
  set.seed(411)
  d <- test_dgp(300, "continuous")
  f <- suppressWarnings(LSIV(
    d$Y, d$A, d$Z, d$X, y_type = "continuous",
    Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 1L,
    y_num.grid = 21L, outcome_include.observed = FALSE,
    outcome_density = "gaussian", ratio_method = "ulsif",
    ratio_bandwidth = "median_separate", ratio_kernels = 20L
  ))
  expect_identical(f$ratio_method, "ulsif")
  expect_identical(f$ratio_bandwidth, "median_separate")
  expect_identical(f$nuisance$density_ratio_ensemble$candidates, "ulsif")
})

ratio_direction_data <- function() {
  d <- expand.grid(replicate = seq_len(240L), Z = 0:1, X1 = c(-1, 1))
  d$A <- as.integer(d$replicate <= ifelse(d$Z == as.integer(d$X1 > 0), 192, 48))
  d$Y <- seq(-0.9, 0.9, length.out = nrow(d))
  list(Y = d$Y, A = d$A, Z = d$Z, X = cbind(X1 = d$X1),
       folds = list(seq(1L, nrow(d), by = 2L), seq(2L, nrow(d), by = 2L)))
}

test_that("direct ratio fits follow the covariate-specific first-stage direction", {
  d <- ratio_direction_data()
  baseline_calls <- ratio_calls <- list()
  baseline_fit <- function(Y, A, Z, X, grid, weights, folds, ...) {
    direction <- if (all(Z == d$Z)) "1" else "0"
    baseline_calls[[direction]] <<- list(Z = Z, folds = folds)
    base <- if (direction == "1") c(0.5, 0.5) else c(0.25, 0.75)
    list(density = matrix(base, length(Y), 2L, byrow = TRUE),
         method = "controlled baseline", candidates = "np")
  }
  ratio_fit <- function(Y, A, Z, X, grid, weights, folds,
                        baseline_density, fold = 5L, ...) {
    direction <- if (all(Z == d$Z)) "1" else "0"
    ratio_calls[[direction]] <<- list(Z = Z, folds = folds,
                                     base = baseline_density)
    ratio <- if (direction == "1") c(0.4, 1.6) else c(1.6, 0.8)
    list(ratio = matrix(ratio, length(Y), 2L, byrow = TRUE),
         method = "controlled direct ratio", candidates = c("kliep", "ulsif"),
         weights = matrix(0.5, length(folds), 2L), ratio_method = "kliep_ulsif")
  }
  testthat::local_mocked_bindings(
    lsiv_baseline_density_np_normal = baseline_fit,
    lsiv_kliep_ulsif_ratio_grid = ratio_fit,
    lsiv_kliep_ulsif_median_ratio_grid = ratio_fit,
    lsiv_single_density_ratio_grid = ratio_fit,
    .package = "LSIV"
  )
  settings <- list(
    c("kliep_ulsif", "cv_separate"), c("kliep_ulsif", "cv_common"),
    c("kliep_ulsif", "median_common"), c("kliep_ulsif", "median_separate"),
    c("kliep", "cv_separate"), c("ulsif", "median_separate")
  )
  for (setting in settings) {
    fit <- suppressWarnings(lsiv_estimate(
      d$Y, d$A, d$Z, d$X, type = "continuous", grid = c(-1, 1),
      weights = c(1, 1), y_ref = -1, folds = d$folds,
      learner = "glm", treatment_formula = A ~ Z * X1,
      outcome_density = "np_normal", ratio_method = setting[1L],
      ratio_bandwidth = setting[2L]
    ))
    forward <- d$X[, 1L] > 0
    reverse <- !forward
    expect_identical(fit$zM, as.integer(forward))
    expect_setequal(names(baseline_calls), c("0", "1"))
    expect_setequal(names(ratio_calls), c("0", "1"))
    expect_equal(ratio_calls[["1"]]$Z, d$Z)
    expect_equal(ratio_calls[["0"]]$Z, 1 - d$Z)
    expect_identical(ratio_calls[["1"]]$folds, d$folds)
    expect_identical(ratio_calls[["0"]]$folds, d$folds)
    expect_equal(ratio_calls[["1"]]$base[1L, ], c(0.5, 0.5))
    expect_equal(ratio_calls[["0"]]$base[1L, ], c(0.25, 0.75))
    expect_equal(fit$nuisance$fY_z0[forward, ],
                 matrix(c(0.5, 0.5), sum(forward), 2L, byrow = TRUE))
    expect_equal(fit$nuisance$fY_z1[forward, ],
                 matrix(c(0.2, 0.8), sum(forward), 2L, byrow = TRUE))
    expect_equal(fit$nuisance$fY_z0[reverse, ],
                 matrix(c(0.4, 0.6), sum(reverse), 2L, byrow = TRUE))
    expect_equal(fit$nuisance$fY_z1[reverse, ],
                 matrix(c(0.25, 0.75), sum(reverse), 2L, byrow = TRUE))
    expect_equal(fit$nuisance$RY[reverse, ],
                 matrix(c(0.625, 1.25), sum(reverse), 2L, byrow = TRUE))
    expect_equal(fit$nuisance$RY,
                 fit$nuisance$fY_z1 / fit$nuisance$fY_z0)
    expect_equal(rowSums(fit$estimation_nuisance$fY_z0), rep(1, length(d$Y)))
    expect_equal(rowSums(fit$estimation_nuisance$fY_z1), rep(1, length(d$Y)))
    expect_identical(fit$nuisance$density_ratio_ensemble$by_zM[["0"]]$numerator_z,
                     0L)
    expect_identical(fit$nuisance$baseline_density$by_zM[["0"]]$denominator_z,
                     1L)
    expect_true(all(is.finite(fit$eif)))
  }
})

test_that("a failed reverse ratio fit preserves the forward direct fit", {
  d <- ratio_direction_data()
  testthat::local_mocked_bindings(
    lsiv_baseline_density_np_normal = function(Y, ...) {
      list(density = matrix(0.5, length(Y), 2L),
           method = "controlled baseline", candidates = "np")
    },
    lsiv_kliep_ulsif_ratio_grid = function(Y, A, Z, ...) {
      if (all(Z == 1 - d$Z)) return(NULL)
      list(ratio = matrix(c(0.4, 1.6), length(Y), 2L, byrow = TRUE),
           method = "controlled direct ratio", candidates = c("kliep", "ulsif"),
           ratio_method = "kliep_ulsif")
    },
    .package = "LSIV"
  )
  expect_warning(fit <- lsiv_estimate(
    d$Y, d$A, d$Z, d$X, type = "continuous", grid = c(-1, 1),
    weights = c(1, 1), y_ref = -1, folds = d$folds,
    learner = "glm", treatment_formula = A ~ Z * X1,
    outcome_density = "np_normal", ratio_method = "kliep_ulsif"
  ), "falling back to the logistic ratio")
  expect_identical(fit$ratio_method, "kliep_ulsif + logistic (partial fallback)")
  expect_identical(fit$nuisance$ratio_method_by_zM[["1"]], "kliep_ulsif")
  expect_identical(fit$nuisance$ratio_method_by_zM[["0"]],
                   "logistic (kliep_ulsif fallback)")
  expect_equal(fit$nuisance$RY[d$X[, 1L] > 0, ],
               matrix(c(0.4, 1.6), sum(d$X[, 1L] > 0), 2L, byrow = TRUE))
  expect_equal(rowSums(fit$nuisance$fY_z0), rep(1, length(d$Y)))
  expect_equal(rowSums(fit$nuisance$fY_z1), rep(1, length(d$Y)))
})

test_that("uLSIF fits both directions when the conditional first stage changes sign", {
  skip_if_not_installed("densityratio")
  set.seed(20261112)
  n <- 480L
  X <- cbind(X1 = stats::runif(n, -1, 1))
  Z <- stats::rbinom(n, 1, 0.5)
  A <- stats::rbinom(n, 1, stats::plogis(-0.3 + 5 * (2 * Z - 1) * X[, 1L]))
  Y <- 0.5 * X[, 1L] + A + stats::rnorm(n)
  fit <- suppressWarnings(lsiv_estimate(
    Y, A, Z, X, type = "continuous", learner = "glm", crossfit = 2L,
    treatment_formula = A ~ Z * X1, outcome_density = "gaussian",
    ratio_method = "ulsif", ratio_bandwidth = "median_separate",
    ratio_kernels = 15L, ratio_folds = 2L, n_grid = 21L,
    include_observed = FALSE
  ))
  expect_setequal(unique(fit$zM), 0:1)
  expect_identical(fit$ratio_method, "ulsif")
  expect_setequal(names(fit$nuisance$density_ratio_ensemble$by_zM), c("0", "1"))
  expect_true(all(is.finite(fit$nuisance$RY)))
  expect_true(all(fit$nuisance$RY > 0))
  expect_equal(fit$nuisance$RY,
               fit$nuisance$fY_z1 / fit$nuisance$fY_z0)
  expect_equal(as.vector(fit$nuisance$fY_z0 %*% fit$weights), rep(1, n),
               tolerance = 1e-10)
  expect_equal(as.vector(fit$nuisance$fY_z1 %*% fit$weights), rep(1, n),
               tolerance = 1e-10)
  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(is.finite(fit$vcov)))
})

test_that("direct estimation is invariant to reversing instrument labels", {
  skip_if_not_installed("densityratio")
  set.seed(20261113)
  n <- 360L
  X <- cbind(X1 = stats::rnorm(n))
  Z <- stats::rbinom(n, 1, 0.5)
  A <- stats::rbinom(n, 1, stats::plogis(-0.7 + 2 * Z + 0.2 * X[, 1L]))
  Y <- 0.4 * X[, 1L] + A + stats::rnorm(n)
  folds <- lsiv_make_folds(n, A, Z, 2L)
  run <- function(instrument) {
    set.seed(20261114)
    suppressWarnings(lsiv_estimate(
      Y, A, instrument, X, type = "continuous", learner = "glm", folds = folds,
      outcome_density = "gaussian", ratio_method = "ulsif",
      ratio_bandwidth = "median_separate", ratio_kernels = 15L,
      ratio_folds = 2L, n_grid = 21L, include_observed = FALSE
    ))
  }
  forward <- run(Z)
  reverse <- run(1 - Z)
  expect_true(all(forward$zM == 1L))
  expect_true(all(reverse$zM == 0L))
  expect_identical(names(forward$nuisance$density_ratio_ensemble$by_zM), "1")
  expect_identical(names(reverse$nuisance$density_ratio_ensemble$by_zM), "0")
  expect_equal(forward$nuisance$fY_z0, reverse$nuisance$fY_z1,
               tolerance = 1e-8)
  expect_equal(forward$nuisance$fY_z1, reverse$nuisance$fY_z0,
               tolerance = 1e-8)
  expect_equal(forward$estimate, reverse$estimate, tolerance = 1e-6)
  expect_equal(forward$vcov, reverse$vcov, tolerance = 1e-6)
})

test_that("reciprocal fallback uses training counts below 100 in either direction", {
  for (flip in c(FALSE, TRUE)) {
    n <- 600L
    X <- cbind(X1 = c(rep(-1, 100L), rep(1, 200L),
                     rep(-1, 99L), rep(1, 201L)))
    Z <- rep(0:1, length.out = n)
    A <- integer(n)
    for (xx in c(-1, 1)) {
      for (zz in 0:1) {
        for (fold in 0:1) {
          ii <- which(X[, 1L] == xx & Z == zz &
                        as.integer(seq_len(n) > 300L) == fold)
          rate <- if (zz == as.integer(xx > 0)) 0.8 else 0.2
          A[ii[seq_len(round(length(ii) * rate))]] <- 1L
        }
      }
    }
    if (flip) Z <- 1L - Z
    observed_Z <- Z
    Y <- seq(-0.9, 0.9, length.out = n)
    folds <- list(seq_len(300L), 301:600)
    testthat::local_mocked_bindings(
      lsiv_baseline_density_np_normal = function(Y, A, Z, ...) {
        direction <- if (all(Z == observed_Z)) {
          1L
        } else {
          0L
        }
        base <- if (direction == 1L) c(0.5, 0.5) else c(0.25, 0.75)
        list(density = matrix(base, length(Y), 2L, byrow = TRUE),
             method = "controlled baseline", candidates = "np")
      },
      lsiv_kliep_ulsif_ratio_grid = function(Y, A, Z, ...) {
        direction <- if (all(Z == observed_Z)) {
          1L
        } else {
          0L
        }
        ratio <- if (direction == 1L) c(0.4, 1.6) else c(1.6, 0.8)
        list(ratio = matrix(ratio, length(Y), 2L, byrow = TRUE),
             method = "controlled direct ratio", candidates = c("kliep", "ulsif"),
             ratio_method = "kliep_ulsif")
      },
      .package = "LSIV"
    )
    fit <- suppressWarnings(lsiv_estimate(
      Y, A, Z, X, type = "continuous", grid = c(-1, 1), weights = c(1, 1),
      y_ref = -1, folds = folds, learner = "glm", treatment_formula = A ~ Z * X1,
      outcome_density = "np_normal", ratio_method = "kliep_ulsif"
    ))
    small_direction <- as.character(as.integer(flip))
    direction <- fit$nuisance$ratio_direction
    expect_identical(direction$min_group, 100L)
    expect_equal(unname(direction$training_counts[, small_direction]), c(99, 100))
    small <- X[, 1L] < 0
    fallback <- small & seq_len(n) <= 300L
    expect_identical(direction$reciprocal_fallback, fallback)
    expect_identical(direction$fit_direction[fallback],
                     rep.int(1L - as.integer(flip), sum(fallback)))
    expect_identical(direction$fit_direction[!fallback], fit$zM[!fallback])
    expected <- if (flip) c(0.625, 1.25) else c(0.4, 1.6)
    expect_equal(fit$nuisance$RY[fallback, ],
                 matrix(expected, sum(fallback), 2L, byrow = TRUE))
    expected_direct <- if (flip) c(0.4, 1.6) else c(0.625, 1.25)
    direct <- small & seq_len(n) > 300L
    expect_equal(fit$nuisance$RY[direct, ],
                 matrix(expected_direct, sum(direct), 2L, byrow = TRUE))
  }
})
