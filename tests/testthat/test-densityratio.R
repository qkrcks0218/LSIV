## ---------------------------------------------------------------------------
## The KLIEP density-ratio engine ported from the UDID package, and its use as
## the R_Y estimator of eq. (15) inside SIV().
## ---------------------------------------------------------------------------

test_that("KLIEP and KLIEP/uLSIF use the requested numerical floor", {
  expect_identical(formals(siv_kliep_ratio_grid)$dr_floor, 1e-4)
  expect_identical(formals(siv_kliep_ulsif_median_ratio_grid)$dr_floor, 1e-4)
  expect_identical(formals(siv_kliep_ulsif_ratio_grid)$dr_floor, 1e-4)
  expect_identical(
    formals(siv_kliep_ulsif_median_ratio_grid)$density_floor, 1e-4
  )
  expect_identical(formals(siv_kliep_ulsif_ratio_grid)$density_floor, 1e-4)
})

test_that("siv_kliep() recovers a known Gaussian density ratio", {
  ## p1 = N(mu, 1), p2 = N(0, 1) => log r(y) = mu*y - mu^2/2
  set.seed(401)
  mu <- 0.5
  x1 <- matrix(stats::rnorm(4000, mu)); x2 <- matrix(stats::rnorm(4000, 0))
  fit <- siv_kliep(x1, x2, n_y_cols = 1, sigma = "median")
  yy <- seq(-1.5, 1.5, by = 0.5)
  est <- siv_kliep_predict(matrix(yy), fit)
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
    f <- siv_kliep(x1, x2, n_y_cols = 1, sigma = sg)
    max(abs(log(siv_kliep_predict(matrix(yy), f)) - truth))
  }
  ## the cross-validated search is the hyper-parameter tuning; it must not be
  ## worse than the median heuristic on a problem where the heuristic is poor
  expect_lt(err("auto"), err("median"))
})

test_that("the bandwidth search is reproducible under set.seed()", {
  ## The CV fold labels are drawn with R's RNG, not arma::randperm(), precisely
  ## so that this holds.
  set.seed(404); a <- matrix(stats::rnorm(300, 0.5)); b <- matrix(stats::rnorm(300))
  set.seed(405); f1 <- siv_kliep(a, b, sigma = "auto")
  set.seed(405); f2 <- siv_kliep(a, b, sigma = "auto")
  expect_identical(f1$kernel_info$sigma, f2$kernel_info$sigma)
  expect_equal(f1$kernel_weights, f2$kernel_weights)
  ## and the median path too (kernel centres are sampled with R's RNG)
  set.seed(406); g1 <- siv_kliep(a, b, sigma = "median")
  set.seed(406); g2 <- siv_kliep(a, b, sigma = "median")
  expect_equal(g1$kernel_weights, g2$kernel_weights)
})

test_that("siv_kliep() validates its inputs", {
  a <- matrix(stats::rnorm(50)); b <- matrix(stats::rnorm(50))
  expect_error(siv_kliep(a, cbind(b, b)), "same number of columns")
  expect_error(siv_kliep(a[0, , drop = FALSE], b), "non-empty")
  expect_error(siv_kliep(a, b, n_y_cols = 2), "n_y_cols")
})

test_that("the per-dimension kernel helpers behave as documented", {
  set.seed(407)
  x <- matrix(stats::rnorm(20), 10, 2)
  ctr <- matrix(stats::rnorm(8), 4, 2)
  sg <- c(0.7, 1.3)
  K <- SIV:::siv_gaussian_kernel_matrix_bw(x, ctr, sg)
  expect_equal(dim(K), c(10L, 4L))
  ## reference computation
  ref <- outer(seq_len(10), seq_len(4), Vectorize(function(i, c)
    exp(-sum((x[i, ] - ctr[c, ])^2 / (2 * sg^2)))))
  expect_equal(K, ref, tolerance = 1e-12)
  ## the ratio is the kernel matrix times the weights
  w <- c(0.1, 0.4, 0.2, 0.3)
  expect_equal(as.numeric(SIV:::siv_gaussian_kernel_ratio_bw(x, ctr, sg, w)),
               as.numeric(K %*% w), tolerance = 1e-12)
  ## scalar-bandwidth path agrees with the per-dimension one at equal sigmas
  expect_equal(as.numeric(SIV:::siv_gaussian_kernel_ratio(x, ctr, 0.9, w)),
               as.numeric(SIV:::siv_gaussian_kernel_ratio_bw(x, ctr, c(0.9, 0.9), w)),
               tolerance = 1e-12)
})

test_that("SIV(ratio_method = 'kliep') runs and reports what it used", {
  skip_on_cran()
  set.seed(408)
  d <- test_dgp(1500, "continuous")
  f <- suppressWarnings(SIV(
    d$Y, d$A, d$Z, d$X, y_type = "continuous",
    y_num.grid = 61L, outcome_include.observed = FALSE,
    Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 2L, outcome_density = "gaussian",
    ratio_method = "kliep"
  ))
  expect_true(is.finite(coef(f)))
  expect_true(is.finite(siv_se(f)))
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
  gw <- siv_grid(d$Y, d$A, "continuous", n_grid = 41L,
                 include_observed = FALSE)
  folds <- siv_make_folds(nrow(d$X), d$A, d$Z, 2L)
  oracle <- test_oracle_nuisance(d$X, gw$grid, gw$weights, "continuous")
  fit <- siv_kliep_ulsif_ratio_grid(
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
  tr <- siv_kliep_ulsif_transform_info(x1, x0)
  z <- siv_kliep_ulsif_transform(
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
  tr <- siv_kliep_ulsif_common_transform_info(x1, x0)
  z <- siv_kliep_ulsif_common_transform(rbind(x1, x0), tr)
  expect_equal(unname(colMeans(z)), rep(0, 3L), tolerance = 1e-12)
  expect_equal(unname(apply(z, 2L, stats::sd)), rep(1, 3L),
               tolerance = 1e-12)

  cv_fit <- siv_kliep_ulsif_YXcommonbw_fit(
    x1, x0, kernel_num = 20L, fold = 2L, methods = "ulsif"
  )
  median_fit <- siv_kliep_ulsif_median_fit(
    x1, x0, kernel_num = 20L, methods = "ulsif"
  )
  expect_s3_class(cv_fit$fits$ulsif, "siv_kliep_ulsif_common")
  expect_s3_class(median_fit$fits$ulsif, "siv_kliep_ulsif_common")
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
  gw <- siv_grid(d$Y, d$A, "continuous", n_grid = 31L,
                 include_observed = FALSE)
  folds <- siv_make_folds(nrow(d$X), d$A, d$Z, 2L)
  oracle <- test_oracle_nuisance(d$X, gw$grid, gw$weights, "continuous")
  fit <- siv_kliep_ulsif_median_ratio_grid(
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

test_that("the simplified ratio interface is validated", {
  set.seed(409)
  d <- test_dgp(300, "continuous")
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "continuous",
                   ratio_method = "nonesuch"), "should be one of")
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "continuous",
                   ratio_bandwidth = "nonesuch"), "should be one of")
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "continuous",
                   ratio_ensemble.loss = "ise"), "unused argument")
  expect_error(SIV(
    d$Y, d$A, d$Z, d$X, y_type = "continuous",
    ratio_bandwidth.cv.grid = c(0, 0.5)
  ), "strictly between 0 and 1")
  expect_false("estimate_ratio" %in% names(formals(SIV)))
  expect_false("plugin" %in% eval(formals(SIV)$ratio_method))

  ## Binary outcomes always use the automatic plug-in ratio.
  db <- test_dgp(300, "binary")
  f <- SIV(db$Y, db$A, db$Z, db$X, y_type = "binary",
           ratio_method = "ulsif")
  expect_identical(f$ratio_method, "plugin (automatic)")
})

test_that("ratio ensemble minimizes MSPE without a loss-selection option", {
  H <- diag(2)
  h <- c(0.8, 0.2)
  fit <- SIV:::siv_ratio_ensemble_weights(H, h)

  expect_equal(fit$weight, c(0.8, 0.2), tolerance = 1e-12)
  expect_equal(fit$weight, fit$mspe)
  expect_equal(fit$mse, fit$mspe)
  expect_identical(fit$loss, "mspe")
  expect_error(SIV:::siv_ratio_ensemble_weights(H, h, loss = "logscore"),
               "unused argument")
})

test_that("standalone uLSIF uses the selected bandwidth rule", {
  skip_if_not_installed("densityratio")
  set.seed(411)
  d <- test_dgp(300, "continuous")
  f <- suppressWarnings(SIV(
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
