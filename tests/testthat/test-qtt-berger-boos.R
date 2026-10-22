qtt_binary_oracle <- function(seed = 2L, n = 3000L) {
  set.seed(seed)
  Z <- rbinom(n, 1, .5)
  A <- rbinom(n, 1, ifelse(Z, .8, .2))
  Y0 <- rbinom(n, 1, .51)
  Y1 <- rbinom(n, 1, .8)
  f <- matrix(rep(c(.49, .51), each = n), n, 2L)
  list(Y = ifelse(A, Y1, Y0), A = A, Z = Z,
       fp = list(delta = matrix(1, n, 2L), fY_z0 = f, fY_z1 = f),
       nuisance = list(pZ1 = rep(.5, n), pA1_z0 = rep(.2, n),
                       pA1_z1 = rep(.8, n), fY_z0 = f, fY_z1 = f))
}

qtt_oracle_fit <- function(d, ...) {
  siv_qtt(d$Y, d$A, d$Z, "binary", 0:1, c(1, 1),
           d$fp, d$nuisance, quantiles = .5, ...)
}

test_that("jump-aware binary inversion includes the true quantile difference", {
  d <- qtt_binary_oracle()
  fit <- siv_estimate(d$Y, d$A, d$Z, matrix(0, length(d$Y), 1L),
                      type = "binary", grid = 0:1, weights = c(1, 1),
                      nuisance = d$nuisance, estimand = "QTT", QTT.quantile = .5)
  # Both population medians are 1, so the true QTT is zero. The old equality
  # test returned [1, 1] on this sample despite the correctly specified nuisances.
  expect_equal(unname(fit$ci_matrix[1L, ]), c(0, 1))
  expect_true(all(fit$qtt$adjusted_pvalue[, 1L] > .05))
  expect_equal(fit$qtt$cdf_left, c(0, fit$qtt$cdf[1L]), tolerance = 1e-12)
  expect_equal(unname(fit$qtt$treated_quantile), 1)
})

test_that("QTT c1 tracks the nominal error and can be overridden", {
  d <- qtt_binary_oracle(n = 400L)
  for (level in c(.95, .90, .99)) {
    fit <- qtt_oracle_fit(d, level = level)
    expect_equal(fit$c1, .02 * (1 - level))
  }
  expect_equal(qtt_oracle_fit(d, level = .90, c1 = .005)$c1, .005)
  for (bad in list(0, .10, NA_real_, c(.001, .002), "0.001")) {
    expect_error(qtt_oracle_fit(d, level = .90, c1 = bad), "QTT_c1")
  }
  expect_error(qtt_oracle_fit(d, num_bootstrap = 1.5), "positive integer")

  set.seed(911)
  f <- suppressWarnings(SIV(
    d$Y, d$A, d$Z, estimand = "QTT", QTT_quantile = .5, level = .90,
    Overall.Estimation_num.crossfit = 1L, Overall.Estimation_num.fold = 2L
  ))
  expect_equal(f$qtt$c1, .002)
})

test_that("binomial rank intervals cover Bernoulli quantiles including jumps", {
  # Enumerate all possible sufficient statistics instead of relying on a Monte
  # Carlo coverage assertion or reproducing the rank-selection formula.
  for (n in c(5L, 25L, 100L)) {
    for (p in c(.2, .49, .5, .51, .8)) {
      for (q in c(.25, .5, .75)) {
        true_quantile <- as.numeric(q > 1 - p)
        intervals <- t(vapply(0:n, function(k) {
          as.numeric(siv_quantile_ci(c(rep(0, n - k), rep(1, k)), q, .001))
        }, numeric(2L)))
        coverage <- sum(dbinom(0:n, n, p) *
                          (intervals[, 1L] <= true_quantile &
                             intervals[, 2L] >= true_quantile))
        expect_gte(coverage, 1 - .001 - 1e-12)
      }
    }
  }
  expect_equal(unname(siv_quantile_ci(c(2, 4), .5, .001)[1L, ]),
               c(-Inf, Inf))
})

test_that("treated quantile inference uses no bootstrap randomness", {
  d <- qtt_binary_oracle(n = 400L)
  seed_before <- .Random.seed
  one <- qtt_oracle_fit(d, num_bootstrap = 1L)
  expect_identical(.Random.seed, seed_before)
  many <- qtt_oracle_fit(d, num_bootstrap = 10000L)
  expect_equal(one$ci, many$ci)
  expect_equal(one$treated_quantile_ci, many$treated_quantile_ci)
  expect_identical(one$num_bootstrap, 0L)
  expect_match(one$treated_quantile_ci_method, "order-statistic")
})

test_that("known discrete QTT support is independent of the density model", {
  d <- qtt_binary_oracle()
  f <- suppressWarnings(siv_estimate(
    d$Y, d$A, d$Z, matrix(0, length(d$Y), 1L), type = "continuous",
    grid = 0:1, weights = c(1, 1), nuisance = d$nuisance,
    estimand = "QTT", QTT.quantile = .5, QTT.support = 0:1
  ))
  expect_true(f$qtt$discrete)
  expect_equal(f$qtt$support, 0:1)
  expect_equal(unname(f$ci_matrix[1L, ]), c(0, 1))
  expect_error(qtt_oracle_fit(d, support = c(1, 0)), "strictly increasing")
  expect_error(qtt_oracle_fit(d, support = 0), "every untreated outcome")
})

test_that("continuous inversion includes accepted cells between integration nodes", {
  set.seed(2)
  n <- 600L
  Z <- rbinom(n, 1, .5)
  A <- rbinom(n, 1, ifelse(Z, .8, .2))
  Y0 <- runif(n, -.95, 1.05)
  Y1 <- rnorm(n, 1.05, .002)
  Y <- ifelse(A, Y1, Y0)
  grid <- c(-.95, seq(-.8, 1, by = .2), 1.05)
  mids <- (grid[-1L] + grid[-length(grid)]) / 2
  weights <- c(mids[1L] - grid[1L], diff(mids),
                tail(grid, 1L) - tail(mids, 1L))
  density <- matrix(.5, n, length(grid))
  fp <- list(delta = matrix(1, n, length(grid)),
             fY_z0 = density, fY_z1 = density)
  nuisance <- list(pZ1 = rep(.5, n), pA1_z0 = rep(.2, n), pA1_z1 = rep(.8, n))
  f <- siv_qtt(Y, A, Z, "continuous", grid, weights, fp, nuisance,
                quantiles = .5)
  # The true effect is 1 and its counterfactual median .05 lies off the grid.
  # Previously this sample returned [1.0496, 1.2505], missing an accepted value.
  expect_lte(f$ci[1L, "lower"], 1)
  expect_gte(f$ci[1L, "upper"], 1)
  accepted <- f$accepted_intervals[[1L]]
  expect_true(any(accepted[, "lower"] <= .05 & accepted[, "upper"] >= .05))
  expect_true(all(Y[A == 0] %in% f$threshold))
  expect_identical(f$threshold[1L], -Inf)
  expect_identical(tail(f$threshold_upper, 1L), Inf)

  # Direct evaluation away from all breakpoints must match its containing cell.
  t <- .05
  cf <- siv_cf_functionals(Y, A, Z, grid, weights, fp$delta,
                           density, density, nuisance$pZ1,
                           nuisance$pA1_z0, nuisance$pA1_z1,
                           outer(grid, t, "<=") * 1, outer(Y, t, "<=") * 1)
  index <- max(which(f$threshold <= t))
  direct_p <- pmin(1, siv_qtt_pvalue(cf$num, A, .5)$p + f$c1)
  expect_equal(unname(f$adjusted_pvalue[index, 1L]), unname(direct_p),
               tolerance = 1e-12)
  expect_equal(f$cdf[index], unname(cf$parameter), tolerance = 1e-12)

  bounded <- siv_qtt(Y, A, Z, "continuous", grid, weights, fp, nuisance,
                      quantiles = .5, bounds = c(-1, 2))
  expect_lte(bounded$ci[1L, "lower"], 1)
  expect_gte(bounded$ci[1L, "upper"], 1)
  expect_true(all(bounded$accepted_intervals[[1L]][, "lower"] >= -1))
  expect_true(all(bounded$accepted_intervals[[1L]][, "upper"] <= 2))
})

test_that("projection keeps complete cells, disconnected sets, and unbounded tails", {
  C1 <- c(lower = 3, upper = 4)
  expect_equal(siv_qtt_ci(2L, c(0, .5, 1), C1,
                          threshold_upper = c(.5, 1, Inf)),
               c(lower = 2, upper = 3.5))
  expect_equal(siv_qtt_ci(c(1L, 3L), c(-Inf, 0, 2), C1,
                          threshold_upper = c(0, 2, Inf)),
               c(lower = -Inf, upper = Inf))
  expect_equal(siv_qtt_ci(integer(), 0:1, C1, support = 0:1),
               c(lower = 2, upper = 4))
})

test_that("repeated QTT intervals respect known binary bounds", {
  expect_identical(siv_qtt_endpoint_median(c(-Inf, -Inf, 1)), -Inf)
  expect_identical(siv_qtt_endpoint_median(c(Inf, Inf, 1, NA)), Inf)
  set.seed(912)
  d <- test_dgp(400L, "binary")
  f <- suppressWarnings(SIV(
    d$Y, d$A, d$Z, d$X, estimand = "QTT", QTT_quantile = .001,
    Overall.Estimation_num.crossfit = 3L, Overall.Estimation_num.fold = 2L
  ))
  expect_true(all(is.finite(f$ci_matrix)))
  expect_gte(f$ci_matrix[1L, "lower"], -1)
  expect_lte(f$ci_matrix[1L, "upper"], 1)
  expect_false(anyNA(f$ci_matrix))
})

test_that("declared bounds constrain exterior ranks and empty inversions", {
  expect_equal(unname(siv_quantile_ci(c(2, 4), .5, .001,
                                       bounds = c(0, 52))[1L, ]), c(0, 52))
  expect_equal(siv_qtt_ci(integer(), 0:52, c(lower = -Inf, upper = Inf),
                          bounds = c(0, 52)), c(lower = -52, upper = 52))
  # The Job Corps treated-quantile interval is narrower than the full domain.
  expect_equal(siv_qtt_ci(integer(), 0:52, c(lower = 11, upper = 16),
                          bounds = c(0, 52)), c(lower = -41, upper = 16))
  expect_equal(siv_qtt_ci(c(1L, 3L), c(-Inf, 0, 52),
                          c(lower = 11, upper = 16),
                          threshold_upper = c(0, 52, Inf), bounds = c(0, 52)),
               c(lower = -41, upper = 16))
  # Counterfactual support alone imposes no bound on the treated outcome.
  expect_equal(siv_qtt_ci(integer(), 0:1, c(lower = 3, upper = 4),
                          support = 0:1), c(lower = 2, upper = 4))
})

test_that("bounded jump outcomes have finite CIs even after empty inversion", {
  set.seed(447)
  n <- 3000L
  Z <- rbinom(n, 1, .5)
  A <- rbinom(n, 1, ifelse(Z, .8, .2))
  Y0 <- sample(c(0, 52), n, replace = TRUE, prob = c(.8, .2))
  Y <- ifelse(A, 14, Y0)
  density <- matrix(rep(c(.8, .2), each = n), n, 2L)
  fp <- list(delta = matrix(1, n, 2L), fY_z0 = density, fY_z1 = density)
  nuisance <- list(pZ1 = rep(.5, n), pA1_z0 = rep(.2, n),
                   pA1_z1 = rep(.8, n), fY_z0 = density, fY_z1 = density)
  fit <- function(...) siv_qtt(Y, A, Z, "continuous", c(0, 52), c(1, 1),
                               fp, nuisance, quantiles = .25, ...)
  unknown <- fit()
  expect_true(unknown$empty_inversion[1L])
  expect_equal(unname(unknown$ci[1L, ]), c(-Inf, Inf))
  bounded <- fit(bounds = c(0, 52))
  expect_true(bounded$empty_inversion[1L])
  expect_equal(unname(bounded$ci[1L, ]), c(-38, 14))
  expect_equal(bounded$ci, bounded$fallback_ci)
  expect_true(all(bounded$threshold >= 0 & bounded$threshold_upper <= 52))
  discrete <- fit(support = 0:52, bounds = c(0, 52))
  expect_false(discrete$empty_inversion[1L])
  expect_equal(unname(discrete$ci[1L, ]), c(14, 14))
  core <- suppressWarnings(siv_estimate(
    Y, A, Z, matrix(0, n, 1L), type = "continuous", grid = c(0, 52),
    weights = c(1, 1), nuisance = nuisance, estimand = "QTT",
    QTT.quantile = .25, QTT.support = 0:52, QTT.bounds = c(0, 52)
  ))
  expect_equal(core$qtt$bounds, c(0, 52))
  expect_equal(unname(core$ci_matrix[1L, ]), c(14, 14))
})

test_that("QTT bounds are validated and passed through the public interface", {
  d <- qtt_binary_oracle(n = 400L)
  for (bad in list(c(1, 0), 0, c(0, 1, 2), c(NA, 1), c(0, Inf), "0,1")) {
    expect_error(qtt_oracle_fit(d, bounds = bad), "QTT_bounds")
  }
  expect_error(qtt_oracle_fit(d, bounds = c(0, .5)), "every observed outcome")
  expect_error(siv_qtt(d$Y, d$A, d$Z, "continuous", 0:1, c(1, 1),
                       d$fp, d$nuisance, support = 0:2, bounds = c(0, 1)),
                "QTT_support.*within.*QTT_bounds")
  set.seed(913)
  f <- suppressWarnings(SIV(
    d$Y * 52, d$A, d$Z, y_type = "continuous", y_grid = c(0, 52),
    y_weights = c(1, 1), y_ref = 0, outcome_density = "gaussian",
    ratio_method = "logistic", estimand = "QTT", QTT_quantile = .001,
    QTT_support = c(0, 52), QTT_bounds = c(0, 52),
    Overall.Estimation_num.crossfit = 1L, Overall.Estimation_num.fold = 2L
  ))
  expect_equal(f$qtt$bounds, c(0, 52))
  expect_true(all(is.finite(f$ci_matrix)))
  expect_gte(f$ci_matrix[1L, "lower"], -52)
  expect_lte(f$ci_matrix[1L, "upper"], 52)
})
