test_that("ratio projection preserves normalization and the IV bound together", {
  nu <- list(pA1_z0 = .25, pA1_z1 = .5,
             fY_z0 = matrix(c(.5, .5), 1L),
             fY_z1 = matrix(c(.05, .95), 1L))
  raw <- nu
  fit <- SIV:::siv_constrain_density_ratio(nu, c(1, 1))
  adjusted <- fit$nuisance

  ## Clipping followed by ordinary renormalization violates the upper bound.
  clipped <- pmin(nu$fY_z1 / nu$fY_z0, 1.5)
  renormalized <- clipped / sum(clipped * nu$fY_z0)
  expect_gt(max(renormalized), 1.5)
  expect_equal(adjusted$RY, matrix(c(.5, 1.5), 1L), tolerance = 1e-10)
  expect_equal(adjusted$fY_z1, matrix(c(.25, .75), 1L), tolerance = 1e-10)
  expect_equal(rowSums(adjusted$fY_z1), 1, tolerance = 1e-12)
  expect_lte(max(adjusted$RY), 1.5)
  expect_identical(nu, raw)
  expect_equal(fit$diagnostics$projected_units, 1)
  expect_lte(fit$diagnostics$max_normalization_error, 1e-12)
  expect_gte(fit$diagnostics$min_falsification_slack, -1e-12)
})

test_that("ratio projection handles both instrument orderings and direct ratios", {
  nu <- list(pA1_z0 = c(.25, .5), pA1_z1 = c(.5, .25),
             fY_z0 = rbind(c(.5, .5), c(.05, .95)),
             fY_z1 = rbind(c(.05, .95), c(.5, .5)))
  fit <- SIV:::siv_constrain_density_ratio(nu, c(1, 1))
  out <- fit$nuisance
  ordered <- out$RY
  ordered[2L, ] <- 1 / ordered[2L, ]
  expect_equal(ordered, matrix(c(.5, 1.5), 2L, 2L, byrow = TRUE),
               tolerance = 1e-10)
  expect_equal(rowSums(out$fY_z0), c(1, 1), tolerance = 1e-12)
  expect_equal(rowSums(out$fY_z1), c(1, 1), tolerance = 1e-12)
  expect_equal(out$RY, out$fY_z1 / out$fY_z0)

  ## The supplied targeted ratio, not an independently fitted numerator,
  ## determines the reconstructed density. Its orientation is always Z=1/Z=0.
  targeted <- nu
  targeted$RY <- nu$fY_z1 / nu$fY_z0
  targeted$fY_z0[2L, ] <- c(.5, .5)
  targeted$fY_z1[1L, ] <- c(.5, .5)
  direct <- SIV:::siv_constrain_density_ratio(targeted, c(1, 1))
  expect_equal(direct$nuisance$fY_z0, out$fY_z0, tolerance = 1e-10)
  expect_equal(direct$nuisance$fY_z1, out$fY_z1, tolerance = 1e-10)
})

test_that("continuous-grid weights and numerical ratio bounds are respected", {
  weights <- c(.2, .5, 1.3)
  nu <- list(pA1_z0 = c(.1, .7), pA1_z1 = c(.7, .1),
             fY_z0 = rbind(c(.4, .3, .2), c(.6, .05, .7)),
             fY_z1 = rbind(c(.6, .05, .7), c(.4, .3, .2)),
             RY = rbind(c(.01, 100, .5), c(100, .01, 2)))
  fit <- SIV:::siv_constrain_density_ratio(nu, weights, dr_floor = .25)
  out <- fit$nuisance
  expect_equal(as.vector(out$fY_z0 %*% weights), c(1, 1), tolerance = 1e-12)
  expect_equal(as.vector(out$fY_z1 %*% weights), c(1, 1), tolerance = 1e-12)
  expect_gte(min(out$RY), .25 - 1e-12)
  expect_lte(max(out$RY), 4 + 1e-12)
  d <- out$pA1_z0 - out$pA1_z1
  implication <- ((1 - out$pA1_z1) * out$fY_z1 -
                    (1 - out$pA1_z0) * out$fY_z0) / d
  expect_gte(min(implication), -1e-12)
})

test_that("compatible oracle inputs are unchanged and degenerate inputs fail clearly", {
  set.seed(880)
  d <- test_dgp(50, "binary")
  nu <- test_oracle_nuisance(d$X, c(0, 1), c(1, 1), "binary")
  fit <- SIV:::siv_constrain_density_ratio(nu, c(1, 1),
                                          density_floor = 0, dr_floor = 0)
  expect_equal(fit$nuisance$fY_z0, nu$fY_z0, tolerance = 1e-12)
  expect_equal(fit$nuisance$fY_z1, nu$fY_z1, tolerance = 1e-12)
  expect_equal(fit$diagnostics$projected_units, 0)
  expect_error(SIV:::siv_constrain_density_ratio(nu, c(1, 1), dr_floor = 1),
               "dr_floor")
  bad <- nu
  bad$pA1_z1 <- bad$pA1_z0
  expect_error(SIV:::siv_constrain_density_ratio(bad, c(1, 1)),
               "nonzero conditional first-stage")
  bad <- nu
  bad$fY_z0[1, 1] <- NA_real_
  expect_error(SIV:::siv_constrain_density_ratio(bad, c(1, 1)),
               "finite and nonnegative")
})

test_that("large guarded ratios retain weighted normalization to roundoff", {
  set.seed(883)
  n <- 200L
  M <- 61L
  weights <- exp(runif(M, -3, 2))
  nu <- list(
    pA1_z0 = runif(n, .01, .99), pA1_z1 = runif(n, .01, .99),
    fY_z0 = matrix(exp(rnorm(n * M, sd = 4)), n),
    fY_z1 = matrix(exp(rnorm(n * M, sd = 4)), n),
    RY = matrix(exp(rnorm(n * M, sd = 10)), n)
  )
  fit <- SIV:::siv_constrain_density_ratio(nu, weights)
  out <- fit$nuisance
  expect_equal(as.vector(out$fY_z0 %*% weights), rep(1, n), tolerance = 1e-12)
  expect_equal(as.vector(out$fY_z1 %*% weights), rep(1, n), tolerance = 1e-12)
  expect_gte(fit$diagnostics$min_falsification_slack, -1e-12)
  expect_lte(fit$diagnostics$max_normalization_error, 1e-12)
})

test_that("SIV estimates are constrained but falsification remains unrestricted", {
  set.seed(881)
  n <- 400L
  Z <- rep(c(0, 1), length.out = n)
  A <- rbinom(n, 1, ifelse(Z == 0, .2, .6))
  Y <- rbinom(n, 1, ifelse(Z == 0, .1, .8))
  X <- matrix(numeric(0), n, 0L)
  nu <- list(pZ1 = rep(.5, n), pA1_z0 = rep(.2, n),
             pA1_z1 = rep(.6, n),
             fY_z0 = matrix(c(.9, .1), n, 2L, byrow = TRUE),
             fY_z1 = matrix(c(.2, .8), n, 2L, byrow = TRUE))
  original <- nu
  folds <- list(seq_len(n))
  set.seed(882)
  expected_test <- SIV:::siv_falsification(
    Y, A, Z, X, c(0, 1), c(1, 1), nu,
    num_bootstrap = 39L, num_bins = 2L
  )
  set.seed(882)
  fit <- siv_estimate(Y, A, Z, X = NULL, type = "binary", nuisance = nu,
                      grid = c(0, 1), weights = c(1, 1), folds = folds,
                      falsification = TRUE, num.false.boot = 39L,
                      num.false.bin = 2L)
  expect_identical(fit$falsification, expected_test)
  expect_true(fit$falsification$reject)
  expect_identical(fit$nuisance, original)
  expect_identical(nu, original)
  expect_gt(fit$ratio_constraint$projected_fraction, 0)
  expect_gte(fit$ratio_constraint$min_falsification_slack, -1e-12)
  expect_equal(fit$fixed_point$fY_z0, fit$estimation_nuisance$fY_z0)
  expect_equal(fit$fixed_point$fY_z1, fit$estimation_nuisance$fY_z1)
  expect_equal(fit$fixed_point$RY,
               fit$estimation_nuisance$fY_z1 / fit$estimation_nuisance$fY_z0)
  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(is.finite(fit$eif)))

  without_test <- siv_estimate(
    Y, A, Z, X = NULL, type = "binary", nuisance = nu,
    grid = c(0, 1), weights = c(1, 1), folds = folds, falsification = FALSE
  )
  expect_null(without_test$falsification)
  expect_identical(fit$estimate, without_test$estimate)
  expect_identical(fit$estimation_nuisance, without_test$estimation_nuisance)
})

test_that("the public SIV interface constrains every repeated cross-fit", {
  set.seed(884)
  d <- test_dgp(300, "binary", p = 1L)
  fit <- suppressWarnings(SIV(
    d$Y, d$A, d$Z, d$X,
    Overall.Estimation_num.crossfit = 2L,
    Overall.Estimation_num.fold = 2L,
    falsification = TRUE, falsification_num.boot = 9L,
    falsification_num.bin = c(2L, 2L)
  ))
  expect_length(fit$median_adjustment$ratio_constraints, 2L)
  for (constraint in fit$median_adjustment$ratio_constraints) {
    expect_true(constraint$enabled)
    expect_gte(constraint$min_falsification_slack, -1e-12)
    expect_lte(constraint$max_normalization_error, 1e-12)
  }
  expect_length(fit$falsification$median_adjustment$p_values, 2L)
})
