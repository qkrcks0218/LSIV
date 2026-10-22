test_that("siv_grid uses the observed support for binary/discrete Y", {
  set.seed(21)
  A <- rbinom(40, 1, 0.5); A[1:5] <- 0
  Yb <- rbinom(40, 1, 0.5); Yb[A == 0][1:2] <- c(0, 1)
  g <- siv_grid(Yb, A, "binary")
  expect_equal(g$grid, c(0, 1))
  expect_equal(g$weights, c(1, 1))
  Yd <- sample(c(0, 1, 2, 5), 40, TRUE)
  gd <- siv_grid(Yd, A, "discrete")
  expect_equal(gd$grid, sort(unique(c(Yd[A == 0], 0))))
  expect_true(all(gd$weights == 1))
  ## a non-0/1 support must be rejected for type = "binary"
  expect_error(siv_grid(Yd, A, "binary"), "observed support")
})

test_that("siv_grid builds a fine grid with midpoint weights for continuous Y", {
  set.seed(22)
  A <- rbinom(200, 1, 0.5); Y <- rnorm(200)
  g <- siv_grid(Y, A, "continuous", n_grid = 51, include_observed = FALSE)
  expect_equal(length(g$grid), length(g$weights))
  expect_true(0 %in% g$grid)                     # reference value is on the grid
  expect_equal(sum(g$weights), diff(range(g$grid)), tolerance = 1e-8)
  expect_true(all(g$weights > 0))
})

test_that("Theorem 3.1 recovers the true delta and gamma from the true theta (binary)", {
  set.seed(23)
  d <- test_dgp(3000, "binary")
  g <- siv_grid(d$Y, d$A, "binary")
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "binary")
  fp <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1,
                        g$grid, g$weights, solver = "closed-form")
  expect_equal(fp$delta, nu$delta, tolerance = 1e-8)
  expect_equal(fp$gamma_z0, nu$gamma_z0, tolerance = 1e-8)
  expect_equal(fp$gamma_z1, nu$gamma_z1, tolerance = 1e-8)
  ## the normalisation delta(y_ref) = 1 holds exactly
  expect_equal(unname(fp$delta[, 1]), rep(1, nrow(fp$delta)), tolerance = 1e-12)
})

test_that("the closed-form and contraction solvers agree (binary)", {
  set.seed(24)
  d <- test_dgp(500, "binary")
  g <- siv_grid(d$Y, d$A, "binary")
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "binary")
  a <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1, g$grid, g$weights,
                       solver = "closed-form")
  b <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1, g$grid, g$weights,
                       solver = "contraction", maxit = 2000)
  expect_equal(a$delta, b$delta, tolerance = 1e-6)
  expect_true(all(b$converged))
})

test_that("the contraction converges from arbitrary starting values (Theorem 3.2)", {
  set.seed(25)
  d <- test_dgp(300, "continuous")
  g <- siv_grid(d$Y, d$A, "continuous", n_grid = 81, include_observed = FALSE)
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "continuous")
  M <- length(g$grid); n <- nrow(d$X)
  g1 <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1, g$grid, g$weights,
                        maxit = 3000)
  set.seed(99)
  init <- matrix(runif(n * M, 0.2, 3), n, M)
  g2 <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1, g$grid, g$weights,
                        maxit = 3000, g_init = init)
  expect_equal(g1$delta, g2$delta, tolerance = 1e-6)
})

test_that("Theorem 3.1 recovers the true delta on a fine grid (continuous)", {
  set.seed(26)
  d <- test_dgp(400, "continuous")
  g <- siv_grid(d$Y, d$A, "continuous", n_grid = 401, include_observed = FALSE)
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "continuous")
  fp <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1, g$grid, g$weights,
                        maxit = 3000)
  expect_lt(max(abs(log(fp$delta) - log(nu$delta))), 1e-4)
  expect_lt(max(abs(fp$gamma_z1 / nu$gamma_z1 - 1)), 1e-6)
})

test_that("the fixed point handles either ordering of z_m and z_M", {
  ## Flipping the labels of Z must leave delta unchanged.
  set.seed(27)
  d <- test_dgp(300, "binary")
  g <- siv_grid(d$Y, d$A, "binary")
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "binary")
  a <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1, g$grid, g$weights,
                       solver = "contraction")
  b <- siv_fixed_point(nu$pA1_z1, nu$pA1_z0, nu$fY_z1, nu$fY_z0, g$grid, g$weights,
                       solver = "contraction")
  expect_equal(a$delta, b$delta, tolerance = 1e-6)
  expect_true(all(a$zM == 1)); expect_true(all(b$zM == 0))
})

test_that("siv_fixed_point validates its inputs", {
  expect_error(
    siv_fixed_point(rep(.2, 3), rep(.6, 3), matrix(.5, 3, 2), matrix(.5, 3, 2),
                    c(0, 1), c(1, 1), y_ref = 7),
    "grid points")
  expect_error(
    siv_fixed_point(rep(.2, 3), rep(.6, 3), matrix(.5, 3, 3), matrix(.5, 3, 3),
                    c(0, 1, 2), c(1, 1, 1), solver = "closed-form"),
    "two-point")
})

## ---------------------------------------------------------------------------
## The influence function must follow the manuscript's replication code
## (Simulation_F_Both/1.Sim_Conti_SepIV.R, "Calibration: beta" onwards).
##
## In that code  gamma.z0.CF = Pr(A=1|z0,X)/Pr(A=0|z0,X)  is the RAW odds and
##               beta.z0.CF  = gamma.z0.CF / denom.z0.
## siv_fixed_point() returns gamma already divided by denom once (eq. 18), so
## siv_eif() must NOT divide by denom a second time.  Dividing twice left the
## weight alpha*beta too small by a factor of denom (range 1.04-2.31) and broke
## the identity checked below; it changed every reported standard error.
## ---------------------------------------------------------------------------

test_that("siv_eif() uses beta = odds / E{alpha | A=0,z,X}", {
  set.seed(501)
  d <- test_dgp(600, "continuous")
  g <- siv_grid(d$Y, d$A, "continuous", n_grid = 61L, include_observed = FALSE)
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "continuous")
  fp <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1,
                        g$grid, g$weights, y_ref = 0, solver = "contraction",
                        dr_floor = 0)
  ef <- siv_eif(d$Y, d$A, d$Z, g$grid, g$weights, fp$delta,
                fp$gamma_z0, fp$gamma_z1, fp$fY_z0, fp$fY_z1,
                nu$pZ1, nu$pA1_z0, nu$pA1_z1)
  n <- length(d$Y); M <- length(g$grid)
  W <- matrix(g$weights, n, M, byrow = TRUE)
  denom_z0 <- rowSums(fp$delta * nu$fY_z0 * W)
  denom_z1 <- rowSums(fp$delta * nu$fY_z1 * W)
  expect_equal(ef$beta_z0, (nu$pA1_z0 / (1 - nu$pA1_z0)) / denom_z0)
  expect_equal(ef$beta_z1, (nu$pA1_z1 / (1 - nu$pA1_z1)) / denom_z1)
  ## `denom` is genuinely non-constant and away from 1 for some observations,
  ## so the two conventions really differ.  Do not use a mean threshold here:
  ## its value depends strongly on the manuscript's non-monotone shape a(y).
  expect_gt(max(abs(denom_z0 - 1)), 0.05)
})

test_that("the reconstructed f(Y0 | X) agrees across the two arms of Z", {
  ## This is the replication code's own sanity check
  ## `mean((fy_from_z0.CF - fy_from_z1.CF)^2)  # should be small`.
  ## It is ~1e-25 with the correct beta and ~9e-04 with beta/denom.
  set.seed(502)
  for (s in 1:2) {
    d <- test_dgp(500, "continuous")
    g <- siv_grid(d$Y, d$A, "continuous", n_grid = 61L, include_observed = FALSE)
    nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "continuous")
    fp <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1,
                          g$grid, g$weights, y_ref = 0,
                          solver = "contraction", dr_floor = 0)
    ef <- siv_eif(d$Y, d$A, d$Z, g$grid, g$weights, fp$delta,
                  fp$gamma_z0, fp$gamma_z1, fp$fY_z0, fp$fY_z1,
                  nu$pZ1, nu$pA1_z0, nu$pA1_z1)
    w0 <- fp$delta * ef$beta_z0; w1 <- fp$delta * ef$beta_z1
    f00 <- nu$fY_z0 * ((1 - nu$pA1_z0) * (1 - nu$pZ1))
    f01 <- nu$fY_z1 * ((1 - nu$pA1_z1) * nu$pZ1)
    resid <- mean(((f00 + w0 * f00) / (1 - nu$pZ1) -
                   (f01 + w1 * f01) / nu$pZ1)^2)
    expect_lt(resid, 1e-18)
  }
})
