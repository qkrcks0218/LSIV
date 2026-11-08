test_that("expit and logit are inverses and numerically stable", {
  v <- c(-800, -50, -1, 0, 1, 50, 800)
  p <- expit(v)
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(expit(0), 0.5)
  pp <- c(1e-8, 0.1, 0.5, 0.9, 1 - 1e-8)
  expect_equal(expit(logit(pp)), pp, tolerance = 1e-10)
})

test_that("lsiv_basis_spec validates its arguments and prints", {
  expect_s3_class(lsiv_basis_spec("linear"), "lsiv_basis_spec")
  expect_error(lsiv_basis_spec("linear", y_ref = c(0, 1)), "single number")
  expect_error(lsiv_basis_spec("ns", df = 0), "positive integer")
  expect_output(print(lsiv_basis_spec("ns", df = 3, interact = "all")), "ns")
})

test_that("the selection-bias basis vanishes at the reference value", {
  ## This is the odds-ratio normalisation delta(y_ref, X) = 1 imposed by all
  ## four papers.
  set.seed(1)
  X <- cbind(A = rnorm(20), B = rnorm(20))
  for (tp in c("linear", "poly", "ns", "bs")) {
    for (yr in c(0, 0.7)) {
      sp <- LSIV:::lsiv_basis_fit(lsiv_basis_spec(tp, df = 3, y_ref = yr,
                                                  interact = "all"),
                                  y_obs = seq(-2, 3, length.out = 50), X = X)
      B <- lsiv_basis_eval(sp, rep(yr, 20), X)
      expect_equal(max(abs(B)), 0, tolerance = 1e-10,
                   info = paste(tp, yr))
      expect_equal(ncol(B), LSIV:::lsiv_basis_dim(sp), info = tp)
    }
  }
})

test_that("spline settings are configurable and change the basis dimension", {
  set.seed(2)
  X <- cbind(A = rnorm(10))
  y <- seq(-1, 2, length.out = 30)
  d3 <- LSIV:::lsiv_basis_dim(LSIV:::lsiv_basis_fit(lsiv_basis_spec("ns", df = 3), y, X))
  d5 <- LSIV:::lsiv_basis_dim(LSIV:::lsiv_basis_fit(lsiv_basis_spec("ns", df = 5), y, X))
  expect_gt(d5, d3)
  ## explicit knots / boundary knots are honoured
  sp <- LSIV:::lsiv_basis_fit(
    lsiv_basis_spec("bs", degree = 2, knots = c(0, 1), Boundary.knots = c(-2, 3)), y, X)
  expect_equal(sp$knots, c(0, 1))
  expect_equal(sp$Boundary.knots, c(-2, 3))
  expect_equal(LSIV:::lsiv_basis_dim(sp), 4L)   # length(knots) + degree
  ## interactions multiply the dimension
  spi <- LSIV:::lsiv_basis_fit(lsiv_basis_spec("ns", df = 3, interact = "all"), y, X)
  expect_equal(LSIV:::lsiv_basis_dim(spi), 2L * LSIV:::lsiv_basis_dim(
    LSIV:::lsiv_basis_fit(lsiv_basis_spec("ns", df = 3), y, X)))
})

test_that("Gauss-Hermite quadrature integrates normal moments exactly", {
  gh <- LSIV:::lsiv_gh_normal(mu = c(0, 1.5), sd = 2, n_nodes = 20)
  m1 <- as.vector(gh$nodes %*% gh$weights)
  m2 <- as.vector((gh$nodes^2) %*% gh$weights)
  expect_equal(m1, c(0, 1.5), tolerance = 1e-10)
  expect_equal(m2, c(4, 1.5^2 + 4), tolerance = 1e-9)
  expect_equal(sum(gh$weights), 1, tolerance = 1e-12)
})
