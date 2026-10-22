test_that("expit and logit are inverses and numerically stable", {
  v <- c(-800, -50, -1, 0, 1, 50, 800)
  p <- expit(v)
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(expit(0), 0.5)
  pp <- c(1e-8, 0.1, 0.5, 0.9, 1 - 1e-8)
  expect_equal(expit(logit(pp)), pp, tolerance = 1e-10)
})

test_that("siv_basis_spec validates its arguments and prints", {
  expect_s3_class(siv_basis_spec("linear"), "siv_basis_spec")
  expect_error(siv_basis_spec("linear", y_ref = c(0, 1)), "single number")
  expect_error(siv_basis_spec("ns", df = 0), "positive integer")
  expect_output(print(siv_basis_spec("ns", df = 3, interact = "all")), "ns")
})

test_that("the selection-bias basis vanishes at the reference value", {
  ## This is the odds-ratio normalisation delta(y_ref, X) = 1 imposed by all
  ## four papers.
  set.seed(1)
  X <- cbind(A = rnorm(20), B = rnorm(20))
  for (tp in c("linear", "poly", "ns", "bs")) {
    for (yr in c(0, 0.7)) {
      sp <- SIV:::siv_basis_fit(siv_basis_spec(tp, df = 3, y_ref = yr,
                                               interact = "all"),
                                y_obs = seq(-2, 3, length.out = 50), X = X)
      B <- siv_basis_eval(sp, rep(yr, 20), X)
      expect_equal(max(abs(B)), 0, tolerance = 1e-10,
                   info = paste(tp, yr))
      expect_equal(ncol(B), SIV:::siv_basis_dim(sp), info = tp)
    }
  }
})

test_that("spline settings are configurable and change the basis dimension", {
  set.seed(2)
  X <- cbind(A = rnorm(10))
  y <- seq(-1, 2, length.out = 30)
  d3 <- SIV:::siv_basis_dim(SIV:::siv_basis_fit(siv_basis_spec("ns", df = 3), y, X))
  d5 <- SIV:::siv_basis_dim(SIV:::siv_basis_fit(siv_basis_spec("ns", df = 5), y, X))
  expect_gt(d5, d3)
  ## explicit knots / boundary knots are honoured
  sp <- SIV:::siv_basis_fit(
    siv_basis_spec("bs", degree = 2, knots = c(0, 1), Boundary.knots = c(-2, 3)), y, X)
  expect_equal(sp$knots, c(0, 1))
  expect_equal(sp$Boundary.knots, c(-2, 3))
  expect_equal(SIV:::siv_basis_dim(sp), 4L)   # length(knots) + degree
  ## interactions multiply the dimension
  spi <- SIV:::siv_basis_fit(siv_basis_spec("ns", df = 3, interact = "all"), y, X)
  expect_equal(SIV:::siv_basis_dim(spi), 2L * SIV:::siv_basis_dim(
    SIV:::siv_basis_fit(siv_basis_spec("ns", df = 3), y, X)))
})

test_that("Gauss-Hermite quadrature integrates normal moments exactly", {
  gh <- SIV:::siv_gh_normal(mu = c(0, 1.5), sd = 2, n_nodes = 20)
  m1 <- as.vector(gh$nodes %*% gh$weights)
  m2 <- as.vector((gh$nodes^2) %*% gh$weights)
  expect_equal(m1, c(0, 1.5), tolerance = 1e-10)
  expect_equal(m2, c(4, 1.5^2 + 4), tolerance = 1e-9)
  expect_equal(sum(gh$weights), 1, tolerance = 1e-12)
})
