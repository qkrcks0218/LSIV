test_that("siv_sandwich reproduces the closed-form sandwich for a GLM", {
  ## For a logistic regression the stacked-EE sandwich must coincide with the
  ## usual Huber-White estimator.
  set.seed(11)
  n <- 400
  x <- rnorm(n)
  y <- rbinom(n, 1, expit(0.4 + 0.8 * x))
  D <- cbind(1, x)
  fit <- glm(y ~ x, family = binomial())
  ef <- function(th) D * as.vector(y - expit(D %*% th))
  sw <- siv_sandwich(ef, coef(fit))
  ## closed form: (X'WX)^{-1} X' diag(r^2) X (X'WX)^{-1}
  p <- fitted(fit); r <- y - p
  bread <- solve(t(D) %*% (D * p * (1 - p)))
  meat <- t(D) %*% (D * r^2)
  V <- bread %*% meat %*% bread
  expect_equal(unname(sw$vcov), unname(V), tolerance = 1e-6)
  expect_true(all(sw$checks))
})

test_that("siv_sandwich scaling is 1/n and diagnostics behave", {
  set.seed(12)
  make <- function(n) {
    y <- rnorm(n, 2, 3)
    list(ef = function(th) cbind(y - th[1], (y - th[1])^2 - th[2]),
         par = c(mean(y), var(y) * (n - 1) / n))
  }
  s1 <- make(1000); s2 <- make(4000)
  v1 <- siv_sandwich(s1$ef, s1$par)$vcov[1, 1]
  v2 <- siv_sandwich(s2$ef, s2$par)$vcov[1, 1]
  ## variance of a sample mean: exactly var(y)/n
  expect_equal(v1 * 1000 / 9, 1, tolerance = 0.2)
  expect_lt(v2, v1)
  chk <- siv_check_vcov(diag(3), p = 3)
  expect_true(all(chk))
  expect_false(siv_check_vcov(matrix(c(1, 2, 0, 1), 2), p = 2)["symmetric"])
  expect_false(siv_check_vcov(diag(c(1, -1)), p = 2)["diag_nonneg"])
  expect_false(siv_check_vcov(diag(2), p = 3)["dim_ok"])
})

test_that("siv_sandwich errors when estfun and theta are non-conformable", {
  expect_error(siv_sandwich(function(th) matrix(0, 10, 3), c(1, 2)),
               "returned 3 columns")
})

test_that("the sandwich accounts for nuisance estimation (bread/meat/vcov)", {
  skip_on_cran()
  ## Fit the Liu (2021) IPW estimator, then recompute the variance of `theta`
  ## while PRETENDING the instrument model f_{Z|C} were known.  The two must
  ## differ, which demonstrates that the nuisance score equations really do
  ## enter the reported sandwich covariance.
  set.seed(61)
  d <- test_dgp_liu2021(3000)
  f <- SIV_Liu(d$Y, d$A, d$Z, d$X, estimator = "ipw", outcome_type = "binary",
               pi_formula = ~ Z + C1 + Z:C1)
  nm <- names(coef(f))
  ith <- match("theta", nm)

  ## the reported (correct) variance
  v_full <- vcov(f)["theta", "theta"]

  ## drop the alpha block from BOTH bread and meat -> treats f_{Z|C} as known
  keep <- !grepl("^alpha", colnames(f$estfun))
  A2 <- f$jacobian[keep, keep, drop = FALSE]
  B2 <- f$meat[keep, keep, drop = FALSE]
  Ai <- solve(A2)
  V2 <- Ai %*% B2 %*% t(Ai) / f$n
  j <- match("theta", colnames(A2))
  v_naive <- V2[j, j]

  expect_true(is.finite(v_full) && is.finite(v_naive))
  expect_gt(v_full, 0)
  expect_gt(v_naive, 0)
  expect_false(isTRUE(all.equal(v_full, v_naive, tolerance = 1e-3)))

  ## bread and meat have the documented shapes and the sandwich reassembles
  p <- length(coef(f)) - 1L                      # ETT is a derived contrast
  expect_equal(dim(f$jacobian), c(p, p))
  expect_equal(dim(f$meat), c(p, p))
  expect_equal(dim(f$estfun), c(f$n, p))
  ## meat = n^{-1} sum_i M_i M_i'
  expect_equal(f$meat, crossprod(f$estfun) / f$n, tolerance = 1e-10,
               ignore_attr = TRUE)
  ## vcov = A^{-1} B A^{-T} / n on the stacked block
  V <- solve(f$jacobian) %*% f$meat %*% t(solve(f$jacobian)) / f$n
  expect_equal(unname(vcov(f)[seq_len(p), seq_len(p)]), unname((V + t(V)) / 2),
               tolerance = 1e-8)
  ## estimating functions are (numerically) centred at the solution
  expect_lt(max(abs(colMeans(f$estfun))), 1e-3)
})
