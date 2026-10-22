## ---------------------------------------------------------------------------
## SIV_LATE(): the efficient influence function estimator of the LATE.
##
## The checks here are run on a CONTROLLED DGP whose nuisance functions and
## whose true LATE are both available in closed form, so that the algebra of
## the influence function is verified independently of any estimation step.
##
##   X ~ N(0,1);  pi(X) = expit(b X);  Z ~ Ber(pi(X))
##   U ~ Unif(0,1);  p_z(X) = expit(a_z + 0.4 X)  with a_1 > a_0
##   A = 1{U <= p_Z(X)}                 (latent index => A^(1) >= A^(0), i.e.
##                                       monotonicity holds by construction)
##   Y = 0.5 X + eps + A * tau(X, U),   tau(X, U) = 1 + X + 2U
##
## so that, in closed form,
##   m_A(z, X) = p_z(X)
##   m_Y(z, X) = 0.5 X + int_0^{p_z(X)} (1 + X + 2u) du
##             = 0.5 X + p_z(X)(1 + X) + p_z(X)^2
## and the complier average causal effect is E{tau(X,U) | p_0(X) < U <= p_1(X)}.
## ---------------------------------------------------------------------------

late_dgp <- function(n, b = 0.5) {
  X <- stats::rnorm(n)
  pi <- expit(b * X)
  Z <- stats::rbinom(n, 1, pi)
  U <- stats::runif(n)
  p0 <- expit(-0.5 + 0.4 * X)
  p1 <- expit(1.0 + 0.4 * X)
  A <- as.numeric(U <= ifelse(Z == 1, p1, p0))
  tau <- 1 + X + 2 * U
  Y <- 0.5 * X + stats::rnorm(n, 0, 0.5) + A * tau
  list(Y = Y, A = A, Z = Z, X = cbind(X = X), U = U, tau = tau, pi = pi,
       mA_z0 = p0, mA_z1 = p1,
       mY_z0 = 0.5 * X + p0 * (1 + X) + p0^2,
       mY_z1 = 0.5 * X + p1 * (1 + X) + p1^2)
}

## The reference value used by the tests below, computed once at high precision.
late_truth <- function(seed = 20240101, n = 2e6) {
  set.seed(seed)
  d <- late_dgp(n)
  c(wald = mean(d$mY_z1 - d$mY_z0) / mean(d$mA_z1 - d$mA_z0),
    complier = mean(d$tau[d$U > d$mA_z0 & d$U <= d$mA_z1]))
}

test_that("the Wald estimand equals the complier average causal effect", {
  skip_on_cran()
  tr <- late_truth()
  ## Under monotonicity these are the same number; agreement to Monte Carlo
  ## error confirms that the estimand implemented is the LATE.
  expect_equal(unname(tr["wald"]), unname(tr["complier"]), tolerance = 1e-3)
})

test_that("the LATE influence function is exactly mean zero", {
  set.seed(101)
  d <- late_dgp(2000)
  ef <- SIV:::siv_late_eif(d$Y, d$A, d$Z, d$pi, d$mY_z0, d$mY_z1,
                           d$mA_z0, d$mA_z1)
  ## By construction P_n[psi_Y - theta psi_A] = 0, so the EIF has mean 0.
  expect_lt(abs(mean(ef$eif)), 1e-10)
  expect_true(all(is.finite(ef$eif)))
  ## theta is exactly the ratio of the two uncentred score means
  expect_equal(ef$late, mean(ef$psi_Y) / mean(ef$psi_A))
})

test_that("the residual correction divides by Pr(Z = observed z | X)", {
  ## This is the one step of the algebra that is easy to get wrong: the
  ## denominator must be pi(X)^Z {1-pi(X)}^{1-Z}, not pi(X).  With a strongly
  ## varying propensity the wrong denominator is a visible bias.
  skip_on_cran()
  set.seed(102)
  d <- late_dgp(4e5, b = 1.8)
  truth <- mean(d$mY_z1 - d$mY_z0) / mean(d$mA_z1 - d$mA_z0)

  right <- SIV:::siv_late_eif(d$Y, d$A, d$Z, d$pi, d$mY_z0, d$mY_z1,
                              d$mA_z0, d$mA_z1)$late
  ## the same expression with pi(X) substituted into both arms
  wrong <- local({
    sgn <- 2 * d$Z - 1
    mY <- ifelse(d$Z == 1, d$mY_z1, d$mY_z0)
    mA <- ifelse(d$Z == 1, d$mA_z1, d$mA_z0)
    mean(sgn * (d$Y - mY) / d$pi + (d$mY_z1 - d$mY_z0)) /
      mean(sgn * (d$A - mA) / d$pi + (d$mA_z1 - d$mA_z0))
  })
  expect_lt(abs(right - truth), 0.02)
  expect_gt(abs(wrong - truth), 0.05)
  expect_gt(abs(wrong - truth), 5 * abs(right - truth))
})

test_that("SIV_LATE() recovers a known LATE and its SE is calibrated", {
  skip_on_cran()
  truth <- unname(late_truth()["wald"])
  set.seed(103)
  B <- 200; n <- 2000
  est <- se <- numeric(B)
  for (b in seq_len(B)) {
    d <- late_dgp(n)
    ## oracle nuisances: this isolates the estimator from nuisance estimation
    f <- SIV_LATE(d$Y, d$A, d$Z, d$X,
                  nuisance = list(pZ1 = d$pi, mA_z0 = d$mA_z0, mA_z1 = d$mA_z1,
                                  mY_z0 = d$mY_z0, mY_z1 = d$mY_z1))
    est[b] <- coef(f); se[b] <- siv_se(f)
  }
  ## unbiased
  expect_lt(abs(mean(est) - truth), 4 * sd(est) / sqrt(B))
  ## the influence-function SE matches the Monte Carlo spread
  expect_equal(mean(se) / sd(est), 1, tolerance = 0.1)
  ## and delivers nominal coverage
  cov95 <- mean(abs(est - truth) <= stats::qnorm(0.975) * se)
  expect_gt(cov95, 0.90)
  expect_lt(cov95, 0.99)
})

test_that("the LATE score is Neyman orthogonal", {
  skip_on_cran()
  set.seed(104)
  d <- late_dgp(5e5)
  truth <- mean(d$mY_z1 - d$mY_z0) / mean(d$mA_z1 - d$mA_z0)
  bias <- function(eps) {
    f <- SIV_LATE(d$Y, d$A, d$Z, d$X, nuisance = list(
      pZ1   = pmin(pmax(d$pi + eps * sin(3 * d$X), 0.02), 0.98),
      mA_z0 = pmin(pmax(d$mA_z0 + eps * sin(d$X), 0.01), 0.99),
      mA_z1 = pmin(pmax(d$mA_z1 + eps * sin(d$X + 1), 0.01), 0.99),
      mY_z0 = d$mY_z0 + eps * cos(2 * d$X),
      mY_z1 = d$mY_z1 + eps * cos(2 * d$X + 1)))
    abs(coef(f) - truth)
  }
  ## Halving the nuisance error should cut the bias by roughly four, not two.
  b1 <- bias(0.20); b2 <- bias(0.10); b3 <- bias(0.05)
  expect_gt(b1 / b2, 2.5)
  expect_gt(b2 / b3, 2.0)
})

test_that("SIV_LATE() runs end to end and returns a well-formed fit", {
  set.seed(105)
  d <- test_dgp(600, "continuous")
  f <- SIV_LATE(d$Y, d$A, d$Z, d$X)
  expect_s3_class(f, "SIV_LATE")
  expect_s3_class(f, "siv_fit")
  expect_named(coef(f), "LATE")
  expect_true(is.finite(coef(f)))
  expect_true(is.finite(siv_se(f)))
  expect_gt(siv_se(f), 0)
  expect_identical(dim(vcov(f)), c(1L, 1L))
  expect_true(all(siv_check_vcov(vcov(f), p = 1)))
  expect_equal(unname(confint(f)[1, ]),
               unname(coef(f) + c(-1, 1) * stats::qnorm(0.975) * siv_se(f)))
  expect_true(all(is.finite(f$eif)))
  expect_lt(abs(mean(f$eif)), 1e-8)
  ## cross-fitting and the no-covariate case both work
  expect_true(is.finite(coef(SIV_LATE(d$Y, d$A, d$Z, d$X, crossfit = 2))))
  expect_true(is.finite(coef(SIV_LATE(d$Y, d$A, d$Z, NULL))))
})

test_that("SIV_LATE() bounds the propensity score and reports the trimming", {
  set.seed(106)
  d <- late_dgp(400)
  f <- SIV_LATE(d$Y, d$A, d$Z, d$X, trim = 0.05)
  ## every estimated probability is a valid probability, and pi is bounded
  expect_true(all(f$nuisance$pZ1 >= 0.05 & f$nuisance$pZ1 <= 0.95))
  expect_true(all(f$nuisance$mA_z0 >= 0 & f$nuisance$mA_z0 <= 1))
  expect_true(all(f$nuisance$mA_z1 >= 0 & f$nuisance$mA_z1 <= 1))
  expect_true(all(is.finite(unlist(f$nuisance[c("mY_z0", "mY_z1")]))))
  expect_error(SIV_LATE(d$Y, d$A, d$Z, d$X, trim = 0.7), "trim")
  ## an aggressive bound must actually bind and be reported
  fw <- suppressWarnings(SIV_LATE(d$Y, d$A, d$Z, d$X, trim = 0.45))
  expect_gt(fw$diagnostics[["n_trimmed"]], 0)
})

test_that("SIV_LATE() validates its inputs and its nuisance argument", {
  set.seed(107)
  d <- test_dgp(300, "continuous")
  expect_error(SIV_LATE(d$Y, d$A, d$Z, d$X, nuisance = list(pZ1 = d$Z)),
               "missing element")
  expect_error(SIV_LATE(d$Y, d$A, d$Z, d$X,
                        nuisance = list(pZ1 = 0.5, mA_z0 = 0.2, mA_z1 = 0.6,
                                        mY_z0 = 0, mY_z1 = 1)),
               "must have length")
  ## LATE needs Y for treated units too, unlike the ATT estimators
  Yna <- d$Y; Yna[d$A == 1][1] <- NA
  expect_error(SIV_LATE(Yna, d$A, d$Z, d$X), "missing values")
})

test_that("a degenerate first stage is warned about, not silently returned", {
  set.seed(108)
  d <- late_dgp(300)
  ## Synthetic degeneracy that drives the influence-function denominator to
  ## EXACTLY zero: with m_A(0,X) = m_A(1,X) = A both the augmentation term
  ## Delta_A(X) and the residual A - m_A(Z,X) vanish, so psi_A == 0.  (Merely
  ## setting m_A(1,X) = m_A(0,X) is not enough: the inverse-probability
  ## residual then still contributes sampling noise, and tau_A is a small
  ## random number rather than zero.)
  expect_warning(
    SIV_LATE(d$Y, d$A, d$Z, d$X,
             nuisance = list(pZ1 = d$pi, mA_z0 = d$A, mA_z1 = d$A,
                             mY_z0 = d$mY_z0, mY_z1 = d$mY_z1)),
    "first stage")
})

test_that("an irrelevant instrument yields a large standard error", {
  ## A null instrument does not trip the degeneracy guard -- at n = 800 the
  ## estimated first stage is a small random number, not numerically zero --
  ## but it must be reflected in the reported uncertainty.
  set.seed(1081)
  n <- 800
  X <- stats::rnorm(n)
  Z <- stats::rbinom(n, 1, 0.5)
  A <- stats::rbinom(n, 1, expit(0.3 * X))   # A does not depend on Z
  Y <- stats::rnorm(n) + A
  f <- SIV_LATE(Y, A, Z, cbind(X = X))
  expect_true(is.finite(coef(f)))
  expect_true(is.finite(siv_se(f)))
  expect_gt(siv_se(f), 0.3)                  # cf. ~0.06 with a strong instrument
})

test_that("SIV_LATE() accepts the Super Learner nuisance estimator", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  set.seed(109)
  d <- test_dgp(400, "continuous")
  fs <- SIV_LATE(d$Y, d$A, d$Z, d$X, learner = "SL", SL.list = 1)
  fg <- SIV_LATE(d$Y, d$A, d$Z, d$X, learner = "glm")
  expect_identical(fs$learner, "SL")
  expect_true(is.finite(coef(fs)))
  expect_true(is.finite(siv_se(fs)))
  ## With a GLM-only library the ensemble is the GLM, so the two agree closely.
  expect_equal(unname(coef(fs)), unname(coef(fg)), tolerance = 1e-6)
})
