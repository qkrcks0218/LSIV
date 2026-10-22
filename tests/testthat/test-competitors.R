## ---------------------------------------------------------------------------
## Competing estimators: Liu et al. (2021), Sun et al. (2021),
## Sun/Miao/Wickramarachchi (2025).  These tests check (i) that the estimating
## functions have mean zero at the true parameter values, (ii) that the
## sandwich covariance is well formed, and (iii) that each estimator recovers
## the parameters of the DGP of its own paper.
## ---------------------------------------------------------------------------

test_that("Liu (2021) estimating functions are unbiased at the true parameters", {
  skip_on_cran()
  set.seed(41)
  d <- test_dgp_liu2021(60000)
  for (est in c("ipw", "or", "dr")) {
    su <- siv_liu_setup(d$Y, d$A, d$Z, d$X, estimator = est,
                        outcome_type = "binary",
                        pi_formula = ~ Z + C1 + Z:C1,
                        fy_formula = ~ C1 + C2 + Z + Z:C1)
    p <- su$par0; bl <- su$blocks
    p[bl$alpha] <- c(0.2, 0.4, -0.5)                    # logit P(Z=1|C)
    p[bl$mh2]   <- mean(d$X[, "C1"])
    if (est != "or") p[bl$pi] <- c(0.4, 2, 0.8, -1.6)   # (Int, Z, C1, Z:C1)
    p[bl$psi]   <- -0.6
    p[bl$theta] <- mean(d$Y0[d$A == 1])
    p[bl$mu1]   <- mean(d$Y[d$A == 1])
    E <- su$estfun(p)[, su$free, drop = FALSE]
    z <- colMeans(E) / (apply(E, 2, stats::sd) / sqrt(nrow(E)))
    expect_lt(max(abs(z)), 4.5)
  }
})

test_that("Liu (2021) estimators recover their own DGP with a valid sandwich", {
  skip_on_cran()
  set.seed(42)
  d <- test_dgp_liu2021(8000)
  for (est in c("ipw", "or", "dr")) {
    f <- SIV_Liu(d$Y, d$A, d$Z, d$X, estimator = est, outcome_type = "binary",
                 pi_formula = ~ Z + C1 + Z:C1, fy_formula = ~ C1 + C2 + Z + Z:C1)
    expect_true(f$converged)
    expect_true(all(f$checks))
    expect_s3_class(f, "SIV_Liu")
    ## psi within 4 sandwich SEs of the true value -0.6
    expect_lt(abs(coef(f)[["psi.y"]] + 0.6), 4 * siv_se(f)[["psi.y"]])
    ## ETT = mu1 - theta must hold exactly
    expect_equal(unname(coef(f)["ETT"]),
                 unname(coef(f)["mu1"] - coef(f)["theta"]), tolerance = 1e-10)
    expect_equal(f$se_type, "sandwich (stacked estimating equations)")
  }
})

test_that("Liu (2021) works with a continuous outcome and a spline basis", {
  skip_on_cran()
  set.seed(43)
  d <- test_dgp_liu2021_cont(3000)
  f <- SIV_Liu(d$Y, d$A, d$Z, d$X, estimator = "ipw", outcome_type = "continuous",
               basis = siv_basis_spec("linear"), pi_formula = ~ Z + C1 + Z:C1)
  expect_true(all(f$checks))
  expect_lt(abs(coef(f)[["psi.y"]] - 0.3), 4 * siv_se(f)[["psi.y"]])
  ## a spline basis is accepted and increases the number of psi parameters
  fs <- SIV_Liu(d$Y, d$A, d$Z, d$X, estimator = "ipw", outcome_type = "continuous",
                basis = siv_basis_spec("ns", df = 3), pi_formula = ~ Z + C1 + Z:C1)
  expect_equal(sum(grepl("^psi", names(coef(fs)))), 3L)
  expect_true(all(fs$checks))
})

test_that("Sun (2021) MNAR estimators recover their own DGP", {
  skip_on_cran()
  set.seed(44)
  d <- test_dgp_sun2021(20000)
  for (est in c("ipw", "or", "dr")) {
    f <- SIV_Sun2021(d$Yobs, d$R, d$Z, d$X, estimator = est,
                     outcome_type = "binary",
                     eta_formula = ~ Z + X1 + X2,
                     fz_formula  = ~ X1 + X2 + X1:X2,
                     fy_formula  = ~ X1 * X2 * Z,
                     h1_formula  = ~ Z + X1 + X2)
    expect_true(f$converged)
    expect_true(all(f$checks))
    ## true gamma = 1.8 and true mu = 0.769 (paper, Section 5)
    expect_lt(abs(coef(f)[["gamma.y"]] - 1.8), 4 * siv_se(f)[["gamma.y"]])
    expect_lt(abs(coef(f)[["mu"]] - 0.769), 4 * siv_se(f)[["mu"]])
    expect_true(grepl("sandwich", f$se_type))
  }
})

test_that("Sun (2021) handles missing outcome values correctly", {
  set.seed(45)
  d <- test_dgp_sun2021(2000)
  expect_true(anyNA(d$Yobs))
  f <- SIV_Sun2021(d$Yobs, d$R, d$Z, d$X, estimator = "ipw",
                   eta_formula = ~ Z + X1 + X2, fz_formula = ~ X1 + X2 + X1:X2)
  expect_true(is.finite(coef(f)[["mu"]]))
  ## must reject an outcome that is missing where R == 1
  bad <- d$Yobs; bad[which(d$R == 1)[1]] <- NA
  expect_error(SIV_Sun2021(bad, d$R, d$Z, d$X), "observed whenever")
})

test_that("Sun (2025) estimating functions are unbiased at the true parameters", {
  skip_on_cran()
  set.seed(46)
  d <- test_dgp_sun2025(120000, "binary")
  Yc <- ifelse(d$R == 1, d$Y, 0); Ac <- 1 - d$R
  for (aug in c(FALSE, TRUE)) {
    su <- siv_sun2025_setup(Yc, Ac, d$Z, d$X, outcome_type = "binary",
                            augmented = aug,
                            alpha_formula = ~ Z + U1 + U2,
                            fz_formula = ~ U1 + U2 + U1:U2,
                            fy_formula = ~ U1 + U2)
    p <- su$par0; bl <- su$blocks
    p[bl$alpha] <- c(2, -3, 0.8, 1)
    p[bl$zeta]  <- c(1, 2, -1, -0.8)
    p[bl$eta]   <- c(0.5, -2, 1)
    p[bl$gamma] <- 2
    p[bl$mu]    <- mean(d$Y)
    E <- su$estfun(p)[, seq_len(bl$mu), drop = FALSE]
    z <- colMeans(E) / (apply(E, 2, stats::sd) / sqrt(nrow(E)))
    expect_lt(max(abs(z)), 4.5)
  }
})

test_that("Sun (2025) is doubly robust and returns a valid sandwich", {
  skip_on_cran()
  set.seed(47)
  d <- test_dgp_sun2025(15000, "binary")
  Yc <- ifelse(d$R == 1, d$Y, 0); Ac <- 1 - d$R
  mu_true <- mean(d$Y)
  ## C1: all correct; C2: p(z|u) wrong; C3: p(y|u) wrong. mu survives all three.
  specs <- list(
    C1 = list(fz = ~ U1 + U2 + U1:U2, fy = ~ U1 + U2),
    C2 = list(fz = ~ U1 + I(U1^2),    fy = ~ U1 + U2),
    C3 = list(fz = ~ U1 + U2 + U1:U2, fy = ~ U1 + I(U1^2))
  )
  for (nm in names(specs)) {
    f <- SIV_Sun(Yc, Ac, d$Z, d$X, outcome_type = "binary",
                 alpha_formula = ~ Z + U1 + U2,
                 fz_formula = specs[[nm]]$fz, fy_formula = specs[[nm]]$fy)
    expect_true(f$converged)
    expect_true(all(f$checks))
    expect_lt(abs(coef(f)[["mu"]] - mu_true), 5 * siv_se(f)[["mu"]])
  }
})

test_that("Sun (2025) refuses an outcome model that involves the instrument", {
  set.seed(48)
  d <- test_dgp(300, "binary")
  expect_error(SIV_Sun(d$Y, d$A, d$Z, d$X, outcome_type = "binary",
                       fy_formula = ~ Z + X01),
               "must not involve")
})

test_that("Sun (2025) maps E(Y0) to the ATT exactly", {
  set.seed(49)
  d <- test_dgp(1500, "binary")
  f <- SIV_Sun(d$Y, d$A, d$Z, d$X, outcome_type = "binary")
  cf <- coef(f)
  expect_equal(unname(cf["theta"]),
               unname((cf["mu"] - cf["c0"]) / cf["pA"]), tolerance = 1e-10)
  expect_equal(unname(cf["ATT"]), unname(cf["mu1"] - cf["theta"]), tolerance = 1e-10)
  expect_true(all(siv_check_vcov(vcov(f), length(cf))))
  expect_gt(vcov(f)["ATT", "ATT"], 0)
})

test_that("Sun (2025) works with a continuous outcome (supplement S5.2 design)", {
  skip_on_cran()
  set.seed(50)
  d <- test_dgp_sun2025(6000, "continuous")
  Yc <- ifelse(d$R == 1, d$Y, 0); Ac <- 1 - d$R
  f <- SIV_Sun(Yc, Ac, d$Z, d$X, outcome_type = "continuous",
               alpha_formula = ~ Z + U1 + U2,
               fz_formula = ~ U1 + U2 + U1:U2,
               fy_formula = ~ U1 + U2)
  expect_true(f$converged)
  expect_true(all(f$checks))
  ## true tilting parameter gamma = 2 and the realised full-data mean
  expect_lt(abs(coef(f)[["gamma.y"]] - 2), 4 * siv_se(f)[["gamma.y"]])
  expect_lt(abs(coef(f)[["mu"]] - mean(d$Y)), 4 * siv_se(f)[["mu"]])
})
