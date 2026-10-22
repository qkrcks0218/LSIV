## ---------------------------------------------------------------------------
## MySL(), ported from the UDID package, and the learner abstraction that lets
## SIV() and SIV_LATE() use it interchangeably with the default GLM regression
## learner.  SIV()'s continuous baseline density is handled separately.
##
## Every machine-learning package is a *suggested* dependency, so all tests
## that actually fit an ensemble are skipped when SuperLearner is absent.
## ---------------------------------------------------------------------------

test_that("siv_bound_prob() maps probabilities into [tol, 1 - tol]", {
  p <- c(-1, 0, 1e-9, 0.5, 1 - 1e-9, 1, 2, NA, NaN, Inf)
  b <- siv_bound_prob(p, tol = 0.01)
  expect_true(all(is.finite(b)))
  expect_true(all(b >= 0.01 & b <= 0.99))
  expect_equal(b[4], 0.5)
  ## non-finite input becomes 0.5 rather than propagating
  expect_equal(unname(b[8:10]), rep(0.5, 3))
})

test_that("MySL tuning matches the Version1 replication", {
  small <- SIV:::siv_mysl_tuning(100)
  large <- SIV:::siv_mysl_tuning(9240)
  expect_equal(small$MLPL, c(1, 2))
  expect_equal(small$MTRY, c(1, 2))
  expect_equal(small$NMN, 10)
  expect_equal(large$NMN, 25)
  expect_equal(small$MLPdecay, 10^c(-1, -3))

  grids <- SIV:::siv_mysl_library_tuning(9240)
  expect_identical(grids$glmnet,
                   list(alpha = c(1, 0.5, 0), useMin = c(TRUE, FALSE)))
  expect_identical(grids$earth, list(degree = 1:5))
  expect_identical(grids$gam, list(deg.gam = 1:5))
  expect_identical(grids$xgboost,
                   list(n.trees = c(100, 300, 500), max_depth = 1:4))
  expect_identical(grids$polymars, list(knots = 2:4))
  expect_identical(grids$ranger,
                   list(num.trees = c(500, 1000, 1500), mtry = c(1, 2)))
  expect_identical(grids$gbm,
                   list(ntree = c(100, 300, 500), intdepth = 1:3,
                        sh = c(0.1, 0.01), nmn = 25))
  expect_identical(grids$mlp,
                   list(L1 = c(1, 2), L2 = 0, L3 = 0,
                        decay = 10^c(-1, -3)))
  expect_identical(SIV:::siv_mysl_preprocess(),
                   c("center", "scale", "pca"))
})

test_that("siv_sl_available() drops unavailable groups and validates input", {
  expect_error(SIV:::siv_sl_available(0), "1:9")
  expect_error(SIV:::siv_sl_available(10), "1:9")
  ## group 1 needs SuperLearner only; it is either kept or the call errors
  if (requireNamespace("SuperLearner", quietly = TRUE)) {
    expect_identical(SIV:::siv_sl_available(1L), 1L)
  }
})

test_that("siv_make_folds() partitions the sample and stratifies on (A, Z)", {
  set.seed(301)
  n <- 400
  A <- stats::rbinom(n, 1, 0.4); Z <- stats::rbinom(n, 1, 0.5)
  expect_identical(SIV:::siv_make_folds(n, A, Z, 1), list(seq_len(n)))
  f <- SIV:::siv_make_folds(n, A, Z, 4)
  expect_length(f, 4L)
  expect_setequal(sort(unlist(f)), seq_len(n))
  ## every fold sees all four (A, Z) cells
  for (k in f) expect_equal(length(unique(interaction(A[k], Z[k]))), 4L)
})

test_that("two-fold splitting exactly matches Data_Version1", {
  set.seed(305)
  A <- rep(c(0, 0, 1, 1), c(11, 12, 13, 14))
  Z <- rep(c(0, 1, 0, 1), c(11, 12, 13, 14))
  cells <- list(
    which(A == 0 & Z == 0), which(A == 0 & Z == 1),
    which(A == 1 & Z == 0), which(A == 1 & Z == 1)
  )
  reference <- vector("list", 2L)
  for (j in seq_along(cells)) {
    permuted <- sample(cells[[j]], length(cells[[j]]))
    n_first <- if (j %in% c(1L, 3L)) {
      floor(length(permuted) / 2)
    } else {
      ceiling(length(permuted) / 2)
    }
    reference[[1L]] <- c(reference[[1L]], permuted[seq_len(n_first)])
    reference[[2L]] <- c(reference[[2L]],
                         setdiff(permuted, reference[[1L]][
                           reference[[1L]] %in% permuted
                         ]))
  }
  set.seed(305)
  expect_identical(SIV:::siv_make_folds(length(A), A, Z, 2L), reference)
})

test_that("the glm learner reproduces a direct glm/lm fit", {
  set.seed(302)
  n <- 300
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  df$b <- stats::rbinom(n, 1, expit(0.5 * df$x1))
  df$g <- 1 + df$x1 - 0.5 * df$x2 + stats::rnorm(n)

  fb <- SIV:::siv_learn_fit(df, "b", c("x1", "x2"), stats::binomial(), "glm")
  expect_equal(SIV:::siv_learn_predict(fb, df),
               unname(stats::predict(stats::glm(b ~ x1 + x2, df, family = stats::binomial()),
                                     type = "response")))
  fg <- SIV:::siv_learn_fit(df, "g", c("x1", "x2"), stats::gaussian(), "glm")
  expect_equal(SIV:::siv_learn_predict(fg, df),
               unname(stats::predict(stats::lm(g ~ x1 + x2, df))))
  expect_equal(SIV:::siv_learn_sigma(fg, df, "g"),
               stats::sigma(stats::lm(g ~ x1 + x2, df)))
})

test_that("a constant response degrades to a constant prediction", {
  df <- data.frame(y = rep(1, 50), x = stats::rnorm(50))
  f <- SIV:::siv_learn_fit(df, "y", "x", stats::binomial(), "glm")
  expect_s3_class(f, "siv_learner_const")
  expect_equal(SIV:::siv_learn_predict(f, df), rep(1, 50))
})

test_that("MySL() fits an ensemble and predicts valid probabilities", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  set.seed(303)
  n <- 300
  X <- matrix(stats::rnorm(n * 2), n, 2, dimnames = list(NULL, c("x1", "x2")))
  d <- data.frame(y = stats::rbinom(n, 1, expit(0.8 * X[, 1])), X)
  fit <- MySL(d, locY = 1, locX = 2:3, Ydist = stats::binomial(), SL.list = 1)
  expect_s3_class(fit, "SuperLearner")
  p <- stats::predict(fit, newdata = d[, 2:3], onlySL = TRUE)$pred
  expect_length(as.numeric(p), n)
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("MySL() errors informatively when no learner group is available", {
  skip_if_not_installed("SuperLearner")
  ## group 4 needs 'gam'; if it is installed this check is vacuous, so pick a
  ## group whose package is genuinely absent, and skip when all are present.
  missing_groups <- Filter(function(g) {
    pk <- SIV:::siv_sl_requires[[as.character(g)]]
    length(pk) && !all(vapply(pk, requireNamespace, logical(1), quietly = TRUE))
  }, 2:9)
  skip_if(length(missing_groups) == 0,
          "all suggested learner packages are installed")
  expect_error(suppressWarnings(SIV:::siv_sl_available(missing_groups[1])),
               "none of the requested learner groups")
})

test_that("SIV() uses the requested Super Learner library", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  set.seed(304)
  d <- test_dgp(400, "binary")
  fs <- suppressWarnings(SIV(d$Y, d$A, d$Z, d$X, y_type = "binary",
                             Overall.Estimation_SL.list = 1,
                             Overall.Estimation_num.crossfit = 1L,
                             Overall.Estimation_num.fold = 1L))
  fg <- suppressWarnings(siv_estimate(
    d$Y, d$A, d$Z, d$X, type = "binary", learner = "glm"
  ))
  expect_true(is.finite(coef(fs)))
  expect_true(is.finite(siv_se(fs)))
  ## a GLM-only Super Learner library IS the GLM, so the binary-outcome
  ## conditional pmf -- and hence the estimate -- agrees to numerical precision
  expect_equal(unname(coef(fs)), unname(coef(fg)), tolerance = 1e-6)
})

test_that("SIV() fixes Super Learner while SIV_LATE() retains its learner", {
  expect_false("learner" %in% names(formals(SIV)))
  expect_identical(formals(SIV)$Overall.Estimation_SL.list, 1L)
  expect_identical(eval(formals(SIV_LATE)$learner)[1], "glm")
})
