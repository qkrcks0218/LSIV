test_that("MySL XGBoost candidates fit their requested boosting rounds", {
  skip_if_not_installed("xgboost")
  set.seed(715)
  n <- 80L
  d <- data.frame(y = rbinom(n, 1, .5), x1 = rnorm(n), x2 = rnorm(n))
  fit <- MySL(d, 1L, 2:3, Ydist = binomial(), SL.list = 5L)
  expected <- expand.grid(SIV:::siv_mysl_library_tuning(n)$xgboost)$n.trees
  actual <- vapply(fit$fitLibrary, function(z) {
    length(unique(xgboost::xgb.model.dt.tree(model = z$object)$Tree))
  }, integer(1L))
  expect_length(actual, length(expected))
  expect_equal(unname(actual), expected)
  expect_true(all(is.finite(fit$SL.predict)))
})

test_that("MySL binomial neural-network candidates use their decay settings", {
  skip_if_not_installed("RSNNS")
  set.seed(716)
  n <- 80L
  d <- data.frame(y = rbinom(n, 1, .5), x1 = rnorm(n), x2 = rnorm(n))
  ## Zero-sized optional layers intentionally generate RSNNS warnings.
  fit <- suppressWarnings(MySL(d, 1L, 2:3, Ydist = binomial(), SL.list = 9L))
  expected <- expand.grid(SIV:::siv_mysl_library_tuning(n)$mlp)$decay
  actual <- vapply(fit$fitLibrary, function(z) {
    z$object$finalModel$learnFuncParams[1L]
  }, numeric(1L))
  expect_equal(unname(actual), expected)
  expect_true(all(is.finite(fit$SL.predict)))
})

test_that("two-fold splitting retains singleton and empty strata correctly", {
  for (sizes in list(c(5L, 5L, 10L, 1L), c(1L, 1L, 1L, 1L),
                    c(8L, 0L, 1L, 9L))) {
    A <- rep(c(0, 0, 1, 1), sizes)
    Z <- rep(c(0, 1, 0, 1), sizes)
    for (seed in 1:4) {
      set.seed(seed)
      folds <- SIV:::siv_make_folds(length(A), A, Z, 2L)
      expect_identical(sort(unlist(folds)), seq_along(A))
      expect_length(intersect(folds[[1L]], folds[[2L]]), 0L)
    }
  }
})

test_that("SIV rejects missing and infinite outcomes before fitting", {
  set.seed(717)
  for (type in c("binary", "continuous")) {
    d <- test_dgp(100L, type, p = 1L)
    for (a in 0:1) for (value in c(NA_real_, NaN, Inf, -Inf)) {
      y <- d$Y
      y[which(d$A == a)[1L]] <- value
      expect_error(SIV(y, d$A, d$Z, d$X, y_type = type,
        Overall.Estimation_num.crossfit = 1L,
        Overall.Estimation_num.fold = 2L), "`Y`.*(missing|finite|observed)")
    }
  }
})
