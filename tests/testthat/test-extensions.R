## Distributional estimands and the falsification test share the nuisance
## estimates and fixed point used by the original ATT implementation.

test_that("the generalized counterfactual EIF reproduces the ATT EIF algebra", {
  set.seed(81)
  d <- test_dgp(350, "continuous")
  g <- siv_grid(d$Y, d$A, "continuous", n_grid = 81,
                include_observed = FALSE)
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "continuous")
  fp <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1,
                        g$grid, g$weights, dr_floor = 0)
  ef <- siv_eif(d$Y, d$A, d$Z, g$grid, g$weights, fp$delta,
                fp$gamma_z0, fp$gamma_z1, fp$fY_z0, fp$fY_z1,
                nu$pZ1, nu$pA1_z0, nu$pA1_z1)
  cf <- siv_cf_functionals(
    d$Y, d$A, d$Z, g$grid, g$weights, fp$delta, fp$fY_z0, fp$fY_z1,
    nu$pZ1, nu$pA1_z0, nu$pA1_z1,
    G_grid = matrix(g$grid, ncol = 1L),
    G_obs = matrix(d$Y, ncol = 1L)
  )

  expect_equal(d$A * d$Y - cf$num[, 1L], ef$num, tolerance = 1e-10)
  expect_equal(unname(cf$parameter),
               mean(d$A * d$Y - ef$num) / mean(d$A), tolerance = 1e-10)
})

test_that("SIV returns selected ATT and QTT estimates with Berger-Boos CIs", {
  set.seed(82)
  d <- test_dgp(700, "continuous")
  g <- siv_grid(d$Y, d$A, "continuous", n_grid = 101,
                include_observed = FALSE)
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "continuous")
  f <- siv_estimate(
    d$Y, d$A, d$Z, d$X, type = "continuous", nuisance = nu,
    grid = g$grid, weights = g$weights,
    estimand = c("ATT", "QTT"), QTT.quantile = c(0.25, 0.5, 0.75),
    QTT.c1 = 0.005, num.QTT = 99
  )

  expect_named(f$estimate, c("ATT", "QTT_0.25", "QTT_0.5", "QTT_0.75"))
  expect_equal(dim(f$vcov), c(4L, 4L))
  expect_true(is.finite(f$vcov["ATT", "ATT"]))
  expect_true(all(is.na(diag(f$vcov)[-1L])))
  expect_true(all(is.finite(f$qtt$estimate)))
  expect_identical(f$qtt$threshold,
                   c(-Inf, sort(unique(c(g$grid, d$Y[d$A == 0])))))
  expect_equal(confint(f), f$ci_matrix, ignore_attr = TRUE)
  expect_true(all(f$qtt$treated_quantile_ci[, 1L] <=
                  f$qtt$treated_quantile_ci[, 2L]))
  expect_output(print(f), "Berger")

  qtt_only <- siv_estimate(
    d$Y, d$A, d$Z, d$X, type = "continuous", nuisance = nu,
    grid = g$grid, weights = g$weights, estimand = "QTT",
    QTT.quantile = 0.5, QTT.c1 = 0.005, num.QTT = 39
  )
  expect_named(qtt_only$estimate, "QTT_0.5")
  expect_false("ATT" %in% names(qtt_only$estimate))
})

test_that("empty QTT inversion uses specified support or an unbounded interval", {
  Y <- c(-2, 1, 5, NA_real_)
  threshold <- 0:5
  treated_ci <- c(lower = 3, upper = 4)

  empty_ci <- siv_qtt_ci(integer(0), threshold, treated_ci, Y)
  expect_equal(empty_ci, c(lower = -Inf, upper = Inf))
  expect_equal(siv_qtt_ci(integer(0), threshold, treated_ci, Y,
                          support = threshold),
               c(lower = -2, upper = 4))

  accepted_ci <- siv_qtt_ci(c(2L, 4L), threshold, treated_ci, Y)
  expect_equal(accepted_ci, c(lower = 0, upper = 3))
})

test_that("SIV median-combines repeated QTT estimates and Berger-Boos intervals", {
  set.seed(821)
  d <- test_dgp(500, "binary")
  f <- suppressWarnings(SIV(
    d$Y, d$A, d$Z, d$X, y_type = "binary",
    estimand = c("ATT", "QTT"), QTT_quantile = 0.5,
    QTT_c1 = 0.005, QTT_numBB = 19L,
    Overall.Estimation_num.crossfit = 3L,
    Overall.Estimation_num.fold = 2L
  ))
  ma <- f$median_adjustment
  qname <- "QTT_0.5"
  lower <- vapply(ma$confidence_intervals, function(ci) {
    unname(ci[qname, "lower"])
  }, numeric(1L))
  upper <- vapply(ma$confidence_intervals, function(ci) {
    unname(ci[qname, "upper"])
  }, numeric(1L))

  expect_equal(unname(f$estimate[qname]),
               stats::median(ma$estimates[, qname]))
  expect_equal(unname(f$ci_matrix[qname, "lower"]), stats::median(lower))
  expect_equal(unname(f$ci_matrix[qname, "upper"]), stats::median(upper))
  expect_equal(f$qtt$estimate, f$estimate[qname])
  expect_equal(f$qtt$ci, f$ci_matrix[qname, , drop = FALSE])
  expect_equal(f$qtt$empty_inversion_count,
               colSums(f$qtt$empty_inversion_by_repeat))
  expect_equal(f$qtt$empty_inversion,
               apply(f$qtt$empty_inversion_by_repeat, 2L, any))
  expect_equal(confint(f)[qname, ], f$ci_matrix[qname, ], ignore_attr = TRUE)
  expect_match(f$se_type, "median of repeated Berger-Boos intervals")
})

test_that("the ATT-only API remains backward compatible", {
  set.seed(83)
  d <- test_dgp(350, "binary")
  f <- SIV(d$Y, d$A, d$Z, d$X, y_type = "binary")
  expect_named(f$estimate, "ATT")
  expect_named(f$ci, c("lower", "upper"))
  expect_equal(dim(f$vcov), c(1L, 1L))
  expect_null(f$qtt)
  expect_null(f$falsification)
})

test_that("the falsification test builds the requested marginal dictionary", {
  set.seed(84)
  d <- test_dgp(500, "binary", p = 2L)
  g <- siv_grid(d$Y, d$A, "binary")
  nu <- test_oracle_nuisance(d$X, g$grid, g$weights, "binary")
  f <- siv_estimate(
    d$Y, d$A, d$Z, d$X, type = "binary", nuisance = nu,
    grid = g$grid, weights = g$weights, falsification = TRUE,
    num.false.boot = 59
  )
  ft <- f$falsification

  expect_true(is.list(ft))
  expect_equal(ft$num_bootstrap, 59L)
  expect_length(ft$bootstrap_statistics, 59L)
  expect_equal(ft$n_dictionary, ft$n_y_bins * ft$n_x_cells_observed)
  expect_equal(nrow(ft$moments), ft$n_dictionary)
  expect_equal(unname(ft$num_bins_requested), c(20L, 5L, 5L))
  expect_equal(ft$floor, 0.5^(0.5 + ncol(d$X) / 2))
  expect_true(ft$floor_default)
  expect_true(all(vapply(ft$x_bins, `[[`, character(1L), "kind") ==
                  "quantile"))
  expect_length(ft$categorical_x, 0L)
  expect_true(ft$p_value >= 0 && ft$p_value <= 1)
  expect_identical(ft$reject, ft$statistic > ft$critical_value)

  ## The default cells are marginal: one covariate at a time, so the number of
  ## cells adds over covariates instead of multiplying, and every row of the
  ## dictionary names the single covariate that defines it.
  expect_identical(ft$x_cells, "marginal")
  expect_equal(ft$n_x_cells_observed,
               sum(vapply(ft$x_bins, `[[`, integer(1L), "n_bins")))
  expect_setequal(unique(ft$moments$x_variable), colnames(d$X))
})

test_that("the falsification dictionary is always marginal", {
  expect_false("x.cells.falsification" %in% names(formals(SIV)))
})

test_that("bin width and categorical-variable status are separate choices", {
  x <- rep(1:3, each = 10)
  quantile_spec <- siv_make_empirical_bins(
    x, n_bins = 2, name = "X1", categorical = FALSE
  )
  categorical_spec <- siv_make_empirical_bins(
    x, n_bins = 2, name = "X1", categorical = TRUE
  )

  expect_equal(quantile_spec$kind, "quantile")
  expect_equal(quantile_spec$n_bins, 2L)
  expect_equal(categorical_spec$kind, "levels")
  expect_equal(categorical_spec$n_bins, 3L)
})

test_that("new SIV arguments are validated", {
  set.seed(85)
  d <- test_dgp(250, "binary")
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "binary",
                   estimand = "ATE"), "should be one of")
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "binary",
                   falsification = NA), "non-missing logical")
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "binary",
                   estimand = "QTT", QTT_quantile = 0), "strictly between")
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "binary",
                   falsification = TRUE,
                   falsification_num.bin = c(20, 5)),
               "ncol")
  expect_error(SIV(d$Y, d$A, d$Z, d$X, y_type = "binary",
                   falsification = TRUE, falsification_num.boot = 2,
                   falsification_categorical = "not_a_column"),
               "Unknown")
})

test_that("QTT and falsification warn when cross-fitting is disabled", {
  set.seed(86)
  d <- test_dgp(300, "binary")
  f <- suppressWarnings(SIV(
    d$Y, d$A, d$Z, d$X, y_type = "binary",
    falsification = TRUE, falsification_num.boot = 9,
    falsification_num.bin = c(2, 2, 2, 2)
  ))
  expect_equal(f$crossfit, 2L)
  expect_true(f$falsification$median_adjustment$applied)
  expect_identical(f$falsification$median_adjustment$num_crossfit, 3L)
  expect_length(f$falsification$median_adjustment$p_values, 3L)
  expect_equal(
    f$falsification$p_value,
    stats::median(f$falsification$median_adjustment$p_values)
  )
  expect_identical(
    f$falsification$reject,
    f$falsification$p_value <= 1 - f$level
  )
  expect_equal(
    f$median_adjustment$falsification_p_values,
    f$falsification$median_adjustment$p_values
  )
  expect_match(f$falsification$method, "median p-value")

  f_single <- suppressWarnings(SIV(
    d$Y, d$A, d$Z, d$X, y_type = "binary",
    Overall.Estimation_num.crossfit = 1L,
    Overall.Estimation_num.fold = 2L,
    falsification = TRUE, falsification_num.boot = 9,
    falsification_num.bin = c(2, 2, 2, 2)
  ))
  expect_false(f_single$falsification$median_adjustment$applied)
  expect_equal(
    f_single$falsification$p_value,
    f_single$falsification$median_adjustment$p_values[[1L]]
  )
  expect_identical(
    f_single$falsification$reject,
    f_single$falsification$statistic >
      f_single$falsification$critical_value
  )

  expect_warning(
    SIV(d$Y, d$A, d$Z, d$X, y_type = "binary", estimand = "QTT",
        QTT_quantile = 0.5, QTT_c1 = 0.005, QTT_numBB = 9,
        Overall.Estimation_num.crossfit = 1L,
        Overall.Estimation_num.fold = 1L),
    "require cross-fitted"
  )
})
