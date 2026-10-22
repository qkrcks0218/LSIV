## ---------------------------------------------------------------------------
## The compiled fast paths in src/siv_fastpath.cpp must reproduce the R
## reference implementations in R/siv_core.R exactly (up to floating-point
## summation order).  The R versions remain authoritative: they are what runs
## when the package sources are sourced rather than installed, and they are the
## definition the C++ is checked against here.
##
## `options(SIV.use_cpp = FALSE)` forces the R path; the cached probe in
## siv_cpp_ok() has to be reset when the option is flipped.
## ---------------------------------------------------------------------------

with_cpp <- function(use, code) {
  old <- getOption("SIV.use_cpp", TRUE)
  on.exit({
    options(SIV.use_cpp = old)
    assign("cache", NULL, envir = environment(SIV:::siv_cpp_ok))
  }, add = TRUE)
  options(SIV.use_cpp = use)
  assign("cache", NULL, envir = environment(SIV:::siv_cpp_ok))
  force(code)
}

test_that("the compiled fast path is available and can be switched off", {
  expect_true(with_cpp(TRUE, SIV:::siv_cpp_ok()))
  expect_false(with_cpp(FALSE, SIV:::siv_cpp_ok()))
})

test_that("siv_delta_moments_cpp() matches its R reference", {
  set.seed(201)
  for (dim in list(c(50, 3), c(200, 40), c(37, 501))) {
    n <- dim[1]; M <- dim[2]
    delta <- matrix(stats::runif(n * M, 0.05, 20), n, M)
    fY <- matrix(stats::runif(n * M, 1e-6, 3), n, M)
    grid <- sort(stats::rnorm(M, 0, 2))
    w <- stats::runif(M, 0.01, 0.5)
    a <- SIV:::siv_delta_moments_cpp(delta, grid, fY, w)
    b <- SIV:::siv_delta_moments_r(delta, grid, fY, w)
    expect_equal(dim(a), c(n, 4L))
    expect_identical(colnames(a), c("denom", "y", "alpha", "alphay"))
    expect_equal(a, b, tolerance = 1e-10)
  }
})

test_that("siv_contraction_cpp() matches its R reference step for step", {
  set.seed(202)
  for (dim in list(c(60, 2), c(120, 25), c(80, 200))) {
    n <- dim[1]; M <- dim[2]
    w <- rep(1 / M, M)
    fzm <- matrix(stats::runif(n * M, 0.05, 1), n, M)
    fzM <- matrix(stats::runif(n * M, 0.05, 1), n, M)
    fzm <- fzm / as.vector(fzm %*% w)
    fzM <- fzM / as.vector(fzM %*% w)
    RY <- fzM / fzm
    p_zm <- stats::runif(n, 0.05, 0.4)
    p_zM <- stats::runif(n, 0.5, 0.9)
    g0 <- matrix(1 / sum(w), n, M)
    g0 <- g0 / as.vector(g0 %*% w)
    a <- SIV:::siv_contraction_cpp(g0, p_zm, p_zM, fzm, fzM, RY, w, 500L, 1e-10)
    b <- SIV:::siv_contraction_r(g0, p_zm, p_zM, fzm, fzM, RY, w, 500L, 1e-10)
    ## identical iteration count and identical convergence flags
    expect_identical(a$iter, b$iter)
    expect_identical(as.logical(a$converged), as.logical(b$converged))
    expect_equal(a$g, b$g, tolerance = 1e-10)
    ## the C++ must not modify its input in place
    expect_equal(g0, matrix(1 / sum(w), n, M) / as.vector(matrix(1 / sum(w), n, M) %*% w))
  }
})

test_that("siv_fixed_point() and siv_eif() agree with and without C++", {
  skip_on_cran()
  for (type in c("binary", "continuous")) {
    for (seed in 1:3) {
      set.seed(seed)
      d <- test_dgp(500, type)
      g <- siv_grid(d$Y, d$A, type, n_grid = 121L, include_observed = FALSE)
      nu <- test_oracle_nuisance(d$X, g$grid, g$weights, type)
      run <- function(use) with_cpp(use, {
        fp <- siv_fixed_point(nu$pA1_z0, nu$pA1_z1, nu$fY_z0, nu$fY_z1,
                              g$grid, g$weights, solver = "contraction",
                              dr_floor = 1e-2)
        ef <- siv_eif(d$Y, d$A, d$Z, g$grid, g$weights, fp$delta,
                      fp$gamma_z0, fp$gamma_z1, fp$fY_z0, fp$fY_z1,
                      nu$pZ1, nu$pA1_z0, nu$pA1_z1)
        list(fp = fp, ef = ef)
      })
      a <- run(TRUE); b <- run(FALSE)
      expect_identical(a$fp$iter, b$fp$iter)
      expect_equal(a$fp$g, b$fp$g, tolerance = 1e-10)
      expect_equal(a$fp$delta, b$fp$delta, tolerance = 1e-10)
      expect_equal(a$fp$gamma_z0, b$fp$gamma_z0, tolerance = 1e-10)
      expect_equal(a$ef$att, b$ef$att, tolerance = 1e-10)
      expect_equal(a$ef$var, b$ef$var, tolerance = 1e-10)
    }
  }
})

test_that("SIV() gives the same answer with and without the compiled path", {
  skip_on_cran()
  for (type in c("binary", "continuous")) {
    set.seed(203)
    d <- test_dgp(500, type)
    set.seed(204)
    f1 <- with_cpp(TRUE,  suppressWarnings(
      SIV(d$Y, d$A, d$Z, d$X, y_type = type, y_num.grid = 121L,
          Overall.Estimation_num.crossfit = 1L,
          Overall.Estimation_num.fold = 1L,
          outcome_include.observed = FALSE, outcome_density = "gaussian",
          ratio_method = "logistic")))
    set.seed(204)
    f2 <- with_cpp(FALSE, suppressWarnings(
      SIV(d$Y, d$A, d$Z, d$X, y_type = type, y_num.grid = 121L,
          Overall.Estimation_num.crossfit = 1L,
          Overall.Estimation_num.fold = 1L,
          outcome_include.observed = FALSE, outcome_density = "gaussian",
          ratio_method = "logistic")))
    expect_equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-10)
    expect_equal(unname(siv_se(f1)), unname(siv_se(f2)), tolerance = 1e-10)
  }
})

test_that("the compiled contraction rejects non-conformable input", {
  m <- matrix(1, 4L, 3L)
  expect_error(SIV:::siv_contraction_cpp(m, rep(1, 4), rep(1, 4), m, m,
                                         matrix(1, 4L, 2L), rep(1, 3), 10L, 1e-8),
               "non-conformable")
  expect_error(SIV:::siv_contraction_cpp(m, rep(1, 4), rep(1, 4), m, m, m,
                                         rep(1, 2), 10L, 1e-8),
               "weights")
  expect_error(SIV:::siv_delta_moments_cpp(m, rep(0, 3), matrix(1, 4L, 2L), rep(1, 3)),
               "non-conformable")
})
