## ----------------------------------------------------------------------------
## Nonparametric density-ratio estimation for R_Y (eq. 15), ported from the
## UDID package.
##
## WHY THIS EXISTS
## ---------------
## The fixed point of Theorem 3.1 needs
##     R_Y(y, x) = f(y | A=0, Z=1, x) / f(y | A=0, Z=0, x)          (eq. 15)
## and Section 5.1 of the manuscript recommends estimating that RATIO directly
## rather than dividing two separately estimated conditional densities.  The
## manuscript's own replication code does exactly that: for a continuous
## outcome `Simulation_Version1/1.Sim_Conti_SepIV.R` fits `densityratio::kliep()`
## (and uLSIF) with a cross-validated Gaussian-kernel bandwidth, cross-fitted,
## and ensembles the candidates.
##
## The compiled KLIEP engine comes from `UDID_Package/src/udid_kliep.cpp` and
## remains available through `siv_kliep()`. Public `SIV()` ratio methods use the
## KLIEP/uLSIF implementation in `densityratio_ensemble.R`.
##
## RELATION TO THE UDID CODE
## -------------------------
##   UDID                              SIV
##   ------------------------------    ------------------------------------
##   KLIEP_bw()                        siv_kliep()
##   .median_heuristic_bw()            siv_kliep_median_bw()
##   compute_density_ratio_Kernel()    siv_kliep_predict()
##   kliep_bw_cpp()                    siv_kliep_bw_cpp()      (src/siv_kliep.cpp)
##   gaussian_kernel_*()               siv_gaussian_kernel_*()
##
## Two deliberate changes:
##   1. the CV fold labels are drawn in R, not by `arma::randperm()`, so that
##      the bandwidth search obeys `set.seed()`;
##   2. the median-heuristic path uses the compiled `siv_kliep_optimize_cpp()`
##      instead of `densratio:::KLIEP_optimize_alpha`; this compiled alternative
##      does not itself depend on `densityratio`.
## ----------------------------------------------------------------------------

## Internal: per-dimension median heuristic (Garreau et al., 2017).
## sigma_Y from the pairwise distances among the Y columns, sigma_X from the X
## columns.  This is UDID's `.median_heuristic_bw()`.
siv_kliep_median_bw <- function(x1, n_y_cols) {
  p <- ncol(x1)
  idx_Y <- seq_len(n_y_cols)
  idx_X <- if (n_y_cols < p) (n_y_cols + 1L):p else integer(0)
  sig <- function(idx) {
    if (!length(idx)) return(1)
    v <- stats::median(stats::dist(x1[, idx, drop = FALSE]))
    if (!is.finite(v) || v < 1e-6) 1 else as.numeric(v)
  }
  sv <- numeric(p)
  sv[idx_Y] <- sig(idx_Y)
  if (length(idx_X)) sv[idx_X] <- sig(idx_X)
  sv
}

#' KLIEP density-ratio estimator with per-dimension bandwidths
#'
#' Estimates \eqn{r(w) = p_1(w)/p_2(w)} from two samples by Kullback-Leibler
#' Importance Estimation (KLIEP), with separate Gaussian-kernel bandwidths for
#' the outcome columns and the covariate columns.
#'
#' This is the engine of the \pkg{UDID} package (`KLIEP_bw()`), ported so that
#' [SIV()] can estimate the density ratio \eqn{R_Y} of eq. (15) *directly*, the
#' way the manuscript's replication code does, instead of dividing two
#' separately estimated conditional densities.
#'
#' @section Hyper-parameter tuning:
#' `sigma = "auto"` runs the compiled coordinate-descent search: it alternates
#' between optimising \eqn{\sigma_Y} and \eqn{\sigma_X}, at each step maximising
#' the `fold`-fold cross-validated KLIEP log-likelihood
#' \eqn{n^{-1}\sum_i \log\{\phi(w_i)^\top\alpha\}} over a decreasing ladder of
#' step sizes (10, 1, 0.1, ..., 1e-5).  `sigma = "median"` skips the search and
#' uses the median heuristic, optimising only the kernel weights -- much faster,
#' and the default inside [SIV()] for that reason.
#' The fold labels are drawn with R's RNG, so the search is reproducible under
#' [set.seed()].
#'
#' @param x1 Numerator sample, an `n1 x p` matrix.
#' @param x2 Denominator sample, an `n2 x p` matrix with the same columns.
#' @param n_y_cols Number of leading columns treated as "outcome" columns, which
#'   get their own bandwidth.  Defaults to 1.
#' @param sigma `"median"` (median heuristic; weights only) or `"auto"`
#'   (cross-validated bandwidth search).  A numeric vector of length `p` fixes
#'   the bandwidths directly.
#' @param kernel_num Number of Gaussian kernel centres, sampled from `x1`.
#' @param fold Number of cross-validation folds used when `sigma = "auto"`.
#' @param n_cd_iter Coordinate-descent sweeps when `sigma = "auto"`.
#' @return A list with `kernel_weights` and `kernel_info` (`kernel`,
#'   `kernel_num`, `sigma`, `centers`), in the format expected by
#'   [siv_kliep_predict()].
#' @references
#' Sugiyama, M., Suzuki, T., Nakajima, S., Kashima, H., von Bunau, P. and
#' Kawanabe, M. (2008). Direct importance estimation for covariate shift
#' adaptation. *Annals of the Institute of Statistical Mathematics* 60, 699-746.
#'
#' Garreau, D., Jitkrittum, W. and Kanagawa, M. (2017). Large sample analysis of
#' the median heuristic. arXiv:1707.07269.
#' @seealso [siv_kliep_predict()], and the `ratio_method` argument of [SIV()].
#' @examples
#' set.seed(1)
#' x1 <- cbind(rnorm(200, 0.5), rnorm(200))
#' x2 <- cbind(rnorm(200, 0.0), rnorm(200))
#' fit <- siv_kliep(x1, x2, n_y_cols = 1)
#' head(siv_kliep_predict(x1, fit))
#' @export
siv_kliep <- function(x1, x2, n_y_cols = 1L, sigma = c("median", "auto"),
                      kernel_num = 100L, fold = 5L, n_cd_iter = 3L) {
  if (is.vector(x1)) x1 <- matrix(x1)
  if (is.vector(x2)) x2 <- matrix(x2)
  x1 <- as.matrix(x1); x2 <- as.matrix(x2)
  storage.mode(x1) <- "double"; storage.mode(x2) <- "double"
  if (ncol(x1) != ncol(x2)) {
    stop("`x1` and `x2` must have the same number of columns.", call. = FALSE)
  }
  if (!nrow(x1) || !nrow(x2)) {
    stop("`x1` and `x2` must both be non-empty.", call. = FALSE)
  }
  n_y_cols <- as.integer(n_y_cols)
  if (n_y_cols < 1L || n_y_cols > ncol(x1)) {
    stop("`n_y_cols` must be between 1 and ncol(x1).", call. = FALSE)
  }
  numeric_sigma <- is.numeric(sigma)
  if (!numeric_sigma) sigma <- match.arg(sigma)

  nx1 <- nrow(x1)
  kernel_num <- min(as.integer(kernel_num), nx1)
  centers <- x1[sample.int(nx1, kernel_num), , drop = FALSE]

  sigma_init <- siv_kliep_median_bw(x1, n_y_cols)

  if (numeric_sigma) {
    sigma_vec <- rep_len(as.numeric(sigma), ncol(x1))
    phi1 <- siv_gaussian_kernel_matrix_bw(x1, centers, sigma_vec)
    phi2 <- siv_gaussian_kernel_matrix_bw(x2, centers, sigma_vec)
    kw <- as.numeric(siv_kliep_optimize_cpp(phi1, phi2))
  } else if (identical(sigma, "median")) {
    ## median heuristic: no bandwidth search, optimise the weights only
    sigma_vec <- sigma_init
    phi1 <- siv_gaussian_kernel_matrix_bw(x1, centers, sigma_vec)
    phi2 <- siv_gaussian_kernel_matrix_bw(x2, centers, sigma_vec)
    kw <- as.numeric(siv_kliep_optimize_cpp(phi1, phi2))
  } else {
    ## "auto": cross-validated coordinate-descent bandwidth search.
    ## Fold labels come from R's RNG so that set.seed() controls the search.
    fold <- max(2L, as.integer(fold))
    cv_split <- integer(nx1)
    cv_split[sample.int(nx1)] <- (seq_len(nx1) - 1L) %% fold + 1L
    res <- siv_kliep_bw_cpp(x1, x2, centers, sigma_init * 10, n_y_cols,
                            cv_split, fold, as.integer(n_cd_iter))
    sigma_vec <- as.numeric(res$sigma)
    kw <- as.numeric(res$kernel_weights)
    centers <- as.matrix(res$centers)
  }

  list(kernel_weights = kw,
       kernel_info = list(kernel = "Gaussian", kernel_num = kernel_num,
                          sigma = sigma_vec, centers = centers),
       fold = fold)
}

#' Evaluate a fitted KLIEP density ratio
#'
#' @param x Matrix of evaluation points, with the same columns as the samples
#'   used to fit `object`.
#' @param object A fit returned by [siv_kliep()].
#' @return Numeric vector of estimated density ratios, one per row of `x`.
#' @seealso [siv_kliep()]
#' @examples
#' set.seed(1)
#' fit <- siv_kliep(matrix(rnorm(200, 0.5)), matrix(rnorm(200)), n_y_cols = 1)
#' head(siv_kliep_predict(matrix(seq(-2, 2, length.out = 5)), fit))
#' @export
siv_kliep_predict <- function(x, object) {
  if (is.vector(x)) x <- matrix(x)
  x <- as.matrix(x); storage.mode(x) <- "double"
  sig <- object$kernel_info$sigma
  ctr <- as.matrix(object$kernel_info$centers)
  kw <- as.numeric(object$kernel_weights)
  out <- if (length(sig) > 1L) {
    siv_gaussian_kernel_ratio_bw(x, ctr, as.numeric(sig), kw)
  } else {
    siv_gaussian_kernel_ratio(x, ctr, as.numeric(sig), kw)
  }
  as.numeric(out)
}

## Internal: cross-fitted KLIEP estimate of R_Y on the outcome grid.
##
## Numerator sample  : (Y, X) among  A == 0, Z == 1
## Denominator sample: (Y, X) among  A == 0, Z == 0
## so the fitted ratio evaluated at (y, x_i) is exactly
##     R_Y(y, x_i) = f(y | A=0, Z=1, x_i) / f(y | A=0, Z=0, x_i).
## This mirrors `Ry_x.kliep.estimator` in 1.Sim_Conti_SepIV.R.
##
## Returns an n x M matrix, or NULL if the ratio could not be fitted on some
## fold (the caller then falls back to the plug-in ratio).
siv_kliep_ratio_grid <- function(Y, A, Z, X, grid, folds, sigma = "median",
                                 kernel_num = 100L, fold = 5L,
                                 dr_floor = 1e-4, min_arm = 10L) {
  n <- length(Y); M <- length(grid); p <- ncol(X)
  RY <- matrix(NA_real_, n, M)
  ## Note on sample sizes.  The numerator sample is (Y, X) among A == 0 AND
  ## Z == 1.  A strong instrument makes that arm SMALL by construction -- under
  ## the manuscript's DGP, log gamma = -2 + 4Z, so almost everyone with Z = 1 is
  ## treated and only ~15% of the untreated have Z = 1.  KLIEP is therefore
  ## fitted on a few dozen points even at n = 1000; `kernel_num` is capped at
  ## that sample size inside siv_kliep().  `min_arm` is a floor below which the
  ## fit is not attempted at all and the caller falls back to the plug-in ratio.
  arm <- integer(0)
  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (!length(tr)) tr <- te
    i1 <- intersect(tr, which(A == 0 & Z == 1))
    i0 <- intersect(tr, which(A == 0 & Z == 0))
    arm <- c(arm, length(i1), length(i0))
    if (length(i1) < min_arm || length(i0) < min_arm) {
      return(NULL)
    }
    x1 <- cbind(Y[i1], X[i1, , drop = FALSE])
    x2 <- cbind(Y[i0], X[i0, , drop = FALSE])
    fit <- tryCatch(siv_kliep(x1, x2, n_y_cols = 1L, sigma = sigma,
                              kernel_num = kernel_num, fold = fold),
                    error = function(e) NULL)
    if (is.null(fit) || !all(is.finite(fit$kernel_weights)) ||
        sum(fit$kernel_weights) <= 0) return(NULL)
    ## evaluate on the full (unit, grid) product
    nte <- length(te)
    big <- cbind(rep(grid, each = nte),
                 X[rep(te, times = M), , drop = FALSE])
    r <- siv_kliep_predict(big, fit)
    RY[te, ] <- matrix(r, nte, M)
  }
  if (anyNA(RY) || !all(is.finite(RY))) return(NULL)
  ## the `dr.LB` guard of the replication code
  if (dr_floor > 0) RY <- pmin(pmax(RY, dr_floor), 1 / dr_floor)
  RY
}
