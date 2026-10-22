## ----------------------------------------------------------------------------
## Core building blocks of the PROPOSED estimator of Park and Tchetgen Tchetgen,
## "Nonparametric inference with an instrumental variable under a separable
## binary treatment choice model".
##
## Everything in this file is a direct transcription of the numbered results of
## that manuscript; the corresponding equation numbers are given inline.
##
##   Theorem 3.1, eq. (14)  -> siv_Psi()            (the mapping Psi)
##   Theorem 3.1, eq. (15)  -> R_Y (density ratio), computed by the caller
##   Theorem 3.1, eq. (16)  -> siv_fixed_point()    (the fixed point g)
##   Theorem 3.1, eq. (17)  -> siv_fixed_point()$delta   (odds ratio function)
##   Theorem 3.1, eq. (18)  -> siv_fixed_point()$gamma   (baseline odds)
##   Theorem 3.2, eq. (19)  -> siv_fixed_point(solver = "contraction")
##   Theorem 4.2 / Supp. eqs. (S.129)-(S.134) -> siv_eif()
##   Theorem 5.1, eq. (22)  -> siv_att_from_eif()
##
## Grid conventions
## ----------------
## The outcome support is represented by a grid `grid` of length M together
## with integration weights `weights` of length M, so that
##      \int h(y) dy   is approximated by   sum_m h(grid[m]) * weights[m].
## For a binary or discrete outcome the grid is the observed support and the
## weights are all 1 (counting measure); for a continuous outcome the grid is a
## fine grid and the weights are the mid-point cell widths.  This is exactly the
## convention used in the replication code of the manuscript (`Y.grid`, `Y.bw`).
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
## Compiled fast paths
##
## `src/siv_fastpath.cpp` provides compiled equivalents of the contraction of
## Theorem 3.2 and of the delta-weighted grid integrals of siv_eif().  The R
## implementations below remain the reference: they are always defined, they
## are what `tests/testthat/test-rcpp.R` checks the compiled code against, and
## they are used automatically whenever the compiled code is unavailable --
## which is the case when the package's R files are sourced directly rather
## than installed (as `0.SeparableIV.R` does in development mode).
##
## Set `options(SIV.use_cpp = FALSE)` to force the R path, e.g. to compare the
## two implementations.
## ----------------------------------------------------------------------------

## Internal: is the compiled fast path usable?  Probed once and cached.
siv_cpp_ok <- local({
  cache <- NULL
  function() {
    if (!isTRUE(getOption("SIV.use_cpp", TRUE))) return(FALSE)
    if (is.null(cache)) {
      cache <<- isTRUE(tryCatch({
        m <- matrix(1, 1L, 1L)
        siv_contraction_cpp(m, 1, 1, m, m, m, 1, 1L, 1e-10)
        siv_delta_moments_cpp(m, 0, m, 1)
        TRUE
      }, error = function(e) FALSE))
    }
    cache
  }
})

## Internal: R reference implementation of the contraction of Theorem 3.2,
## eq. (19).  `siv_contraction_cpp()` reproduces this exactly.
siv_contraction_r <- function(g, pA1_zm, pA1_zM, fY_zm, fY_zM, RY, weights,
                              maxit, tol) {
  n <- nrow(g); M <- ncol(g)
  Wm <- matrix(weights, n, M, byrow = TRUE)
  converged <- rep(FALSE, n)
  iter <- 0L
  for (j in seq_len(maxit)) {
    iter <- j
    gn <- siv_Psi(g, pA1_zm, pA1_zM, fY_zm, fY_zM, RY, weights)
    gn <- pmax(gn, 0)
    den <- as.vector(gn %*% weights)
    ## Units whose update degenerated keep their previous (already
    ## normalised) value; they are reported as non-converged.
    bad <- !is.finite(den) | den <= 0
    gn[bad, ] <- g[bad, , drop = FALSE]
    den[bad] <- 1
    gn <- gn / den
    ## convergence in the (weighted) L1 distance used by the reference code
    dif <- rowSums(abs(gn - g) * Wm)
    g <- gn
    converged <- (dif < tol) & !bad
    if (all(converged)) break
  }
  list(g = g, iter = iter, converged = converged)
}

## Internal: the four delta-weighted grid integrals used by siv_eif(), for one
## level of Z.  Returns an n x 4 matrix (denom, y, alpha, alphay).
siv_delta_moments_r <- function(delta, grid, fY, weights) {
  n <- nrow(delta); M <- ncol(delta)
  W <- matrix(weights, n, M, byrow = TRUE)
  Ygrid <- matrix(grid, n, M, byrow = TRUE)
  base <- delta * fY * W
  out <- cbind(rowSums(base), rowSums(Ygrid * base),
               rowSums(delta * base), rowSums(Ygrid * delta * base))
  colnames(out) <- c("denom", "y", "alpha", "alphay")
  out
}

## Internal: dispatch to the compiled implementation when available.
siv_delta_moments <- function(delta, grid, fY, weights) {
  if (siv_cpp_ok()) siv_delta_moments_cpp(delta, grid, fY, weights)
  else siv_delta_moments_r(delta, grid, fY, weights)
}

## Internal: one application of the mapping Psi of Theorem 3.1, eq. (14),
## vectorised over units.  `zm` / `zM` are length-n 0/1 vectors identifying, for
## each unit, the level of Z minimising / maximising f(A = 1 | Z, X) (see the
## paragraph preceding Theorem 3.1).
##
##   Psi{g(y,X); theta} =
##     [ f(A=1|Z=zm,X) g(y,X)
##       + f(A=0|Z=zm,X) f(y|A=0,Z=zm,X) - f(A=0|Z=zM,X) f(y|A=0,Z=zM,X) ]
##     * \int g(t,X) R_Y(t,X) dt   /   [ f(A=1|Z=zM,X) R_Y(y,X) ]
##
## with  R_Y(y,X) = f(y|A=0,Z=zM,X) / f(y|A=0,Z=zm,X)   (eq. 15).
siv_Psi <- function(g, pA1_zm, pA1_zM, fY_zm, fY_zM, RY, weights) {
  pA0_zm <- 1 - pA1_zm
  pA0_zM <- 1 - pA1_zM
  ## \int g(t,X) R_Y(t,X) dt , one number per unit
  S <- as.vector((g * RY) %*% weights)
  num <- pA1_zm * g + pA0_zm * fY_zm - pA0_zM * fY_zM
  (num * S) / (pA1_zM * RY)
}

#' Solve the fixed-point equation of Theorem 3.1
#'
#' Given the observed-data nuisance functions
#' \eqn{\theta = \{f_Z, f_A, f_Y\}} evaluated on an outcome grid, this function
#' returns the fixed point \eqn{g(\cdot, X;\theta)} of the mapping \eqn{\Psi}
#' defined in eq. (14) of Park and Tchetgen Tchetgen, together with the induced
#' odds-ratio function \eqn{\delta} (eq. 17) and baseline odds function
#' \eqn{\gamma} (eq. 18).
#'
#' Two solvers are available.  `"contraction"` implements the normalised
#' iterative update of Theorem 3.2, eq. (19),
#' \eqn{h^{(j+1)} = \tilde\Psi(h^{(j)})} with
#' \eqn{\tilde\Psi(g) = \Psi(g) / \int \Psi(g)}, which Theorem 3.2 shows
#' converges to the fixed point from *any* non-negative starting function
#' integrating to one.  `"closed-form"` is available only for a binary outcome
#' (`M = 2`), where the fixed-point equation reduces to a scalar quadratic and
#' can be solved exactly; this reproduces the `g0y1.binary()` routine of the
#' manuscript's replication code.
#'
#' @param pA1_z0,pA1_z1 Numeric vectors of length `n` holding
#'   \eqn{f(A = 1 \mid Z = 0, X_i)} and \eqn{f(A = 1 \mid Z = 1, X_i)}.
#' @param fY_z0,fY_z1 Numeric `n x M` matrices holding
#'   \eqn{f(Y = y_m \mid A = 0, Z = 0, X_i)} and
#'   \eqn{f(Y = y_m \mid A = 0, Z = 1, X_i)}, normalised to integrate to one
#'   against `weights`.
#' @param grid Numeric vector of length `M` giving the outcome grid.
#' @param weights Numeric vector of length `M` of integration weights.
#' @param y_ref Reference outcome value \eqn{y_R} at which \eqn{\delta} is
#'   normalised to one (eq. 2).  Must belong to `grid`.
#' @param solver Either `"contraction"` (Theorem 3.2) or `"closed-form"`
#'   (binary outcome only).
#' @param maxit,tol Maximum number of iterations and convergence tolerance for
#'   the contraction.
#' @param g_init Optional `n x M` matrix of starting values.  Defaults to the
#'   uniform density on the grid, which is admissible by Theorem 3.2.
#' @param density_floor Small positive number used to floor densities and
#'   density ratios for numerical stability (Assumption (A3) requires them to be
#'   bounded away from zero).
#' @param g_floor Lower bound on the fixed point `g`, expressed as a fraction of
#'   each unit's maximum, enforcing the part of Assumption (A3) that requires
#'   `g` to be bounded away from zero.  The returned `g_floor_frac` reports how
#'   often the bound binds; a non-negligible value indicates that the supplied
#'   nuisance estimates violate (A3).
#' @param dr_floor Lower bound (and, reciprocally, upper bound) applied to the
#'   estimated density ratio \eqn{R_Y} of eq. (15).  A ratio of two separately
#'   estimated conditional densities is unstable in the tails; the replication
#'   code of the manuscript guards against this with `dr.LB = 1e-2`.  [SIV()]
#'   uses `1e-4` for internally estimated nuisance functions.  `0`
#'   disables the guard, which is appropriate when the densities are known
#'   exactly.
#' @param RY Optional `n x M` matrix supplying \eqn{R_Y} directly, e.g. from a
#'   targeted density-ratio estimator (Section 5.1 of the manuscript recommends
#'   estimating the ratio directly rather than as a ratio of density estimates).
#'   It must be given as \eqn{f(y \mid A=0,Z=1,X)/f(y \mid A=0,Z=0,X)}; the
#'   function inverts it for units whose \eqn{z_M = 0}.
#'
#' @return A list with elements `g` (`n x M` fixed point, normalised to
#'   integrate to one), `delta` (`n x M` odds-ratio function, eq. 17),
#'   `gamma_z0`, `gamma_z1` (baseline odds function, eq. 18), `iter`,
#'   `converged` and `zM` (the level of `Z` maximising the treatment
#'   probability, per unit), together with `contraction_modulus`.
#'
#' @examples
#' ## Two-point outcome support: the closed-form and iterative solvers agree.
#' set.seed(1)
#' n <- 5
#' pA1_z0 <- rep(0.2, n); pA1_z1 <- rep(0.6, n)
#' p0 <- runif(n, 0.3, 0.7); p1 <- runif(n, 0.3, 0.7)
#' fY_z0 <- cbind(1 - p0, p0); fY_z1 <- cbind(1 - p1, p1)
#' a <- siv_fixed_point(pA1_z0, pA1_z1, fY_z0, fY_z1, c(0, 1), c(1, 1),
#'                      solver = "closed-form")
#' b <- siv_fixed_point(pA1_z0, pA1_z1, fY_z0, fY_z1, c(0, 1), c(1, 1),
#'                      solver = "contraction")
#' max(abs(a$delta - b$delta))
#' @export
siv_fixed_point <- function(pA1_z0, pA1_z1, fY_z0, fY_z1, grid, weights,
                            y_ref = 0,
                            solver = c("contraction", "closed-form"),
                            maxit = 1000L, tol = 1e-10,
                            g_init = NULL, density_floor = 1e-4,
                            dr_floor = 1e-4, RY = NULL, g_floor = 1e-8) {
  solver <- match.arg(solver)
  n <- length(pA1_z0)
  M <- length(grid)
  stopifnot(length(weights) == M, is.matrix(fY_z0), is.matrix(fY_z1))
  stopifnot(nrow(fY_z0) == n, ncol(fY_z0) == M, nrow(fY_z1) == n, ncol(fY_z1) == M)
  iref <- which.min(abs(grid - y_ref))
  if (abs(grid[iref] - y_ref) > 1e-8) {
    stop("`y_ref` must be one of the grid points.", call. = FALSE)
  }

  fY_z0 <- pmax(fY_z0, density_floor)
  fY_z1 <- pmax(fY_z1, density_floor)

  ## zM = argmax_z f(A = 1 | Z = z, X); zm is the other level (Section 3.2).
  zM <- as.integer(pA1_z1 >= pA1_z0)
  pA1_zM <- ifelse(zM == 1L, pA1_z1, pA1_z0)
  pA1_zm <- ifelse(zM == 1L, pA1_z0, pA1_z1)
  fY_zM <- fY_z1; fY_zM[zM == 0L, ] <- fY_z0[zM == 0L, , drop = FALSE]
  fY_zm <- fY_z0; fY_zm[zM == 0L, ] <- fY_z1[zM == 0L, , drop = FALSE]
  ## eq. (15).  A user-supplied density ratio (e.g. from a targeted
  ## density-ratio estimator, which Section 5.1 of the manuscript recommends
  ## over the ratio of two separate density estimates) takes precedence.
  if (is.null(RY)) {
    RY <- fY_zM / fY_zm
  } else {
    RY <- as.matrix(RY)
    stopifnot(nrow(RY) == n, ncol(RY) == M)
    RY[zM == 0L, ] <- 1 / RY[zM == 0L, , drop = FALSE]
  }
  if (dr_floor > 0) {
    ## Floor and cap the estimated density ratio, then rebuild and renormalise
    ## the numerator density so that the two stay mutually consistent.  This
    ## mirrors the `dr.LB` guard in the replication code
    ## (Simulation_F_Both/1.Sim_Conti_SepIV.R), and is needed because a ratio of
    ## two independently estimated densities is unstable in the tails.
    RY <- pmin(pmax(RY, dr_floor), 1 / dr_floor)
    fY_zM <- RY * fY_zm
    fY_zM <- fY_zM / as.vector(fY_zM %*% weights)
    RY <- fY_zM / fY_zm
    ## keep fY_z0 / fY_z1 consistent with the adjusted numerator
    fY_z1[zM == 1L, ] <- fY_zM[zM == 1L, , drop = FALSE]
    fY_z0[zM == 0L, ] <- fY_zM[zM == 0L, , drop = FALSE]
  }
  RY <- pmax(RY, density_floor)

  if (solver == "closed-form") {
    if (M != 2L) {
      stop("`solver = \"closed-form\"` is only available for a two-point ",
           "outcome support (M = 2).", call. = FALSE)
    }
    g <- siv_fixed_point_binary(pA1_zm, pA1_zM, fY_zm, fY_zM, RY, weights)
    iter <- 1L; converged <- rep(TRUE, n)
  } else {
    if (is.null(g_init)) {
      g_init <- matrix(1 / sum(weights), n, M)
    }
    g <- g_init / as.vector(g_init %*% weights)
    ct <- if (siv_cpp_ok()) {
      siv_contraction_cpp(g, pA1_zm, pA1_zM, fY_zm, fY_zM, RY, weights,
                          as.integer(maxit), tol)
    } else {
      siv_contraction_r(g, pA1_zm, pA1_zM, fY_zm, fY_zM, RY, weights,
                        maxit, tol)
    }
    g <- ct$g; iter <- ct$iter; converged <- ct$converged
  }

  ## Assumption (A3) requires g(y, X) to be bounded away from zero.  Enforce it
  ## relative to each unit's own maximum, so that delta = g / f(y|A=0,Z=zm,X)
  ## cannot be exactly zero (which would make log delta = -Inf and degrade every
  ## downstream integral).  `g_floor_frac` records how often the floor binds:
  ## a non-negligible value means the nuisance estimates violate (A3).
  gmax <- apply(g, 1L, max)
  gmin <- g_floor * gmax
  n_floored <- sum(g < gmin)
  if (n_floored > 0L) {
    g <- pmax(g, gmin)
    g <- g / as.vector(g %*% weights)
  }
  g_floor_frac <- n_floored / length(g)

  ## eq. (17): delta(y, X) = {g(y,X)/g(y_R,X)} * {f(y_R|A=0,Z=zm,X)/f(y|A=0,Z=zm,X)}
  dimnames(g) <- NULL
  delta <- g / fY_zm
  delta <- delta / delta[, iref]
  dimnames(delta) <- NULL

  ## eq. (18): gamma(z, X) = {f(A=1|Z=z,X)/f(A=0|Z=z,X)} / \int delta(y,X) f(y|A=0,Z=z,X) dy
  denom_z0 <- as.vector((delta * fY_z0) %*% weights)
  denom_z1 <- as.vector((delta * fY_z1) %*% weights)
  gamma_z0 <- (pA1_z0 / (1 - pA1_z0)) / denom_z0
  gamma_z1 <- (pA1_z1 / (1 - pA1_z1)) / denom_z1

  ## At the fitted fixed point, this is the largest value of lambda(y, X) in the
  ## contraction proof.  Values close to one indicate a weak contraction margin
  ## and potentially slow or numerically sensitive iteration.  The integral S
  ## depends on the fitted g, so this quantity does not define an exogenous
  ## pointwise lower bound on R_Y.
  S_scale <- as.vector((g * RY) %*% weights)
  contraction_modulus <- max((pA1_zm * S_scale / pA1_zM) / RY)

  list(g = g, delta = delta, gamma_z0 = gamma_z0, gamma_z1 = gamma_z1,
       fY_z0 = fY_z0, fY_z1 = fY_z1,
       denom_z0 = denom_z0, denom_z1 = denom_z1,
       iter = iter, converged = converged, zM = zM, RY = RY,
       g_floor_frac = g_floor_frac,
       contraction_modulus = contraction_modulus,
       fY_zm = fY_zm, fY_zM = fY_zM, y_ref_index = iref)
}

## Internal: exact solution of the fixed-point equation when the outcome support
## has two points.  Writing g = (g_1, g_2) with g_1 + g_2 = 1 (unit weights) and
## D(y) = f(A=0|Z=zm,X) f(y|A=0,Z=zm,X) - f(A=0|Z=zM,X) f(y|A=0,Z=zM,X), the
## fixed-point equation Psi(g) = g evaluated at the second support point reduces
## to the quadratic
##   f(A=1|Z=zm) {R_Y(y_2) - R_Y(y_1)} g_2^2
## + [ D(y_2){R_Y(y_2) - R_Y(y_1)} + f(A=1|Z=zm) R_Y(y_1)
##     - f(A=1|Z=zM) R_Y(y_2) ] g_2
## + D(y_2) R_Y(y_1) = 0 ,
## which is exactly the routine `g0y1.binary()` of the replication code.
siv_fixed_point_binary <- function(pA1_zm, pA1_zM, fY_zm, fY_zM, RY, weights) {
  n <- length(pA1_zm)
  w1 <- weights[1L]; w2 <- weights[2L]
  D2 <- (1 - pA1_zm) * fY_zm[, 2L] - (1 - pA1_zM) * fY_zM[, 2L]
  r1 <- RY[, 1L]; r2 <- RY[, 2L]
  ## Coefficients of  c1 g2^2 + c2 g2 + c3 = 0 .  With unit weights these are
  ## exactly `eq.c1`, `eq.c2`, `eq.c3` of `g0y1.binary()` in the replication code.
  c1 <- pA1_zm * w2 * (r2 - r1)
  c2 <- D2 * w2 * (r2 - r1) + pA1_zm * r1 - pA1_zM * r2
  c3 <- D2 * r1
  ub <- 1 / w2                       ## g2 <= 1/w2 because g1 w1 + g2 w2 = 1
  mid <- ub / 2
  g2 <- numeric(n)
  quad <- abs(c1) > 1e-6
  ## degenerate (linear) case
  g2[!quad] <- -c3[!quad] / c2[!quad]
  ## quadratic case: take the root lying in the admissible interval [0, 1/w2]
  disc <- pmax(c2^2 - 4 * c1 * c3, 0)
  s1 <- (-c2 + sqrt(disc)) / (2 * c1)
  s2 <- (-c2 - sqrt(disc)) / (2 * c1)
  pick <- abs(s1 - mid) <= mid
  g2[quad] <- ifelse(pick[quad], s1[quad], s2[quad])
  g2 <- pmin(pmax(g2, 1e-12), ub - 1e-12)
  g <- cbind((1 - g2 * w2) / w1, g2)
  g / as.vector(g %*% weights)
}

#' Efficient influence function for the ATT under the logit-separable IV model
#'
#' Implements the efficient influence function of Theorem 4.2 of Park and
#' Tchetgen Tchetgen -- equivalently eq. (21) of the main paper evaluated at the
#' unique \eqn{\omega} of Theorem 4.2, whose closed form is given by
#' eqs. (S.129)-(S.134) of the Supplementary Material.
#'
#' The linear system solved for \eqn{(B(X), C(X))} is eq. (S.131);
#' \eqn{L(Y^{(0)},X)} is the display below eq. (S.133); and
#' \eqn{\omega = \{L - E(L\mid X)\}\{Z - E(Z\mid X)\}} is eq. (S.133).
#'
#' @param Y,A,Z Observed outcome, treatment and instrument (length `n`).
#' @param grid,weights Outcome grid and integration weights (length `M`).
#' @param delta `n x M` matrix of odds-ratio values \eqn{\delta(y_m, X_i)}.
#' @param gamma_z0,gamma_z1 Baseline odds \eqn{\gamma(z, X_i)}.
#' @param fY_z0,fY_z1 `n x M` matrices of \eqn{f(y \mid A = 0, Z = z, X)}.
#' @param pZ1 Vector of \eqn{f(Z = 1 \mid X_i)}.
#' @param pA1_z0,pA1_z1 Vectors of \eqn{f(A = 1 \mid Z = z, X_i)}.
#'
#' @return A list with the per-unit numerator of the influence function
#'   (`num`), the estimated ATT (`att`), the centred influence function
#'   (`eif`), the variance estimate (`var`) and several intermediate
#'   quantities.
#' @keywords internal
#' @export
siv_eif <- function(Y, A, Z, grid, weights, delta, gamma_z0, gamma_z1,
                    fY_z0, fY_z1, pZ1, pA1_z0, pA1_z1) {
  n <- length(Y)
  M <- length(grid)
  W <- matrix(weights, n, M, byrow = TRUE)
  Ygrid <- matrix(grid, n, M, byrow = TRUE)
  pZ0 <- 1 - pZ1

  ## Indicator matrix picking out the observed outcome among the untreated.
  pos <- matrix(0, n, M)
  idx <- which(A == 0)
  jj <- vapply(Y[idx], function(y) which.min(abs(grid - y)), integer(1L))
  pos[cbind(idx, jj)] <- 1

  ## Conditional expectations mu(z, X; g) = E{g alpha | A=0,z,X} / E{alpha | A=0,z,X}
  ## which is eq. (4) of the main paper.
  int_z <- function(V, z) {
    f <- if (z == 0) fY_z0 else fY_z1
    as.vector(rowSums(V * f * W))
  }
  ## The four delta-weighted integrals per level of Z are computed in a single
  ## sweep over the grid (compiled when available); see `siv_delta_moments()`.
  m0 <- siv_delta_moments(delta, grid, fY_z0, weights)
  m1 <- siv_delta_moments(delta, grid, fY_z1, weights)
  denom_z0 <- as.vector(m0[, "denom"])
  denom_z1 <- as.vector(m1[, "denom"])
  mu_z0 <- as.vector(m0[, "y"]) / denom_z0
  mu_z1 <- as.vector(m1[, "y"]) / denom_z1
  mu_alpha_z0 <- as.vector(m0[, "alpha"]) / denom_z0
  mu_alpha_z1 <- as.vector(m1[, "alpha"]) / denom_z1
  mu_alphay_z0 <- as.vector(m0[, "alphay"]) / denom_z0
  mu_alphay_z1 <- as.vector(m1[, "alphay"]) / denom_z1

  ## beta(z, X) = {f(A=1|z,X)/f(A=0|z,X)} / E{alpha | A = 0, z, X}
  ##
  ## This is the `beta.z0.CF <- gamma.z0.CF/denom.z0.CF` of the manuscript's
  ## replication code (Simulation_F_Both/1.Sim_Conti_SepIV.R, "Calibration:
  ## beta"), in which `gamma.z0.CF` is the RAW odds Pr(A=1|z,X)/Pr(A=0|z,X).
  ##
  ## It is computed here from `pA1_z*` rather than from the `gamma_z*`
  ## arguments so that the scaling is unambiguous.  Note that it coincides with
  ## the `gamma_z0`/`gamma_z1` returned by siv_fixed_point() (eq. 18), which are
  ## already divided by `denom` once.  Dividing those by `denom` a SECOND time
  ## -- which is what this function used to do -- makes the weight
  ## alpha * beta too small by a factor of `denom` (observed range 1.04-2.31),
  ## which breaks the identity that f(Y0=y|X) reconstructed from the z = 0 and
  ## the z = 1 arm must agree: that residual was 9e-04 before this correction
  ## and is 1e-25 after it.  The estimator stayed consistent either way
  ## (E[(1-A) * weight * (Y - mu)] = 0 for any scalar multiple of the weight),
  ## so only the influence function -- and hence every reported standard error
  ## -- was affected.
  beta_z0 <- (pA1_z0 / (1 - pA1_z0)) / denom_z0
  beta_z1 <- (pA1_z1 / (1 - pA1_z1)) / denom_z1

  ## weight = alpha * beta = delta(y,X) gamma(z,X) (the object multiplying
  ## (1 - A) in eq. 21).
  wgt_z0 <- delta * beta_z0
  wgt_z1 <- delta * beta_z1
  wgt_obs_z0 <- rowSums(wgt_z0 * pos)
  wgt_obs_z1 <- rowSums(wgt_z1 * pos)

  ## Full-data density f(Y0 = y | X) reconstructed from theta via eq. (9).
  fya0z0 <- fY_z0 * ((1 - pA1_z0) * pZ0)
  fya0z1 <- fY_z1 * ((1 - pA1_z1) * pZ1)
  fya1z0 <- wgt_z0 * fya0z0
  fya1z1 <- wgt_z1 * fya0z1
  fy_from_z0 <- (fya0z0 + fya1z0) / pZ0
  fy_from_z1 <- (fya0z1 + fya1z1) / pZ1
  fY <- (fy_from_z0 + fy_from_z1) / 2
  mean_y <- function(V) as.vector(rowSums(V * fY * W))

  ## b(z = 1, y, X) = Pr(Z = 1 | Y0 = y, A = 1, X)
  b_z1 <- (1 + wgt_z0) * (beta_z1 * pZ1) /
    (beta_z1 * pZ1 + beta_z0 * pZ0 + delta * (beta_z1 * beta_z0))

  ## (B(X), C(X)) from the 2 x 2 linear system (S.131).  The second row is
  ## scaled by beta(0, X) on both sides, which leaves the solution unchanged.
  A11 <- beta_z1 * pZ1 + beta_z0 * beta_z1 * pZ1 * mu_alpha_z0
  A12 <- beta_z1 * pZ1 + 2 * beta_z0 * pZ0 - beta_z1 * pZ0 + beta_z0 * beta_z1 * mu_alpha_z0
  A21 <- beta_z0 * pZ1 + beta_z0 * beta_z1 * pZ1 * mu_alpha_z1
  A22 <- beta_z0 + beta_z0 * beta_z1 * mu_alpha_z1
  b1 <- beta_z1 * pZ1 * mu_z0 + 2 * beta_z0 * pZ0 * mu_z0 - beta_z1 * pZ0 * mu_z0 +
    beta_z0 * beta_z1 * mu_alphay_z0
  b2 <- beta_z0 * mu_z1 + beta_z0 * beta_z1 * mu_alphay_z1
  det <- A11 * A22 - A12 * A21
  Bx <- (-b1 * A22 + b2 * A12) / det
  Cx <- (-b2 * A11 + b1 * A21) / det

  ## L(Y0, X) below eq. (S.133)
  coef1 <- 1 / (b_z1 - pZ1)
  coef2 <- b_z1 / (b_z1 - pZ1)
  L <- coef1 * Ygrid + coef2 * Bx + coef1 * Cx

  mu_L <- mean_y(L)
  ## omega(Y0, Z, X) = {L - E(L | X)} {Z - E(Z | X)}   (eq. S.133)
  w_z0 <- (L - mu_L) * (0 - pZ1)
  w_z1 <- (L - mu_L) * (1 - pZ1)
  w_obs_z0 <- rowSums(w_z0 * pos)
  w_obs_z1 <- rowSums(w_z1 * pos)
  mu_w_z0 <- int_z(w_z0 * delta, 0) / denom_z0
  mu_w_z1 <- int_z(w_z1 * delta, 1) / denom_z1

  mu   <- ifelse(Z == 0, mu_z0, mu_z1)
  wgt  <- ifelse(Z == 0, wgt_obs_z0, wgt_obs_z1)
  om   <- ifelse(Z == 0, w_obs_z0, w_obs_z1)
  mu_w <- ifelse(Z == 0, mu_w_z0, mu_w_z1)

  ## Numerator of eq. (21) / eq. (S.134)
  num <- A * Y -
    ((1 - A) * wgt * (Y - mu) + A * mu -
       (1 - A) * wgt * (om - mu_w) - A * mu_w - (1 - A) * om)

  pA <- mean(A)
  att <- mean(num) / pA
  eif <- (num - att * A) / pA
  ## Assumption (A5) (strong relevance) is what keeps `b_z1 - pZ1` -- the
  ## denominator of L -- away from zero; report the realised minimum.
  min_relevance <- min(abs(b_z1 - pZ1))

  list(num = num, att = att, eif = eif, var = stats::var(eif) / n,
       min_relevance = min_relevance,
       delta = delta, beta_z0 = beta_z0, beta_z1 = beta_z1,
       mu_z0 = mu_z0, mu_z1 = mu_z1, B = Bx, C = Cx, L = L, omega_obs = om,
       fY = fY)
}
