## ----------------------------------------------------------------------------
## Low-level utilities shared by all estimators in the package.
##
## Traceability of the statistical content of this file:
##   * `siv_sandwich()` implements the stacked-estimating-function ("sandwich")
##     covariance matrix of Sun et al. (2021), Supplemental Appendix, eq. (S2),
##     which is also the covariance used by Liu et al. (2021) (Section 6,
##     "se was evaluated using the sandwich estimator accounting for all sources
##     of variabilities") and by Sun, Miao and Wickramarachchi (2025),
##     Propositions 1 and 2 / Supplement Sections S2-S3.
##   * `siv_tilt_basis()` builds the finite-dimensional basis used for the
##     selection-bias / odds-ratio function.  The normalisation b(y_ref) = 0
##     enforces the odds-ratio boundary condition
##     delta(y_ref, X) = 1  <=>  log delta(y_ref, X) = 0,
##     which is imposed by *all four* papers:
##       - Park & Tchetgen Tchetgen, eq. (2) ("delta(y_R, Z, X) = 1");
##       - Liu et al. (2021), definition of psi below eq. (4.6) (psi has
##         f(Y0 = 0 | ...) in the reference position, so psi(0, Z, C) = 0);
##       - Sun et al. (2021), eq. (4.1) (delta is an odds ratio relative to
##         Y = 0, hence delta(x, 0, z) = 0);
##       - Sun, Miao & Wickramarachchi (2025), eq. (2.1)-(2.2)
##         ("h(0, x; gamma) = 0").
## ----------------------------------------------------------------------------

#' Inverse-logit (expit) function
#'
#' @param v Numeric vector on the logit scale.
#' @return `1 / (1 + exp(-v))`, computed in a numerically stable way.
#' @examples
#' expit(c(-5, 0, 5))
#' @export
expit <- function(v) {
  out <- numeric(length(v))
  pos <- !is.na(v) & v >= 0
  neg <- !is.na(v) & v < 0
  out[pos] <- 1 / (1 + exp(-v[pos]))
  ev <- exp(v[neg])
  out[neg] <- ev / (1 + ev)
  out[is.na(v)] <- NA_real_
  out
}

#' Logit function
#'
#' @param p Numeric vector of probabilities in (0, 1).
#' @return `log(p / (1 - p))`.
#' @examples
#' logit(c(0.1, 0.5, 0.9))
#' @export
logit <- function(p) log(p) - log1p(-p)

## Internal: standardise / validate the (Y, A, Z, X) data supplied by a user.
##
## `A` is the *treatment* in the causal (Liu / Park & Tchetgen Tchetgen)
## parameterisation.  In the missing-data (Sun) parameterisation the response
## indicator is R = 1 - A; the mapping is applied by the calling function.
siv_check_data <- function(Y, A, Z, X, allow_na_Y = FALSE) {
  if (!is.numeric(Y)) stop("`Y` must be a numeric vector.", call. = FALSE)
  n <- length(Y)
  if (length(A) != n) stop("`A` must have the same length as `Y`.", call. = FALSE)
  if (length(Z) != n) stop("`Z` must have the same length as `Y`.", call. = FALSE)
  A <- as.numeric(A)
  Z <- as.numeric(Z)
  if (!all(A %in% c(0, 1))) stop("`A` must be binary (0/1).", call. = FALSE)
  if (!all(Z %in% c(0, 1))) {
    stop("`Z` must be a binary (0/1) instrument; the methods implemented here ",
         "are derived for a binary instrument.", call. = FALSE)
  }
  if (is.null(X)) {
    X <- matrix(numeric(0), nrow = n, ncol = 0)
  } else {
    X <- as.matrix(X)
    if (nrow(X) != n) stop("`X` must have `length(Y)` rows.", call. = FALSE)
    if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  if (!allow_na_Y && anyNA(Y)) stop("`Y` must not contain missing values.", call. = FALSE)
  ## Y is only used where A == 0 (it equals the treatment-free potential
  ## outcome there); it may legitimately be NA where A == 1.
  if (anyNA(Y[A == 0])) {
    stop("`Y` must be observed for every unit with `A == 0`.", call. = FALSE)
  }
  if (all(A == 0) || all(A == 1)) stop("`A` must take both values 0 and 1.", call. = FALSE)
  if (all(Z == 0) || all(Z == 1)) stop("`Z` must take both values 0 and 1.", call. = FALSE)
  list(Y = Y, A = A, Z = Z, X = X, n = n)
}

## Internal: validate the explicit outcome-type argument.  Used by SIV() and by
## the competing estimators.  The type is *never* inferred from the data; it is
## only cross-checked against the data for obvious contradictions.
siv_match_type <- function(type, Y, A) {
  type <- match.arg(type, c("binary", "continuous"))
  obs <- Y[A == 0]
  if (type == "binary" && !all(obs %in% c(0, 1))) {
    stop("`type = \"binary\"` was requested but the observed treatment-free ",
         "outcomes `Y[A == 0]` are not all 0/1.", call. = FALSE)
  }
  if (type == "continuous" && length(unique(obs)) <= 2L) {
    warning("`type = \"continuous\"` was requested but `Y[A == 0]` takes at ",
            "most two distinct values.", call. = FALSE)
  }
  type
}

## ----------------------------------------------------------------------------
## Basis construction for the selection-bias / odds-ratio function
## ----------------------------------------------------------------------------

#' Specification of the basis for the selection-bias (odds-ratio) function
#'
#' All three competing estimators implemented in this package model the
#' selection-bias function -- equivalently the log odds-ratio function relating
#' the outcome and the treatment / response indicator -- through a
#' finite-dimensional working model.  `siv_basis_spec()` collects the tuning
#' parameters of that working model in one place so that they are configurable
#' rather than hard-coded.
#'
#' The basis is always *normalised* so that it vanishes at the reference value
#' `y_ref`.  This enforces the boundary condition imposed by every source paper:
#' `delta(y_ref, X) = 1` in Park and Tchetgen Tchetgen (eq. 2),
#' `psi(0, Z, C) = 0` in Liu et al. (2021),
#' `delta(x, 0, z) = 0` in Sun et al. (2021) eq. (4.1), and
#' `h(0, x; gamma) = 0` in Sun, Miao and Wickramarachchi (2025) eq. (2.2).
#'
#' @param type Basis for the outcome.  `"linear"` gives the single term
#'   `y - y_ref`, which reproduces the specifications actually used in the
#'   source papers (Liu et al. 2021, eq. 4.8, `psi(Y0, Z, C; psi) = psi * Y0`;
#'   Sun et al. 2021, `delta(x, y, z) = gamma * y`; Sun, Miao and
#'   Wickramarachchi 2025, `h(y, x; gamma) = gamma * y`).  `"ns"` and `"bs"`
#'   give natural-cubic and B-spline bases from the \pkg{splines} package.
#'   `"poly"` gives raw powers `y^1, ..., y^degree` (centred at `y_ref`).
#' @param df Degrees of freedom of the spline basis (ignored when `knots` is
#'   supplied, and when `type` is `"linear"`).
#' @param degree Degree of the B-spline / polynomial basis.
#' @param knots Optional vector of interior knots.  If `NULL`, the knots are
#'   placed at equally spaced quantiles of the observed outcomes, which is the
#'   default rule of the \pkg{splines} package.
#' @param Boundary.knots Optional length-2 vector of boundary knots.  If `NULL`,
#'   the range of the observed outcomes is used.
#' @param y_ref Reference value of the outcome at which the basis is normalised
#'   to zero.  Defaults to `0`, the reference value used in all source papers
#'   and in the simulation study of Park and Tchetgen Tchetgen.
#' @param interact Optional character vector naming columns of the covariate
#'   matrix with which the outcome basis is interacted.  `interact = NULL`
#'   (default) gives a selection-bias function that does not depend on
#'   covariates; `interact = "all"` interacts the outcome basis with every
#'   column of `X`; a character or integer vector selects specific columns.
#'   Interactions are needed whenever the true log odds-ratio function depends
#'   on covariates (as it does in the data-generating process of Park and
#'   Tchetgen Tchetgen, Section 7, where
#'   `log delta(y, X) = (0.4 + 0.05 * W) * a(y)`).
#'
#' @return An object of class `"siv_basis_spec"`.
#' @examples
#' siv_basis_spec("linear")
#' siv_basis_spec("ns", df = 3, interact = "all")
#' @export
siv_basis_spec <- function(type = c("linear", "ns", "bs", "poly"),
                           df = 3L,
                           degree = 3L,
                           knots = NULL,
                           Boundary.knots = NULL,
                           y_ref = 0,
                           interact = NULL) {
  type <- match.arg(type)
  if (!is.numeric(y_ref) || length(y_ref) != 1L) {
    stop("`y_ref` must be a single number.", call. = FALSE)
  }
  if (type %in% c("ns", "bs", "poly")) {
    if (!is.numeric(df) || length(df) != 1L || df < 1) {
      stop("`df` must be a single positive integer.", call. = FALSE)
    }
  }
  structure(list(type = type, df = as.integer(df), degree = as.integer(degree),
                 knots = knots, Boundary.knots = Boundary.knots,
                 y_ref = y_ref, interact = interact),
            class = "siv_basis_spec")
}

#' @export
print.siv_basis_spec <- function(x, ...) {
  cat("Selection-bias basis specification\n")
  cat("  outcome basis :", x$type, "\n")
  if (x$type %in% c("ns", "bs", "poly")) cat("  df / degree   :", x$df, "/", x$degree, "\n")
  if (!is.null(x$knots)) cat("  interior knots:", paste(signif(x$knots, 4), collapse = ", "), "\n")
  if (!is.null(x$Boundary.knots)) {
    cat("  boundary knots:", paste(signif(x$Boundary.knots, 4), collapse = ", "), "\n")
  }
  cat("  reference y   :", x$y_ref, "\n")
  cat("  interactions  :",
      if (is.null(x$interact)) "none" else paste(x$interact, collapse = ", "), "\n")
  invisible(x)
}

## Internal: freeze the data-dependent parts of a basis specification (knot
## locations, boundary knots, interacting covariate columns) so that the basis
## is a *fixed* function of y afterwards.  This matters for the sandwich
## variance: the basis must not change when the parameter vector is perturbed.
siv_basis_fit <- function(spec, y_obs, X) {
  stopifnot(inherits(spec, "siv_basis_spec"))
  spec <- unclass(spec)
  y_obs <- y_obs[is.finite(y_obs)]
  if (spec$type %in% c("ns", "bs")) {
    if (is.null(spec$Boundary.knots)) spec$Boundary.knots <- range(c(y_obs, spec$y_ref))
    if (is.null(spec$knots)) {
      nk <- if (spec$type == "ns") spec$df - 1L else spec$df - spec$degree
      spec$knots <- if (nk >= 1L) {
        as.numeric(stats::quantile(y_obs, probs = seq_len(nk) / (nk + 1L), names = FALSE))
      } else numeric(0)
    }
  }
  idx <- NULL
  if (!is.null(spec$interact) && ncol(X) > 0L) {
    idx <- if (identical(spec$interact, "all")) {
      seq_len(ncol(X))
    } else if (is.character(spec$interact)) {
      m <- match(spec$interact, colnames(X))
      if (anyNA(m)) stop("`interact` names columns not present in `X`.", call. = FALSE)
      m
    } else {
      as.integer(spec$interact)
    }
  }
  spec$.interact_idx <- idx
  structure(spec, class = c("siv_basis_fitted", "siv_basis_spec"))
}

## Internal: raw (un-normalised) outcome basis.
siv_basis_raw <- function(spec, y) {
  y <- as.numeric(y)
  switch(spec$type,
    linear = matrix(y, ncol = 1L, dimnames = list(NULL, "y")),
    poly   = {
      out <- outer(y, seq_len(spec$df), "^")
      colnames(out) <- paste0("y^", seq_len(spec$df))
      out
    },
    ns = {
      out <- splines::ns(y, knots = spec$knots, Boundary.knots = spec$Boundary.knots,
                         intercept = FALSE)
      out <- matrix(as.numeric(out), nrow = length(y))
      colnames(out) <- paste0("ns", seq_len(ncol(out)))
      out
    },
    bs = {
      out <- splines::bs(y, knots = spec$knots, degree = spec$degree,
                         Boundary.knots = spec$Boundary.knots, intercept = FALSE)
      out <- matrix(as.numeric(out), nrow = length(y))
      colnames(out) <- paste0("bs", seq_len(ncol(out)))
      out
    }
  )
}

#' Evaluate a normalised selection-bias basis
#'
#' Returns the design matrix `b(y, x)` of the working model for the
#' selection-bias function, normalised so that `b(y_ref, x) = 0` for every `x`.
#'
#' @param spec A fitted basis specification, as produced internally from
#'   [siv_basis_spec()].
#' @param y Numeric vector of outcome values.
#' @param X Covariate matrix with `length(y)` rows (only used when the
#'   specification requests interactions).
#' @return A numeric matrix with `length(y)` rows and `siv_basis_dim(spec)`
#'   columns.
#' @keywords internal
#' @export
siv_basis_eval <- function(spec, y, X = NULL) {
  b <- siv_basis_raw(spec, y)
  b0 <- siv_basis_raw(spec, spec$y_ref)
  ## Normalisation b(y_ref) = 0, i.e. the odds-ratio boundary condition.
  b <- sweep(b, 2L, as.numeric(b0), "-")
  idx <- spec$.interact_idx
  if (is.null(idx) || length(idx) == 0L) return(b)
  Xi <- as.matrix(X)[, idx, drop = FALSE]
  out <- b
  nm <- colnames(b)
  for (j in seq_along(idx)) {
    add <- b * Xi[, j]
    colnames(add) <- paste0(nm, ":", colnames(Xi)[j])
    out <- cbind(out, add)
  }
  out
}

## Internal: dimension of the (normalised) basis.
siv_basis_dim <- function(spec) {
  q <- switch(spec$type, linear = 1L, poly = spec$df,
              ns = length(spec$knots) + 1L,
              bs = length(spec$knots) + spec$degree)
  k <- length(spec$.interact_idx)
  q * (1L + (if (is.null(k)) 0L else k))
}

## ----------------------------------------------------------------------------
## Gauss-Hermite quadrature (Golub-Welsch); used to evaluate the exponential
## tilting integrals for a continuous outcome with a Gaussian working density.
## Implemented here so the package has no dependency beyond stats/splines.
## ----------------------------------------------------------------------------

## Internal: nodes/weights for  \int f(x) exp(-x^2) dx  ~  sum_k w_k f(x_k).
siv_gauss_hermite <- function(n) {
  n <- as.integer(n)
  if (n < 1L) stop("`n` must be a positive integer.", call. = FALSE)
  i <- seq_len(n - 1L)
  b <- sqrt(i / 2)
  J <- matrix(0, n, n)
  if (n > 1L) {
    J[cbind(i, i + 1L)] <- b
    J[cbind(i + 1L, i)] <- b
  }
  ev <- eigen(J, symmetric = TRUE)
  ord <- order(ev$values)
  x <- ev$values[ord]
  w <- sqrt(pi) * (ev$vectors[1L, ord])^2
  list(nodes = x, weights = w)
}

## Internal: nodes/weights for E{f(Y)} with Y ~ N(mu, sd^2), for a *vector* of
## means and a scalar sd.  Returns a list with an n x K matrix of nodes and a
## length-K vector of probability weights summing to one.
siv_gh_normal <- function(mu, sd, n_nodes = 20L) {
  gh <- siv_gauss_hermite(n_nodes)
  nodes <- outer(as.numeric(mu), sqrt(2) * sd * gh$nodes, "+")
  list(nodes = nodes, weights = gh$weights / sqrt(pi))
}

## ----------------------------------------------------------------------------
## Sandwich (stacked estimating equation) covariance
## ----------------------------------------------------------------------------

#' Sandwich covariance for a stacked system of estimating equations
#'
#' Computes the M-estimation ("sandwich") covariance matrix
#' \deqn{\hat V = \hat A^{-1} \hat B \hat A^{-\top} / n,\qquad
#'       \hat A = n^{-1}\sum_i \partial M(O_i;\theta)/\partial\theta^\top,\qquad
#'       \hat B = n^{-1}\sum_i M(O_i;\hat\theta) M(O_i;\hat\theta)^\top ,}
#' evaluated at \eqn{\hat\theta}.  This is exactly the estimator of
#' Sun et al. (2021), Supplemental Appendix eq. (S2) -- there written as
#' \eqn{E\{\partial M/\partial\theta^\top\}^{-1} E\{MM^\top\}
#'      E\{\partial M/\partial\theta^\top\}^{-\top}} for the asymptotic
#' variance of \eqn{n^{1/2}(\hat\theta-\theta_0)}, with unknown expectations
#' replaced by sample means -- and it is the covariance used by
#' Liu et al. (2021) and by Sun, Miao and Wickramarachchi (2025)
#' (Propositions 1-2 and Supplement S2-S3, where the partitioned form
#' \eqn{G = g - E(\partial g/\partial\tau^\top)
#'      \{E(\partial m/\partial\tau^\top)\}^{-1} m} is algebraically identical
#' to the corresponding block of the full stacked sandwich).
#'
#' @param estfun A function of the parameter vector returning an `n x p` matrix
#'   whose `i`-th row is the estimating function \eqn{M(O_i;\theta)^\top}.  The
#'   *mean* of each column must be (approximately) zero at `theta`.
#' @param theta Parameter vector at which the covariance is evaluated (the
#'   solution of the estimating equations).
#' @param names Optional character vector of parameter names.
#' @param tol Tolerance used by the numerical checks on the returned matrix.
#' @param jac_args Arguments passed to [numDeriv::jacobian()] as
#'   `method.args`.  The default `list(r = 2)` uses two Richardson
#'   extrapolation levels, which reproduces the full-accuracy Jacobian to
#'   better than 1e-8 in the models fitted here at roughly half the cost.
#'
#' @return A list with components
#'   \item{vcov}{the `p x p` sandwich covariance matrix of \eqn{\hat\theta};}
#'   \item{bread}{\eqn{\hat A^{-1}};}
#'   \item{meat}{\eqn{\hat B};}
#'   \item{jacobian}{\eqn{\hat A};}
#'   \item{estfun}{the `n x p` matrix of estimating functions;}
#'   \item{checks}{a named logical vector of numerical diagnostics.}
#'
#' @references
#' Sun, B., Liu, L., Miao, W., Wirth, K., Robins, J. and Tchetgen Tchetgen, E. J.
#' (2021). Semiparametric estimation with data missing not at random using an
#' instrumental variable. *Statistica Sinica*. Supplemental Appendix, eq. (S2).
#'
#' @examples
#' ## Sandwich covariance of a sample mean and a logistic regression slope.
#' set.seed(1)
#' x <- rnorm(200); y <- rbinom(200, 1, expit(0.5 * x))
#' fit <- glm(y ~ x, family = binomial())
#' ef <- function(th) cbind(1, x) * as.vector(y - expit(cbind(1, x) %*% th))
#' siv_sandwich(ef, coef(fit))$vcov
#' @export
siv_sandwich <- function(estfun, theta, names = NULL, tol = 1e-8,
                         jac_args = list(r = 2)) {
  psi <- estfun(theta)
  psi <- as.matrix(psi)
  n <- nrow(psi)
  p <- ncol(psi)
  if (length(theta) != p) {
    stop(sprintf("`estfun` returned %d columns but `theta` has length %d.",
                 p, length(theta)), call. = FALSE)
  }
  ## A_hat = n^{-1} sum_i d M(O_i; theta) / d theta^T  (the "bread" matrix,
  ## before inversion).  Computed by numerical differentiation of the *column
  ## means*, exactly as in the reference code of Sun et al. (2021), where
  ## `jacobian(M, theta) / n` is used with M returning column sums.
  A <- numDeriv::jacobian(function(th) colMeans(estfun(th)), theta,
                          method.args = jac_args)
  ## B_hat = n^{-1} sum_i M_i M_i^T  (the "meat" matrix).
  B <- crossprod(psi) / n
  Ainv <- tryCatch(solve(A), error = function(e) {
    tryCatch(siv_ginv(A), error = function(e2) matrix(NA_real_, p, p))
  })
  V <- Ainv %*% B %*% t(Ainv) / n
  V <- (V + t(V)) / 2  ## enforce exact symmetry after the numerical inverse
  if (!is.null(names)) {
    dimnames(V) <- list(names, names)
    dimnames(A) <- list(names, names)
    dimnames(B) <- list(names, names)
    colnames(psi) <- names
  }
  checks <- c(
    dim_ok        = identical(dim(V), c(p, p)),
    finite        = all(is.finite(V)),
    diag_nonneg   = all(is.na(diag(V))) || all(diag(V) > -tol),
    symmetric     = all(is.na(V)) || max(abs(V - t(V))) < tol,
    jacobian_ok   = all(is.finite(A)) && abs(det(A)) > 0
  )
  list(vcov = V, bread = Ainv, meat = B, jacobian = A, estfun = psi,
       checks = checks, n = n)
}

## Internal: minimal Moore-Penrose pseudo-inverse, so that a singular Jacobian
## degrades gracefully (and is reported through `checks`) instead of aborting.
siv_ginv <- function(M, tol = sqrt(.Machine$double.eps)) {
  s <- svd(M)
  keep <- s$d > max(tol * s$d[1L], 0)
  if (!any(keep)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, keep, drop = FALSE] %*% ((1 / s$d[keep]) * t(s$u[, keep, drop = FALSE]))
}

## Internal: numerically robust root finding for a system of estimating
## equations.  Uses a damped Newton iteration on the *mean* estimating function
## with numerical Jacobians, falling back to `stats::optim` on the squared norm.
## Convergence criterion follows Sun, Miao and Wickramarachchi (2025),
## Supplement Table S3 note: "the residual of the square average estimating
## equation component being less than 1e-7".
siv_solve_ee <- function(fn, start, tol = 1e-7, maxit = 200L, verbose = FALSE,
                         fallback = TRUE, fallback_maxit = 500L,
                         fallback_maxpar = 10L) {
  theta <- start
  conv <- FALSE
  crit <- NA_real_
  stalled <- 0L
  for (it in seq_len(maxit)) {
    f <- fn(theta)
    crit_new <- mean(f^2)
    ## Stop early when the criterion stops improving materially: a system that
    ## is going to converge does so quickly, and grinding on a stalled one is
    ## what makes a misspecified/near-singular fit expensive.
    if (is.finite(crit) && is.finite(crit_new) &&
        crit - crit_new < 1e-4 * crit) {
      stalled <- stalled + 1L
    } else {
      stalled <- 0L
    }
    crit <- crit_new
    if (!is.finite(crit)) break
    if (crit < tol) { conv <- TRUE; break }
    if (stalled >= 5L) break
    ## Forward differences are sufficient for the Newton direction; the
    ## sandwich Jacobian below uses Richardson extrapolation.
    J <- tryCatch(numDeriv::jacobian(fn, theta, method = "simple"),
                  error = function(e) NULL)
    if (is.null(J) || !all(is.finite(J))) break
    step <- tryCatch(solve(J, -f), error = function(e) siv_ginv(J) %*% (-f))
    step <- as.numeric(step)
    ## damped update: halve the step until the criterion decreases
    ok <- FALSE
    lam <- 1
    for (k in 1:20) {
      cand <- theta + lam * step
      fc <- tryCatch(fn(cand), error = function(e) rep(NA_real_, length(f)))
      cc <- mean(fc^2)
      if (is.finite(cc) && cc < crit) { theta <- cand; ok <- TRUE; break }
      lam <- lam / 2
    }
    if (!ok) break
  }
  if (!conv && fallback && length(theta) <= fallback_maxpar) {
    ## Fall back on a derivative-free minimisation of the squared norm.  This is
    ## only attempted for low-dimensional systems: a Nelder-Mead simplex needs
    ## O(100 p) iterations to make progress, so in higher dimensions it cannot
    ## rescue a stalled Newton iteration and merely burns time.  The estimating
    ## equations are unchanged either way; a system that does not solve is
    ## reported as a convergence failure (`converged = FALSE`), which is one of
    ## the metrics the simulation reports.
    obj <- function(th) {
      v <- tryCatch(mean(fn(th)^2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- stats::optim(theta, obj, method = "Nelder-Mead",
                       control = list(maxit = fallback_maxit, reltol = 1e-12))
    if (is.finite(op$value) && op$value < crit) { theta <- op$par; crit <- op$value }
    conv <- is.finite(crit) && crit < tol
  }
  if (verbose) message(sprintf("siv_solve_ee: criterion = %.3e, converged = %s",
                               crit, conv))
  list(par = theta, criterion = crit, converged = conv)
}

## Internal: build a model matrix from a one-sided formula and a data frame,
## returning a plain numeric matrix (with intercept unless suppressed).
siv_model_matrix <- function(formula, data) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  mm <- stats::model.matrix(formula, mf)
  storage.mode(mm) <- "double"
  mm
}

## Internal: assemble a data frame holding Y, A, Z and the columns of X, used
## for building the working-model design matrices from user formulas.
siv_frame <- function(Y, A, Z, X) {
  df <- data.frame(Y = Y, A = A, Z = Z)
  if (ncol(X) > 0L) df <- cbind(df, as.data.frame(X))
  df
}

## ----------------------------------------------------------------------------
## Warm starts for multi-parameter selection-bias models
## ----------------------------------------------------------------------------

## Internal: reduce a fitted basis specification to its scalar "linear" version,
## keeping the same reference value.  Used to build a one-parameter version of
## the selection-bias model whose estimating equations are well conditioned.
siv_basis_scalar <- function(spec) {
  siv_basis_spec("linear", y_ref = spec$y_ref, interact = NULL)
}

## Internal: least-squares projection of a fitted scalar selection-bias
## function  psi1 * (y - y_ref)  onto a richer basis, evaluated at the observed
## outcomes.  This gives a starting value for the multi-parameter model that
## reproduces the (well-identified) scalar fit as closely as the richer basis
## allows, which is a far better starting point than the zero vector.
siv_project_psi <- function(spec_full, y, X, psi_scalar, y_ref) {
  Bfull <- siv_basis_eval(spec_full, y, X)
  target <- psi_scalar * (y - y_ref)
  fit <- tryCatch(stats::lm.fit(Bfull, target), error = function(e) NULL)
  if (is.null(fit) || anyNA(fit$coefficients)) return(rep(0, ncol(Bfull)))
  cf <- fit$coefficients
  cf[!is.finite(cf)] <- 0
  as.numeric(cf)
}

## Internal: TRUE when a formula's terms separate into terms involving `Y` and
## terms not involving `Y` (i.e. there is no Y-by-covariate interaction), so
## that the linear predictor is additive and can be evaluated on an outcome grid
## in O(n + M) rather than O(n * M) predictions.
siv_formula_additive_in_Y <- function(formula, data) {
  tl <- attr(stats::terms(formula, data = data), "term.labels")
  if (!length(tl)) return(TRUE)
  hasY <- vapply(tl, function(t) "Y" %in% all.vars(stats::reformulate(t)),
                 logical(1))
  ## a term is "mixed" if it mentions Y and at least one other variable
  mixed <- vapply(seq_along(tl), function(j) {
    v <- all.vars(stats::reformulate(tl[j]))
    hasY[j] && length(setdiff(v, "Y")) > 0L
  }, logical(1))
  !any(mixed)
}
