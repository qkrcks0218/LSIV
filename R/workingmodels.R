## ----------------------------------------------------------------------------
## Finite-dimensional working models shared by the three competing estimators.
##
## All three competing papers replace one or more infinite-dimensional nuisance
## functions by a parametric working model and solve a system of estimating
## equations.  This file collects
##   * the outcome-density working model  f(y | .; eta)   (`siv_ydens_*`), and
##   * the exponential-tilting integrals
##       E[ exp{s * psi(Y)} v(Y) | . ]
##     that all three papers require (Liu et al. 2021 eq. (4.6);
##     Sun et al. 2021 eq. (4.5); Sun, Miao & Wickramarachchi 2025 eqs. (2.4)
##     and (3.11)).
##
## Two families are supported, matching the specifications actually used in the
## source papers:
##   "binomial"   -- a logistic model over a finite outcome support.  Used by
##                   Liu et al. (2021) eq. (5.2), Sun et al. (2021) (saturated
##                   logistic for P(Y = 1 | R = 1, X, Z)) and by Sun, Miao and
##                   Wickramarachchi (2025) Section 4 (p(Y = 1 | u; eta) =
##                   expit{h3(u) eta}).
##   "gaussian"   -- a homoscedastic normal model.  Used by Sun, Miao and
##                   Wickramarachchi (2025), Supplement Section S5.2:
##                   p(y | u; eta) = (2 pi sigma^2)^{-1/2}
##                                   exp[-(y - h1(u) eta)^2 / (2 sigma^2)],
##                   with eta = (sigma, eta^T)^T.
##
## For "binomial" the tilting integrals are exact finite sums over the support;
## for "gaussian" they are evaluated by Gauss-Hermite quadrature.
## ----------------------------------------------------------------------------

## Internal: number of free parameters of the outcome-density working model.
siv_ydens_npar <- function(family, D, support) {
  switch(family,
         binomial = ncol(D) * (length(support) - 1L),
         gaussian = ncol(D) + 1L,
         stop("unknown family"))
}

## Internal: starting values for the outcome-density working model, obtained by
## a (possibly weighted) fit on the subset for which Y is observed.
siv_ydens_start <- function(family, D, y, subset, support, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(y))
  switch(family,
    binomial = {
      out <- numeric(0)
      for (m in seq_len(length(support) - 1L)) {
        yy <- as.numeric(abs(y - support[m + 1L]) < 1e-8)
        fit <- suppressWarnings(stats::glm.fit(D[subset, , drop = FALSE], yy[subset],
                                               weights = weights[subset],
                                               family = stats::binomial()))
        cf <- fit$coefficients; cf[!is.finite(cf)] <- 0
        out <- c(out, cf)
      }
      out
    },
    gaussian = {
      fit <- stats::lm.wfit(D[subset, , drop = FALSE], y[subset], w = weights[subset])
      cf <- fit$coefficients; cf[!is.finite(cf)] <- 0
      s <- stats::sd(stats::residuals(fit))
      c(cf, max(s, 1e-3))
    })
}

## Internal: probability mass / density of the working model at the observed y.
## Returns a length-n vector.
siv_ydens_pdf <- function(family, par, D, y, support) {
  switch(family,
    binomial = {
      P <- siv_ydens_grid(family, par, D, support)
      idx <- vapply(y, function(v) which.min(abs(support - v)), integer(1L))
      P[cbind(seq_along(y), idx)]
    },
    gaussian = {
      k <- ncol(D)
      m <- as.vector(D %*% par[seq_len(k)])
      stats::dnorm(y, m, abs(par[k + 1L]))
    })
}

## Internal: probability mass function on the support (binomial only), n x M.
siv_ydens_grid <- function(family, par, D, support) {
  if (family != "binomial") stop("`siv_ydens_grid()` requires family = \"binomial\".")
  M <- length(support)
  k <- ncol(D)
  eta <- matrix(0, nrow(D), M)          # first support point is the reference
  for (m in seq_len(M - 1L)) {
    eta[, m + 1L] <- as.vector(D %*% par[(m - 1L) * k + seq_len(k)])
  }
  ## Multinomial-logit style renormalisation; for M = 2 this is exactly the
  ## logistic model logit P(Y = support[2]) = D %*% par.
  e <- exp(eta - apply(eta, 1L, max))
  e / rowSums(e)
}

## Internal: per-observation score of the outcome-density working model,
## restricted to `subset` (an n-vector of 0/1 or logical weights).  Returns an
## n x npar matrix whose column means vanish at the restricted MLE.
siv_ydens_score <- function(family, par, D, y, support, subset_w) {
  n <- nrow(D)
  switch(family,
    binomial = {
      P <- siv_ydens_grid(family, par, D, support)
      k <- ncol(D)
      out <- matrix(0, n, k * (length(support) - 1L))
      for (m in seq_len(length(support) - 1L)) {
        ind <- as.numeric(abs(y - support[m + 1L]) < 1e-8)
        ind[is.na(ind)] <- 0
        out[, (m - 1L) * k + seq_len(k)] <- D * (subset_w * (ind - P[, m + 1L]))
      }
      out
    },
    gaussian = {
      k <- ncol(D)
      s <- abs(par[k + 1L])
      m <- as.vector(D %*% par[seq_len(k)])
      r <- y - m; r[is.na(r)] <- 0
      cbind(D * (subset_w * r / s^2),
            subset_w * (r^2 - s^2) / s^3)
    })
}

## Internal: exponential-tilting moments.
##
##   num_j = E[ exp(sgn * b(Y, .) %*% psi) * v_j(Y) | . ]  and
##   den   = E[ exp(sgn * b(Y, .) %*% psi)            | . ]
##
## under the working density f(y | .; par).  `vfun` is a function of a numeric
## vector of outcome values and an index vector of units, returning a matrix
## with one column per requested moment.
##
## `sgn = +1` gives Liu et al. (2021) eq. (4.6)
##   E{g(Y0,C)|A=1,Z,C} = E[exp(psi) g | A=0,Z,C] / E[exp(psi) | A=0,Z,C];
## `sgn = -1` gives Sun et al. (2021) eq. (4.5) / Sun, Miao and Wickramarachchi
## (2025) eq. (2.4)
##   p(y|R=0,x) = exp(-h) p(y|R=1,x) / E[exp(-h) | R=1,x].
siv_tilt_moments <- function(family, par, D, support, basis_spec, X, psi,
                             vfun, sgn = 1, n_nodes = 20L) {
  n <- nrow(D)
  if (family == "binomial") {
    P <- siv_ydens_grid(family, par, D, support)        # n x M
    M <- length(support)
    den <- numeric(n)
    numl <- NULL
    for (m in seq_len(M)) {
      b <- siv_basis_eval(basis_spec, rep(support[m], n), X)   # n x q
      w <- exp(sgn * as.vector(b %*% psi)) * P[, m]
      den <- den + w
      vm <- vfun(rep(support[m], n))
      vm <- as.matrix(vm)
      if (is.null(numl)) numl <- matrix(0, n, ncol(vm))
      numl <- numl + vm * w
    }
    return(list(num = numl, den = den))
  }
  ## gaussian: Gauss-Hermite quadrature
  k <- ncol(D)
  s <- abs(par[k + 1L])
  mu <- as.vector(D %*% par[seq_len(k)])
  gh <- siv_gh_normal(mu, s, n_nodes)
  den <- numeric(n)
  numl <- NULL
  for (kk in seq_along(gh$weights)) {
    yk <- gh$nodes[, kk]
    b <- siv_basis_eval(basis_spec, yk, X)
    w <- exp(sgn * as.vector(b %*% psi)) * gh$weights[kk]
    den <- den + w
    vk <- as.matrix(vfun(yk))
    if (is.null(numl)) numl <- matrix(0, n, ncol(vk))
    numl <- numl + vk * w
  }
  list(num = numl, den = den)
}

## Internal: split the non-intercept columns of a baseline design matrix into
## those that involve the instrument (h1, in the notation of Liu et al. 2021,
## eqs. 4.3 and 4.5) and those that do not (h2, eq. 4.4).
siv_split_ZC <- function(formula, data) {
  mm <- siv_model_matrix(formula, data)
  asg <- attr(mm, "assign")
  tl <- attr(stats::terms(formula, data = data), "term.labels")
  has_int <- "(Intercept)" %in% colnames(mm)
  involves_Z <- logical(ncol(mm))
  for (j in seq_len(ncol(mm))) {
    if (asg[j] == 0L) { involves_Z[j] <- NA; next }   # intercept
    involves_Z[j] <- grepl("(^|:)Z($|:)", tl[asg[j]])
  }
  list(mm = mm, intercept = has_int,
       idx_h1 = which(!is.na(involves_Z) & involves_Z),
       idx_h2 = which(!is.na(involves_Z) & !involves_Z),
       idx_int = which(is.na(involves_Z)))
}
