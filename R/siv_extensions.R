## ----------------------------------------------------------------------------
## Distributional estimands and falsification diagnostics for SIV().
##
## The QTT implementation follows Supplementary Material A.6.3: a
## Berger-Boos confidence set combines an order-statistic confidence interval for the
## observed treated quantile with EIF tests of candidate counterfactual
## quantiles.  The falsification implementation follows Section 6 and eqs.
## (30)-(31). The unrestricted cross-fitted nuisance estimates, before the
## density-ratio projection used for SIV estimation, enter the uncentered
## influence-function score and a centered Gaussian
## multiplier bootstrap supplies the critical value.
## ----------------------------------------------------------------------------

## Internal: EIF numerators for one or more counterfactual functionals
## E{G(Y(0)) | A = 1}.  `G_grid` is M x J and `G_obs` is n x J.  This is the
## general, uniformly bounded-G version of the calculation in siv_eif(); see
## Supplementary Material C.6--C.7.  Keeping this helper separate makes the
## QTT CDF calculation vectorized over all candidate thresholds.
siv_cf_functionals <- function(Y, A, Z, grid, weights, delta,
                               fY_z0, fY_z1, pZ1,
                               pA1_z0, pA1_z1, G_grid = NULL, G_obs = NULL) {
  n <- length(Y)
  M <- length(grid)
  W <- matrix(weights, n, M, byrow = TRUE)
  pZ1 <- siv_bound_prob(pZ1, 1e-8)
  pZ0 <- 1 - pZ1
  pA1_z0 <- siv_bound_prob(pA1_z0, 1e-8)
  pA1_z1 <- siv_bound_prob(pA1_z1, 1e-8)

  ## Quantities not depending on G.
  base_z0 <- delta * fY_z0 * W
  base_z1 <- delta * fY_z1 * W
  denom_z0 <- rowSums(base_z0)
  denom_z1 <- rowSums(base_z1)
  mu_alpha_z0 <- rowSums(delta * base_z0) / denom_z0
  mu_alpha_z1 <- rowSums(delta * base_z1) / denom_z1
  beta_z0 <- (pA1_z0 / (1 - pA1_z0)) / denom_z0
  beta_z1 <- (pA1_z1 / (1 - pA1_z1)) / denom_z1

  wgt_z0 <- delta * beta_z0
  wgt_z1 <- delta * beta_z1

  ## The observed-outcome lookup is needed only when A = 0.  A finite
  ## placeholder is used elsewhere so that QTT can issue its own, clearer error
  ## when treated outcomes are missing.
  jj <- rep.int(1L, n)
  idx0 <- which(A == 0)
  jj[idx0] <- vapply(Y[idx0], function(y) which.min(abs(grid - y)), integer(1L))
  didx <- cbind(seq_len(n), jj)
  wgt_obs_z0 <- wgt_z0[didx]
  wgt_obs_z1 <- wgt_z1[didx]

  ## Reconstruct f(Y(0) | X), exactly as in siv_eif().
  fya0z0 <- fY_z0 * ((1 - pA1_z0) * pZ0)
  fya0z1 <- fY_z1 * ((1 - pA1_z1) * pZ1)
  fya1z0 <- wgt_z0 * fya0z0
  fya1z1 <- wgt_z1 * fya0z1
  fY <- ((fya0z0 + fya1z0) / pZ0 + (fya0z1 + fya1z1) / pZ1) / 2

  ## b(1,y,X), followed by the coefficient functions in L.
  b_z1 <- (1 + wgt_z0) * (beta_z1 * pZ1) /
    (beta_z1 * pZ1 + beta_z0 * pZ0 + delta * beta_z1 * beta_z0)
  coef1 <- 1 / (b_z1 - pZ1)
  coef2 <- b_z1 / (b_z1 - pZ1)

  A11 <- beta_z1 * pZ1 + beta_z0 * beta_z1 * pZ1 * mu_alpha_z0
  A12 <- beta_z1 * pZ1 + 2 * beta_z0 * pZ0 - beta_z1 * pZ0 +
    beta_z0 * beta_z1 * mu_alpha_z0
  A21 <- beta_z0 * pZ1 + beta_z0 * beta_z1 * pZ1 * mu_alpha_z1
  A22 <- beta_z0 + beta_z0 * beta_z1 * mu_alpha_z1
  det <- A11 * A22 - A12 * A21

  int_c1 <- rowSums(coef1 * fY * W)
  int_c2 <- rowSums(coef2 * fY * W)
  c1_obs <- coef1[didx]
  c2_obs <- coef2[didx]
  dc1_z0 <- delta * coef1 * fY_z0 * W
  dc2_z0 <- delta * coef2 * fY_z0 * W
  dc1_z1 <- delta * coef1 * fY_z1 * W
  dc2_z1 <- delta * coef2 * fY_z1 * W

  evaluate <- function(G_grid, G_obs) {
    G_grid <- as.matrix(G_grid)
    G_obs <- as.matrix(G_obs)
    if (nrow(G_grid) != M) {
      stop("`G_grid` must have `length(grid)` rows.", call. = FALSE)
    }
    if (nrow(G_obs) != n || ncol(G_obs) != ncol(G_grid)) {
      stop("`G_obs` must have `length(Y)` rows and the same number of columns ",
           "as `G_grid`.", call. = FALSE)
    }

    G_lookup <- G_grid[jj, , drop = FALSE]
    ## G-dependent conditional moments.  Matrix multiplication evaluates all J
    ## functionals in one pass over the outcome grid.
    mu_z0 <- (base_z0 %*% G_grid) / denom_z0
    mu_z1 <- (base_z1 %*% G_grid) / denom_z1
    mu_alphaG_z0 <- ((delta * base_z0) %*% G_grid) / denom_z0
    mu_alphaG_z1 <- ((delta * base_z1) %*% G_grid) / denom_z1

    b1 <- beta_z1 * pZ1 * mu_z0 +
      2 * beta_z0 * pZ0 * mu_z0 - beta_z1 * pZ0 * mu_z0 +
      beta_z0 * beta_z1 * mu_alphaG_z0
    b2 <- beta_z0 * mu_z1 + beta_z0 * beta_z1 * mu_alphaG_z1
    Bx <- (-b1 * A22 + b2 * A12) / det
    Cx <- (-b2 * A11 + b1 * A21) / det

    ## L and E(L | X).
    mu_L <- (coef1 * fY * W) %*% G_grid + Bx * int_c2 + Cx * int_c1
    ## As in siv_eif(), L(Y_i, X_i) is evaluated at the closest outcome-grid
    ## point; G_obs itself is retained in the leading residual term below.
    L_obs <- G_lookup * c1_obs + Bx * c2_obs + Cx * c1_obs
    omega_obs <- (L_obs - mu_L) * (Z - pZ1)

    ## E{omega | A=1,z,X}.  Since omega=(L-E(L|X))(z-e), only weighted
    ## integrals of the three terms in L are needed.
    EdL_z0 <- (dc1_z0 %*% G_grid + Bx * rowSums(dc2_z0) +
                 Cx * rowSums(dc1_z0)) / denom_z0
    EdL_z1 <- (dc1_z1 %*% G_grid + Bx * rowSums(dc2_z1) +
                 Cx * rowSums(dc1_z1)) / denom_z1
    mu_w_z0 <- (0 - pZ1) * (EdL_z0 - mu_L)
    mu_w_z1 <- (1 - pZ1) * (EdL_z1 - mu_L)

    mu <- mu_z0
    mu_w <- mu_w_z0
    if (any(Z == 1)) {
      mu[Z == 1, ] <- mu_z1[Z == 1, , drop = FALSE]
      mu_w[Z == 1, ] <- mu_w_z1[Z == 1, , drop = FALSE]
    }
    wgt <- ifelse(Z == 0, wgt_obs_z0, wgt_obs_z1)

    ## Uncentered numerator for E{G(Y(0)) | A=1}; this is the bracketed
    ## counterfactual part of eq. (21), and is the quantity used in A.6.3.
    num <- (1 - A) * wgt * (G_obs - mu) + A * mu -
      (1 - A) * wgt * (omega_obs - mu_w) - A * mu_w -
      (1 - A) * omega_obs
    pA <- mean(A)
    parameter <- colMeans(num) / pA
    eif <- (num - A %o% parameter) / pA

    list(num = num, parameter = parameter, eif = eif,
         var = apply(eif, 2L, stats::var) / n,
         min_relevance = min(abs(b_z1 - pZ1)),
         beta_z0 = beta_z0, beta_z1 = beta_z1, fY = fY,
         B = Bx, C = Cx, mu_z0 = mu_z0, mu_z1 = mu_z1,
         omega_obs = omega_obs)
  }
  if (is.null(G_grid) && is.null(G_obs)) return(evaluate)
  evaluate(G_grid, G_obs)
}

## Internal: generalized inverse of a numerically estimated CDF.
siv_inverse_cdf <- function(threshold, cdf, q, discrete = FALSE) {
  cdf <- cummax(pmin(pmax(as.numeric(cdf), 0), 1))
  if (discrete) {
    hit <- which(cdf >= q)
    return(threshold[if (length(hit)) hit[1L] else length(threshold)])
  }
  keep <- !duplicated(cdf)
  if (sum(keep) < 2L) return(threshold[which.min(abs(cdf - q))])
  as.numeric(stats::approx(cdf[keep], threshold[keep], xout = q,
                           rule = 2, ties = "ordered")$y)
}

## Internal: distribution-free treated-quantile intervals from binomial ranks.
## Exterior ranks use declared outcome bounds, or infinities when unknown.
siv_quantile_ci <- function(y, quantiles, error, bounds = NULL) {
  y <- sort(y)
  n <- length(y)
  lower_rank <- stats::qbinom(error / 2, n, quantiles)
  upper_rank <- stats::qbinom(1 - error / 2, n, quantiles) + 1
  lower <- rep(if (is.null(bounds)) -Inf else bounds[1L], length(quantiles))
  upper <- rep(if (is.null(bounds)) Inf else bounds[2L], length(quantiles))
  lower[lower_rank > 0] <- y[lower_rank[lower_rank > 0]]
  upper[upper_rank <= n] <- y[upper_rank[upper_rank <= n]]
  cbind(lower = lower, upper = upper)
}

## Internal: boundary-normal tests used in quantile inversion.  At a zero-
## variance boundary, a satisfied inequality is never rejected.
siv_qtt_pvalue <- function(num, A, q, alternative = "two.sided") {
  n <- length(A)
  pA <- mean(A)
  rho <- colMeans(num) / pA - q
  score <- sweep((num - q * A) / pA, 2L, rho, "-")
  se <- sqrt(colMeans(score^2) / n)
  if (any(!is.finite(rho)) || any(!is.finite(se))) {
    stop("Non-finite counterfactual quantile scores; check nuisance estimates ",
         "and instrument relevance.", call. = FALSE)
  }
  positive <- se > 0
  p <- rep(1, length(rho))
  if (alternative == "two.sided") {
    p[positive] <- 2 * stats::pnorm(-abs(rho[positive] / se[positive]))
    p[!positive & abs(rho) > sqrt(.Machine$double.eps)] <- 0
  } else if (alternative == "greater") {
    p[positive] <- stats::pnorm(rho[positive] / se[positive], lower.tail = FALSE)
    p[!positive & rho > sqrt(.Machine$double.eps)] <- 0
  } else if (alternative == "less") {
    p[positive] <- stats::pnorm(rho[positive] / se[positive])
    p[!positive & rho < -sqrt(.Machine$double.eps)] <- 0
  } else {
    stop("Unknown quantile-test alternative.", call. = FALSE)
  }
  list(p = p, se = se)
}

## Internal: project the complete accepted set through the treated-quantile CI.
## For continuous outcomes, threshold_upper supplies each accepted cell's right
## edge. An empty inversion uses declared counterfactual support or outcome
## bounds, falling back to all of R only when neither is known.
siv_qtt_ci <- function(accepted, threshold, treated_ci, Y = NULL,
                       threshold_upper = threshold, support = NULL,
                       bounds = NULL) {
  cf_domain <- if (!is.null(support)) range(support) else bounds
  if (!is.null(bounds)) {
    treated_ci <- c(lower = max(treated_ci["lower"], bounds[1L]),
                    upper = min(treated_ci["upper"], bounds[2L]))
  }
  if (!length(accepted)) {
    if (is.null(cf_domain)) return(c(lower = -Inf, upper = Inf))
    cf_limits <- cf_domain
  } else {
    cf_limits <- c(min(threshold[accepted]), max(threshold_upper[accepted]))
    if (!is.null(cf_domain)) {
      cf_limits <- c(max(cf_limits[1L], cf_domain[1L]),
                     min(cf_limits[2L], cf_domain[2L]))
      if (cf_limits[1L] > cf_limits[2L]) cf_limits <- cf_domain
    }
  }
  c(lower = unname(treated_ci["lower"]) - cf_limits[2L],
    upper = unname(treated_ci["upper"]) - cf_limits[1L])
}

## Internal: retain infinite confidence limits when combining repetitions.
siv_qtt_endpoint_median <- function(x) {
  keep <- !is.na(x)
  if (any(keep)) stats::median(x[keep]) else NA_real_
}

## Internal: QTT estimates and Berger-Boos confidence intervals.  The nominal
## error probability is alpha = 1 - level; its default nuisance allocation is
## c1 = 0.02 * alpha.  num_bootstrap is retained only for API compatibility.
siv_qtt <- function(Y, A, Z, type, grid, weights, fp, nuisance,
                    quantiles = c(0.25, 0.5, 0.75), level = 0.95,
                    c1 = 0.02 * (1 - level), num_bootstrap = 1000L,
                    support = NULL, bounds = NULL) {
  type <- match.arg(type, c("binary", "discrete", "continuous"))
  if (is.null(bounds) && type == "binary") bounds <- c(0, 1)
  if (!is.null(bounds)) {
    if (!is.numeric(bounds) || length(bounds) != 2L ||
        any(!is.finite(bounds)) || bounds[1L] > bounds[2L]) {
      stop("`QTT_bounds` must be a finite ordered numeric pair.", call. = FALSE)
    }
    bounds <- as.numeric(bounds)
    if (type == "binary") bounds <- c(max(0, bounds[1L]), min(1, bounds[2L]))
    if (bounds[1L] > bounds[2L] || any(!is.finite(Y)) ||
        any(Y < bounds[1L] | Y > bounds[2L])) {
      stop("`QTT_bounds` must contain every observed outcome.", call. = FALSE)
    }
  }
  quantiles <- sort(unique(as.numeric(quantiles)))
  if (!length(quantiles) || any(!is.finite(quantiles)) ||
      any(quantiles <= 0 | quantiles >= 1)) {
    stop("`QTT_quantile` must contain probabilities strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must lie strictly between 0 and 1.", call. = FALSE)
  }
  alpha <- 1 - level
  if (!is.numeric(c1) || length(c1) != 1L || !is.finite(c1) ||
      c1 <= 0 || c1 >= alpha) {
    stop("`QTT_c1` must lie strictly between 0 and `1 - level`.", call. = FALSE)
  }
  if (!is.numeric(num_bootstrap) || length(num_bootstrap) != 1L ||
      !is.finite(num_bootstrap) || num_bootstrap < 1 ||
      num_bootstrap != floor(num_bootstrap)) {
    stop("`QTT_numBB` must be a positive integer.", call. = FALSE)
  }
  y1 <- Y[A == 1]
  if (length(y1) < 2L || any(!is.finite(y1))) {
    stop("QTT estimation requires at least two finite outcomes among units ",
         "with `A == 1`.", call. = FALSE)
  }
  if (is.null(support) && type != "continuous") {
    support <- if (type == "binary") c(0, 1) else grid
  }
  discrete <- !is.null(support)
  if (discrete) {
    if (!is.numeric(support) || !length(support) || any(!is.finite(support)) ||
        is.unsorted(support, strictly = TRUE)) {
      stop("`QTT_support` must be a finite, strictly increasing numeric vector.",
           call. = FALSE)
    }
    if (any(!Y[A == 0] %in% support)) {
      stop("`QTT_support` must contain every untreated outcome.", call. = FALSE)
    }
    if (!is.null(bounds) && any(support < bounds[1L] | support > bounds[2L])) {
      stop("`QTT_support` must lie within `QTT_bounds`.", call. = FALSE)
    }
    threshold <- as.numeric(support)
    threshold_upper <- threshold
  } else {
    ## Every indicator in the score is constant on these half-open cells.
    ## Include untreated observations even when they are not integration nodes.
    breaks <- sort(unique(c(grid, Y[A == 0])))
    threshold <- c(-Inf, breaks)
    threshold_upper <- c(breaks, Inf)
    if (!is.null(bounds)) {
      threshold <- sort(unique(c(bounds, breaks[breaks >= bounds[1L] &
                                                 breaks <= bounds[2L]])))
      threshold_upper <- c(threshold[-1L], bounds[2L])
    }
  }
  q_label <- vapply(quantiles, function(q) {
    format(q, trim = TRUE, scientific = FALSE, digits = 15)
  }, character(1L))
  nm <- paste0("QTT_", q_label)
  J <- length(threshold)
  K <- length(quantiles)
  cdf_raw <- cdf_left <- numeric(J)
  pvalue <- se_rho <- se_left <- matrix(NA_real_, J, K,
                                        dimnames = list(NULL, nm))

  ## Prepare the G-independent EIF calculations once, then evaluate thresholds
  ## in bounded-size batches rather than allocating n by J score matrices.
  evaluate <- siv_cf_functionals(
    Y, A, Z, grid, weights, fp$delta, fp$fY_z0, fp$fY_z1,
    nuisance$pZ1, nuisance$pA1_z0, nuisance$pA1_z1
  )
  for (first in seq.int(1L, J, by = 64L)) {
    index <- seq.int(first, min(first + 63L, J))
    t <- threshold[index]
    right <- evaluate(outer(grid, t, "<=") * 1, outer(Y, t, "<=") * 1)
    cdf_raw[index] <- right$parameter
    if (discrete) {
      left <- evaluate(outer(grid, t, "<") * 1, outer(Y, t, "<") * 1)
      cdf_left[index] <- left$parameter
    }
    for (j in seq_along(quantiles)) {
      if (discrete) {
        ## H0 requires F(t-) <= q and F(t) >= q.  Bonferroni combines
        ## the two one-sided boundary tests into one valid quantile test.
        lp <- siv_qtt_pvalue(left$num, A, quantiles[j], "greater")
        rp <- siv_qtt_pvalue(right$num, A, quantiles[j], "less")
        raw_p <- pmin(1, 2 * pmin(lp$p, rp$p))
        se_left[index, j] <- lp$se
      } else {
        rp <- siv_qtt_pvalue(right$num, A, quantiles[j])
        raw_p <- rp$p
      }
      pvalue[index, j] <- pmin(1, raw_p + c1)
      se_rho[index, j] <- rp$se
    }
  }
  cdf_monotone <- cummax(pmin(pmax(cdf_raw, 0), 1))
  finite <- is.finite(threshold)
  tau0 <- vapply(quantiles, function(q) {
    siv_inverse_cdf(threshold[finite], cdf_monotone[finite], q,
                    discrete = discrete)
  }, numeric(1L))
  tau1 <- as.numeric(stats::quantile(y1, probs = quantiles, names = FALSE,
                                     type = if (discrete) 1 else 7))
  estimate <- stats::setNames(tau1 - tau0, nm)
  C1 <- siv_quantile_ci(y1, quantiles, c1, bounds = bounds)
  rownames(C1) <- nm
  ci <- matrix(NA_real_, K, 2L,
               dimnames = list(nm, c("lower", "upper")))
  empty_inversion <- stats::setNames(logical(K), nm)
  accepted_count <- stats::setNames(integer(K), nm)
  accepted_intervals <- stats::setNames(vector("list", K), nm)
  fallback_ci <- ci
  for (j in seq_along(quantiles)) {
    accepted <- which(pvalue[, j] > alpha)
    empty_inversion[j] <- !length(accepted)
    accepted_count[j] <- length(accepted)
    accepted_intervals[[j]] <- cbind(lower = threshold[accepted],
                                     upper = threshold_upper[accepted])
    ci[j, ] <- siv_qtt_ci(accepted, threshold, C1[j, ],
                          threshold_upper = threshold_upper, support = support,
                          bounds = bounds)
    fallback_ci[j, ] <- siv_qtt_ci(integer(), threshold, C1[j, ],
                                   support = support, bounds = bounds)
  }
  list(estimate = estimate, ci = ci, quantile = quantiles,
       treated_quantile = stats::setNames(tau1, nm),
       counterfactual_quantile = stats::setNames(tau0, nm),
       treated_quantile_ci = C1,
       treated_quantile_ci_method = "Binomial order-statistic confidence interval",
       c1 = c1, nominal_error = alpha, level = level,
       num_bootstrap = 0L, num_bootstrap_requested = num_bootstrap,
       threshold = threshold, threshold_upper = threshold_upper,
       support = support, bounds = bounds, discrete = discrete,
       cdf = cdf_raw, cdf_left = if (discrete) cdf_left else NULL,
       cdf_monotone = cdf_monotone,
       adjusted_pvalue = pvalue, rho_se = se_rho,
       rho_se_left = if (discrete) se_left else NULL,
       accepted_intervals = accepted_intervals,
       empty_inversion = empty_inversion, accepted_count = accepted_count,
       fallback_ci = fallback_ci,
       uninformative_ci = if (is.null(bounds)) c(lower = -Inf, upper = Inf) else
         c(lower = bounds[1L] - bounds[2L], upper = bounds[2L] - bounds[1L]),
       method = "Berger-Boos inversion with order-statistic treated-quantile CI")
}

## Internal: construct either exact categorical cells or disjoint empirical-
## quantile bins.  The requested count determines only the quantile grid; it is not
## used to guess whether a variable is categorical.
siv_make_empirical_bins <- function(x, n_bins, name, categorical = FALSE) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop("Falsification binning requires finite values of `", name, "`.",
         call. = FALSE)
  }
  n_bins_raw <- as.numeric(n_bins)[1L]
  if (!is.finite(n_bins_raw) || n_bins_raw < 1 ||
      abs(n_bins_raw - round(n_bins_raw)) > sqrt(.Machine$double.eps)) {
    stop("Requested bin counts must be positive integers.", call. = FALSE)
  }
  n_bins <- as.integer(round(n_bins_raw))
  lev <- sort(unique(x))
  if (isTRUE(categorical) || length(lev) == 1L) {
    return(list(kind = "levels", levels = lev, breaks = NULL,
                n_bins = length(lev), name = name))
  }
  probs <- seq(0, 1, length.out = n_bins + 1L)
  cutoffs <- unique(as.numeric(stats::quantile(
    x, probs = probs, names = FALSE, type = 7
  )))
  ## Replacing the empirical endpoints by infinities gives one interval for
  ## each adjacent pair of distinct quantile cutoffs and avoids empty exterior
  ## bins.  Repeated quantiles are necessarily collapsed.
  if (length(cutoffs) < 2L) {
    return(list(kind = "levels", levels = lev, breaks = NULL,
                n_bins = length(lev), name = name))
  }
  breaks <- cutoffs
  breaks[1L] <- -Inf
  breaks[length(breaks)] <- Inf
  list(kind = "quantile", levels = NULL, breaks = breaks,
       n_bins = length(breaks) - 1L, name = name)
}

siv_apply_empirical_bins <- function(spec, x) {
  if (spec$kind == "levels") return(match(as.numeric(x), spec$levels))
  findInterval(as.numeric(x), spec$breaks, all.inside = TRUE,
               rightmost.closed = TRUE)
}

siv_bin_labels <- function(spec) {
  if (spec$kind == "levels") return(paste0(spec$name, "=", spec$levels))
  lo <- spec$breaks[-length(spec$breaks)]
  hi <- spec$breaks[-1L]
  paste0(spec$name, " in ", ifelse(is.infinite(lo), "(-Inf", paste0("[", lo)),
         ", ", ifelse(is.infinite(hi), "Inf)", paste0(hi, "]")))
}

## Internal: influence-function falsification test of eqs. (30)--(31).
siv_falsification <- function(Y, A, Z, X, grid, weights, nuisance,
                              num_bootstrap = 1000L,
                              num_bins = NULL, variance_floor = NULL,
                              categorical = NULL, level = 0.95,
                              bootstrap_chunk = 100L) {
  n <- length(Y)
  X <- as.matrix(X)
  dx <- ncol(X)
  if (is.null(num_bins)) num_bins <- c(20L, rep.int(5L, dx))
  num_bins <- as.numeric(num_bins)
  if (length(num_bins) != dx + 1L || any(!is.finite(num_bins)) ||
      any(num_bins < 1) ||
      any(abs(num_bins - round(num_bins)) > sqrt(.Machine$double.eps))) {
    stop("`falsification_num.bin` must contain `1 + ncol(X)` positive integers. ",
         "one bin count for Y followed by one for each X column.",
         call. = FALSE)
  }
  num_bins <- as.integer(round(num_bins))
  names(num_bins) <- c("Y", colnames(X))
  num_bootstrap <- as.integer(num_bootstrap)[1L]
  if (!is.finite(num_bootstrap) || num_bootstrap < 1L) {
    stop("`falsification_num.boot` must be a positive integer.", call. = FALSE)
  }
  floor_was_null <- is.null(variance_floor)
  if (floor_was_null) variance_floor <- 0.5^(0.5 + dx / 2)
  variance_floor <- as.numeric(variance_floor)[1L]
  if (!is.finite(variance_floor) || variance_floor <= 0) {
    stop("`falsification_floor` must be `NULL` or one positive finite number.",
         call. = FALSE)
  }

  if (is.null(categorical)) {
    categorical <- integer(0)
  } else if (is.character(categorical)) {
    bad <- setdiff(categorical, colnames(X))
    if (length(bad)) {
      stop("Unknown `falsification_categorical` column(s). ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
    categorical <- match(categorical, colnames(X))
  } else if (is.numeric(categorical)) {
    categorical <- as.integer(categorical)
    if (anyNA(categorical) || any(categorical < 1L | categorical > dx)) {
      stop("Numeric `falsification_categorical` entries must index columns ",
           "of `X`.", call. = FALSE)
    }
  } else {
    stop("`falsification_categorical` must be `NULL`, column names, or ",
         "column indices.", call. = FALSE)
  }
  binary_x <- if (dx) vapply(seq_len(dx), function(j) {
    values <- sort(unique(X[, j]))
    length(values) <= 2L && all(values %in% c(0, 1))
  }, logical(1L)) else logical(0)
  categorical <- sort(unique(c(categorical, which(binary_x))))

  y_for_bins <- Y[is.finite(Y)]
  y_is_binary <- length(unique(y_for_bins)) <= 2L &&
    all(unique(y_for_bins) %in% c(0, 1))
  yspec <- siv_make_empirical_bins(y_for_bins, num_bins[1L], "Y",
                                   categorical = y_is_binary)
  y_grid_code <- siv_apply_empirical_bins(yspec, grid)
  y_obs_code <- rep.int(NA_integer_, n)
  y_obs_code[is.finite(Y)] <- siv_apply_empirical_bins(yspec, Y[is.finite(Y)])
  By <- yspec$n_bins
  Ymem <- outer(y_grid_code, seq_len(By), "==") * 1

  if (dx) {
    xspec <- lapply(seq_len(dx), function(j) {
      siv_make_empirical_bins(X[, j], num_bins[j + 1L], colnames(X)[j],
                              categorical = j %in% categorical)
    })
    xcode <- vapply(seq_len(dx), function(j) {
      siv_apply_empirical_bins(xspec[[j]], X[, j])
    }, integer(n))
    xcode <- matrix(xcode, nrow = n, ncol = dx,
                    dimnames = list(NULL, colnames(X)))
  } else {
    xspec <- list()
    xcode <- matrix(integer(0), n, 0L)
  }

  ## The dictionary in eq. (31) always uses marginal X cells.  It is the union
  ## over covariates of the cells formed by one covariate at a time.
  if (!dx) {
    cells <- list(list(index = seq_len(n), variable = NA_character_,
                       definition = "(no covariates)"))
  } else {
    cells <- unlist(lapply(seq_len(dx), function(j) {
      labels <- siv_bin_labels(xspec[[j]])
      lapply(seq_len(xspec[[j]]$n_bins), function(b) {
        list(index = which(xcode[, j] == b), variable = colnames(X)[j],
             definition = labels[b])
      })
    }), recursive = FALSE)
    ## An empirical bin can still be empty when a covariate is categorical with
    ## levels that do not all occur; such a cell carries no moment.
    cells <- cells[vapply(cells, function(z) length(z$index) > 0L, logical(1L))]
  }
  ncell <- length(cells)

  e <- siv_bound_prob(nuisance$pZ1, 1e-8)
  q0 <- 1 - siv_bound_prob(nuisance$pA1_z0, 1e-8)
  q1 <- 1 - siv_bound_prob(nuisance$pA1_z1, 1e-8)
  d <- q1 - q0
  d_floor <- sqrt(.Machine$double.eps)
  if (any(abs(d) < d_floor)) {
    warning("The estimated conditional first-stage difference is numerically ",
            "zero for some observations; it was truncated for the ",
            "falsification score.", call. = FALSE)
    d <- ifelse(d < 0, -1, 1) * pmax(abs(d), d_floor)
  }
  qZ <- ifelse(Z == 1, q1, q0)
  hz <- Z / e - (1 - Z) / (1 - e)
  W <- matrix(weights, n, length(grid), byrow = TRUE)

  fits <- vector("list", ncell)
  rows <- vector("list", ncell)
  sqrt_n <- sqrt(n)
  observed <- 0
  counter <- 0L
  for (cc in seq_len(ncell)) {
    ii <- cells[[cc]]$index
    pbin0 <- (nuisance$fY_z0[ii, , drop = FALSE] * W[ii, , drop = FALSE]) %*% Ymem
    pbin1 <- (nuisance$fY_z1[ii, , drop = FALSE] * W[ii, , drop = FALSE]) %*% Ymem
    mu0 <- q0[ii] * pbin0
    mu1 <- q1[ii] * pbin1
    ng <- mu1 - mu0
    muZ <- mu0
    if (any(Z[ii] == 1)) muZ[Z[ii] == 1, ] <- mu1[Z[ii] == 1, , drop = FALSE]

    gobs <- matrix(0, length(ii), By)
    ok <- which(A[ii] == 0 & !is.na(y_obs_code[ii]))
    if (length(ok)) gobs[cbind(ok, y_obs_code[ii][ok])] <- 1
    phi <- ng / d[ii] + hz[ii] *
      (((1 - A[ii]) * gobs - muZ) / d[ii] -
         ng * ((1 - A[ii]) - qZ[ii]) / d[ii]^2)
    m <- colSums(phi) / n
    sigma <- sqrt(pmax(colSums(phi^2) / n - m^2, 0))
    denom <- pmax(variance_floor, sigma)
    component <- -sqrt_n * m / denom
    observed <- max(observed, component)
    fits[[cc]] <- list(index = ii, phi = phi, m = m,
                       sigma = sigma, denom = denom, component = component)

    rows[[cc]] <- data.frame(
      index = counter + seq_len(By), x_cell = cc,
      y_bin = seq_len(By), n_x_cell = length(ii),
      y_definition = siv_bin_labels(yspec),
      x_variable = rep.int(cells[[cc]]$variable, By),
      x_definition = rep.int(cells[[cc]]$definition, By),
      moment = m, sigma = sigma, denominator = denom,
      statistic = component, stringsAsFactors = FALSE
    )
    counter <- counter + By
  }
  moments <- do.call(rbind, rows)

  ## Centered Gaussian multiplier bootstrap, evaluated in chunks to keep the
  ## memory requirement O(N * chunk + N * number_of_Y_bins).
  boot <- numeric(num_bootstrap)
  chunk <- max(1L, as.integer(bootstrap_chunk))
  starts <- seq.int(1L, num_bootstrap, by = chunk)
  for (first in starts) {
    last <- min(num_bootstrap, first + chunk - 1L)
    bb <- last - first + 1L
    mult <- matrix(stats::rnorm(n * bb), n, bb)
    sum_mult <- colSums(mult)
    mx <- rep.int(0, bb)
    for (ff in fits) {
      process <- crossprod(ff$phi, mult[ff$index, , drop = FALSE]) -
        ff$m %o% sum_mult
      process <- -process / sqrt_n
      process <- process / ff$denom
      mx <- pmax(mx, apply(process, 2L, max))
    }
    boot[first:last] <- mx
  }
  alpha <- 1 - level
  critical <- as.numeric(stats::quantile(boot, 1 - alpha, names = FALSE,
                                         type = 7))
  pvalue <- (1 + sum(boot >= observed)) / (num_bootstrap + 1)

  list(statistic = observed, critical_value = critical,
       p_value = pvalue, reject = observed > critical,
       level = level, floor = variance_floor,
       floor_default = floor_was_null,
       num_bootstrap = num_bootstrap, bootstrap_statistics = boot,
       num_bins_requested = num_bins,
       num_bins_realized = c(Y = By, stats::setNames(
         vapply(xspec, `[[`, integer(1L), "n_bins"), colnames(X)
       )),
       x_cells = "marginal",
       n_y_bins = By, n_x_cells_observed = ncell,
       n_dictionary = nrow(moments), moments = moments,
       y_bins = yspec, x_bins = xspec,
       categorical_x = colnames(X)[categorical],
       method = "cross-fitted IF test with centered Gaussian multiplier bootstrap")
}
