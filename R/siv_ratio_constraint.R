## Internal: prepare SIV estimation inputs without modifying the unrestricted
## nuisance estimates used by the falsification test. All constraints below
## hold on the supplied outcome grid, with its integration weights.
##
## In the zM/zm ordering of eq. (15), the falsification implication is
##   R_Y <= Pr(A=0 | zm, X) / Pr(A=0 | zM, X).
## Project the ratio onto this bound and its numerical bounds, subject to
##   sum(R_Y * fY_zm * weights) = 1.
## The weighted least-squares projection of the numerically guarded ratio is
##   clip(R_Y + lambda, lower, upper).
## A separate scalar lambda per unit is found by vectorized bisection. Never
## renormalize after clipping, since that can undo the falsification bound.
siv_constrain_density_ratio <- function(nuisance, weights,
                                        density_floor = 1e-4,
                                        dr_floor = 1e-4, tol = 1e-12) {
  f0 <- nuisance$fY_z0
  f1 <- nuisance$fY_z1
  if (!is.matrix(f0) || !is.matrix(f1) ||
      !identical(dim(f0), dim(f1)) || !nrow(f0) ||
      length(weights) != ncol(f0) || any(!is.finite(weights)) ||
      any(weights <= 0)) {
    stop("Density matrices and positive integration weights must have ",
         "matching dimensions.", call. = FALSE)
  }
  n <- nrow(f0)
  p0 <- nuisance$pA1_z0
  p1 <- nuisance$pA1_z1
  if (length(p0) != n || length(p1) != n ||
      any(!is.finite(c(p0, p1))) || any(c(p0, p1) < 0 | c(p0, p1) > 1)) {
    stop("Treatment probabilities must be finite values in [0, 1], ",
         "one per density row.", call. = FALSE)
  }
  if (length(density_floor) != 1L || !is.finite(density_floor) ||
      density_floor < 0 || length(dr_floor) != 1L ||
      !is.finite(dr_floor) || dr_floor < 0 || dr_floor >= 1) {
    stop("`density_floor` must be nonnegative and `dr_floor` must be ",
         "in [0, 1).", call. = FALSE)
  }
  if (any(!is.finite(f0)) || any(!is.finite(f1)) ||
      any(f0 < 0) || any(f1 < 0)) {
    stop("Conditional densities must be finite and nonnegative.", call. = FALSE)
  }
  ## Match the numerical probability guard used by the EIF, but retain the
  ## original probabilities in the unrestricted nuisance object.
  p0 <- siv_bound_prob(p0, 1e-8)
  p1 <- siv_bound_prob(p1, 1e-8)
  if (any(p0 == p1)) {
    stop("SIV estimation requires a nonzero conditional first-stage ",
         "difference; equal treatment probabilities cannot identify the ",
         "counterfactual distribution.", call. = FALSE)
  }
  f0 <- pmax(f0, density_floor)
  f1 <- pmax(f1, density_floor)
  mass0 <- as.vector(f0 %*% weights)
  mass1 <- as.vector(f1 %*% weights)
  if (any(mass0 <= 0) || any(mass1 <= 0) ||
      any(!is.finite(c(mass0, mass1)))) {
    stop("Each conditional density must have positive finite grid mass.",
         call. = FALSE)
  }
  f0 <- f0 / mass0
  f1 <- f1 / mass1
  if (any(f0 <= 0) || any(f1 <= 0)) {
    stop("SIV density ratios require positive densities on the grid; ",
         "use a positive `density_floor`.", call. = FALSE)
  }

  zM <- as.integer(p1 >= p0)
  reverse <- zM == 0L
  base <- f0
  base[reverse, ] <- f1[reverse, , drop = FALSE]
  if (is.null(nuisance$RY)) {
    ratio <- f1 / f0
  } else {
    ratio <- as.matrix(nuisance$RY)
    if (!identical(dim(ratio), dim(f0)) || anyNA(ratio)) {
      stop("`nuisance$RY` must match the density dimensions and contain ",
           "no missing values.", call. = FALSE)
    }
  }
  lower <- max(dr_floor, .Machine$double.eps)
  guard_upper <- if (dr_floor > 0) 1 / dr_floor else Inf
  ratio <- pmin(pmax(ratio, lower), guard_upper)
  if (any(!is.finite(ratio))) {
    stop("The estimated density ratio must be finite after numerical ",
         "guarding.", call. = FALSE)
  }
  ratio[reverse, ] <- 1 / ratio[reverse, , drop = FALSE]
  pA1_zM <- ifelse(reverse, p0, p1)
  pA1_zm <- ifelse(reverse, p1, p0)
  iv_upper <- (1 - pA1_zm) / (1 - pA1_zM)
  upper <- pmin(iv_upper, guard_upper)
  original_ratio <- ratio
  mass <- sweep(base, 2L, weights, "*")
  needs_projection <- rowSums(ratio < lower | ratio > upper) > 0L |
    abs(rowSums(ratio * mass) - 1) > tol
  shift <- numeric(n)
  iterations <- 0L
  if (any(needs_projection)) {
    ii <- which(needs_projection)
    rr <- ratio[ii, , drop = FALSE]
    mm <- mass[ii, , drop = FALSE]
    uu <- upper[ii]
    lo <- lower - apply(rr, 1L, max)
    hi <- uu - apply(rr, 1L, min)
    for (j in seq_len(80L)) {
      mid <- lo + (hi - lo) / 2
      projected <- pmin(pmax(rr + mid, lower), uu)
      integral <- rowSums(projected * mm)
      iterations <- j
      if (all(abs(integral - 1) <= tol)) break
      below <- integral < 1
      lo[below] <- mid[below]
      hi[!below] <- mid[!below]
    }
    if (any(abs(integral - 1) > 1e-10)) {
      stop("Could not normalize the bounded density ratio on the grid.",
           call. = FALSE)
    }
    ## Large guarded ratios can cause cancellation in rr + mid. Repair only
    ## the resulting roundoff by convex mixing toward the appropriate bound.
    ## Unlike division by the integral, this cannot leave the allowed interval.
    residual <- 1 - integral
    below <- which(residual > .Machine$double.eps)
    above <- which(residual < -.Machine$double.eps)
    if (length(below)) {
      room <- uu[below] - projected[below, , drop = FALSE]
      capacity <- rowSums(room * mm[below, , drop = FALSE])
      fraction <- pmin(residual[below] / capacity, 1)
      projected[below, ] <- (1 - fraction) * projected[below, , drop = FALSE] +
        fraction * uu[below]
    }
    if (length(above)) {
      room <- projected[above, , drop = FALSE] - lower
      capacity <- rowSums(room * mm[above, , drop = FALSE])
      fraction <- pmin(-residual[above] / capacity, 1)
      projected[above, ] <- (1 - fraction) * projected[above, , drop = FALSE] +
        fraction * lower
    }
    ratio[ii, ] <- projected
    shift[ii] <- mid
  }
  numerator <- ratio * base
  f0[reverse, ] <- numerator[reverse, , drop = FALSE]
  f1[!reverse, ] <- numerator[!reverse, , drop = FALSE]
  slack <- (1 - pA1_zm) * base - (1 - pA1_zM) * numerator
  change <- abs(ratio - original_ratio)
  out <- nuisance
  out$pA1_z0 <- p0
  out$pA1_z1 <- p1
  out$fY_z0 <- f0
  out$fY_z1 <- f1
  ## Public nuisance ratios always retain the fixed Z=1/Z=0 orientation.
  out$RY <- f1 / f0
  list(nuisance = out, diagnostics = list(
    enabled = TRUE,
    projected_units = sum(needs_projection),
    projected_fraction = mean(needs_projection),
    adjusted_grid_fraction = mean(change > tol),
    max_ratio_adjustment = max(change),
    lower_bound = lower, upper_bound = upper, iv_upper_bound = iv_upper,
    shift = shift, iterations = iterations,
    max_normalization_error = max(abs(as.vector(numerator %*% weights) - 1)),
    min_falsification_slack = min(slack)
  ))
}
