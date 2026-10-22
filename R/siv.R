## ----------------------------------------------------------------------------
## SIV(): the PROPOSED estimator of Park and Tchetgen Tchetgen.
##
## Sections 3 and 5 of the manuscript.  The estimator proceeds in three steps
## (Section 5.1):
##   Step 1  estimate theta = {f_Z, f_A, f_Y} from the observed data;
##   Step 2  obtain delta and gamma from theta via Theorems 3.1-3.2;
##   Step 3  plug into the efficient influence function of Theorem 4.2 and
##           average, eq. (22).
## Inference uses the variance estimator of Theorem 5.1.
## ----------------------------------------------------------------------------

#' Build an outcome grid and integration weights
#'
#' For a binary or discrete outcome the grid is the *observed support* of the
#' treatment-free outcome, `sort(unique(Y[A == 0]))`, with counting-measure
#' weights.  For a continuous outcome a fine grid over the observed range is
#' constructed (optionally augmented with the observed values, as in the
#' replication code of the manuscript) and mid-point cell widths are used as
#' integration weights.
#'
#' @param Y,A Outcome and treatment vectors.
#' @param type One of `"binary"`, `"discrete"`, `"continuous"`.
#' @param y_ref Reference outcome value; the grid is forced to contain it.
#' @param n_grid Number of equally spaced points for a continuous outcome.
#' @param include_observed For a continuous outcome, whether to add the observed
#'   values `Y[A == 0]` to the grid (the default, matching the replication code).
#' @param expand Multiple of `sd(Y[A == 0])` by which the grid range is expanded
#'   beyond the observed range of the untreated outcomes, for a continuous
#'   outcome.  Making the expansion proportional to the outcome's own scale
#'   keeps the grid invariant to the units in which `Y` is recorded.  If the
#'   standard deviation is not a positive finite number (a degenerate `Y[A ==
#'   0]`), `expand` is used as an absolute amount instead.
#' @return A list with components `grid` and `weights`.
#' @examples
#' set.seed(1); A <- rbinom(20, 1, 0.5); Y <- rbinom(20, 1, 0.5)
#' siv_grid(Y, A, "binary")
#' @export
siv_grid <- function(Y, A, type, y_ref = 0, n_grid = 101L,
                     include_observed = TRUE, expand = 0.1) {
  yo <- Y[A == 0]
  if (type %in% c("binary", "discrete")) {
    grid <- sort(unique(c(yo, y_ref)))
    if (type == "binary") {
      if (!setequal(grid, c(0, 1))) {
        stop("`type = \"binary\"` requires the observed support of `Y[A == 0]` ",
             "to be {0, 1}.", call. = FALSE)
      }
      grid <- c(0, 1)
    }
    return(list(grid = grid, weights = rep(1, length(grid))))
  }
  ## The expansion is measured in standard deviations of the untreated
  ## outcomes, so that the grid does not depend on the units of `Y`.
  sd_yo <- stats::sd(yo)
  pad <- expand * if (is.finite(sd_yo) && sd_yo > 0) sd_yo else 1
  rng <- range(yo) + c(-pad, pad)
  grid <- seq(rng[1L], rng[2L], length.out = n_grid)
  if (include_observed) grid <- c(grid, yo)
  grid <- sort(unique(grid))
  ## force one grid point to sit exactly at the reference value
  grid[which.min(abs(grid - y_ref))] <- y_ref
  grid <- sort(unique(grid))
  M <- length(grid)
  ctr <- (grid[-M] + grid[-1L]) / 2
  weights <- c(ctr[1L] - grid[1L], diff(ctr), grid[M] - ctr[M - 1L])
  list(grid = grid, weights = weights)
}

## Internal: default nuisance estimation with K-fold cross-fitting.
## Returns cross-fitted values of f(Z=1|X), f(A=1|Z=z,X) and f(y|A=0,Z=z,X).
siv_nuisance_default <- function(Y, A, Z, X, type, grid, weights, folds,
                                 instrument_formula, treatment_formula,
                                 outcome_formula, density_floor,
                                 outcome_density = "np_gss_cv", n_bins = 12L,
                                 density_cv_folds = 5L, gss_formula = NULL,
                                 ratio_formula = NULL,
                                 learner = "SL", SL.list = 1L,
                                 ratio_method = "kliep_ulsif",
                                 ratio_bandwidth = "cv_separate",
                                 ratio_bandwidth.cv.grid =
                                   c(0.5, 0.75, 0.9),
                                 ratio_kernels = 100L,
                                 ratio_folds = 5L, dr_floor = 1e-4,
                                 Aarg = NULL, Zarg = NULL) {
  n <- length(Y)
  M <- length(grid)
  df <- siv_frame(Y, A, Z, X)
  pZ1 <- numeric(n)
  pA1_z0 <- pA1_z1 <- numeric(n)
  fY_z0 <- fY_z1 <- matrix(0, n, M)
  RY <- if (type == "continuous") matrix(1, n, M) else NULL
  ratio_ensemble_used <- FALSE
  ratio_ensemble_fallback <- NULL
  ratio_ensemble_fit <- NULL

  Xnm <- colnames(X)
  rhsX <- if (length(Xnm)) paste(Xnm, collapse = " + ") else "1"
  outcome_formula_user <- !is.null(outcome_formula)
  if (is.null(instrument_formula)) instrument_formula <- stats::as.formula(paste("Z ~", rhsX))
  if (is.null(treatment_formula)) treatment_formula <- stats::as.formula(paste("A ~ Z +", rhsX))
  if (is.null(outcome_formula))  outcome_formula  <- stats::as.formula(paste("Y ~ Z +", rhsX))
  ## Targeted density-ratio model.  With a binary instrument, Bayes' rule gives
  ##   R_Y(y,X) = f(y|A=0,Z=1,X)/f(y|A=0,Z=0,X)
  ##            = { P(Z=1|Y=y,A=0,X)/P(Z=0|Y=y,A=0,X) }
  ##              * { P(Z=0|A=0,X)  /P(Z=1|A=0,X)   } ,
  ## so the density RATIO of eq. (15) -- the quantity the fixed point of
  ## Theorem 3.1 actually needs -- can be estimated from a single logistic
  ## regression of Z on (Y, X) among the untreated, instead of dividing two
  ## separately estimated densities.  Section 5.1 of the manuscript recommends
  ## exactly this: "one can employ strategies that target the density ratio
  ## directly, yielding more stable and reliable estimates".
  if (is.null(ratio_formula)) {
    ybas <- if (type == "continuous") "splines::ns(Y, df = 3)" else "factor(Y)"
    ratio_formula <- stats::as.formula(paste("Z ~", ybas, "+", rhsX))
  }
  outcome_learner <- if (type == "continuous" &&
                         outcome_formula_user) "glm" else learner

  baseline_density_fit <- NULL
  if (type == "continuous" && outcome_density %in%
      c("np_normal", "np_cv", "np_gss_cv")) {
    baseline_density_fit <- switch(
      outcome_density,
      np_normal = siv_baseline_density_np_normal(
        Y, A, Z, X, grid, weights, folds,
        density_floor = density_floor
      ),
      np_cv = siv_baseline_density_np_cv(
        Y, A, Z, X, grid, weights, folds,
        density_floor = density_floor
      ),
      np_gss_cv = siv_baseline_density_np_gss_cv(
        Y, A, Z, X, grid, weights, folds,
        density_floor = density_floor,
        density_cv_folds = density_cv_folds,
        gss_formula = gss_formula
      )
    )
    fY_z0 <- baseline_density_fit$density
  }

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (length(tr) == 0L) tr <- te            # crossfit = 1 -> in-sample

    ## f(Z = 1 | X)
    mz <- siv_learn_fit(df[tr, , drop = FALSE], "Z", Xnm, stats::binomial(),
                        learner, SL.list, formula = instrument_formula,
                        SL.N = n)
    pZ1[te] <- siv_learn_predict(mz, df[te, , drop = FALSE])

    ## f(A = 1 | Z, X)
    ma <- siv_learn_fit(df[tr, , drop = FALSE], "A", c("Z", Xnm),
                        stats::binomial(), learner, SL.list,
                        formula = treatment_formula, SL.N = n)
    d0 <- df[te, , drop = FALSE]; d0$Z <- 0
    d1 <- df[te, , drop = FALSE]; d1$Z <- 1
    pA1_z0[te] <- siv_learn_predict(ma, d0)
    pA1_z1[te] <- siv_learn_predict(ma, d1)

    ## f(Y | A = 0, Z, X)
    tr0 <- intersect(tr, which(A == 0))
    if (type == "continuous" && outcome_density %in%
        c("np_normal", "np_cv", "np_gss_cv")) {
      ## These nonparametric paths estimate only the baseline arm.  The Z = 1
      ## density is reconstructed below from the directly estimated ratio.
      fY_z1[te, ] <- fY_z0[te, ]
    } else if (type == "continuous" && outcome_density == "gaussian") {
      dtr0 <- df[tr0, , drop = FALSE]
      my <- siv_learn_fit(dtr0, "Y", c("Z", Xnm), stats::gaussian(),
                          outcome_learner, SL.list, formula = outcome_formula,
                          SL.N = n)
      s  <- siv_learn_sigma(my, dtr0, "Y")
      m0 <- siv_learn_predict(my, d0)
      m1 <- siv_learn_predict(my, d1)
      fY_z0[te, ] <- stats::dnorm(matrix(grid, length(te), M, byrow = TRUE), m0, s)
      fY_z1[te, ] <- stats::dnorm(matrix(grid, length(te), M, byrow = TRUE), m1, s)
    } else if (type == "continuous") {
      ## "binned": discretise the outcome into `n_bins` quantile bins among the
      ## untreated, fit a one-vs-rest logistic model for each bin, and convert
      ## the bin probabilities to a density on `grid` by dividing by the bin
      ## width.  This leaves the SHAPE of f(Y | A = 0, Z, X) unrestricted (up to
      ## the binning) while keeping the covariate dependence parametric, which
      ## matters because the true conditional density under a logit-separable
      ## mechanism is an exponentially tilted -- hence non-Gaussian -- density.
      brks <- unique(stats::quantile(Y[A == 0],
                                     probs = seq(0, 1, length.out = n_bins + 1L),
                                     names = FALSE))
      brks[1L] <- min(grid) - 1e-8
      brks[length(brks)] <- max(grid) + 1e-8
      nb <- length(brks) - 1L
      binw <- diff(brks)
      gbin <- cut(grid, breaks = brks, labels = FALSE, include.lowest = TRUE)
      P0 <- P1 <- matrix(0, length(te), nb)
      for (b in seq_len(nb)) {
        dtr <- df[tr0, , drop = FALSE]
        dtr$Y <- as.numeric(cut(dtr$Y, breaks = brks, labels = FALSE,
                                include.lowest = TRUE) == b)
        if (all(dtr$Y == 0) || all(dtr$Y == 1)) {
          P0[, b] <- P1[, b] <- mean(dtr$Y)
          next
        }
        mm <- siv_learn_fit(dtr, "Y", c("Z", Xnm), stats::binomial(),
                            outcome_learner, SL.list, formula = outcome_formula,
                            SL.N = n)
        P0[, b] <- siv_learn_predict(mm, d0)
        P1[, b] <- siv_learn_predict(mm, d1)
      }
      P0 <- pmax(P0, density_floor); P1 <- pmax(P1, density_floor)
      P0 <- P0 / rowSums(P0); P1 <- P1 / rowSums(P1)
      dens0 <- sweep(P0, 2L, binw, "/")
      dens1 <- sweep(P1, 2L, binw, "/")
      fY_z0[te, ] <- dens0[, gbin, drop = FALSE]
      fY_z1[te, ] <- dens1[, gbin, drop = FALSE]
    } else if (type == "binary") {
      ## Binary: pool all untreated observations and fit one model for
      ## Pr(Y = 1 | A = 0, Z, X). Predictions at Z = 0 and Z = 1 from this
      ## common fit preserve the full untreated sample in each outer fold and
      ## avoid instability from fitting the two instrument arms separately.
      dtr0 <- df[tr0, , drop = FALSE]
      my <- siv_learn_fit(dtr0, "Y", c("Z", Xnm), stats::binomial(),
                          learner, SL.list, formula = outcome_formula,
                          SL.N = n)
      py1_z0 <- siv_bound_prob(siv_learn_predict(my, d0), density_floor)
      py1_z1 <- siv_bound_prob(siv_learn_predict(my, d1), density_floor)
      fY_z0[te, ] <- cbind(1 - py1_z0, py1_z0)
      fY_z1[te, ] <- cbind(1 - py1_z1, py1_z1)
    } else {
      ## Discrete: one logistic model per support point (one-vs-rest),
      ## renormalised to a proper probability mass function.
      P0 <- P1 <- matrix(0, length(te), M)
      for (m in seq_len(M)) {
        dtr <- df[tr0, , drop = FALSE]
        dtr$Y <- as.numeric(abs(dtr$Y - grid[m]) < 1e-8)
        if (all(dtr$Y == 0) || all(dtr$Y == 1)) {
          P0[, m] <- P1[, m] <- mean(dtr$Y)
          next
        }
        mm <- siv_learn_fit(dtr, "Y", c("Z", Xnm), stats::binomial(),
                            learner, SL.list, formula = outcome_formula,
                            SL.N = n)
        P0[, m] <- siv_learn_predict(mm, d0)
        P1[, m] <- siv_learn_predict(mm, d1)
      }
      P0 <- pmax(P0, density_floor); P1 <- pmax(P1, density_floor)
      fY_z0[te, ] <- P0 / rowSums(P0 * matrix(weights, length(te), M, byrow = TRUE))
      fY_z1[te, ] <- P1 / rowSums(P1 * matrix(weights, length(te), M, byrow = TRUE))
    }

    ## targeted density ratio
    if (type == "continuous") {
      mr <- suppressWarnings(stats::glm(ratio_formula,
                                        data = df[tr0, , drop = FALSE],
                                        family = stats::binomial()))
      mq <- suppressWarnings(stats::glm(
        stats::as.formula(paste("Z ~", rhsX)),
        data = df[tr0, , drop = FALSE], family = stats::binomial()))
      qx <- stats::predict(mq, newdata = df[te, , drop = FALSE], type = "link")
      nte <- length(te)
      dte <- df[te, , drop = FALSE]
      if (siv_formula_additive_in_Y(ratio_formula, df)) {
        ## Additive in the Y-terms and the X-terms, so the linear predictor
        ## factorises:  lp(y, x) = lp(y, x_ref) + lp(y_1, x) - lp(y_1, x_ref).
        ## This costs M + n predictions instead of n * M.
        ref <- dte[1L, , drop = FALSE]
        gref <- ref[rep(1L, M), , drop = FALSE]; gref$Y <- grid
        a <- stats::predict(mr, newdata = gref, type = "link")     # length M
        d1 <- dte; d1$Y <- grid[1L]
        cvec <- stats::predict(mr, newdata = d1, type = "link")    # length n_te
        lp <- outer(as.numeric(cvec), as.numeric(a - a[1L]), "+")
      } else {
        big <- dte[rep(seq_len(nte), times = M), , drop = FALSE]
        big$Y <- rep(grid, each = nte)
        lp <- tryCatch(stats::predict(mr, newdata = big, type = "link"),
                       error = function(e) rep(0, nrow(big)))
        lp <- matrix(lp, nte, M)
      }
      lp[!is.finite(lp)] <- 0
      RY[te, ] <- exp(lp - qx)
    }
  }
  ## ---- continuous-outcome density ratio ------------------------------------
  ## The logistic fit above is retained as a stable fallback for every direct
  ## nonparametric ratio method.
  if (type == "continuous" && ratio_method != "logistic") {
    ratio_ensemble_fit <- if (identical(ratio_method, "kliep_ulsif")) {
      if (identical(ratio_bandwidth, "median_common")) {
        siv_kliep_ulsif_median_ratio_grid(
          Y, Aarg, Zarg, X, grid, weights, folds, fY_z0,
          kernel_num = ratio_kernels, fold = ratio_folds,
          dr_floor = dr_floor, density_floor = density_floor
        )
      } else {
        siv_kliep_ulsif_ratio_grid(
          Y, Aarg, Zarg, X, grid, weights, folds, fY_z0,
          kernel_num = ratio_kernels, fold = ratio_folds,
          dr_floor = dr_floor, density_floor = density_floor,
          bandwidth_mode = ratio_bandwidth,
          cv_grid = ratio_bandwidth.cv.grid
        )
      }
    } else {
      siv_single_density_ratio_grid(
        Y, Aarg, Zarg, X, grid, weights, folds, fY_z0,
        method = ratio_method, bandwidth = ratio_bandwidth,
        cv_grid = ratio_bandwidth.cv.grid,
        kernel_num = ratio_kernels, fold = ratio_folds,
        dr_floor = dr_floor, density_floor = density_floor
      )
    }
    if (is.null(ratio_ensemble_fit)) {
      arm <- attr(ratio_ensemble_fit, "arm_sizes")
      ratio_ensemble_fallback <- "logistic"
      warning("`ratio_method = \"", ratio_method, "\"`: the direct ",
              "density-ratio estimator could not be fitted on every fold; ",
              "falling back to the logistic ratio",
              if (!is.null(arm))
                sprintf(" (per-fold sizes of the (A=0, Z=1) / (A=0, Z=0) arms: %s)",
                        paste(apply(arm, 1L, paste, collapse = "/"),
                              collapse = ", ")) else "",
              ".", call. = FALSE)
    } else {
      RY <- ratio_ensemble_fit$ratio
      ratio_ensemble_used <- TRUE
    }
  }

  ## normalise the conditional densities against the integration weights
  Wm <- matrix(weights, n, M, byrow = TRUE)
  fY_z0 <- pmax(fY_z0, density_floor); fY_z1 <- pmax(fY_z1, density_floor)
  fY_z0 <- fY_z0 / rowSums(fY_z0 * Wm)
  fY_z1 <- fY_z1 / rowSums(fY_z1 * Wm)
  if (type == "binary") {
    ## Preserve the defining Bernoulli complement after numerical normalisation.
    fY_z0[, 1L] <- 1 - fY_z0[, 2L]
    fY_z1[, 1L] <- 1 - fY_z1[, 2L]
  } else if (type == "continuous") {
    ## rebuild the Z = 1 density from the targeted ratio, then renormalise
    fY_z1 <- pmax(RY, dr_floor) * fY_z0
    fY_z1 <- fY_z1 / rowSums(fY_z1 * Wm)
    RY <- fY_z1 / fY_z0
  }
  list(pZ1 = pZ1, pA1_z0 = pA1_z0, pA1_z1 = pA1_z1,
       fY_z0 = fY_z0, fY_z1 = fY_z1, RY = RY,
       ratio_ensemble_used = ratio_ensemble_used,
       ratio_ensemble_fallback = ratio_ensemble_fallback,
       density_ratio_ensemble = if (is.null(ratio_ensemble_fit)) NULL else list(
         method = ratio_ensemble_fit$method,
         candidates = ratio_ensemble_fit$candidates,
         weights = ratio_ensemble_fit$weights,
         mspe_weights = ratio_ensemble_fit$mspe_weights,
         mse_weights = ratio_ensemble_fit$mse_weights,
         kl_weights = ratio_ensemble_fit$kl_weights,
         sigma_candidates = ratio_ensemble_fit$sigma_candidates,
         sigma_factor = ratio_ensemble_fit$sigma_factor,
         selected_sigma = ratio_ensemble_fit$selected_sigma,
         selected_bandwidth = ratio_ensemble_fit$selected_bandwidth,
         bandwidth_range = ratio_ensemble_fit$bandwidth_range,
         arm_sizes = ratio_ensemble_fit$arm_sizes,
         ratio_method = ratio_ensemble_fit$ratio_method,
         ratio_bandwidth = ratio_ensemble_fit$ratio_bandwidth,
         ensemble_loss = ratio_ensemble_fit$ensemble_loss,
         cv_grid = ratio_ensemble_fit$cv_grid
       ),
       baseline_density = if (is.null(baseline_density_fit)) NULL else list(
         method = baseline_density_fit$method,
         candidates = baseline_density_fit$candidates,
         weights = baseline_density_fit$weights,
         mspe_weights = baseline_density_fit$mspe_weights,
         mse_weights = baseline_density_fit$mse_weights,
         kl_weights = baseline_density_fit$kl_weights,
         ensemble_loss = baseline_density_fit$ensemble_loss,
         bandwidth = baseline_density_fit$bandwidth
       ))
}

siv_estimate <- function(Y, A, Z, X = NULL,
                type = c("binary", "continuous", "discrete"),
                nuisance = NULL,
                grid = NULL, weights = NULL,
                y_ref = 0,
                crossfit = 1L, folds = NULL,
                solver = NULL,
                instrument_formula = NULL,
                treatment_formula = NULL,
                outcome_formula = NULL,
                outcome_density = c("np_gss_cv", "np_normal", "np_cv",
                                    "gaussian", "binned"),
                n_bins = 12L,
                density_cv_folds = 5L,
                gss_formula = NULL,
                learner = "SL", SL.list = 1L,
                ratio_method = c("kliep_ulsif", "kliep", "ulsif", "logistic"),
                ratio_bandwidth = c("cv_separate", "cv_common",
                                    "median_separate", "median_common"),
                ratio_bandwidth.cv.grid = c(0.5, 0.75, 0.9),
                ratio_kernels = 100L, ratio_folds = 5L,
                ratio_formula = NULL,
                n_grid = 101L, include_observed = TRUE,
                density_floor = 1e-4, dr_floor = 1e-4,
                maxit = 1000L, tol = 1e-10,
                level = 0.95,
                estimand = "ATT",
                QTT.quantile = c(0.25, 0.5, 0.75),
                QTT.c1 = 0.02 * (1 - level),
                num.QTT = 1000L,
                diagnostic_warnings = FALSE,
                falsification = FALSE,
                num.false.boot = 1000L,
                num.false.bin = c(20L, rep.int(5L, ncol(X))),
                falsification.floor = NULL,
                categorical.falsification = NULL,
                QTT.support = NULL, QTT.bounds = NULL) {
  crossfit_missing <- missing(crossfit)
  cl <- match.call()
  type <- match.arg(type)
  if (length(outcome_density) == 1L && identical(outcome_density, "np_gss")) {
    outcome_density <- "np_gss_cv"
  }
  outcome_density <- match.arg(outcome_density)
  ratio_method <- match.arg(ratio_method)
  ratio_bandwidth <- match.arg(ratio_bandwidth)
  if (!is.numeric(ratio_bandwidth.cv.grid) ||
      !length(ratio_bandwidth.cv.grid) ||
      any(!is.finite(ratio_bandwidth.cv.grid)) ||
      any(ratio_bandwidth.cv.grid <= 0 | ratio_bandwidth.cv.grid >= 1)) {
    stop("`ratio_bandwidth.cv.grid` must contain finite probabilities ",
         "strictly between 0 and 1.", call. = FALSE)
  }
  ratio_bandwidth.cv.grid <- unique(as.numeric(ratio_bandwidth.cv.grid))
  estimand <- unique(match.arg(estimand, c("ATT", "QTT"), several.ok = TRUE))
  if (!is.logical(falsification) || length(falsification) != 1L ||
      is.na(falsification)) {
    stop("`falsification` must be one non-missing logical value.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between 0 and 1.", call. = FALSE)
  }
  needs_crossfit <- "QTT" %in% estimand || isTRUE(falsification)
  if (needs_crossfit && is.null(nuisance)) {
    if (crossfit_missing) {
      crossfit <- 5L
    } else if (as.integer(crossfit)[1L] < 2L) {
      warning("QTT and falsification inference require cross-fitted nuisance ",
              "estimates; `crossfit = 1` was explicitly requested.",
              call. = FALSE)
    }
  }
  dat <- siv_check_data(Y, A, Z, X, allow_na_Y = FALSE)
  Y <- dat$Y; A <- dat$A; Z <- dat$Z; X <- dat$X; n <- dat$n
  if (any(!is.finite(Y))) {
    stop("`Y` must contain finite outcomes for ATT and QTT estimation.",
         call. = FALSE)
  }
  if (type != "discrete") siv_match_type(type, Y, A)
  requested_ratio_method <- ratio_method
  if (type != "continuous") ratio_method <- "plugin"
  if (is.null(solver)) solver <- if (type == "binary") "closed-form" else "contraction"
  solver <- match.arg(solver, c("closed-form", "contraction"))
  if (solver == "closed-form" && type != "binary") {
    stop("`solver = \"closed-form\"` is only available for `type = \"binary\"`.",
         call. = FALSE)
  }

  if (is.null(grid)) {
    gw <- siv_grid(Y, A, type, y_ref = y_ref, n_grid = n_grid,
                   include_observed = include_observed)
    grid <- gw$grid
    if (is.null(weights)) weights <- gw$weights
  } else {
    if (!is.numeric(grid) || any(!is.finite(grid)) ||
        length(grid) < if (type == "continuous") 2L else 1L ||
        is.unsorted(grid, strictly = TRUE)) {
      stop("`grid` must be a finite, strictly increasing numeric vector",
           if (type == "continuous") " with at least two values." else ".",
           call. = FALSE)
    }
    if (is.null(weights)) {
      if (type %in% c("binary", "discrete")) {
        weights <- rep(1, length(grid))
      } else {
        M_grid <- length(grid)
        centres <- (grid[-M_grid] + grid[-1L]) / 2
        weights <- c(centres[1L] - grid[1L], diff(centres),
                     grid[M_grid] - centres[M_grid - 1L])
      }
    }
  }
  if (!is.numeric(weights) || length(weights) != length(grid) ||
      any(!is.finite(weights)) || any(weights <= 0)) {
    stop("`weights` must contain one positive finite value per grid point.",
         call. = FALSE)
  }
  M <- length(grid)

  ## Folds are stratified on (A, Z) so that every fold contains all four
  ## cells; `siv_make_folds()` is shared with SIV_LATE().
  if (is.null(folds)) folds <- siv_make_folds(n, A, Z, crossfit)

  user_nuisance <- !is.null(nuisance)
  if (is.null(nuisance)) {
    nuisance <- siv_nuisance_default(Y, A, Z, X, type, grid, weights, folds,
                                     instrument_formula, treatment_formula,
                                     outcome_formula, density_floor,
                                     outcome_density = outcome_density,
                                     n_bins = n_bins,
                                     density_cv_folds = density_cv_folds,
                                     gss_formula = gss_formula,
                                     ratio_formula = ratio_formula,
                                     learner = learner, SL.list = SL.list,
                                     ratio_method = ratio_method,
                                     ratio_bandwidth = ratio_bandwidth,
                                     ratio_bandwidth.cv.grid =
                                       ratio_bandwidth.cv.grid,
                                     ratio_kernels = ratio_kernels,
                                     ratio_folds = ratio_folds,
                                     dr_floor = dr_floor,
                                     Aarg = A, Zarg = Z)
  } else {
    need <- c("pZ1", "pA1_z0", "pA1_z1", "fY_z0", "fY_z1")
    miss <- setdiff(need, names(nuisance))
    if (length(miss)) {
      stop("`nuisance` is missing element(s): ", paste(miss, collapse = ", "),
           ".", call. = FALSE)
    }
    if (!identical(dim(nuisance$fY_z0), c(n, M)) ||
        !identical(dim(nuisance$fY_z1), c(n, M))) {
      stop(sprintf("`nuisance$fY_z0` and `nuisance$fY_z1` must be %d x %d matrices.",
                   n, M), call. = FALSE)
    }
  }

  ## Keep `nuisance` unrestricted for the falsification test. Estimation uses
  ## a separate, normalized density-ratio projection satisfying its observable
  ## inequality. Do not apply another floor or renormalization in the solver,
  ## since that could undo the bound or make the ratio and densities disagree.
  ratio_fit <- siv_constrain_density_ratio(nuisance, weights,
                                           density_floor = density_floor,
                                           dr_floor = dr_floor)
  estimation_nuisance <- ratio_fit$nuisance
  fp <- siv_fixed_point(estimation_nuisance$pA1_z0,
                        estimation_nuisance$pA1_z1,
                        estimation_nuisance$fY_z0, estimation_nuisance$fY_z1,
                        grid, weights, y_ref = y_ref, solver = solver,
                        maxit = maxit, tol = tol, density_floor = 0,
                        dr_floor = 0)

  ef <- siv_eif(Y, A, Z, grid, weights, fp$delta, fp$gamma_z0, fp$gamma_z1,
                fp$fY_z0, fp$fY_z1, estimation_nuisance$pZ1,
                estimation_nuisance$pA1_z0, estimation_nuisance$pA1_z1)

  ## Diagnostics for Assumptions (A3) and (A5).  These are the two conditions
  ## that can fail silently when the nuisance functions are poorly estimated, so
  ## they are surfaced rather than left to produce an implausible estimate.
  ## The (A3)/(A5) diagnostics are always RETURNED in `$diagnostics`; whether
  ## they are also signalled as warnings is controlled by `diagnostic_warnings`,
  ## which is off by default because the thresholds fire routinely in simulation
  ## studies and drown the console.  Inspect `fit$diagnostics` instead.
  if (isTRUE(diagnostic_warnings) && fp$g_floor_frac > 1e-3) {
    warning(sprintf(paste0("The fixed point hit its lower bound for %.1f%% of ",
                           "(unit, outcome) pairs: the estimated nuisance ",
                           "functions violate Assumption (A3). Consider a ",
                           "coarser outcome grid or supplying better nuisance ",
                           "estimates via `nuisance =`."),
                    100 * fp$g_floor_frac), call. = FALSE)
  }
  if (isTRUE(diagnostic_warnings) &&
      is.finite(fp$contraction_modulus) && fp$contraction_modulus > 0.99) {
    warning(sprintf(paste0("The modulus of the fixed-point mapping reached ",
                           "%.6f. This indicates a weak contraction margin; ",
                           "check convergence and the sensitivity of the ",
                           "result to the outcome grid and density-ratio ",
                           "estimator."),
                    fp$contraction_modulus), call. = FALSE)
  }
  if (isTRUE(diagnostic_warnings) &&
      is.finite(ef$min_relevance) && ef$min_relevance < 1e-3) {
    warning(sprintf(paste0("min |Pr(Z=1|Y(0),A=1,X) - Pr(Z=1|X)| = %.2e is ",
                           "close to zero, so Assumption (A5) (strong ",
                           "relevance) is nearly violated and the influence ",
                           "function is numerically unstable."),
                    ef$min_relevance), call. = FALSE)
  }

  qtt <- if ("QTT" %in% estimand) {
    siv_qtt(Y, A, Z, type, grid, weights, fp, estimation_nuisance,
            quantiles = QTT.quantile, level = level, c1 = QTT.c1,
            num_bootstrap = num.QTT, support = QTT.support, bounds = QTT.bounds)
  } else NULL
  falsification_result <- if (isTRUE(falsification)) {
    siv_falsification(Y, A, Z, X, grid, weights, nuisance,
                      num_bootstrap = num.false.boot,
                      num_bins = num.false.bin,
                      variance_floor = falsification.floor,
                      categorical = categorical.falsification, level = level)
  } else NULL

  est <- numeric(0)
  if ("ATT" %in% estimand) est <- c(est, ATT = ef$att)
  if (!is.null(qtt)) est <- c(est, qtt$estimate)
  V <- matrix(NA_real_, length(est), length(est),
              dimnames = list(names(est), names(est)))
  if ("ATT" %in% names(est)) V["ATT", "ATT"] <- ef$var
  se <- sqrt(diag(V))
  names(se) <- names(est)
  zcrit <- stats::qnorm(1 - (1 - level) / 2)
  ci_matrix <- matrix(NA_real_, length(est), 2L,
                      dimnames = list(names(est), c("lower", "upper")))
  if ("ATT" %in% names(est)) {
    ci_matrix["ATT", ] <- ef$att + c(-1, 1) * zcrit * sqrt(ef$var)
  }
  if (!is.null(qtt)) ci_matrix[names(qtt$estimate), ] <- qtt$ci
  ## Preserve the original length-two `$ci` component for an ATT-only fit.
  ci <- if (identical(names(est), "ATT")) {
    stats::setNames(ci_matrix["ATT", ], c("lower", "upper"))
  } else ci_matrix
  se_type <- if (is.null(qtt)) {
    "influence function (Theorem 5.1)"
  } else if ("ATT" %in% names(est)) {
    paste("ATT: influence function (Theorem 5.1);",
          "QTT: Berger-Boos test inversion (Supplement A.6.3)")
  } else {
    "Berger-Boos test inversion (Supplement A.6.3); QTT has no regular Wald SE"
  }
  out <- list(
    estimate  = est,
    vcov      = V,
    se        = se,
    ci        = ci,
    ci_matrix = ci_matrix,
    level     = level,
    estimand  = estimand,
    qtt       = qtt,
    falsification = falsification_result,
    eif       = ef$eif,
    nuisance  = nuisance,
    estimation_nuisance = estimation_nuisance,
    ratio_constraint = ratio_fit$diagnostics,
    delta     = fp$delta,
    gamma_z0  = fp$gamma_z0,
    gamma_z1  = fp$gamma_z1,
    zM        = fp$zM,
    fixed_point = fp,
    influence_detail = ef,
    converged = all(fp$converged),
    prop_converged = mean(fp$converged),
    diagnostics = c(g_floor_frac = fp$g_floor_frac,
                    min_relevance = ef$min_relevance,
                    contraction_modulus = fp$contraction_modulus),
    n         = n,
    crossfit  = length(folds),
    folds     = folds,
    type      = type,
    solver    = solver,
    outcome_density = outcome_density,
    ratio_method = if (user_nuisance) {
      "supplied"
    } else if (type != "continuous") {
      "plugin (automatic)"
    } else if (isTRUE(nuisance$ratio_ensemble_used)) {
      requested_ratio_method
    } else if (requested_ratio_method != "logistic") {
      paste0("logistic (", requested_ratio_method, " fallback)")
    } else {
      "logistic"
    },
    ratio_bandwidth = if (type == "continuous") ratio_bandwidth else NULL,
    ratio_bandwidth.cv.grid = if (type == "continuous") {
      ratio_bandwidth.cv.grid
    } else NULL,
    grid      = grid,
    weights   = weights,
    se_type   = se_type,
    method    = "Park & Tchetgen Tchetgen separable IV (proposed)",
    call      = cl
  )
  class(out) <- c("SIV", "siv_fit")
  out
}

siv_finite_median <- function(x) {
  keep <- is.finite(x)
  if (any(keep)) stats::median(x[keep]) else NA_real_
}

siv_repeated_median_adjustment <- function(estimates, variances) {
  estimates <- as.matrix(estimates)
  variances <- as.matrix(variances)
  if (!identical(dim(estimates), dim(variances))) {
    stop("`estimates` and `variances` must have identical dimensions.",
         call. = FALSE)
  }
  valid_estimate <- is.finite(estimates)
  median_estimate <- apply(estimates, 2L, siv_finite_median)
  names(median_estimate) <- colnames(estimates)
  valid_variance <- valid_estimate & is.finite(variances)
  adjustable <- colSums(valid_variance) > 0L
  adjusted_variance <- rep(NA_real_, ncol(variances))
  adjusted_variance[adjustable] <- vapply(which(adjustable), function(j) {
    keep <- valid_variance[, j]
    stats::median(variances[keep, j] +
                    (estimates[keep, j] - median_estimate[j])^2)
  }, numeric(1L))
  names(adjusted_variance) <- colnames(variances)
  list(
    estimate = median_estimate,
    variance = adjusted_variance,
    adjustable = adjustable,
    valid_estimate = valid_estimate,
    valid_variance = valid_variance
  )
}

#' Separable instrumental-variable estimator
#'
#' Estimates the ATT or QTT under the separable-IV model and optionally runs
#' the marginal-cell falsification test. Instrument and treatment nuisance
#' functions are estimated with [MySL()]. For a binary outcome,
#' `Pr(Y = 1 | A = 0, Z, X)` is fitted once per cross-fitting fold using all
#' untreated observations, with `Z` included as a predictor. Predictions at
#' `Z = 0` and `Z = 1` give the two conditional mass functions, and the
#' probability for `Y = 0` is computed as the exact complement. For binary and
#' discrete outcomes, the density ratio is formed automatically from the
#' fitted mass functions. The `ratio_method` and `ratio_bandwidth` arguments
#' are used only for continuous outcomes.
#'
#' \NTC{For a continuous outcome, the default baseline density combines
#' \code{np} and \code{gss} candidates using inner-validation MSPE alone.
#' MSPE denotes mean squared prediction error, also known as integrated
#' squared error (ISE). The default KLIEP/uLSIF density-ratio ensemble uses
#' the same criterion, evaluated on the reconstructed conditional densities.
#' Both ensembles minimize the quadratic score directly over nonnegative
#' weights summing to one, including the endpoints. Their loss criterion is
#' fixed and has no user option.}
#'
#' For repeated cross-fitting estimates \eqn{\hat\tau_s}, `SIV()` reports
#' \deqn{\hat\tau_{median} = median_s(\hat\tau_s)}
#' For the ATT, whose split-specific variance estimates are
#' \eqn{\hat\sigma_s^2}, it reports
#' \deqn{\hat\sigma_{median}^2 = median_s\{\hat\sigma_s^2 +
#' (\hat\tau_s - \hat\tau_{median})^2\}.}
#' For each QTT, the Berger-Boos interval is calculated within every
#' cross-fitting repetition. The final lower and upper endpoints are the
#' respective medians of the split-specific endpoints.
#' \NTC{Treated-quantile intervals use binomial order-statistic ranks, with no
#' bootstrap simulation. Discrete quantiles are tested using both CDF limits
#' at a candidate support point. Continuous inversion includes complete cells
#' between every integration node and observed untreated outcome. An empty
#' inversion uses the full specified counterfactual support or known outcome
#' bounds. Known outcome bounds also constrain the treated-quantile interval
#' and the accepted counterfactual cells. The counterfactual fallback is
#' unbounded when neither support nor bounds are known. Unknown treated-outcome
#' bounds can also produce infinite endpoints through exterior binomial ranks.
#' Infinite endpoints are retained in the endpoint medians.}
#' When the falsification test is requested, its final p-value is the median of
#' the repetition-specific multiplier-bootstrap p-values. The final rejection
#' decision compares this median p-value with `1 - level`.
#'
#' @param Y Numeric outcome vector.
#' @param A Binary treatment indicator.
#' @param Z Binary instrument.
#' @param X Optional matrix or data frame of covariates.
#' @param estimand Any of `"ATT"` and `"QTT"`.
#' @param level \NTC{Confidence level, defaulting to 0.95. The nominal error
#'   probability is \code{1 - level}.}
#' @param y_type Outcome type.
#' @param Overall.Estimation_SL.list Learner groups passed to [MySL()].
#' @param Overall.Estimation_num.crossfit Number of independently randomized
#'   cross-fitting estimates combined by median adjustment.
#' @param Overall.Estimation_num.fold Number of folds used within each
#'   cross-fitting estimate.
#' @param y_grid,y_weights Optional outcome grid and integration weights. When
#'   `y_grid` is supplied without `y_weights`, counting-measure weights are used
#'   for binary or discrete outcomes and midpoint-cell weights are computed for
#'   a continuous outcome.
#' @param y_num.grid Number of base points in an automatically created
#'   continuous outcome grid.
#' @param y_ref Reference outcome value. When `NULL`, the untreated-outcome
#'   mean is used for a continuous outcome and zero is used for a binary or
#'   discrete outcome.
#' @param outcome_density \NTC{Continuous-outcome baseline-density estimator.
#'   The default \code{"np_gss_cv"} combines likelihood-CV \code{np} and
#'   \code{gss} density candidates using inner-validation MSPE alone.
#'   \code{"np_normal"} uses only \code{np} with a normal-reference bandwidth,
#'   and \code{"np_cv"} uses only \code{np} with \code{bwmethod = "cv.ml"}.
#'   The former name \code{"np_gss"} remains an alias for \code{"np_gss_cv"}.}
#' @param outcome_formula Optional formula for the Gaussian or binned outcome
#'   working model. Supplying it selects the formula-based GLM for that model.
#' @param outcome_num.bins Number of bins for the binned density estimator.
#' @param outcome_cv.folds Number of validation folds for the np/gss ensemble.
#' @param outcome_gss.formula Optional gss density formula. When `NULL`, the
#'   default contains the main effects of `Y` and every covariate, without
#'   interactions.
#' @param outcome_include.observed Whether the continuous grid includes the
#'   observed untreated outcomes.
#' @param outcome_floor Positive numerical floor for outcome densities.
#' @param ratio_method Continuous-outcome density-ratio method.
#' @param ratio_bandwidth \NTC{Bandwidth selection rule. The default
#'   \code{"cv_separate"} cross-validates separate outcome and standardized-
#'   covariate bandwidths. Other choices use common bandwidths or median
#'   heuristics.}
#' @param ratio_bandwidth.cv.grid Distance-quantile probabilities used to
#'   construct the candidate bandwidths for `"cv_common"` and
#'   `"cv_separate"`.
#' @param ratio_kernels Number of kernel centers.
#' @param ratio_folds Number of density-ratio validation folds.
#' @param ratio_formula Formula used by the logistic density-ratio method.
#' @param ratio_floor Positive numerical floor for density ratios.
#' @param contraction_solver Fixed-point solver.
#' @param contraction_maxit,contraction_tol Contraction controls.
#' @param QTT_quantile Quantile indices for QTT.
#' @param QTT_c1 \NTC{Berger-Boos adjustment constant. Defaults to
#'   \code{0.02 * (1 - level)}, giving 0.001 for a 95 percent interval.
#'   A user-specified value must lie strictly between zero and
#'   \code{1 - level}.}
#' @param QTT_numBB \NTC{Legacy bootstrap-count argument retained for call
#'   compatibility. It has no effect on QTT inference, which now uses an
#'   order-statistic treated-quantile interval.}
#' @param QTT_support \NTC{Optional complete finite support of the treatment-free
#'   outcome, as a strictly increasing numeric vector. Supplying it selects
#'   discrete QTT inference independently of the density-estimation
#'   \code{y_type}. For example, use \code{0:52} for integer weeks worked.
#'   When omitted, binary outcomes use \code{c(0, 1)}, discrete outcomes use
#'   \code{y_grid}, and continuous outcomes use interval inversion.}
#' @param QTT_bounds \NTC{Optional known lower and upper bounds for both
#'   potential outcomes, supplied as a finite ordered numeric pair. For example,
#'   use \code{c(0, 52)} for weeks worked. These bounds constrain treated-quantile
#'   intervals, accepted counterfactual cells, and the fallback for an empty
#'   inversion, so the QTT interval lies within the implied difference range.
#'   Binary outcomes default to \code{c(0, 1)}. Otherwise the default is
#'   \code{NULL}, and neither integration-grid endpoints nor sample extrema
#'   are interpreted as population bounds.}
#' @param falsification Whether to run the falsification test.
#' @param falsification_num.boot Number of multiplier-bootstrap draws.
#' @param falsification_num.bin Requested Y and marginal-X bin counts.
#' @param falsification_floor Positive variance floor. A data-dimension default
#'   is used when this is `NULL`.
#' @param falsification_categorical Names or indices of categorical X columns.
#' @param diagnostic_warnings Whether numerical diagnostics signal warnings.
#' @section Default options:
#' Character-vector arguments use their first value unless another value is
#' supplied. The defaults are summarized below.
#' \itemize{
#'   \item \NTC{Target and outcome.} The default target is `estimand = "ATT"` with
#'   `level = 0.95`. The outcome defaults to `y_type = "binary"`, an automatic
#'   grid and weights, 101 base grid points for a continuous outcome, and an
#'   automatic reference outcome. The reference is the untreated-outcome mean
#'   for a continuous outcome and zero otherwise.
#'   \item \NTC{Overall estimation.} `Overall.Estimation_SL.list = 1L` selects the
#'   first [MySL()] learner group. Three independently randomized estimates are
#'   computed by default, each using two folds. Point estimates use their
#'   repeated-cross-fitting medians. The ATT variance uses the median adjustment,
#'   while each QTT interval uses the medians of its Berger-Boos endpoints. A
#'   requested falsification test reports the median of its repeated p-values.
#'   \item \NTC{Continuous outcome density.} \NTC{The default is
#'   \code{outcome_density = "np_gss_cv"}, with likelihood-CV \code{np} and
#'   \code{gss} candidates, 5 validation folds, and MSPE-only weights.
#'   The binned alternative uses 12 bins. Observed untreated outcomes are
#'   included in the grid, \code{outcome_floor = 1e-4}, and formula arguments
#'   are \code{NULL} by default.}
#'   \item \NTC{Continuous density ratio.} \NTC{The defaults are
#'   \code{ratio_method = "kliep_ulsif"}, \code{ratio_bandwidth = "cv_separate"},
#'   and \code{ratio_bandwidth.cv.grid = c(0.5, 0.75, 0.9)}, with 100 kernel
#'   centers, 5 validation folds, a fixed MSPE criterion, no logistic formula,
#'   and \code{ratio_floor = 1e-4}. Binary and discrete outcomes instead form
#'   the ratio directly from their fitted mass functions.}
#'   \item \NTC{Contraction.} `contraction_solver = NULL` selects the closed-form
#'   binary solver for a binary outcome and the contraction solver otherwise.
#'   The iteration limit is 1000 and the tolerance is `1e-10`.
#'   \item \NTC{QTT.} \NTC{The default quantile indices are 0.25, 0.5, and 0.75.
#'   \code{QTT_c1 = 0.02 * (1 - level)} allocates two percent of the nominal
#'   error probability, giving 0.001 at the default \code{level = 0.95}.
#'   \code{QTT_support = NULL} uses the outcome type to choose the inversion.
#'   \code{QTT_bounds = NULL} leaves outcome bounds unspecified except for
#'   binary outcomes, which use \code{c(0, 1)}.
#'   The legacy \code{QTT_numBB = 1000L} argument is unused.}
#'   \item \NTC{Falsification and diagnostics.} The falsification test and diagnostic
#'   warnings are off by default. When enabled, the falsification test uses
#'   1000 bootstrap draws and default marginal bin counts. A `NULL`
#'   `falsification_floor` invokes the data-dimension default.
#' }
#' @return An object of class `"SIV"`. Its `median_adjustment` component stores
#'   every split-specific estimate and variance together with the folds and
#'   numerical diagnostics used in the adjustment. Non-finite split-specific
#'   estimates are excluded parameter by parameter, and the corresponding
#'   counts and indicators are recorded in that component.
#' @export
SIV <- function(
    Y,
    A,
    Z,
    X = NULL,
    estimand = "ATT",
    level = 0.95,
    Overall.Estimation_SL.list = 1L,
    Overall.Estimation_num.crossfit = 3L,
    Overall.Estimation_num.fold = 2L,
    y_type = c("binary", "continuous", "discrete"),
    y_grid = NULL,
    y_weights = NULL,
    y_num.grid = 101L,
    y_ref = NULL,
    outcome_density = c("np_gss_cv", "np_normal", "np_cv", "gaussian",
                        "binned"),
    outcome_formula = NULL,
    outcome_num.bins = 12L,
    outcome_cv.folds = 5L,
    outcome_gss.formula = NULL,
    outcome_include.observed = TRUE,
    outcome_floor = 1e-4,
    ratio_method = c("kliep_ulsif", "kliep", "ulsif", "logistic"),
    ratio_bandwidth = c("cv_separate", "cv_common", "median_separate",
                        "median_common"),
    ratio_bandwidth.cv.grid = c(0.5, 0.75, 0.9),
    ratio_kernels = 100L,
    ratio_folds = 5L,
    ratio_formula = NULL,
    ratio_floor = 1e-4,
    contraction_solver = NULL,
    contraction_maxit = 1000L,
    contraction_tol = 1e-10,
    QTT_quantile = c(0.25, 0.5, 0.75),
    QTT_c1 = 0.02 * (1 - level),
    QTT_numBB = 1000L,
    falsification = FALSE,
    falsification_num.boot = 1000L,
    falsification_num.bin = c(20L, rep.int(5L, ncol(X))),
    falsification_floor = NULL,
    falsification_categorical = NULL,
    diagnostic_warnings = FALSE,
    QTT_support = NULL,
    QTT_bounds = NULL) {
  cl <- match.call()
  y_type <- match.arg(y_type)
  validate_count <- function(x, name) {
    value <- as.numeric(x)
    if (length(value) != 1L || !is.finite(value) || value < 1 ||
        abs(value - round(value)) > sqrt(.Machine$double.eps)) {
      stop("`", name, "` must be one positive integer.", call. = FALSE)
    }
    as.integer(round(value))
  }
  num_crossfit <- validate_count(
    Overall.Estimation_num.crossfit, "Overall.Estimation_num.crossfit"
  )
  num_fold <- validate_count(
    Overall.Estimation_num.fold, "Overall.Estimation_num.fold"
  )
  y_num.grid <- validate_count(y_num.grid, "y_num.grid")
  automatic_y_ref <- is.null(y_ref)
  if (automatic_y_ref) {
    y_ref <- if (identical(y_type, "continuous")) {
      reference_data <- siv_check_data(Y, A, Z, X, allow_na_Y = TRUE)
      mean(reference_data$Y[reference_data$A == 0])
    } else {
      0
    }
  }
  if (automatic_y_ref && !is.null(y_grid) &&
      !any(abs(as.numeric(y_grid) - as.numeric(y_ref)) < 1e-12)) {
    if (!is.null(y_weights)) {
      stop("An automatic `y_ref` not already in `y_grid` cannot be inserted ",
           "when `y_weights` is supplied. Add the reference value to both ",
           "arguments or set `y_ref` explicitly.", call. = FALSE)
    }
    y_grid <- sort(unique(c(as.numeric(y_grid), as.numeric(y_ref))))
  }
  if (!is.numeric(outcome_floor) || length(outcome_floor) != 1L ||
      !is.finite(outcome_floor) || outcome_floor <= 0) {
    stop("`outcome_floor` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(ratio_floor) || length(ratio_floor) != 1L ||
      !is.finite(ratio_floor) || ratio_floor <= 0 || ratio_floor >= 1) {
    stop("`ratio_floor` must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  estimates <- variances <- NULL
  repeat_folds <- vector("list", num_crossfit)
  repeat_diagnostics <- vector("list", num_crossfit)
  repeat_ratio_constraints <- vector("list", num_crossfit)
  repeat_confidence_intervals <- vector("list", num_crossfit)
  repeat_qtt <- vector("list", num_crossfit)
  repeat_falsification <- vector("list", num_crossfit)
  repeat_converged <- logical(num_crossfit)
  repeat_prop_converged <- numeric(num_crossfit)
  for (s in seq_len(num_crossfit)) {
    fit <- siv_estimate(
      Y = Y, A = A, Z = Z, X = X,
      estimand = estimand, level = level,
      type = y_type, grid = y_grid, weights = y_weights, y_ref = y_ref,
      crossfit = num_fold, folds = NULL,
      instrument_formula = NULL, treatment_formula = NULL,
      outcome_density = outcome_density, outcome_formula = outcome_formula,
      n_bins = outcome_num.bins, density_cv_folds = outcome_cv.folds,
      gss_formula = outcome_gss.formula,
      include_observed = outcome_include.observed,
      density_floor = outcome_floor,
      learner = "SL", SL.list = Overall.Estimation_SL.list,
      ratio_method = ratio_method, ratio_bandwidth = ratio_bandwidth,
      ratio_bandwidth.cv.grid = ratio_bandwidth.cv.grid,
      ratio_kernels = ratio_kernels, ratio_folds = ratio_folds,
      n_grid = y_num.grid, ratio_formula = ratio_formula,
      dr_floor = ratio_floor,
      solver = contraction_solver, maxit = contraction_maxit,
      tol = contraction_tol,
      QTT.quantile = QTT_quantile, QTT.c1 = QTT_c1, num.QTT = QTT_numBB,
      QTT.support = QTT_support, QTT.bounds = QTT_bounds,
      falsification = falsification,
      num.false.boot = falsification_num.boot,
      num.false.bin = falsification_num.bin,
      falsification.floor = falsification_floor,
      categorical.falsification = falsification_categorical,
      diagnostic_warnings = diagnostic_warnings
    )
    estimates <- rbind(estimates, fit$estimate)
    variances <- rbind(variances, diag(fit$vcov))
    repeat_folds[[s]] <- fit$folds
    repeat_diagnostics[[s]] <- fit$diagnostics
    repeat_ratio_constraints[[s]] <- fit$ratio_constraint
    repeat_confidence_intervals[[s]] <- fit$ci_matrix
    repeat_qtt[[s]] <- fit$qtt
    repeat_falsification[[s]] <- fit$falsification
    repeat_converged[s] <- isTRUE(fit$converged)
    repeat_prop_converged[s] <- fit$prop_converged
  }

  colnames(estimates) <- names(fit$estimate)
  colnames(variances) <- names(fit$estimate)
  adjustment <- siv_repeated_median_adjustment(estimates, variances)
  median_estimate <- adjustment$estimate
  adjusted_variance <- adjustment$variance
  adjustable <- adjustment$adjustable
  valid_estimate <- adjustment$valid_estimate
  valid_variance <- adjustment$valid_variance

  fit$estimate[] <- median_estimate
  diag(fit$vcov)[adjustable] <- adjusted_variance[adjustable]
  fit$se <- sqrt(diag(fit$vcov))
  names(fit$se) <- names(fit$estimate)
  qtt_names <- grep("^QTT_", names(fit$estimate), value = TRUE)
  if (length(qtt_names)) {
    for (nm in qtt_names) {
      fit$ci_matrix[nm, "lower"] <- siv_qtt_endpoint_median(vapply(
        repeat_confidence_intervals,
        function(ci) unname(ci[nm, "lower"]), numeric(1L)
      ))
      fit$ci_matrix[nm, "upper"] <- siv_qtt_endpoint_median(vapply(
        repeat_confidence_intervals,
        function(ci) unname(ci[nm, "upper"]), numeric(1L)
      ))
    }
    fit$qtt$estimate <- fit$estimate[qtt_names]
    fit$qtt$ci <- fit$ci_matrix[qtt_names, , drop = FALSE]
    fit$qtt$counterfactual_quantile <-
      fit$qtt$treated_quantile - fit$qtt$estimate
    empty_by_repeat <- do.call(rbind, lapply(repeat_qtt, function(qtt) {
      unname(qtt$empty_inversion[qtt_names])
    }))
    colnames(empty_by_repeat) <- qtt_names
    fit$qtt$empty_inversion <- stats::setNames(
      apply(empty_by_repeat, 2L, any), qtt_names
    )
    fit$qtt$empty_inversion_count <- stats::setNames(
      colSums(empty_by_repeat), qtt_names
    )
    fit$qtt$empty_inversion_by_repeat <- empty_by_repeat
  }
  if ("ATT" %in% names(fit$estimate) && is.finite(fit$vcov["ATT", "ATT"])) {
    zcrit <- stats::qnorm(1 - (1 - level) / 2)
    fit$ci_matrix["ATT", ] <- fit$estimate["ATT"] +
      c(-1, 1) * zcrit * sqrt(fit$vcov["ATT", "ATT"])
  }
  fit$ci <- if (identical(names(fit$estimate), "ATT")) {
    stats::setNames(fit$ci_matrix["ATT", ], c("lower", "upper"))
  } else {
    fit$ci_matrix
  }
  if (isTRUE(falsification)) {
    p_values <- vapply(repeat_falsification, `[[`, numeric(1L), "p_value")
    statistics <- vapply(repeat_falsification, `[[`, numeric(1L), "statistic")
    critical_values <- vapply(
      repeat_falsification, `[[`, numeric(1L), "critical_value"
    )
    rejections <- vapply(repeat_falsification, `[[`, logical(1L), "reject")
    median_p_value <- stats::median(p_values)
    representative <- which.min(abs(p_values - median_p_value))[1L]
    fit$falsification <- repeat_falsification[[representative]]
    fit$falsification$p_value <- median_p_value
    if (num_crossfit > 1L) {
      fit$falsification$reject <- median_p_value <= 1 - level
    }
    fit$falsification$median_adjustment <- list(
      applied = num_crossfit > 1L,
      num_crossfit = num_crossfit,
      p_values = p_values,
      median_p_value = median_p_value,
      statistics = statistics,
      critical_values = critical_values,
      rejections = rejections,
      representative = representative
    )
    if (num_crossfit > 1L) {
      fit$falsification$method <- paste0(
        fit$falsification$method, "; median p-value across ",
        num_crossfit, " repeated cross-fitting estimates"
      )
    }
  }
  diagnostics_matrix <- do.call(rbind, repeat_diagnostics)
  fit$diagnostics <- apply(diagnostics_matrix, 2L, stats::median)
  fit$converged <- all(repeat_converged)
  fit$prop_converged <- mean(repeat_prop_converged)
  fit$crossfit <- num_fold
  fit$num_crossfit <- num_crossfit
  fit$num_fold <- num_fold
  fit$y_ref <- y_ref
  fit$median_adjustment <- list(
    applied = num_crossfit > 1L,
    num_crossfit = num_crossfit,
    num_fold = num_fold,
    estimates = estimates,
    variances = variances,
    median_estimate = median_estimate,
    confidence_intervals = repeat_confidence_intervals,
    combined_confidence_interval = fit$ci_matrix,
    falsification_p_values = if (isTRUE(falsification)) {
      fit$falsification$median_adjustment$p_values
    } else NULL,
    adjusted_parameters = colnames(estimates)[adjustable],
    adjusted_variance = adjusted_variance,
    valid_estimate = valid_estimate,
    valid_variance = valid_variance,
    num_valid_estimate = colSums(valid_estimate),
    num_valid_variance = colSums(valid_variance),
    diagnostics = diagnostics_matrix,
    ratio_constraints = repeat_ratio_constraints,
    converged = repeat_converged,
    prop_converged = repeat_prop_converged,
    folds = repeat_folds
  )
  if (num_crossfit > 1L && "ATT" %in% names(fit$estimate)) {
    att_se_type <- paste0(
      "median-adjusted repeated cross-fitting (S = ", num_crossfit,
      ", folds = ", num_fold, "); influence function (Theorem 5.1)"
    )
    fit$se_type <- if (is.null(fit$qtt)) {
      att_se_type
    } else {
      paste0("ATT: ", att_se_type,
             "; QTT: median of repeated Berger-Boos intervals ",
             "(Supplement A.6.3 and A.8.4)")
    }
  } else if (num_crossfit > 1L && length(qtt_names)) {
    fit$se_type <- paste0(
      "median of repeated Berger-Boos intervals ",
      "(Supplement A.6.3 and A.8.4); QTT has no regular Wald SE"
    )
  }
  fit$call <- cl
  fit
}
