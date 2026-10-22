## Nonparametric baseline-density estimation for
## f(Y | A = 0, Z = 0, X).  This ports the np/gss candidate library and the
## inner-validation ensemble used in Data_Version1/1.SIV_JC.R and
## Simulation_Version1/1.Sim_Conti_SepIV.R.

## Internal: floor and normalise each row against the outcome-grid weights.
siv_density_normalize <- function(density, weights, density_floor) {
  density <- as.matrix(density)
  density[!is.finite(density)] <- density_floor
  density <- pmax(density, density_floor)
  den <- rowSums(density * matrix(weights, nrow(density), length(weights),
                                  byrow = TRUE))
  bad <- !is.finite(den) | den <= 0
  if (any(bad)) {
    density[bad, ] <- 1
    den[bad] <- sum(weights)
  }
  density / den
}

## Internal: fit the np bandwidth object on the A = 0, Z = 0 training arm.
siv_np_cde_fit <- function(data, xnames,
                           bwmethod = c("cv.ml", "normal-reference")) {
  bwmethod <- match.arg(bwmethod)
  if (length(xnames)) {
    form <- stats::reformulate(xnames, response = "Y")
    bandwidth <- NULL
    invisible(utils::capture.output(
      bandwidth <- do.call(np::npcdensbw,
                           list(formula = form, data = data,
                                bwmethod = bwmethod))
    ))
    return(list(kind = "conditional",
                bandwidth = bandwidth,
                xnames = xnames))
  }
  bandwidth <- NULL
  invisible(utils::capture.output(
    bandwidth <- np::npudensbw(~Y, data = data, bwmethod = bwmethod)
  ))
  list(kind = "unconditional",
       bandwidth = bandwidth,
       xnames = character(0))
}

## Internal: evaluate the np estimate for every new X and outcome-grid pair in
## one call.  np accepts paired evaluation rows, so arranging the data in grid
## blocks avoids the per-observation R loop used by the legacy np/gss path.
siv_np_cde_predict_grid <- function(object, training, newx, grid) {
  nr <- nrow(newx)
  M <- length(grid)
  if (identical(object$kind, "conditional")) {
    fit <- np::npcdens(
      bws = object$bandwidth,
      txdat = training[, object$xnames, drop = FALSE],
      tydat = training$Y,
      exdat = newx[rep(seq_len(nr), times = M), object$xnames, drop = FALSE],
      eydat = rep(grid, each = nr)
    )
  } else {
    fit <- np::npudens(
      bws = object$bandwidth,
      tdat = training["Y"],
      edat = data.frame(Y = grid)
    )
  }
  pred <- as.numeric(stats::fitted(fit))
  if (identical(object$kind, "unconditional")) {
    return(matrix(pred, nr, M, byrow = TRUE))
  }
  matrix(pred, nr, M)
}

## Internal: evaluate the np estimate.  `training` can be an inner-training
## subset while `object$bandwidth` stays fixed at the outer-training value,
## exactly as in the Version1 code.
siv_np_cde_predict <- function(object, training, newx, y) {
  y <- as.numeric(y)
  if (identical(object$kind, "conditional")) {
    fit <- np::npcdens(
      bws = object$bandwidth,
      txdat = training[, object$xnames, drop = FALSE],
      tydat = training$Y,
      exdat = newx[, object$xnames, drop = FALSE],
      eydat = y
    )
  } else {
    fit <- np::npudens(
      bws = object$bandwidth,
      tdat = training["Y"],
      edat = data.frame(Y = y)
    )
  }
  as.numeric(stats::fitted(fit))
}

## Internal: build and fit the gss conditional-density candidate.  With no
## covariates, ssden/dssden supplies the corresponding unconditional estimate.
siv_gss_default_formula <- function(xnames) {
  stats::as.formula(paste("~", paste(c("Y", xnames), collapse = " + ")))
}

siv_gss_cde_fit <- function(data, xnames, ydomain,
                            gss_formula = NULL) {
  ydomain <- range(c(ydomain, data$Y), finite = TRUE)
  if (!length(xnames)) {
    return(list(kind = "unconditional",
                fit = gss::ssden(~Y, data = data["Y"],
                                 domain = list(Y = ydomain)),
                xnames = character(0)))
  }
  if (is.null(gss_formula)) {
    gss_formula <- siv_gss_default_formula(xnames)
  }
  list(kind = "conditional",
       fit = gss::sscden(formula = gss_formula, response = ~Y, data = data,
                         ydomain = data.frame(Y = ydomain)),
       xnames = xnames)
}

## Internal: evaluate one gss density row.  Evaluation failures are handled by
## the caller, which substitutes the np candidate as in the Version1 code.
siv_gss_cde_predict_row <- function(object, newx, y) {
  if (identical(object$kind, "unconditional")) {
    return(as.numeric(gss::dssden(object$fit, data.frame(Y = y))))
  }
  as.numeric(gss::dsscden(object$fit, y = y,
                          x = newx[, object$xnames, drop = FALSE]))
}

## Internal: evaluate both candidates on an observation-by-grid array.
siv_cde_candidate_predict <- function(np_fit, gss_fit, training, newx, grid,
                                      weights, density_floor) {
  nr <- nrow(newx)
  M <- length(grid)
  np_density <- matrix(0, nr, M)
  gss_density <- matrix(0, nr, M)
  for (i in seq_len(nr)) {
    xi <- newx[i, , drop = FALSE]
    xgrid <- xi[rep(1L, M), , drop = FALSE]
    np_row <- siv_np_cde_predict(np_fit, training, xgrid, grid)
    np_density[i, ] <- np_row
    gss_row <- tryCatch(
      siv_gss_cde_predict_row(gss_fit, xi, grid),
      error = function(e) np_row
    )
    if (length(gss_row) != M || any(!is.finite(gss_row))) gss_row <- np_row
    gss_density[i, ] <- gss_row
  }
  list(
    np = siv_density_normalize(np_density, weights, density_floor),
    gss = siv_density_normalize(gss_density, weights, density_floor)
  )
}

## Internal: obtain a candidate density at each observed outcome.  Version1
## put every untreated outcome on the integration grid.  Linear interpolation
## retains the same values in that case and also supports a user-supplied grid.
siv_density_at_observed <- function(density, grid, y) {
  vapply(seq_along(y), function(i) {
    stats::approx(grid, density[i, ], xout = y[i], rule = 2,
                  ties = "ordered")$y
  }, numeric(1L))
}

## \NTC{Choose both np/gss and KLIEP/uLSIF mixture weights by MSPE (also known
## as ISE). Evaluate the quadratic's stationary point and both endpoints on
## the exact simplex. Retain `mse` as an alias for `mspe`.}
siv_density_mspe_weights <- function(H, h) {
  if (!is.numeric(H) || !identical(dim(H), c(2L, 2L)) ||
      !is.numeric(h) || length(h) != 2L ||
      any(!is.finite(H)) || any(!is.finite(h))) {
    stop("MSPE weighting requires two candidate densities and finite scores.",
         call. = FALSE)
  }

  weight <- function(alpha) c(alpha, 1 - alpha)
  mspe_objective <- function(alpha) {
    w <- weight(alpha)
    as.numeric(0.5 * t(w) %*% H %*% w - sum(w * h))
  }
  cross <- (H[1L, 2L] + H[2L, 1L]) / 2
  curvature <- H[1L, 1L] - 2 * cross + H[2L, 2L]
  candidates <- c(0, 1)
  if (curvature > 0) {
    stationary <- (h[1L] - h[2L] - cross + H[2L, 2L]) / curvature
    candidates <- c(candidates, pmin(1, pmax(0, stationary)))
  }
  alpha <- candidates[which.min(vapply(candidates, mspe_objective, numeric(1L)))]
  wmspe <- weight(alpha)
  list(
    weight = wmspe,
    mspe = wmspe,
    mse = wmspe,
    loss = "mspe"
  )
}

## \NTC{Retain the internal ratio helper name for existing callers.}
siv_ratio_ensemble_weights <- function(H, h) {
  siv_density_mspe_weights(H, h)
}

## Internal: split the baseline arm for inner validation.
siv_density_inner_folds <- function(index, nfold) {
  nfold <- max(2L, min(as.integer(nfold), length(index)))
  label <- sample(rep_len(seq_len(nfold), length(index)))
  split(index, label)
}

## Internal: cross-fitted np estimator for f(Y | A = 0, Z = 0, X).
## The grid is evaluated in one vectorised np call per outer fold.
siv_baseline_density_np <- function(
    Y, A, Z, X, grid, weights, folds, density_floor = 1e-4,
    bwmethod = c("normal-reference", "cv.ml"), option = "np_normal") {
  bwmethod <- match.arg(bwmethod)
  if (!requireNamespace("np", quietly = TRUE)) {
    stop("`outcome_density = \"", option, "\"` requires the `np` package.",
         call. = FALSE)
  }
  n <- length(Y)
  M <- length(grid)
  X <- as.matrix(X)
  xnames <- colnames(X)
  df <- data.frame(Y = Y)
  if (length(xnames)) df <- cbind(df, as.data.frame(X))
  baseline <- which(A == 0 & Z == 0)
  out <- matrix(0, n, M)
  fold_weights <- matrix(
    1, length(folds), 1L, dimnames = list(NULL, "np")
  )

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (!length(tr)) tr <- te
    train_index <- intersect(tr, baseline)
    if (length(train_index) < max(10L, ncol(X) + 3L)) {
      stop("Too few A = 0, Z = 0 observations in cross-fitting training fold ",
           k, " for `outcome_density = \"", option, "\"`.",
           call. = FALSE)
    }
    training <- df[train_index, , drop = FALSE]
    np_fit <- siv_np_cde_fit(
      training, xnames, bwmethod = bwmethod
    )
    out[te, ] <- siv_np_cde_predict_grid(
      np_fit, training, df[te, , drop = FALSE], grid
    )
  }

  list(
    density = siv_density_normalize(out, weights, density_floor),
    weights = fold_weights,
    mspe_weights = NULL,
    mse_weights = NULL,
    kl_weights = NULL,
    candidates = "np",
    bandwidth = bwmethod,
    method = paste("cross-fitted np", bwmethod, "estimator")
  )
}

siv_baseline_density_np_normal <- function(
    Y, A, Z, X, grid, weights, folds, density_floor = 1e-4) {
  siv_baseline_density_np(
    Y, A, Z, X, grid, weights, folds, density_floor,
    bwmethod = "normal-reference", option = "np_normal"
  )
}

siv_baseline_density_np_cv <- function(
    Y, A, Z, X, grid, weights, folds, density_floor = 1e-4) {
  siv_baseline_density_np(
    Y, A, Z, X, grid, weights, folds, density_floor,
    bwmethod = "cv.ml", option = "np_cv"
  )
}

## Internal: cross-fitted np/gss ensemble for f(Y | A = 0, Z = 0, X).
siv_baseline_density_np_gss_cv <- function(Y, A, Z, X, grid, weights, folds,
                                           density_floor = 1e-4,
                                           density_cv_folds = 5L,
                                           gss_formula = NULL) {
  if (!requireNamespace("np", quietly = TRUE) ||
      !requireNamespace("gss", quietly = TRUE)) {
    stop("`outcome_density = \"np_gss_cv\"` requires the `np` and `gss` ",
         "packages.", call. = FALSE)
  }
  if (!is.numeric(density_cv_folds) || length(density_cv_folds) != 1L ||
      !is.finite(density_cv_folds) || density_cv_folds < 2) {
    stop("`density_cv_folds` must be an integer of at least 2.",
         call. = FALSE)
  }
  density_cv_folds <- as.integer(density_cv_folds)
  n <- length(Y)
  M <- length(grid)
  X <- as.matrix(X)
  xnames <- colnames(X)
  df <- data.frame(Y = Y)
  if (length(xnames)) df <- cbind(df, as.data.frame(X))
  baseline <- which(A == 0 & Z == 0)
  out <- matrix(0, n, M)
  fold_weights <- matrix(NA_real_, length(folds), 2L,
                         dimnames = list(NULL, c("np", "gss")))
  mspe_weights <- fold_weights

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (!length(tr)) tr <- te
    train_index <- intersect(tr, baseline)
    if (length(train_index) < max(10L, ncol(X) + 3L)) {
      stop("Too few A = 0, Z = 0 observations in cross-fitting training fold ",
           k, " for the np/gss baseline-density estimator.", call. = FALSE)
    }
    training <- df[train_index, , drop = FALSE]
    np_fit <- siv_np_cde_fit(training, xnames)
    gss_fit <- siv_gss_cde_fit(training, xnames, range(grid), gss_formula)

    inner <- siv_density_inner_folds(train_index, density_cv_folds)
    cv_np <- cv_gss <- matrix(NA_real_, length(train_index), M,
                              dimnames = list(as.character(train_index), NULL))
    for (v in seq_along(inner)) {
      va <- inner[[v]]
      inn_tr <- setdiff(train_index, va)
      inner_training <- df[inn_tr, , drop = FALSE]
      inner_gss <- tryCatch(
        siv_gss_cde_fit(inner_training, xnames, range(grid), gss_formula),
        error = function(e) NULL
      )
      if (is.null(inner_gss)) inner_gss <- gss_fit
      cand <- siv_cde_candidate_predict(
        np_fit, inner_gss, inner_training, df[va, , drop = FALSE], grid,
        weights, density_floor
      )
      cv_np[as.character(va), ] <- cand$np
      cv_gss[as.character(va), ] <- cand$gss
    }

    observed <- cbind(
      np = siv_density_at_observed(cv_np, grid, Y[train_index]),
      gss = siv_density_at_observed(cv_gss, grid, Y[train_index])
    )
    H <- matrix(0, 2L, 2L)
    candidates <- list(cv_np, cv_gss)
    for (j in 1:2) {
      for (l in j:2) {
        H[j, l] <- H[l, j] <- mean(rowSums(
          candidates[[j]] * candidates[[l]] *
            matrix(weights, length(train_index), M, byrow = TRUE)
        ))
      }
    }
    h <- colMeans(observed)
    ew <- siv_density_mspe_weights(H, h)
    fold_weights[k, ] <- ew$weight
    mspe_weights[k, ] <- ew$mspe

    pred <- siv_cde_candidate_predict(
      np_fit, gss_fit, training, df[te, , drop = FALSE], grid, weights,
      density_floor
    )
    out[te, ] <- ew$weight[1L] * pred$np + ew$weight[2L] * pred$gss
  }

  list(
    density = siv_density_normalize(out, weights, density_floor),
    weights = fold_weights,
    mspe_weights = mspe_weights,
    mse_weights = mspe_weights,
    kl_weights = NULL,
    ensemble_loss = "mspe",
    candidates = c("np", "gss"),
    bandwidth = "cv.ml",
    method = "cross-fitted np/gss MSPE ensemble"
  )
}

## Backward-compatible internal alias for code written before the option was
## renamed to `np_gss_cv`.
siv_baseline_density_np_gss <- siv_baseline_density_np_gss_cv
