lsiv_nuisance_probabilities <- function(Y, A, Z, X, folds,
                                        instrument_formula, treatment_formula,
                                        learner, SL.list) {
  n <- length(Y)
  df <- lsiv_frame(Y, A, Z, X)
  Xnm <- colnames(X)
  rhsX <- if (length(Xnm)) paste(Xnm, collapse = " + ") else "1"
  if (is.null(instrument_formula)) {
    instrument_formula <- stats::as.formula(paste("Z ~", rhsX))
  }
  if (is.null(treatment_formula)) {
    treatment_formula <- stats::as.formula(paste("A ~ Z +", rhsX))
  }
  pZ1 <- pA1_z0 <- pA1_z1 <- numeric(n)
  direction_counts <- matrix(0L, length(folds), 2L,
                              dimnames = list(NULL, c("0", "1")))
  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (!length(tr)) tr <- te
    mz <- lsiv_learn_fit(df[tr, , drop = FALSE], "Z", Xnm, stats::binomial(),
                         learner, SL.list, formula = instrument_formula,
                         SL.N = n)
    pZ1[te] <- lsiv_learn_predict(mz, df[te, , drop = FALSE])
    ma <- lsiv_learn_fit(df[tr, , drop = FALSE], "A", c("Z", Xnm),
                         stats::binomial(), learner, SL.list,
                         formula = treatment_formula, SL.N = n)
    d0 <- df[te, , drop = FALSE]
    d1 <- d0
    d0$Z <- 0
    d1$Z <- 1
    pA1_z0[te] <- lsiv_learn_predict(ma, d0)
    pA1_z1[te] <- lsiv_learn_predict(ma, d1)
    training0 <- df[tr, , drop = FALSE]
    training1 <- training0
    training0$Z <- 0
    training1$Z <- 1
    training_zM <- as.integer(
      lsiv_learn_predict(ma, training1) >= lsiv_learn_predict(ma, training0)
    )
    direction_counts[k, ] <- c(sum(training_zM == 0L), sum(training_zM == 1L))
  }
  list(pZ1 = pZ1, pA1_z0 = pA1_z0, pA1_z1 = pA1_z1,
       direction_counts = direction_counts)
}

lsiv_nuisance_default <- function(Y, A, Z, X, type, grid, weights, folds,
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
                                  Aarg = NULL, Zarg = NULL,
                                  direction_min = 100L) {
  args <- mget(names(formals(lsiv_nuisance_default)),
               envir = environment(), inherits = FALSE)
  args$direction_min <- NULL
  if (type != "continuous") {
    return(do.call(lsiv_nuisance_direction, args))
  }
  probabilities <- lsiv_nuisance_probabilities(
    Y, A, Z, X, folds, instrument_formula, treatment_formula, learner, SL.list
  )
  if (!is.numeric(direction_min) || length(direction_min) != 1L ||
      !is.finite(direction_min) || direction_min < 0 ||
      direction_min != floor(direction_min)) {
    stop("`direction_min` must be a nonnegative integer.", call. = FALSE)
  }
  zM <- as.integer(probabilities$pA1_z1 >= probabilities$pA1_z0)
  fit_direction <- zM
  for (k in seq_along(folds)) {
    te <- folds[[k]]
    counts <- probabilities$direction_counts[k, ]
    small <- counts[as.character(zM[te])] < direction_min
    primary <- as.integer(counts[["1"]] >= counts[["0"]])
    fit_direction[te[small]] <- primary
  }
  directions <- c(1L, 0L)[c(any(fit_direction == 1L), any(fit_direction == 0L))]
  args$probabilities <- probabilities
  fits <- lapply(directions, function(numerator_z) {
    args$numerator_z <- numerator_z
    do.call(lsiv_nuisance_direction, args)
  })
  names(fits) <- as.character(directions)
  out <- fits[[1L]]
  if (length(fits) > 1L) {
    for (j in seq_along(fits)[-1L]) {
      ii <- which(fit_direction == directions[j])
      for (component in c("fY_z0", "fY_z1", "RY")) {
        out[[component]][ii, ] <- fits[[j]][[component]][ii, , drop = FALSE]
      }
    }
  }
  out$zM <- zM
  out$ratio_direction <- list(
    min_group = direction_min,
    training_counts = probabilities$direction_counts,
    fit_direction = fit_direction,
    reciprocal_fallback = fit_direction != zM
  )
  out$ratio_ensemble_used <- all(vapply(
    fits, function(fit) isTRUE(fit$ratio_ensemble_used), logical(1L)
  ))
  out$ratio_ensemble_fallback <- if (any(vapply(
    fits, function(fit) !is.null(fit$ratio_ensemble_fallback), logical(1L)
  ))) "logistic" else NULL
  out$ratio_method_by_zM <- vapply(fits, function(fit) {
    if (isTRUE(fit$ratio_ensemble_used)) {
      ratio_method
    } else if (ratio_method != "logistic") {
      paste0("logistic (", ratio_method, " fallback)")
    } else {
      "logistic"
    }
  }, character(1L))
  for (component in c("density_ratio_ensemble", "baseline_density")) {
    summaries <- lapply(seq_along(fits), function(j) {
      result <- fits[[j]][[component]]
      if (is.null(result)) return(NULL)
      result$numerator_z <- directions[j]
      result$denominator_z <- 1L - directions[j]
      result
    })
    names(summaries) <- names(fits)
    available <- which(!vapply(summaries, is.null, logical(1L)))
    if (length(available)) {
      out[[component]] <- summaries[[available[1L]]]
      out[[component]]$by_zM <- summaries
    }
  }
  out
}
