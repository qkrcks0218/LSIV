## ----------------------------------------------------------------------------
## Shared nuisance-estimation plumbing: cross-fitting folds and a thin
## fit/predict abstraction over the two supported learners.
##
##   learner = "glm"  the SIV_LATE() default. Parametric working models specified
##                    by the user's formulas (`instrument_formula`,
##                    `treatment_formula`, `outcome_formula`), fitted with
##                    stats::glm / stats::lm.  This is the behaviour the package
##                    has always had and is unchanged.
##   learner = "SL"   the Super Learner ensemble of MySL(), ported from the UDID
##                    package.  Formulas are not used, because MySL() takes
##                    column indices rather than a formula: the covariate set is
##                    the columns of `X`, plus `Z` for the models that condition
##                    on the instrument.
##
## Both `SIV()` and `SIV_LATE()` go through this file. `SIV()` always selects
## Super Learner, while `SIV_LATE()` retains the learner switch.
## ----------------------------------------------------------------------------

## Internal: cross-fitting folds shared by SIV() and SIV_LATE().  Two-fold
## splitting reproduces Data_Version1 exactly; larger fold counts remain
## stratified on (A, Z) so that every fold contains all four cells.
siv_make_folds <- function(n, A, Z, crossfit) {
  crossfit <- max(1L, as.integer(crossfit))
  if (crossfit == 1L) return(list(seq_len(n)))
  if (crossfit == 2L) {
    ## Match the two-fold split in Data_Version1 exactly.  The four cells are
    ## sampled in the order (A,Z) = (0,0), (0,1), (1,0), (1,1).  The first
    ## fold receives floor(n_cell / 2) observations when Z = 0 and
    ## ceiling(n_cell / 2) observations when Z = 1; the complements form the
    ## second fold.
    cells <- list(
      which(A == 0 & Z == 0),
      which(A == 0 & Z == 1),
      which(A == 1 & Z == 0),
      which(A == 1 & Z == 1)
    )
    first <- second <- vector("list", length(cells))
    for (j in seq_along(cells)) {
      ## sample(x) treats a length-one numeric x as the range 1:x.
      permuted <- cells[[j]][sample.int(length(cells[[j]]))]
      n_first <- if (j %in% c(1L, 3L)) {
        floor(length(permuted) / 2)
      } else {
        ceiling(length(permuted) / 2)
      }
      first[[j]] <- permuted[seq_len(n_first)]
      second[[j]] <- setdiff(permuted, first[[j]])
    }
    return(list(unlist(first, use.names = FALSE),
                unlist(second, use.names = FALSE)))
  }
  strat <- interaction(A, Z, drop = TRUE)
  idx <- split(seq_len(n), strat)
  assign <- integer(n)
  for (s in idx) assign[s] <- sample(rep_len(seq_len(crossfit), length(s)))
  split(seq_len(n), assign)
}

## Internal: fit one conditional-mean / conditional-probability model.
##
## `data`     training data frame,
## `yname`    name of the response column,
## `xnames`   names of the predictor columns,
## `family`   stats::gaussian() or stats::binomial(),
## `learner`  "glm" or "SL",
## `formula`  used only by the "glm" learner; if NULL, `yname ~ xnames` is used.
##
## Degenerate responses (constant y) are handled here rather than inside the
## learners: every caller wants the constant back, not an error.
siv_learn_fit <- function(data, yname, xnames, family, learner = c("glm", "SL"),
                          SL.list = 1L, formula = NULL, SL.N = nrow(data)) {
  learner <- match.arg(learner)
  y <- data[[yname]]
  if (!length(xnames) || length(unique(y[is.finite(y)])) <= 1L) {
    return(structure(list(constant = if (length(y)) mean(y, na.rm = TRUE) else 0),
                     class = "siv_learner_const"))
  }
  if (learner == "glm") {
    if (is.null(formula)) {
      rhs <- if (length(xnames)) paste(xnames, collapse = " + ") else "1"
      formula <- stats::as.formula(paste(yname, "~", rhs))
    }
    fit <- if (identical(family$family, "gaussian")) {
      stats::lm(formula, data = data)
    } else {
      suppressWarnings(stats::glm(formula, data = data, family = family))
    }
    return(structure(list(fit = fit, family = family),
                     class = "siv_learner_glm"))
  }
  ## Super Learner.  MySL() addresses columns positionally, so hand it a frame
  ## whose first column is the response and whose remaining columns are the
  ## predictors, in a fixed order.
  d <- data[, c(yname, xnames), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  fit <- MySL(d, locY = 1L, locX = 1L + seq_along(xnames),
              Ydist = family, SL.list = SL.list, SL.N = SL.N)
  structure(list(fit = fit, xnames = xnames, family = family),
            class = "siv_learner_sl")
}

## Internal: predict from a fitted `siv_learn_fit()` object on `newdata`.
## Always returns a plain numeric vector on the response scale.
siv_learn_predict <- function(object, newdata) {
  n <- nrow(newdata)
  out <- switch(
    class(object)[1L],
    siv_learner_const = rep(object$constant, n),
    siv_learner_glm   = if (identical(object$family$family, "gaussian")) {
      as.numeric(stats::predict(object$fit, newdata = newdata))
    } else {
      as.numeric(stats::predict(object$fit, newdata = newdata, type = "response"))
    },
    siv_learner_sl = as.numeric(
      stats::predict(object$fit,
                     newdata = newdata[, object$xnames, drop = FALSE],
                     onlySL = TRUE)$pred),
    stop("unknown learner object", call. = FALSE)
  )
  out
}

## Internal: residual standard deviation of a fitted conditional-mean model,
## used by the Gaussian conditional-density working model.  For a Super Learner
## fit there is no `sigma()` method, so it is computed from the in-sample
## residuals of the ensemble prediction.
siv_learn_sigma <- function(object, data, yname) {
  if (inherits(object, "siv_learner_glm") &&
      identical(object$family$family, "gaussian")) {
    return(stats::sigma(object$fit))
  }
  r <- data[[yname]] - siv_learn_predict(object, data)
  r <- r[is.finite(r)]
  if (!length(r)) return(1)
  max(stats::sd(r), 1e-8)
}
