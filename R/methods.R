## ----------------------------------------------------------------------------
## S3 methods shared by all fitted objects in the package.
##
## Every fitted object carries
##   $estimate  named vector of parameter estimates
##   $vcov      covariance matrix conformable with $estimate
##   $se        standard errors
##   $se_type   a human-readable description of *where the standard errors come
##              from*.  This is deliberately explicit: the competing estimators
##              report a genuine sandwich covariance obtained by stacking their
##              estimating equations, while the proposed estimator reports an
##              influence-function based variance (Theorem 5.1).  The two are
##              never conflated in the printed output.
## ----------------------------------------------------------------------------

#' Extract estimated coefficients from a separable-IV fit
#'
#' @param object A fitted object returned by [SIV()], [SIV_Liu()],
#'   [SIV_Sun()] or [SIV_Sun2021()].
#' @param ... Ignored.
#' @return A named numeric vector of estimates.
#' @examples
#' set.seed(1); n <- 300; X <- cbind(X1 = rnorm(n))
#' Z <- rbinom(n, 1, plogis(0.2 * X[, 1])); Y0 <- rbinom(n, 1, 0.5)
#' A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0)); Y <- ifelse(A, 1, Y0)
#' coef(SIV(Y, A, Z, X, y_type = "binary"))
#' @export
coef.siv_fit <- function(object, ...) object$estimate

#' Extract the estimated covariance matrix from a separable-IV fit
#'
#' For [SIV_Liu()], [SIV_Sun()] and [SIV_Sun2021()] this is the sandwich
#' covariance matrix obtained by stacking *all* estimating equations, including
#' those defining the nuisance parameters (Sun et al. 2021, Supplemental
#' Appendix eq. (S2)).  For [SIV()] it is the influence-function based variance
#' of Theorem 5.1 of Park and Tchetgen Tchetgen.
#'
#' @param object A fitted separable-IV object.
#' @param ... Ignored.
#' @return A covariance matrix.
#' @examples
#' set.seed(1); n <- 300; X <- cbind(X1 = rnorm(n))
#' Z <- rbinom(n, 1, plogis(0.2 * X[, 1])); Y0 <- rbinom(n, 1, 0.5)
#' A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0)); Y <- ifelse(A, 1, Y0)
#' vcov(SIV(Y, A, Z, X, y_type = "binary"))
#' @export
vcov.siv_fit <- function(object, ...) object$vcov

#' Confidence intervals for a separable-IV fit
#'
#' @param object A fitted separable-IV object.
#' @param parm Optional character vector of parameter names.
#' @param level Confidence level.
#' @param ... Ignored.
#' @return A matrix of lower/upper confidence limits.
#' @examples
#' set.seed(1); n <- 300; X <- cbind(X1 = rnorm(n))
#' Z <- rbinom(n, 1, plogis(0.2 * X[, 1])); Y0 <- rbinom(n, 1, 0.5)
#' A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0)); Y <- ifelse(A, 1, Y0)
#' confint(SIV(Y, A, Z, X, y_type = "binary"))
#' @export
confint.siv_fit <- function(object, parm, level = 0.95, ...) {
  est <- object$estimate
  if (missing(parm) || is.null(parm)) parm <- names(est)
  bad <- setdiff(parm, names(est))
  if (length(bad)) {
    stop("Unknown parameter(s): ", paste(bad, collapse = ", "), ".",
         call. = FALSE)
  }
  same_level <- !is.null(object$ci_matrix) && !is.null(object$level) &&
    isTRUE(all.equal(as.numeric(level), as.numeric(object$level)))
  if (same_level) {
    out <- object$ci_matrix[parm, , drop = FALSE]
  } else {
    se <- sqrt(pmax(diag(as.matrix(object$vcov)), 0))
    names(se) <- names(est)
    q <- stats::qnorm(1 - (1 - level) / 2)
    out <- cbind(est[parm] - q * se[parm], est[parm] + q * se[parm])
    if (anyNA(out) && !is.null(object$qtt)) {
      warning("A QTT interval at a new confidence level requires refitting ",
              "`SIV()` with that `level`; returning `NA` for those rows.",
              call. = FALSE)
    }
  }
  colnames(out) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                                 trim = TRUE), " %")
  out
}

#' Standard errors of a separable-IV fit
#'
#' @param object A fitted separable-IV object.
#' @param ... Ignored.
#' @return A named numeric vector of standard errors, with an attribute
#'   `"type"` describing their provenance.
#' @examples
#' set.seed(1); n <- 300; X <- cbind(X1 = rnorm(n))
#' Z <- rbinom(n, 1, plogis(0.2 * X[, 1])); Y0 <- rbinom(n, 1, 0.5)
#' A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0)); Y <- ifelse(A, 1, Y0)
#' siv_se(SIV(Y, A, Z, X, y_type = "binary"))
#' @export
siv_se <- function(object, ...) UseMethod("siv_se")

#' @rdname siv_se
#' @export
siv_se.siv_fit <- function(object, ...) {
  se <- sqrt(pmax(diag(as.matrix(object$vcov)), 0))
  names(se) <- names(object$estimate)
  attr(se, "type") <- object$se_type
  se
}

#' Summarise a separable-IV fit
#'
#' @param object A fitted separable-IV object.
#' @param parm Optional character vector selecting the parameters to report.
#'   Defaults to the target parameters of the method.
#' @param ... Ignored.
#' @return An object of class `"summary.siv_fit"`.
#' @examples
#' set.seed(1); n <- 300; X <- cbind(X1 = rnorm(n))
#' Z <- rbinom(n, 1, plogis(0.2 * X[, 1])); Y0 <- rbinom(n, 1, 0.5)
#' A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0)); Y <- ifelse(A, 1, Y0)
#' summary(SIV(Y, A, Z, X, y_type = "binary"))
#' @export
summary.siv_fit <- function(object, parm = NULL, ...) {
  est <- object$estimate
  se <- sqrt(pmax(diag(as.matrix(object$vcov)), 0))
  names(se) <- names(est)
  if (is.null(parm)) parm <- names(est)
  parm <- intersect(parm, names(est))
  z <- est[parm] / se[parm]
  tab <- cbind(Estimate = est[parm], `Std. Error` = se[parm],
               `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  ci <- confint.siv_fit(object, parm, level = object$level %||% 0.95)
  out <- list(method = object$method, se_type = object$se_type,
              coefficients = tab, ci = ci, n = object$n,
              converged = object$converged,
              checks = object$checks, call = object$call,
              type = object$type, estimator = object$estimator,
              qtt = object$qtt, falsification = object$falsification)
  class(out) <- "summary.siv_fit"
  out
}

#' @export
print.summary.siv_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n", x$method, "\n", sep = "")
  if (!is.null(x$type)) cat("Outcome type: ", x$type, "\n", sep = "")
  if (!is.null(x$estimator)) cat("Estimator:    ", x$estimator, "\n", sep = "")
  cat("Sample size:  ", x$n, "\n", sep = "")
  cat("Std. errors:  ", x$se_type, "\n", sep = "")
  if (!is.null(x$converged)) {
    cat("Converged:    ", isTRUE(x$converged), "\n", sep = "")
  }
  cat("\n")
  stats::printCoefmat(x$coefficients, digits = digits, has.Pvalue = TRUE,
                      P.values = TRUE)
  if (any(!is.finite(x$coefficients[, "Std. Error"])) && !is.null(x$ci)) {
    cat("\nConfidence intervals (QTT uses Berger--Boos inversion):\n")
    print(round(x$ci, digits))
  }
  if (!is.null(x$falsification)) {
    ft <- x$falsification
    cat(sprintf(paste0("\nFalsification test: T = %.*f, critical value = %.*f, ",
                       "p-value = %.*f, reject = %s\n"),
                digits, ft$statistic, digits, ft$critical_value,
                digits, ft$p_value, if (isTRUE(ft$reject)) "yes" else "no"))
    cat(sprintf("  dictionary: %d moments over %d %s X cells x %d Y bins\n",
                ft$n_dictionary, ft$n_x_cells_observed,
                ft$x_cells %||% "joint", ft$n_y_bins))
  }
  if (!is.null(x$checks)) {
    bad <- names(x$checks)[!x$checks]
    if (length(bad)) {
      cat("\nWARNING: failed numerical checks on the covariance matrix: ",
          paste(bad, collapse = ", "), "\n", sep = "")
    }
  }
  cat("\n")
  invisible(x)
}

#' @export
print.siv_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n", x$method, "\n\n", sep = "")
  est <- x$estimate
  se <- sqrt(pmax(diag(as.matrix(x$vcov)), 0))
  keep <- if (!is.null(x$target)) intersect(names(x$target), names(est)) else names(est)
  if (!length(keep)) keep <- names(est)
  tab <- cbind(Estimate = est[keep], SE = se[keep])
  if (!is.null(x$ci_matrix)) {
    tab <- cbind(tab, Lower = x$ci_matrix[keep, "lower"],
                 Upper = x$ci_matrix[keep, "upper"])
  }
  print(round(tab, digits))
  cat("\n(", x$se_type, ")\n", sep = "")
  if (!is.null(x$falsification)) {
    ft <- x$falsification
    cat(sprintf("Falsification: T = %.*f, p-value = %.*f, reject = %s\n",
                digits, ft$statistic, digits, ft$p_value,
                if (isTRUE(ft$reject)) "yes" else "no"))
  }
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Numerical diagnostics of an estimated covariance matrix
#'
#' Verifies that a covariance matrix has the expected dimension, contains only
#' finite entries, has non-negative diagonal entries up to a tolerance and is
#' symmetric up to a tolerance.  These are the checks required of every
#' sandwich covariance produced by this package.
#'
#' @param V A numeric matrix.
#' @param p Expected dimension.
#' @param tol Numerical tolerance.
#' @return A named logical vector.
#' @examples
#' siv_check_vcov(diag(2), p = 2)
#' @export
siv_check_vcov <- function(V, p = NULL, tol = 1e-8) {
  V <- as.matrix(V)
  c(dim_ok      = is.null(p) || identical(dim(V), c(as.integer(p), as.integer(p))),
    finite      = all(is.finite(V)),
    diag_nonneg = all(diag(V) > -tol),
    symmetric   = max(abs(V - t(V))) < tol)
}
