## ----------------------------------------------------------------------------
## SIV_LATE(): the efficient-influence-function (doubly robust / Neyman
## orthogonal) estimator of the LOCAL average treatment effect.
##
## This is a *different estimand* from the rest of the package.  SIV(),
## SIV_Liu() and SIV_Sun() all target the ATT  E(Y^{(1)} - Y^{(0)} | A = 1)
## under the logit-separable treatment choice model.  SIV_LATE() targets the
## complier average causal effect under the standard LATE assumptions
## (Imbens & Angrist 1994; Abadie 2003) -- random assignment of Z given X,
## exclusion, first-stage relevance and monotonicity -- and makes no separability
## assumption at all.  It is included as a reference point, not as a competitor
## for the same target.
##
## Derivation of the influence function
## ------------------------------------
## Write
##   pi(X)     = Pr(Z = 1 | X)
##   m_Y(z, X) = E(Y | Z = z, X)      m_A(z, X) = E(A | Z = z, X)
##   DY(X)     = m_Y(1,X) - m_Y(0,X)  DA(X)     = m_A(1,X) - m_A(0,X)
##   tau_Y     = E{DY(X)}             tau_A     = E{DA(X)}
##   LATE      = tau_Y / tau_A .
##
## tau_Y is an ATE-type functional with Z in the role of the treatment, so its
## efficient influence function is the usual AIPW/doubly robust one:
##
##   phi_Y = Z{Y - m_Y(1,X)}/pi(X) - (1-Z){Y - m_Y(0,X)}/{1 - pi(X)}
##           + DY(X) - tau_Y .
##
## Because Z{Y - m_Y(1,X)} = Z{Y - m_Y(Z,X)} and likewise for the (1-Z) term,
## the two inverse-probability terms collapse to a single expression,
##
##   phi_Y = (2Z - 1){Y - m_Y(Z,X)} / Pr(Z | X) + DY(X) - tau_Y ,          (*)
##
## in which -- and this is the one place the algebra is easy to get wrong --
## `Pr(Z | X)` is the probability of the OBSERVED value of Z,
##
##   Pr(Z | X) = pi(X)^Z {1 - pi(X)}^{1-Z} ,
##
## NOT pi(X).  With Z = 1 the factor is (+1)/pi(X) and with Z = 0 it is
## (-1)/{1 - pi(X)}, which is exactly the pair of AIPW weights.  Using pi(X) in
## both arms is a genuine bias: on a test DGP with pi(X) = expit(1.8 X) it moved
## the estimate by -0.106 against a Monte Carlo error of 0.004.
##
## phi_A is (*) with A in place of Y.  LATE is the ratio tau_Y / tau_A, so by
## the delta method its efficient influence function is
##
##   phi = {phi_Y - LATE * phi_A} / tau_A ,
##
## and since tau_Y - LATE * tau_A = 0 the two centring constants cancel, leaving
## the form implemented below.  This is the score of the "interactive IV" model
## of Chernozhukov et al. (2018, Econometrics Journal, Section 5.2); it is
## Neyman orthogonal, so first-order errors in the estimated nuisance functions
## do not contribute to the first-order bias of the estimator.
##
## Verified numerically (see tests/testthat/test-late.R):
##   * the Wald estimand E{DY}/E{DA} equals the complier average causal effect
##     under a latent-index DGP satisfying monotonicity (agreement 1.3e-05);
##   * with oracle nuisances the influence function has mean exactly 0 and its
##     standard error matches the Monte Carlo standard deviation of the
##     estimator (ratio 1.013, 95% coverage 0.955 over 400 replications);
##   * perturbing every nuisance function by O(eps) moves the estimate by
##     O(eps^2), while the corresponding plug-in ratio moves by O(eps).
## ----------------------------------------------------------------------------

## Internal: assemble the LATE influence function from nuisance estimates.
## Exposed separately so that the algebra can be tested against oracle
## nuisance functions without going through the estimation pipeline.
siv_late_eif <- function(Y, A, Z, pZ1, mY_z0, mY_z1, mA_z0, mA_z1) {
  n <- length(Y)
  ## Pr(Z = observed z | X); see the derivation above.
  prZ <- ifelse(Z == 1, pZ1, 1 - pZ1)
  sgn <- 2 * Z - 1
  mY_obs <- ifelse(Z == 1, mY_z1, mY_z0)
  mA_obs <- ifelse(Z == 1, mA_z1, mA_z0)

  ## Uncentred scores: psi_Y has mean tau_Y and psi_A has mean tau_A.
  psi_Y <- sgn * (Y - mY_obs) / prZ + (mY_z1 - mY_z0)
  psi_A <- sgn * (A - mA_obs) / prZ + (mA_z1 - mA_z0)

  tau_Y <- mean(psi_Y)
  tau_A <- mean(psi_A)
  late <- tau_Y / tau_A
  eif <- (psi_Y - late * psi_A) / tau_A

  list(late = late, eif = eif, var = stats::var(eif) / n,
       psi_Y = psi_Y, psi_A = psi_A, tau_Y = tau_Y, tau_A = tau_A,
       reduced_form = tau_Y, first_stage = tau_A)
}

## Internal: cross-fitted nuisance estimation for SIV_LATE().
siv_late_nuisance <- function(Y, A, Z, X, folds, learner, SL.list,
                              outcome_family,
                              instrument_formula, treatment_formula,
                              outcome_formula, trim) {
  n <- length(Y)
  df <- siv_frame(Y, A, Z, X)
  Xnm <- colnames(X)
  rhsX <- if (length(Xnm)) paste(Xnm, collapse = " + ") else "1"
  if (is.null(instrument_formula)) instrument_formula <- stats::as.formula(paste("Z ~", rhsX))
  if (is.null(treatment_formula)) treatment_formula <- stats::as.formula(paste("A ~ Z +", rhsX))
  if (is.null(outcome_formula))   outcome_formula   <- stats::as.formula(paste("Y ~ Z +", rhsX))

  pZ1 <- mA_z0 <- mA_z1 <- mY_z0 <- mY_z1 <- numeric(n)
  sl_info <- vector("list", length(folds))

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (length(tr) == 0L) tr <- te          # crossfit = 1 -> in-sample
    dtr <- df[tr, , drop = FALSE]
    dte <- df[te, , drop = FALSE]
    d0 <- dte; d0$Z <- 0
    d1 <- dte; d1$Z <- 1

    ## Pr(Z = 1 | X)
    fz <- siv_learn_fit(dtr, "Z", Xnm, stats::binomial(), learner, SL.list,
                        formula = instrument_formula, SL.N = n)
    pZ1[te] <- siv_learn_predict(fz, dte)

    ## E(A | Z, X)
    fa <- siv_learn_fit(dtr, "A", c("Z", Xnm), stats::binomial(), learner,
                        SL.list, formula = treatment_formula, SL.N = n)
    mA_z0[te] <- siv_learn_predict(fa, d0)
    mA_z1[te] <- siv_learn_predict(fa, d1)

    ## E(Y | Z, X)
    fy <- siv_learn_fit(dtr, "Y", c("Z", Xnm), outcome_family, learner,
                        SL.list, formula = outcome_formula, SL.N = n)
    mY_z0[te] <- siv_learn_predict(fy, d0)
    mY_z1[te] <- siv_learn_predict(fy, d1)

    if (learner == "SL") {
      sl_info[[k]] <- list(
        pZ1 = if (inherits(fz, "siv_learner_sl")) fz$fit$coef else NULL,
        mA  = if (inherits(fa, "siv_learner_sl")) fa$fit$coef else NULL,
        mY  = if (inherits(fy, "siv_learner_sl")) fy$fit$coef else NULL)
    }
  }

  ## Numerical safeguards.  The instrument propensity score sits in a
  ## denominator, so it is bounded away from 0 and 1 (this is the `bound_prob()`
  ## guard of the UDID package, exposed here as `trim`).  E(A | Z, X) is a
  ## probability too, but it is never inverted, so it is only clipped to [0, 1]
  ## to keep the first stage interpretable.
  n_trim <- sum(pZ1 < trim | pZ1 > 1 - trim, na.rm = TRUE)
  pZ1 <- siv_bound_prob(pZ1, tol = trim)
  mA_z0 <- pmin(pmax(mA_z0, 0), 1)
  mA_z1 <- pmin(pmax(mA_z1, 0), 1)

  list(pZ1 = pZ1, mA_z0 = mA_z0, mA_z1 = mA_z1,
       mY_z0 = mY_z0, mY_z1 = mY_z1,
       n_trimmed = n_trim, trim = trim,
       sl_info = if (learner == "SL") sl_info else NULL)
}

#' Efficient influence function estimator of the LATE
#'
#' Estimates the local average treatment effect (LATE, the complier average
#' causal effect) with the efficient-influence-function / Neyman-orthogonal
#' estimator, using cross-fitted nuisance functions.
#'
#' @section Estimand:
#' With \eqn{\pi(X) = \Pr(Z = 1 \mid X)},
#' \eqn{m_Y(z,X) = E(Y \mid Z = z, X)} and \eqn{m_A(z,X) = E(A \mid Z = z, X)},
#' write \eqn{\Delta_Y(X) = m_Y(1,X) - m_Y(0,X)} and
#' \eqn{\Delta_A(X) = m_A(1,X) - m_A(0,X)}.  The target is the covariate-adjusted
#' Wald ratio
#' \deqn{\mathrm{LATE} = E\{\Delta_Y(X)\} / E\{\Delta_A(X)\},}
#' which equals \eqn{E(Y^{(1)} - Y^{(0)} \mid \mathrm{complier})} under the
#' assumptions of Imbens and Angrist (1994) and Abadie (2003).
#'
#' \NTC{This is not the estimand of \code{SIV()}, \code{SIV_Liu()} or
#' \code{SIV_Sun()}}, which
#' all target the ATT \eqn{E(Y^{(1)} - Y^{(0)} \mid A = 1)} under the
#' logit-separable treatment choice model.  `SIV_LATE()` assumes monotonicity
#' instead of separability and identifies a complier-weighted effect.  The two
#' coincide only in special cases, so a bias comparison against the ATT is not
#' meaningful; compare each estimator against its own estimand.
#'
#' @section Influence function:
#' The estimator solves \eqn{P_n \psi(O;\theta,\eta) = 0} with
#' \deqn{\psi = \frac{(2Z-1)\{Y - m_Y(Z,X)\}}{\Pr(Z \mid X)} + \Delta_Y(X)
#'       - \theta\left[\frac{(2Z-1)\{A - m_A(Z,X)\}}{\Pr(Z \mid X)}
#'       + \Delta_A(X)\right],}
#' where \eqn{\Pr(Z \mid X) = \pi(X)^Z\{1-\pi(X)\}^{1-Z}} is the probability of
#' the *observed* value of the instrument -- with \eqn{Z = 1} the leading factor
#' is \eqn{1/\pi(X)} and with \eqn{Z = 0} it is \eqn{-1/\{1-\pi(X)\}}, i.e. the
#' two augmented inverse-probability weights.  Substituting \eqn{\pi(X)} into
#' both arms is a common slip and produces a genuinely biased estimator.
#' The efficient influence function is
#' \eqn{\phi = \psi / E\{\Delta_A(X)\}} and the reported standard error is
#' \eqn{\{n^{-1}\widehat{\mathrm{Var}}(\phi)\}^{1/2}}.
#'
#' @param Y Numeric outcome vector, observed for every unit.
#' @param A Binary treatment indicator (1 = treated).
#' @param Z Binary instrument.
#' @param X Matrix or data frame of covariates (may be `NULL`).
#' @param learner Nuisance estimator: `"glm"` (parametric working models, the
#'   package default) or `"SL"` (the Super Learner ensemble of [MySL()]).
#' @param SL.list Integer vector of learner groups passed to [MySL()] when
#'   `learner = "SL"`; see that function for the group numbering.  Groups whose
#'   packages are not installed are dropped with a warning.
#' @param crossfit Number of cross-fitting folds.  `crossfit = 1` (the default,
#'   matching [SIV()]) fits the nuisance models on the full sample.  Folds are
#'   stratified on `(A, Z)`, exactly as in [SIV()].
#' @param folds Optional list of integer vectors giving the cross-fitting folds.
#' @param outcome_family Family for the \eqn{E(Y \mid Z, X)} model.  Defaults to
#'   `binomial()` when `Y` is 0/1 and `gaussian()` otherwise.
#' @param instrument_formula,treatment_formula,outcome_formula Optional formulas
#'   for the `learner = "glm"` working models (`Z ~ .`, `A ~ Z + .`,
#'   `Y ~ Z + .`).  Ignored when `learner = "SL"`, which addresses covariates by
#'   column rather than by formula.
#' @param nuisance Optional named list of pre-computed nuisance estimates with
#'   elements `pZ1`, `mA_z0`, `mA_z1`, `mY_z0`, `mY_z1`, each of length `n`.
#'   Supplying this bypasses estimation entirely, which is how the oracle checks
#'   in the test suite are run.
#' @param trim Lower bound applied to \eqn{\pi(X)}, which is mapped into
#'   `[trim, 1 - trim]` by [siv_bound_prob()].  The instrument propensity score
#'   appears in a denominator, so this guard is what keeps the influence
#'   function finite when the instrument is nearly deterministic given `X`.
#' @param level Confidence level for the reported interval.
#'
#' @return An object of class `"SIV_LATE"` inheriting from `"siv_fit"`, with
#'   components `estimate` (the LATE), `vcov`, `se`, `ci`, `eif`, `nuisance`,
#'   `reduced_form` (\eqn{\hat\tau_Y}), `first_stage` (\eqn{\hat\tau_A}) and
#'   `diagnostics`.  `diagnostics` reports the estimated first stage, the
#'   number of trimmed propensity scores and the range of \eqn{\hat\pi(X)}; a
#'   weak first stage is warned about, because it is the denominator of the
#'   influence function.
#'   Use [coef()], [vcov()], [siv_se()], [confint()] and [summary()] to extract
#'   results, exactly as for the other estimators in the package.
#'
#' @references
#' Imbens, G. W. and Angrist, J. D. (1994). Identification and estimation of
#' local average treatment effects. *Econometrica* 62, 467-475.
#'
#' Abadie, A. (2003). Semiparametric instrumental variable estimation of
#' treatment response models. *Journal of Econometrics* 113, 231-263.
#'
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
#' Newey, W. and Robins, J. (2018). Double/debiased machine learning for
#' treatment and structural parameters. *The Econometrics Journal* 21, C1-C68.
#' Section 5.2 ("interactive IV model").
#'
#' @seealso [SIV()] for the ATT under the logit-separable model; [MySL()] for
#'   the Super Learner ensemble used when `learner = "SL"`.
#'
#' @examples
#' set.seed(2024)
#' n <- 500; X <- cbind(X1 = rnorm(n)); Z <- rbinom(n, 1, plogis(0.2 * X[, 1]))
#' U <- runif(n); p0 <- plogis(-0.5 + 0.2 * X[, 1]); p1 <- plogis(1 + 0.2 * X[, 1])
#' A <- ifelse(Z == 1, U < p1, U < p0); Y <- X[, 1] + A + rnorm(n)
#' fit <- SIV_LATE(Y, A, Z, X)
#' coef(fit)
#' siv_se(fit)
#' @export
SIV_LATE <- function(Y, A, Z, X = NULL,
                     learner = c("glm", "SL"),
                     SL.list = 1L,
                     crossfit = 1L, folds = NULL,
                     outcome_family = NULL,
                     instrument_formula = NULL,
                     treatment_formula = NULL,
                     outcome_formula = NULL,
                     nuisance = NULL,
                     trim = 1e-2,
                     level = 0.95) {
  cl <- match.call()
  learner <- match.arg(learner)
  ## Unlike the ATT estimators, LATE needs Y for treated units as well: it is a
  ## contrast of E(Y | Z, X), not of the treatment-free potential outcome.
  dat <- siv_check_data(Y, A, Z, X, allow_na_Y = FALSE)
  Y <- dat$Y; A <- dat$A; Z <- dat$Z; X <- dat$X; n <- dat$n
  if (!is.numeric(trim) || length(trim) != 1L || trim < 0 || trim >= 0.5) {
    stop("`trim` must be a single number in [0, 0.5).", call. = FALSE)
  }
  if (is.null(outcome_family)) {
    outcome_family <- if (all(Y %in% c(0, 1))) stats::binomial() else stats::gaussian()
  }

  if (is.null(folds)) folds <- siv_make_folds(n, A, Z, crossfit)

  if (is.null(nuisance)) {
    nuisance <- siv_late_nuisance(Y, A, Z, X, folds, learner, SL.list,
                                  outcome_family, instrument_formula,
                                  treatment_formula, outcome_formula, trim)
  } else {
    need <- c("pZ1", "mA_z0", "mA_z1", "mY_z0", "mY_z1")
    miss <- setdiff(need, names(nuisance))
    if (length(miss)) {
      stop("`nuisance` is missing element(s): ", paste(miss, collapse = ", "),
           ".", call. = FALSE)
    }
    bad <- need[vapply(nuisance[need], function(v) length(v) != n, logical(1))]
    if (length(bad)) {
      stop(sprintf("`nuisance` element(s) %s must have length %d.",
                   paste(bad, collapse = ", "), n), call. = FALSE)
    }
    nuisance$n_trimmed <- 0L
    nuisance$trim <- 0
  }

  ef <- siv_late_eif(Y, A, Z, nuisance$pZ1, nuisance$mY_z0, nuisance$mY_z1,
                     nuisance$mA_z0, nuisance$mA_z1)

  ## The first stage is the denominator of the influence function, so a weak one
  ## is surfaced rather than left to produce an implausible estimate.  This is
  ## the LATE analogue of the Assumption (A5) relevance diagnostic of SIV().
  if (!is.finite(ef$tau_A) || abs(ef$tau_A) < 1e-3) {
    warning(sprintf(paste0("The estimated first stage E{E(A|Z=1,X) - ",
                           "E(A|Z=0,X)} = %.2e is close to zero: the ",
                           "instrument is weak and the LATE influence ",
                           "function is numerically unstable."), ef$tau_A),
            call. = FALSE)
  }
  if (nuisance$n_trimmed > 0L) {
    warning(sprintf(paste0("%d of %d estimated instrument propensity scores ",
                           "were outside [%g, %g] and have been trimmed."),
                    nuisance$n_trimmed, n, nuisance$trim, 1 - nuisance$trim),
            call. = FALSE)
  }

  est <- ef$late
  V <- matrix(ef$var, 1L, 1L, dimnames = list("LATE", "LATE"))
  se <- sqrt(V[1L, 1L])
  q <- stats::qnorm(1 - (1 - level) / 2)
  out <- list(
    estimate  = c(LATE = est),
    vcov      = V,
    se        = c(LATE = se),
    ci        = c(lower = est - q * se, upper = est + q * se),
    level     = level,
    eif       = ef$eif,
    nuisance  = nuisance,
    reduced_form = ef$tau_Y,
    first_stage  = ef$tau_A,
    influence_detail = ef,
    converged = all(is.finite(ef$eif)) && is.finite(est),
    diagnostics = c(first_stage = ef$tau_A,
                    reduced_form = ef$tau_Y,
                    n_trimmed = nuisance$n_trimmed,
                    min_pZ1 = min(nuisance$pZ1),
                    max_pZ1 = max(nuisance$pZ1)),
    n         = n,
    learner   = learner,
    crossfit  = length(folds),
    se_type   = "influence function (efficient, Neyman-orthogonal)",
    method    = "Efficient influence function estimator of the LATE",
    call      = cl
  )
  class(out) <- c("SIV_LATE", "siv_fit")
  out
}
