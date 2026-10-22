## ----------------------------------------------------------------------------
## SIV_Sun2021(): the missing-not-at-random (MNAR) estimators of
##   Sun, B., Liu, L., Miao, W., Wirth, K., Robins, J. and Tchetgen Tchetgen,
##   E. J. "Semiparametric estimation with data missing not at random using an
##   instrumental variable."  Statistica Sinica.  (arXiv:1607.03197)
##
## Observed data: (X, Y, R, Z) with Y observed if and only if R = 1.
##   (IV.1) exclusion restriction   P(y | x, z) = P(y | x)
##   (IV.2) IV relevance            pi~(x, z) depends on z
##   eq. (4.1)  selection bias function
##       delta(x,y,z) = log[{P(R=1|x,y,z)/P(R=0|x,y,z)}
##                          / {P(R=1|x,Y=0,z)/P(R=0|x,Y=0,z)}]
##   eq. (4.4)  pi(x,y,z) = 1/[1 + exp{-delta(x,y,z) - eta(x,z)}]
##   eq. (4.5)  f(y|R=0,x,z) = exp{-delta} f(y|R=1,x,z)
##                             / E[exp{-delta} | R=1,x,z]
##   eqs. (4.6)-(4.9)   IPW  (Proposition 1)
##   eqs. (4.10)-(4.11) OR   (Proposition 2)
##   eqs. (4.12)-(4.14) DR   (Proposition 3)
##   eq. (S2) (Supplemental Appendix)  sandwich variance
##
## The implementation mirrors, equation for equation, the reference R code
## printed in the Supplemental Appendix of that paper (pp. S25-S35).
##
## Parameter names used below
##   beta   -> parameters of q(z | x; beta), the IV model            (MLE)
##   alpha  -> parameters of eta(x, z; alpha), the baseline log odds
##   gamma  -> parameters of delta(x, y, z; gamma), the selection bias function
##   eta_y  -> parameters of f(y | R = 1, x, z; eta), complete-case density (MLE)
##   mu     -> the target, E(Y)
## ----------------------------------------------------------------------------

#' Estimating-equation setup for the Sun et al. (2021) MNAR estimators
#'
#' Builds the design matrices, parameter blocks, starting values and the
#' per-observation estimating-function closure used by [SIV_Sun2021()], without
#' solving the equations.  Exported so that eqs. (4.7), (4.8), (4.10) and
#' (4.13) of Sun et al. (2021) can be checked directly.
#'
#' @inheritParams SIV_Sun2021
#' @return A list with the estimating-function closure (`estfun`), starting
#'   values (`par0`), the parameter block index (`blocks`), the parameter
#'   `names` and the design matrices.
#' @examples
#' set.seed(5)
#' n <- 500; X <- cbind(X1 = rnorm(n)); Z <- rbinom(n, 1, plogis(X[, 1]))
#' Y <- rbinom(n, 1, plogis(0.2 * X[, 1])); R <- rbinom(n, 1, plogis(-1 + 2 * Z + Y))
#' Yobs <- ifelse(R == 1, Y, NA_real_)
#' su <- siv_sun2021_setup(Yobs, R, Z, X, estimator = "ipw")
#' dim(su$estfun(su$par0))
#' @export
siv_sun2021_setup <- function(Y, R, Z, X = NULL,
                              estimator = c("dr", "ipw", "or"),
                              outcome_type = c("binary", "continuous"),
                              basis = siv_basis_spec("linear"),
                              eta_formula = NULL, fz_formula = NULL,
                              fy_formula = NULL,
                              h1_formula = NULL,
                              h2_fun = function(Z, X) Z,
                              g_fun = NULL,
                              q1_fun = function(Z, X) Z,
                              v_fun = function(Z, X) Z,
                              support = NULL, n_nodes = 20L) {
  estimator <- match.arg(estimator)
  outcome_type <- match.arg(outcome_type)
  n <- length(Y)
  R <- as.numeric(R); Z <- as.numeric(Z)
  if (!all(R %in% c(0, 1))) stop("`R` must be binary (0/1).", call. = FALSE)
  if (!all(Z %in% c(0, 1))) stop("`Z` must be binary (0/1).", call. = FALSE)
  if (anyNA(Y[R == 1])) stop("`Y` must be observed whenever `R == 1`.", call. = FALSE)
  X <- if (is.null(X)) matrix(numeric(0), n, 0) else as.matrix(X)
  if (is.null(colnames(X)) && ncol(X)) colnames(X) <- paste0("X", seq_len(ncol(X)))
  Yf <- Y; Yf[is.na(Yf)] <- 0          # placeholder; always multiplied by R
  df <- data.frame(Y = Yf, R = R, Z = Z)
  if (ncol(X)) df <- cbind(df, as.data.frame(X))
  Xnm <- colnames(X)
  rhsX <- if (length(Xnm)) paste(Xnm, collapse = " + ") else "1"
  if (is.null(eta_formula)) eta_formula <- stats::as.formula(paste("~ Z +", rhsX))
  if (is.null(fz_formula))  fz_formula  <- stats::as.formula(paste("~", rhsX))
  if (is.null(fy_formula))  fy_formula  <- stats::as.formula(paste("~ Z +", rhsX))
  if (is.null(h1_formula))  h1_formula  <- eta_formula
  if (is.null(support)) support <- sort(unique(Yf[R == 1]))
  family <- if (outcome_type == "binary") "binomial" else "gaussian"

  Deta <- siv_model_matrix(eta_formula, df)     # eta(x, z; alpha)
  Dz   <- siv_model_matrix(fz_formula, df)      # q(z | x; beta)
  Dy   <- siv_model_matrix(fy_formula, df)      # f(y | R = 1, x, z; eta_y)
  H1   <- siv_model_matrix(h1_formula, df)      # eq. (4.7)
  bs   <- siv_basis_fit(basis, Yf[R == 1], X)
  Bobs <- siv_basis_eval(bs, Yf, X)             # only used where R = 1
  q <- ncol(Bobs)
  if (is.null(g_fun)) g_fun <- function(y, X) siv_basis_eval(bs, y, X)

  ## E(h2 | X; beta) and E(q1 | X; beta), E(v | X; beta) with binary Z
  eval_zfun <- function(f) {
    list(obs = as.matrix(f(Z, X)),
         z0  = as.matrix(f(rep(0, n), X)),
         z1  = as.matrix(f(rep(1, n), X)))
  }
  h2 <- eval_zfun(h2_fun); q1 <- eval_zfun(q1_fun); vv <- eval_zfun(v_fun)

  n_be <- ncol(Dz); n_al <- ncol(Deta); n_ey <- siv_ydens_npar(family, Dy, support)
  blocks <- list(beta = seq_len(n_be)); k <- n_be
  blocks$alpha <- k + seq_len(n_al); k <- k + n_al
  blocks$eta_y <- k + seq_len(n_ey); k <- k + n_ey
  blocks$gamma <- k + seq_len(q);    k <- k + q
  blocks$mu    <- k + 1L;            k <- k + 1L
  npar <- k
  pnames <- c(paste0("beta.", colnames(Dz)),
              paste0("alpha.", colnames(Deta)),
              paste0("etaY", seq_len(n_ey)),
              paste0("gamma.", colnames(Bobs)),
              "mu")

  ## eq. (4.4): extended propensity score
  pi_of <- function(par) {
    expit(as.vector(Deta %*% par[blocks$alpha]) + as.vector(Bobs %*% par[blocks$gamma]))
  }
  ## eq. (4.5): E{v(Y) | R = 0, X, Z} under the exponential tilt of the
  ## complete-case density.  `which` selects the integrand.
  m0_of <- function(par, which = c("Y", "basis")) {
    which <- match.arg(which)
    vfun <- if (which == "Y") function(yv) matrix(yv, ncol = 1L) else
      function(yv) g_fun(yv, X)
    tm <- siv_tilt_moments(family, par[blocks$eta_y], Dy, support, bs, X,
                           par[blocks$gamma], vfun, sgn = -1, n_nodes = n_nodes)
    tm$num / tm$den
  }

  estfun <- function(par) {
    pz <- expit(as.vector(Dz %*% par[blocks$beta]))
    sc_beta <- Dz * (Z - pz)                       # MLE score for q(z | x; beta)
    Eh2 <- h2$z1 * pz + h2$z0 * (1 - pz)
    Eq1 <- q1$z1 * pz + q1$z0 * (1 - pz)
    Ev  <- vv$z1 * pz + vv$z0 * (1 - pz)
    pihat <- pi_of(par)
    ## complete-case MLE score for f(y | R = 1, x, z; eta_y)
    sc_ey <- siv_ydens_score(family, par[blocks$eta_y], Dy, Yf, support, R)

    if (estimator == "ipw") {
      ## eq. (4.7)
      e47 <- (R / pihat - 1) * H1
      ## eq. (4.8):  g(X, Y) * {h2(Z, X) - E(h2 | X; beta)}
      gY <- g_fun(Yf, X)
      e48 <- (R / pihat) * gY * as.vector(h2$obs - Eh2)
      ## eq. (4.9)
      e_mu <- R * Yf / pihat - par[blocks$mu]
      cols <- list(sc_beta, e47, sc_ey * 0, e48, cbind(e_mu))
      ## The IPW estimator does not use f(y | R = 1, .); its block is held fixed.
    } else if (estimator == "or") {
      m0Y <- as.vector(m0_of(par, "Y"))
      m0B <- m0_of(par, "basis")
      gY <- g_fun(Yf, X)
      ## eq. (4.10)
      e410 <- as.vector(q1$obs - Eq1) * ((1 - R) * m0B + R * gY)
      ## eq. (4.11)
      e_mu <- R * Yf + (1 - R) * m0Y - par[blocks$mu]
      cols <- list(sc_beta, Deta * 0, sc_ey, e410, cbind(e_mu))
    } else {
      m0Y <- as.vector(m0_of(par, "Y"))
      m0B <- m0_of(par, "basis")
      gY <- g_fun(Yf, X)
      ## eq. (4.7) (shared with IPW: identifies alpha)
      e47 <- (R / pihat - 1) * H1
      ## eq. (4.12): G_DR = R{u - E(u|R=0,X,Z)}/pi + E(u|R=0,X,Z)
      Gb <- (R / pihat) * (gY - m0B) + m0B
      GY <- (R / pihat) * (Yf - m0Y) + m0Y
      ## eq. (4.13)
      e413 <- as.vector(vv$obs - Ev) * Gb
      ## eq. (4.14) with u(X, Y) = Y
      e_mu <- GY - par[blocks$mu]
      cols <- list(sc_beta, e47, sc_ey, e413, cbind(e_mu))
    }
    do.call(cbind, cols)
  }

  free <- rep(TRUE, npar)
  if (estimator == "ipw") free[blocks$eta_y] <- FALSE
  if (estimator == "or")  free[blocks$alpha] <- FALSE

  par0 <- numeric(npar)
  mz <- suppressWarnings(stats::glm.fit(Dz, Z, family = stats::binomial()))
  par0[blocks$beta] <- ifelse(is.finite(mz$coefficients), mz$coefficients, 0)
  mr <- suppressWarnings(stats::glm.fit(Deta, R, family = stats::binomial()))
  par0[blocks$alpha] <- ifelse(is.finite(mr$coefficients), mr$coefficients, 0)
  par0[blocks$eta_y] <- siv_ydens_start(family, Dy, Yf, which(R == 1), support)
  par0[blocks$gamma] <- 0
  par0[blocks$mu] <- mean(Yf[R == 1])
  if (estimator == "or") par0[blocks$alpha] <- 0

  list(estfun = estfun, par0 = par0, blocks = blocks, names = pnames,
       free = free, n = n, estimator = estimator, outcome_type = outcome_type,
       basis = bs, support = support, pi_of = pi_of, m0_of = m0_of,
       Dz = Dz, Deta = Deta, Dy = Dy)
}

#' Sun et al. (2021) MNAR estimators with an instrumental variable
#'
#' Implements the inverse-probability-weighted (IPW), outcome-regression (OR)
#' and doubly robust (DR) estimators of the outcome mean \eqn{\mu_0 = E(Y)}
#' when \eqn{Y} is missing not at random and an instrument \eqn{Z} is available
#' that affects the response mechanism but not the outcome (Sun et al. 2021,
#' Section 4), together with the sandwich covariance of eq. (S2) of that
#' paper's Supplemental Appendix.
#'
#' @section What can be missing, and under what assumptions:
#' Only the outcome \eqn{Y} is subject to missingness; \eqn{(X, Z, R)} are fully
#' observed.  The response indicator is \eqn{R} (`R = 1` when \eqn{Y} is
#' observed).  The identifying assumptions are the exclusion restriction (IV.1),
#' \eqn{P(y \mid x, z) = P(y \mid x)}, IV relevance (IV.2), and a correctly
#' specified selection-bias model.  The missingness is *nonignorable*: the
#' extended propensity score \eqn{\pi(x,y,z) = P(R = 1 \mid x, y, z)} depends on
#' the possibly unobserved \eqn{Y}.
#'
#' @section Estimating equations:
#' With \eqn{\pi(x,y,z;\alpha,\gamma) = \mathrm{expit}\{\eta(x,z;\alpha) +
#' \delta(x,y,z;\gamma)\}} (eq. 4.4) and \eqn{\delta(x,y,z;\gamma) =
#' b(y,x)^\top\gamma} the normalised selection-bias function (eq. 4.1):
#' \describe{
#'   \item{IPW}{eq. (4.7) \eqn{P_n[\{R/\pi - 1\}h_1(X,Z)] = 0} and
#'     eq. (4.8) \eqn{P_n[(R/\pi)g(X,Y)\{h_2(Z,X) - E(h_2\mid X;\hat\beta)\}] = 0};
#'     \eqn{\hat\mu_{IPW} = P_n(RY/\pi)} (eq. 4.9).}
#'   \item{OR}{the complete-case density \eqn{f(y\mid R=1,x,z;\eta)} is fitted by
#'     maximum likelihood on the respondents, the non-respondent density is
#'     obtained by exponential tilting (eq. 4.5), and eq. (4.10)
#'     \eqn{P_n[\{q_1 - E(q_1\mid X;\hat\beta)\}\{(1-R)E(q_2\mid R=0,X,Z) +
#'     Rq_2\}] = 0} is solved for \eqn{\gamma};
#'     \eqn{\hat\mu_{OR} = P_n[RY + (1-R)E(Y\mid R=0,X,Z)]} (eq. 4.11).}
#'   \item{DR}{with \eqn{G_{DR} = (R/\pi)\{u - E(u\mid R=0,X,Z)\} +
#'     E(u\mid R=0,X,Z)} (eq. 4.12), eq. (4.13)
#'     \eqn{P_n[\{v - E(v\mid X;\hat\beta)\}G_{DR}] = 0} is solved jointly with
#'     eq. (4.7); \eqn{\hat\mu_{DR} = P_n G_{DR}} with \eqn{u(X,Y) = Y}
#'     (eq. 4.14).}
#' }
#'
#' @section Variance:
#' All estimating equations -- the MLE scores for \eqn{q(z\mid x;\beta)} and
#' (for OR/DR) \eqn{f(y\mid R=1,x,z;\eta)}, the equations for
#' \eqn{(\alpha,\gamma)} and the equation defining \eqn{\mu} -- are stacked into
#' a single vector \eqn{M(\theta)}, and the covariance is
#' \eqn{\hat A^{-1}\hat B\hat A^{-\top}/n} with
#' \eqn{\hat A = n^{-1}\sum_i \partial M_i/\partial\theta^\top} and
#' \eqn{\hat B = n^{-1}\sum_i M_iM_i^\top}; this is exactly eq. (S2) of the
#' Supplemental Appendix and the `solve(dM) %*% (m %*% t(m)/n) %*% t(solve(dM))/n`
#' construction of the reference code.
#'
#' @param Y Outcome vector; may be `NA` where `R == 0`.
#' @param R Binary response indicator (`1` = observed).
#' @param Z Binary instrument.
#' @param X Matrix or data frame of fully observed covariates.
#' @param estimator One of `"dr"` (default), `"ipw"`, `"or"`.
#' @param outcome_type `"binary"` or `"continuous"`.
#' @param basis A [siv_basis_spec()] describing the selection-bias function
#'   \eqn{\delta}.  The default reproduces \eqn{\delta(x,y,z) = \gamma y}, the
#'   specification used in Section 5 of the paper.
#' @param eta_formula One-sided formula for the baseline log odds
#'   \eqn{\eta(x,z;\alpha)} (default `~ Z + <covariates>`).
#' @param fz_formula One-sided formula for \eqn{q(z\mid x;\beta)}.
#' @param fy_formula One-sided formula for \eqn{f(y\mid R=1,x,z;\eta)}.
#' @param h1_formula Formula giving \eqn{h_1(X,Z)} in eq. (4.7); defaults to
#'   `eta_formula`, which makes eq. (4.7) exactly identify \eqn{\alpha}.
#' @param h2_fun,q1_fun,v_fun Functions of `(Z, X)` giving \eqn{h_2} in
#'   eq. (4.8), \eqn{q_1} in eq. (4.10) and \eqn{v} in eq. (4.13).  All default
#'   to \eqn{Z}, the choice used in the paper's simulation.
#' @param g_fun Function of `(y, X)` giving \eqn{g(X,Y)} in eq. (4.8) and
#'   \eqn{q_2} in eq. (4.10); defaults to the selection-bias basis, which makes
#'   the systems exactly identified.
#' @param support Outcome support for the binomial working model.
#' @param n_nodes Gauss-Hermite nodes for a continuous outcome.
#' @param level Confidence level.
#' @param tol,maxit Control parameters for the estimating-equation solver.
#'
#' @return An object of class `"SIV_Sun2021"` with `estimate` (including `mu`,
#'   the outcome mean, and the selection-bias parameters), the sandwich `vcov`,
#'   `se`, `converged` and fitted nuisance quantities.
#'
#' @references
#' Sun, B., Liu, L., Miao, W., Wirth, K., Robins, J. and Tchetgen Tchetgen,
#' E. J.  Semiparametric estimation with data missing not at random using an
#' instrumental variable.  *Statistica Sinica*, Sections 4.1-4.3 and
#' Supplemental Appendix eq. (S2).
#'
#' @examples
#' set.seed(6)
#' n <- 1000; X <- cbind(X1 = rnorm(n)); Z <- rbinom(n, 1, plogis(X[, 1]))
#' Y <- rbinom(n, 1, plogis(0.2 * X[, 1])); R <- rbinom(n, 1, plogis(-1 + 2 * Z + Y))
#' Yobs <- ifelse(R == 1, Y, NA_real_)
#' fit <- SIV_Sun2021(Yobs, R, Z, X, estimator = "ipw")
#' coef(fit)[c("gamma.y", "mu")]
#' @export
SIV_Sun2021 <- function(Y, R, Z, X = NULL,
                        estimator = c("dr", "ipw", "or"),
                        outcome_type = c("binary", "continuous"),
                        basis = siv_basis_spec("linear"),
                        eta_formula = NULL, fz_formula = NULL, fy_formula = NULL,
                        h1_formula = NULL,
                        h2_fun = function(Z, X) Z,
                        g_fun = NULL,
                        q1_fun = function(Z, X) Z,
                        v_fun = function(Z, X) Z,
                        support = NULL, n_nodes = 20L,
                        level = 0.95, tol = 1e-7, maxit = 200L) {
  cl <- match.call()
  su <- siv_sun2021_setup(Y, R, Z, X, estimator = estimator,
                          outcome_type = outcome_type, basis = basis,
                          eta_formula = eta_formula, fz_formula = fz_formula,
                          fy_formula = fy_formula, h1_formula = h1_formula,
                          h2_fun = h2_fun, g_fun = g_fun, q1_fun = q1_fun,
                          v_fun = v_fun, support = support, n_nodes = n_nodes)
  estfun <- su$estfun; par0 <- su$par0; free <- su$free; pnames <- su$names
  ef_free <- function(pf) { p <- par0; p[free] <- pf; colMeans(estfun(p))[free] }
  sol <- siv_solve_ee(ef_free, par0[free], tol = tol, maxit = maxit)
  par <- par0; par[free] <- sol$par
  sw <- siv_sandwich(function(pf) { p <- par; p[free] <- pf
                                    estfun(p)[, free, drop = FALSE] },
                     par[free], names = pnames[free])
  est <- par[free]; names(est) <- pnames[free]
  V <- sw$vcov
  se <- sqrt(pmax(diag(V), 0))
  q95 <- stats::qnorm(1 - (1 - level) / 2)
  out <- list(
    estimate = est, vcov = V, se = se,
    ci = cbind(lower = est - q95 * se, upper = est + q95 * se),
    level = level, converged = sol$converged, criterion = sol$criterion,
    checks = c(sw$checks, stats::setNames(siv_check_vcov(V, length(est)),
                                          paste0("final_", names(siv_check_vcov(V, length(est)))))),
    bread = sw$bread, meat = sw$meat, jacobian = sw$jacobian,
    estfun = sw$estfun, estfun_fn = estfun, blocks = su$blocks,
    free = free, par_full = par,
    nuisance = list(pi = su$pi_of(par)),
    n = su$n, estimator = su$estimator, outcome_type = su$outcome_type,
    basis = su$basis,
    se_type = "sandwich (stacked estimating equations; Sun et al. 2021, eq. S2)",
    method = sprintf("Sun et al. (2021) MNAR-IV, %s estimator",
                     toupper(su$estimator)),
    target = c(mu = "E(Y)"),
    call = cl
  )
  class(out) <- c("SIV_Sun2021", "siv_fit")
  out
}
