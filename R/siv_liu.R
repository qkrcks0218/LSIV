## ----------------------------------------------------------------------------
## SIV_Liu(): the estimators of
##   Liu, L., Miao, W., Sun, B., Robins, J. and Tchetgen Tchetgen, E. (2021+).
##   "Identification and inference for marginal average treatment effect on the
##    treated with an instrumental variable."  (arXiv:1506.08149)
##
## Section 4 of that paper gives three estimators of  theta = E(Y0 | A = 1):
##   Section 4.1  IPW, eqs. (4.1)-(4.5) and Proposition 2;
##   Section 4.2  OR,  eqs. (4.6)-(4.8) and Proposition 3;
##   Section 4.2  DR,  eqs. (4.11)-(4.12) and Proposition 4.
##
## Notation of the paper -> notation here
##   pi(Y0, Z, C; theta)  extended propensity score Pr(A = 1 | Y0, Z, C)
##   psi(Y0, Z, C; psi)   selection bias / generalized odds ratio function
##   f_{Z|C}(z|c; alpha)  instrument model
##   f(Y | A = 0, Z, C; eta)  outcome density among the untreated
##   theta = E(Y0 | A = 1);   ETT = E(Y | A = 1) - theta
##
## Correspondence with Park & Tchetgen Tchetgen: psi = log delta (odds-ratio
## function, their eq. 2) and the baseline part of logit pi is log gamma.
## ----------------------------------------------------------------------------

#' Estimating-equation setup for the Liu et al. (2021) estimators
#'
#' Builds the design matrices, parameter blocks, starting values and the
#' per-observation estimating-function closure used by [SIV_Liu()], without
#' solving the equations.  Exported so that the estimating equations of
#' Liu et al. (2021), eqs. (4.2)-(4.5), (4.7) and (4.11)-(4.12), can be checked
#' directly -- for instance by verifying that their expectation vanishes at the
#' true parameter values.
#'
#' @inheritParams SIV_Liu
#' @return A list with the estimating-function closure (`estfun`), starting
#'   values (`par0`), the parameter block index (`blocks`), the parameter
#'   `names`, the `free`-parameter indicator, and the design matrices.
#' @examples
#' set.seed(4)
#' n <- 300; X <- cbind(X1 = rnorm(n)); Z <- rbinom(n, 1, plogis(0.2 * X[, 1]))
#' Y0 <- rbinom(n, 1, 0.5); A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0))
#' Y <- ifelse(A, 1, Y0)
#' su <- siv_liu_setup(Y, A, Z, X, estimator = "ipw",
#'                     outcome_type = "binary")
#' dim(su$estfun(su$par0))
#' @export
siv_liu_setup <- function(Y, A, Z, X = NULL,
                    estimator = c("dr", "ipw", "or"),
                    outcome_type = c("binary", "continuous"),
                    basis = siv_basis_spec("linear"),
                    pi_formula = NULL, fz_formula = NULL, fy_formula = NULL,
                    l_fun = function(Z, X) Z,
                    w_fun = function(Z, X) Z,
                    support = NULL, n_nodes = 20L,
                    start = NULL, level = 0.95,
                    tol = 1e-7, maxit = 200L) {
  estimator <- match.arg(estimator)
  outcome_type <- match.arg(outcome_type)
  dat <- siv_check_data(Y, A, Z, X)
  Y <- dat$Y; A <- dat$A; Z <- dat$Z; X <- dat$X; n <- dat$n
  siv_match_type(outcome_type, Y, A)
  df <- siv_frame(Y, A, Z, X)
  Xnm <- colnames(X)
  rhsX <- if (length(Xnm)) paste(Xnm, collapse = " + ") else "1"
  if (is.null(pi_formula)) pi_formula <- stats::as.formula(paste("~ Z +", rhsX))
  if (is.null(fz_formula)) fz_formula <- stats::as.formula(paste("~", rhsX))
  if (is.null(fy_formula)) fy_formula <- stats::as.formula(paste("~ Z +", rhsX))
  if (is.null(support)) support <- sort(unique(Y[A == 0]))
  family <- if (outcome_type == "binary") "binomial" else "gaussian"

  ## --- design matrices ------------------------------------------------------
  sp <- siv_split_ZC(pi_formula, df)
  Dpi <- sp$mm                                  # includes intercept
  Dz  <- siv_model_matrix(fz_formula, df)       # instrument model
  Dy  <- siv_model_matrix(fy_formula, df)       # outcome model
  bs  <- siv_basis_fit(basis, Y[A == 0], X)
  Bobs <- siv_basis_eval(bs, Y, X)              # n x q ; only used where A == 0
  q <- ncol(Bobs)

  ## h1 / h2 columns of the baseline design (eqs. 4.3 and 4.4)
  H1 <- Dpi[, sp$idx_h1, drop = FALSE]
  H2 <- Dpi[, sp$idx_h2, drop = FALSE]
  ## E(h1 | C; alpha) and E(l | C; alpha) with binary Z
  df0 <- df; df0$Z <- 0
  df1 <- df; df1$Z <- 1
  H1_z0 <- siv_model_matrix(pi_formula, df0)[, sp$idx_h1, drop = FALSE]
  H1_z1 <- siv_model_matrix(pi_formula, df1)[, sp$idx_h1, drop = FALSE]
  lv <- as.numeric(l_fun(Z, X)); l0 <- as.numeric(l_fun(rep(0, n), X)); l1 <- as.numeric(l_fun(rep(1, n), X))
  wv <- as.numeric(w_fun(Z, X)); w0 <- as.numeric(w_fun(rep(0, n), X)); w1 <- as.numeric(w_fun(rep(1, n), X))

  ## --- parameter block layout ----------------------------------------------
  n_al <- ncol(Dz)
  n_pi <- ncol(Dpi)                              # intercept + h1 + h2
  n_et <- siv_ydens_npar(family, Dy, support)
  n_m2 <- ncol(H2)                               # empirical means of h2
  blocks <- list(alpha = seq_len(n_al))
  k <- n_al
  blocks$mh2  <- k + seq_len(n_m2); k <- k + n_m2
  if (estimator != "ipw") { blocks$eta <- k + seq_len(n_et); k <- k + n_et }
  blocks$pi   <- k + seq_len(n_pi); k <- k + n_pi
  blocks$psi  <- k + seq_len(q);    k <- k + q
  blocks$theta <- k + 1L; k <- k + 1L
  blocks$mu1   <- k + 1L; k <- k + 1L
  npar <- k

  pnames <- c(paste0("alpha.", colnames(Dz)),
              paste0("mh2.", colnames(H2)),
              if (estimator != "ipw") paste0("eta", seq_len(n_et)) else NULL,
              paste0("pi.", colnames(Dpi)),
              paste0("psi.", colnames(Bobs)),
              "theta", "mu1")

  ## --- pieces of the estimating function ------------------------------------
  ## Extended propensity score pi(Y0, Z, C; theta):  only evaluated where A = 0
  ## (there Y = Y0) or through its value at the observed Y for the DR term.
  pi_of <- function(par) {
    lin <- as.vector(Dpi %*% par[blocks$pi]) + as.vector(Bobs %*% par[blocks$psi])
    expit(lin)
  }
  ## E[exp(psi) v(Y) | A = 0, Z, C; eta] / E[exp(psi) | A = 0, Z, C; eta]
  ## which is eq. (4.6) with g = v.  `Y` and the basis are integrated in a
  ## single quadrature pass, since both are needed by the OR and DR estimators.
  m_both <- function(par) {
    tm <- siv_tilt_moments(family, par[blocks$eta], Dy, support, bs, X,
                           par[blocks$psi],
                           function(yv) cbind(yv, siv_basis_eval(bs, yv, X)),
                           sgn = 1, n_nodes = n_nodes)
    m <- tm$num / tm$den
    list(Y = as.vector(m[, 1L]), basis = m[, -1L, drop = FALSE])
  }
  m_of <- function(par, which = c("Y", "basis")) {
    which <- match.arg(which)
    m_both(par)[[which]]
  }

  estfun <- function(par) {
    pz <- expit(as.vector(Dz %*% par[blocks$alpha]))
    sc_alpha <- Dz * (Z - pz)                                 # MLE score for f_{Z|C}
    Eh1 <- H1_z1 * pz + H1_z0 * (1 - pz)
    El  <- l1 * pz + l0 * (1 - pz)
    Ew  <- w1 * pz + w0 * (1 - pz)
    mh2 <- par[blocks$mh2]
    sc_mh2 <- sweep(H2, 2L, mh2, "-")                         # defines \hat E(h2)

    pihat <- pi_of(par)
    ipw <- (1 - A) / (1 - pihat)                              # (1-A)/{1-pi}

    ## eq. (4.2)
    e42 <- ipw - 1
    ## eq. (4.3)
    e43 <- ipw * (H1 - Eh1)
    ## eq. (4.4)
    e44 <- (1 - ipw) * sweep(H2, 2L, mh2, "-")
    cols <- list(sc_alpha, sc_mh2)

    if (estimator == "ipw") {
      ## eq. (4.5)
      e45 <- (ipw * (lv - El)) * Bobs
      ## Proposition 2: theta_ipw = P_n[ pi Y (1-A) / {Pr(A=1) (1-pi)} ]
      e_th <- pihat * Y * ipw - par[blocks$theta] * A
      cols <- c(cols, list(cbind(e42), e43, e44, e45, cbind(e_th)))
    } else {
      sc_eta <- siv_ydens_score(family, par[blocks$eta], Dy, Y, support, 1 - A)
      mm <- m_both(par)
      mY <- mm$Y                                               # E(Y0 | A=1,Z,C)
      if (estimator == "or") {
        ## eq. (4.7) with g = b (the selection-bias basis) and w scalar
        mB <- mm$basis                                         # n x q
        e47 <- ((wv - Ew)) * (A * mB + (1 - A) * Bobs)
        ## Proposition 3
        e_th <- A * mY - par[blocks$theta] * A
        ## The OR estimator does not use eqs. (4.2)-(4.4); the baseline part of
        ## the extended propensity score is not identified by it, so it is held
        ## at zero and excluded from inference (see `keep` below).
        e_pi <- Dpi * 0
        cols <- c(cols, list(sc_eta, e_pi, e47, cbind(e_th)))
      } else {
        ## DR: eq. (4.11)-(4.12).  The displayed expansion of \tilde Q_g in the
        ## paper,
        ##   (1-A) g /(1-pi) + {A - pi}/(1-pi) * m_g ,
        ## is rewritten in the algebraically identical but numerically stable
        ## form
        ##   (1-A) {g - pi * m_g}/(1-pi) + A * m_g ,
        ## in which pi is only ever evaluated where A = 0 (so that Y = Y0) and
        ## no 0/0 arises as pi -> 1.
        mB <- mm$basis
        Qb <- (1 - A) * (Bobs - pihat * mB) / (1 - pihat) + A * mB
        e412 <- (wv - Ew) * Qb
        ## Proposition 4 uses the *un-tilded* Q_g of eq. (4.11) with
        ## \tilde g(Y,C) = Y:
        ##   Q_g = {(1-A) pi /(1-pi)} {g - m_g} + A m_g ,
        ## whose mean is E(A Y0) = theta * Pr(A = 1).  (The tilded
        ## \tilde Q_g = Q_g + (1-A) g used in eq. (4.12) has mean E(Y0) instead
        ## and is the right object for the covariance restriction, not for the
        ## target parameter.)
        QY <- ((1 - A) * pihat / (1 - pihat)) * (Y - mY) + A * mY
        e_th <- QY - par[blocks$theta] * A
        cols <- c(cols, list(sc_eta, cbind(e42), e43, e44, e412, cbind(e_th)))
      }
    }
    ## mu1 = E(Y | A = 1)
    cols <- c(cols, list(cbind(A * (Y - par[blocks$mu1]))))
    out <- do.call(cbind, cols)
    out
  }

  ## For the OR estimator the baseline propensity block is not identified; hold
  ## it fixed at zero and drop it from the parameter vector.
  free <- rep(TRUE, npar)
  if (estimator == "or") free[blocks$pi] <- FALSE

  ## --- starting values ------------------------------------------------------
  par0 <- numeric(npar)
  mz <- suppressWarnings(stats::glm.fit(Dz, Z, family = stats::binomial()))
  par0[blocks$alpha] <- ifelse(is.finite(mz$coefficients), mz$coefficients, 0)
  par0[blocks$mh2] <- colMeans(H2)
  if (estimator != "ipw") {
    par0[blocks$eta] <- siv_ydens_start(family, Dy, Y, which(A == 0), support)
  }
  ma <- suppressWarnings(stats::glm.fit(Dpi, A, family = stats::binomial()))
  par0[blocks$pi] <- ifelse(is.finite(ma$coefficients), ma$coefficients, 0)
  par0[blocks$psi] <- 0
  par0[blocks$theta] <- mean(Y[A == 0])
  par0[blocks$mu1] <- mean(Y[A == 1])
  if (estimator == "or") par0[blocks$pi] <- 0
  ## user- or warm-start-supplied starting values, given as a named list of
  ## parameter blocks
  if (!is.null(start)) {
    for (nmb in intersect(names(start), names(blocks))) {
      v <- start[[nmb]]
      if (!is.null(v) && length(v) == length(blocks[[nmb]])) par0[blocks[[nmb]]] <- v
    }
  }
  ## the baseline propensity block is not identified by the OR estimating
  ## equations and is held at zero (and excluded from `free`)
  if (estimator == "or") par0[blocks$pi] <- 0

  list(estfun = estfun, par0 = par0, blocks = blocks, names = pnames,
       free = free, n = n, estimator = estimator, outcome_type = outcome_type,
       basis = bs, Dz = Dz, Dpi = Dpi, Dy = Dy, support = support,
       pi_of = pi_of, m_of = m_of)
}


#' Liu et al. (2021) separable-IV estimators of the ATT
#'
#' Implements the inverse-probability-weighted (IPW), outcome-regression (OR)
#' and doubly robust (DR) estimators of \eqn{\theta = E(Y^{(0)} \mid A = 1)}
#' and of the effect of treatment on the treated
#' \eqn{\mathrm{ETT} = E(Y \mid A = 1) - \theta} developed in Section 4 of
#' Liu et al. (2021), together with the sandwich covariance matrix obtained by
#' stacking all estimating equations (as in the data analysis of that paper,
#' Section 6: "se was evaluated using the sandwich estimator accounting for all
#' sources of variabilities").
#'
#' @section Estimating equations:
#' The extended propensity score is modelled as
#' \deqn{\mathrm{logit}\,\pi(Y^{(0)}, Z, C; \theta) =
#'       \theta_0 + h_1(Z,C)^\top\theta_1 + h_2(C)^\top\theta_2
#'       + b(Y^{(0)}, C)^\top\psi ,}
#' where \eqn{h_1} collects the baseline terms that involve the instrument,
#' \eqn{h_2} those that do not, and \eqn{b} is the normalised selection-bias
#' basis (see [siv_basis_spec()]).  Following the recipe below eq. (4.5) of the
#' paper, \eqn{(h_1^\top, h_2^\top, t^\top)^\top =
#' \partial\,\mathrm{logit}\,\pi/\partial\tilde\theta} and \eqn{t = b}, with
#' \eqn{l(Z,C) = Z} a scalar function of \eqn{(Z,C)}, giving exactly
#' \eqn{\dim(\theta)} equations.
#'
#' `estimator = "ipw"` solves eqs. (4.2)-(4.5)
#' \deqn{P_n[(1-A)/\{1-\pi\} - 1] = 0,}
#' \deqn{P_n[\{(1-A)/(1-\pi)\}\{h_1 - E(h_1\mid C;\hat\alpha)\}] = 0,}
#' \deqn{P_n[\{1 - (1-A)/(1-\pi)\}\{h_2 - \hat E(h_2)\}] = 0,}
#' \deqn{P_n[\{(1-A)/(1-\pi)\}\,t(Y,C)\{l(Z,C) - E(l\mid C;\hat\alpha)\}] = 0,}
#' and reports \eqn{\hat\theta_{ipw}} of Proposition 2.
#'
#' `estimator = "or"` estimates \eqn{\eta} by restricted maximum likelihood
#' among the untreated, solves eq. (4.7)
#' \deqn{P_n[\{w(Z,C) - E(w\mid C;\hat\alpha)\}
#'   \{A\,E[g(Y^{(0)},C)\mid A=1,Z,C;\psi,\hat\eta] + (1-A)g(Y,C)\}] = 0}
#' with \eqn{E[g \mid A = 1, Z, C]} given by eq. (4.6), and reports
#' \eqn{\hat\theta_{reg}} of Proposition 3.
#'
#' `estimator = "dr"` solves eq. (4.12)
#' \deqn{P_n[\{w(Z,C) - E(w\mid C;\hat\alpha)\}\,
#'        \tilde Q_g(Y,A,Z,C;\tilde\theta,\hat\eta)] = 0}
#' jointly with eqs. (4.2)-(4.4), where
#' \deqn{\tilde Q_g = \frac{1-A}{1-\pi}g(Y,C)
#'   + \frac{A - \pi}{1-\pi}\,
#'     \frac{E[e^{\psi}g \mid A=0,Z,C;\eta]}{E[e^{\psi}\mid A=0,Z,C;\eta]},}
#' and reports \eqn{\hat\theta_{DR}} of Proposition 4 with
#' \eqn{\tilde g(Y,C) = Y}.
#'
#' @param Y Numeric outcome vector (observed for all units; only its values on
#'   `A == 0` are treated as the treatment-free potential outcome).
#' @param A Binary treatment indicator.
#' @param Z Binary instrument.
#' @param X Matrix or data frame of covariates \eqn{C}.
#' @param estimator One of `"dr"` (default), `"ipw"`, `"or"`.
#' @param outcome_type `"binary"` or `"continuous"`; determines the working
#'   model for \eqn{f(Y \mid A = 0, Z, C;\eta)} (logistic vs. homoscedastic
#'   normal) used by the OR and DR estimators.
#' @param basis A [siv_basis_spec()] describing the selection-bias function
#'   \eqn{\psi}.  The default `siv_basis_spec("linear")` reproduces eq. (4.8)
#'   of the paper, \eqn{\psi(Y^{(0)}, Z, C;\psi) = \psi Y^{(0)}}.
#' @param pi_formula One-sided formula for the *baseline* part of the extended
#'   propensity score (default `~ Z + <all covariates>`).
#' @param fz_formula One-sided formula for the instrument model
#'   \eqn{f_{Z\mid C}} (default `~ <all covariates>`).
#' @param fy_formula One-sided formula for the outcome model
#'   \eqn{f(Y \mid A = 0, Z, C;\eta)} (default `~ Z + <all covariates>`).
#' @param l_fun,w_fun Functions of `(Z, X)` returning the scalar instruments
#'   \eqn{l(Z,C)} of eq. (4.5) and \eqn{w(Z,C)} of eqs. (4.7)/(4.12).  Both
#'   default to \eqn{Z}, the choice used in the paper's simulations.
#' @param support Outcome support used by the binomial working model; defaults
#'   to `sort(unique(Y[A == 0]))`.
#' @param n_nodes Number of Gauss-Hermite nodes for the tilting integrals when
#'   `outcome_type = "continuous"`.
#' @param start Optional named list of starting values, with elements named
#'   after the parameter blocks (`alpha`, `mh2`, `eta`, `pi`, `psi`, `theta`,
#'   `mu1`).  Supplying it disables the automatic warm start.
#' @param warm_start Logical.  When the selection-bias basis has more than one
#'   column, first fit the scalar model of eq. (4.8) and project its solution
#'   onto the richer basis to initialise the full system.  This changes only the
#'   starting value, not the estimating equations.
#' @param level Confidence level.
#' @param tol,maxit Control parameters for the estimating-equation solver.
#'
#' @return An object of class `"SIV_Liu"` (see [summary.siv_fit()]) containing
#'   `estimate` (a named vector with `theta` = \eqn{E(Y^{(0)}\mid A=1)},
#'   `ETT`, the selection-bias parameters and the remaining model parameters),
#'   `vcov` (the sandwich covariance of all of them), `se`, `converged` and the
#'   fitted nuisance quantities.
#'
#' @references
#' Liu, L., Miao, W., Sun, B., Robins, J. and Tchetgen Tchetgen, E. J.
#' Identification and inference for marginal average treatment effect on the
#' treated with an instrumental variable.  Sections 4.1-4.3.
#'
#' @examples
#' set.seed(3)
#' n <- 600; X <- cbind(X1 = rnorm(n)); Z <- rbinom(n, 1, plogis(0.2 * X[, 1]))
#' Y0 <- rbinom(n, 1, 0.5); A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0))
#' Y <- ifelse(A, 1, Y0)
#' fit <- SIV_Liu(Y, A, Z, X, estimator = "dr", outcome_type = "binary")
#' coef(fit)[c("theta", "ETT")]
#' sqrt(diag(vcov(fit)))[c("theta", "ETT")]
#' @export
SIV_Liu <- function(Y, A, Z, X = NULL,
                    estimator = c("dr", "ipw", "or"),
                    outcome_type = c("binary", "continuous"),
                    basis = siv_basis_spec("linear"),
                    pi_formula = NULL, fz_formula = NULL, fy_formula = NULL,
                    l_fun = function(Z, X) Z,
                    w_fun = function(Z, X) Z,
                    support = NULL, n_nodes = 20L,
                    start = NULL, level = 0.95, warm_start = TRUE,
                    tol = 1e-7, maxit = 200L) {
  cl <- match.call()
  ## Two-stage warm start.  When the selection-bias model has more than one
  ## parameter, its estimating equations are poorly conditioned from a zero
  ## start, because a binary instrument supplies only one direction of exogenous
  ## variation.  We therefore first fit the scalar model psi(Y0) = psi * Y0 --
  ## exactly the specification of Liu et al. (2021), eq. (4.8), which is well
  ## conditioned -- and project its fit onto the richer basis to initialise the
  ## full system.  This changes only the starting value, never the estimating
  ## equations, so the estimator and its sandwich are unaffected.
  if (isTRUE(warm_start) && is.null(start)) {
    A0 <- as.numeric(A)
    bs_full <- siv_basis_fit(basis, Y[A0 == 0], if (is.null(X)) NULL else as.matrix(X))
    if (siv_basis_dim(bs_full) > 1L) {
      f0 <- try(SIV_Liu(Y, A, Z, X, estimator = estimator,
                        outcome_type = outcome_type,
                        basis = siv_basis_scalar(bs_full),
                        pi_formula = pi_formula, fz_formula = fz_formula,
                        fy_formula = fy_formula, l_fun = l_fun, w_fun = w_fun,
                        support = support, n_nodes = n_nodes, start = NULL,
                        level = level, warm_start = FALSE,
                        tol = tol, maxit = maxit), silent = TRUE)
      if (!inherits(f0, "try-error") && isTRUE(f0$converged)) {
        b0 <- f0$blocks; p0 <- f0$par_full
        Xm <- if (is.null(X)) matrix(numeric(0), length(Y), 0) else as.matrix(X)
        ii <- which(A0 == 0)         # Y equals Y0 only among the untreated
        start <- list(
          alpha = p0[b0$alpha], mh2 = p0[b0$mh2],
          eta = if (!is.null(b0$eta)) p0[b0$eta] else NULL,
          pi = p0[b0$pi],
          psi = siv_project_psi(bs_full, Y[ii], Xm[ii, , drop = FALSE],
                                p0[b0$psi], bs_full$y_ref),
          theta = p0[b0$theta], mu1 = p0[b0$mu1])
      }
    }
  }
  su <- siv_liu_setup(Y, A, Z, X, estimator = estimator,
                      outcome_type = outcome_type, basis = basis,
                      pi_formula = pi_formula, fz_formula = fz_formula,
                      fy_formula = fy_formula, l_fun = l_fun, w_fun = w_fun,
                      support = support, n_nodes = n_nodes, start = start)
  estfun <- su$estfun; par0 <- su$par0; free <- su$free
  blocks <- su$blocks; pnames <- su$names

  ef_free <- function(pf) {
    p <- par0; p[free] <- pf
    colMeans(estfun(p))[free]
  }
  sol <- siv_solve_ee(ef_free, par0[free], tol = tol, maxit = maxit)
  par <- par0; par[free] <- sol$par

  ## --- sandwich covariance (Sun et al. 2021, Supplemental Appendix eq. S2) ---
  ef_sw <- function(pf) {
    p <- par; p[free] <- pf
    estfun(p)[, free, drop = FALSE]
  }
  sw <- siv_sandwich(ef_sw, par[free], names = pnames[free])

  ## ETT = mu1 - theta : exact linear contrast of the stacked parameters
  ii <- match(c("theta", "mu1"), pnames[free])
  L <- rep(0, sum(free)); L[ii[2L]] <- 1; L[ii[1L]] <- -1
  ett <- sum(L * par[free])
  var_ett <- as.numeric(t(L) %*% sw$vcov %*% L)
  cov_ett <- as.numeric(sw$vcov %*% L)

  V <- rbind(cbind(sw$vcov, cov_ett), c(cov_ett, var_ett))
  est <- c(par[free], ett)
  nm <- c(pnames[free], "ETT")
  names(est) <- nm; dimnames(V) <- list(nm, nm)
  se <- sqrt(pmax(diag(V), 0))
  q95 <- stats::qnorm(1 - (1 - level) / 2)

  out <- list(
    estimate = est, vcov = V, se = se,
    ci = cbind(lower = est - q95 * se, upper = est + q95 * se),
    level = level,
    converged = sol$converged, criterion = sol$criterion,
    checks = c(sw$checks, stats::setNames(siv_check_vcov(V, length(est)),
                                          paste0("final_", names(siv_check_vcov(V, length(est)))))),
    bread = sw$bread, meat = sw$meat, jacobian = sw$jacobian,
    estfun = sw$estfun, estfun_fn = estfun,
    blocks = blocks, free = free, par_full = par,
    nuisance = list(
      pi = su$pi_of(par),
      pZ1 = expit(as.vector(su$Dz %*% par[blocks$alpha])),
      mY = if (estimator != "ipw") as.vector(su$m_of(par, "Y")) else NULL),
    n = su$n, estimator = estimator, outcome_type = outcome_type,
    basis = su$basis,
    se_type = "sandwich (stacked estimating equations)",
    method = sprintf("Liu et al. (2021) separable IV, %s estimator",
                     toupper(estimator)),
    target = c(theta = "E(Y0 | A = 1)", ETT = "E(Y1 - Y0 | A = 1)"),
    call = cl
  )
  class(out) <- c("SIV_Liu", "siv_fit")
  out
}
