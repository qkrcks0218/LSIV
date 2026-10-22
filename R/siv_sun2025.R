## ----------------------------------------------------------------------------
## SIV_Sun(): the estimator of
##   Sun, B., Miao, W. and Wickramarachchi, D. S. (2025).
##   "On doubly robust estimation with nonignorable missing data using
##    instrumental variables."  Statistica Sinica (accepted), SS-2023-0383,
##   together with its Supplementary Material.
##
## The paper is written in the nonignorable-missing-data setting with observed
## data O = (R, RY, X), X = (Z, U), and target mu_0 = E(Y).  The logit-separable
## IV model of Park & Tchetgen Tchetgen is the same model with
##       R = 1 - A ,     Y = Y^{(0)} ,     U = X (the covariates) ,
## i.e. the treatment-free potential outcome plays the role of the outcome and
## is "missing" exactly for the treated.  This function therefore estimates
## mu_0 = E(Y^{(0)}) with the Sun-Miao-Wickramarachchi machinery and then maps
## it to the ATT through the *exact* algebraic identity
##       E(Y^{(0)} | A = 1) = { E(Y^{(0)}) - E[(1-A) Y] } / E(A) ,
##       ATT = E(Y | A = 1) - E(Y^{(0)} | A = 1) ,
## with the extra moments E[(1-A)Y], E(A) and E(Y | A = 1) added to the stacked
## estimating equations so that the sandwich covariance covers the ATT exactly.
##
## Equation map (main paper unless stated otherwise)
##   eq. (2.1)   selection bias function h(y, x)
##   eq. (2.2)   pi(w; alpha, gamma) = expit{ alpha(x) + h(y, x; gamma) }
##   eq. (2.10)  the three moment equations m(O; gamma, tau) for
##               tau = (alpha, zeta, eta):
##                  P_n[{R/pi - 1} d alpha(X;alpha)/d alpha ]        = 0
##                  P_n[ d log p(Z|U; zeta)/d zeta ]                 = 0
##                  P_n[ (R/pi) d log p(Y|U; eta)/d eta ]            = 0
##   Result 1    influence function; q(W; mu, c, d) and g(O; nu, c, d)
##   Example 2   f(W) - f_bar(W) = f_yz(U){Y - p(Y=1|U)}{Z - p(Z=1|U)}
##   Example 3   for continuous Y or Z, use moments v1(Y), v2(Z)
##   eq. (3.11)  p(y | R = 0, x; nu) = {1 - pi(y,x;nu)} p(y|u)
##                                     / [1 - int pi(t,x;nu) p(t|u) dnu(t)]
##   Section 3.4 / Example 4  the computationally simpler estimator, based on
##               g~(O; nu, c, d) = R q(W; mu, c, d) / pi(W; nu), with c(w) = 0,
##               d(w) = yz; it solves (2.10) jointly with
##                  P_n[(R/pi) d_yz(U){Y - p(Y=1|U;eta)}{Z - p(Z=1|U;zeta)}] = 0
##               and takes the Hajek estimator of the outcome mean.
##               This is the estimator used in the paper's own simulation study
##               (Section 4) and application (Section 5), and it is the default
##               here (`augmented = FALSE`).
##   Prop. 1/2 and Supplement S2-S3   sandwich variance
##   Supplement S5.2   the continuous-outcome specification:
##               p(y | u; eta) = (2 pi sigma^2)^{-1/2} exp{-(y - h1(u)'eta)^2
##                               /(2 sigma^2)},  eta = (sigma, eta')' ;
##               convergence criterion "the residual of the square average
##               estimating equation component being less than 1e-7".
## ----------------------------------------------------------------------------

## Internal: E{v(Y) | R = 0, X; nu} under eq. (3.11).
siv_sun2025_m0 <- function(family, eta, Dy, support, X, pi_fun, vfun, n_nodes) {
  n <- nrow(Dy)
  if (family == "binomial") {
    P <- siv_ydens_grid(family, eta, Dy, support)
    den <- numeric(n); numl <- NULL
    for (m in seq_along(support)) {
      yv <- rep(support[m], n)
      w <- (1 - pi_fun(yv)) * P[, m]
      den <- den + w
      vm <- as.matrix(vfun(yv))
      if (is.null(numl)) numl <- matrix(0, n, ncol(vm))
      numl <- numl + vm * w
    }
    return(numl / den)
  }
  k <- ncol(Dy); s <- abs(eta[k + 1L]); mu <- as.vector(Dy %*% eta[seq_len(k)])
  gh <- siv_gh_normal(mu, s, n_nodes)
  den <- numeric(n); numl <- NULL
  for (kk in seq_along(gh$weights)) {
    yv <- gh$nodes[, kk]
    w <- (1 - pi_fun(yv)) * gh$weights[kk]
    den <- den + w
    vk <- as.matrix(vfun(yv))
    if (is.null(numl)) numl <- matrix(0, n, ncol(vk))
    numl <- numl + vk * w
  }
  numl / den
}

#' Estimating-equation setup for the Sun, Miao and Wickramarachchi (2025) estimator
#'
#' Builds the design matrices, parameter blocks, starting values and the
#' per-observation estimating-function closure used by [SIV_Sun()], without
#' solving the equations.  Exported so that eqs. (2.10) and the moment
#' condition of Example 4 can be checked directly.
#'
#' @inheritParams SIV_Sun
#' @return A list with the estimating-function closure (`estfun`), starting
#'   values (`par0`), the parameter block index (`blocks`), the parameter
#'   `names` and the design matrices.
#' @examples
#' set.seed(7)
#' n <- 300; X <- cbind(X1 = rnorm(n)); Z <- rbinom(n, 1, plogis(0.2 * X[, 1]))
#' Y0 <- rbinom(n, 1, 0.5); A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0))
#' Y <- ifelse(A, 1, Y0)
#' su <- siv_sun2025_setup(Y, A, Z, X, outcome_type = "binary")
#' dim(su$estfun(su$par0))
#' @export
siv_sun2025_setup <- function(Y, A, Z, X = NULL,
                              outcome_type = c("binary", "continuous"),
                              basis = siv_basis_spec("linear"),
                              alpha_formula = NULL,
                              fz_formula = NULL,
                              fy_formula = NULL,
                              augmented = FALSE,
                              support = NULL, n_nodes = 20L) {
  outcome_type <- match.arg(outcome_type)
  dat <- siv_check_data(Y, A, Z, X)
  Y <- dat$Y; A <- dat$A; Z <- dat$Z; X <- dat$X; n <- dat$n
  siv_match_type(outcome_type, Y, A)
  R <- 1 - A                                  # R = 1 <=> Y^{(0)} is observed
  ## Y is the treatment-free potential outcome only where A == 0; elsewhere it
  ## is replaced by a placeholder that is always multiplied by R.
  Y0 <- ifelse(A == 0, Y, 0)
  df <- siv_frame(Y0, A, Z, X)
  Xnm <- colnames(X)
  rhsX <- if (length(Xnm)) paste(Xnm, collapse = " + ") else "1"
  ## alpha(x) is a function of x = (z, u), so it may involve Z;
  ## p(z | u) and p(y | u) are functions of u only -- the exclusion restriction
  ## of Assumption 2 is what makes them variation independent (Section 2.3).
  if (is.null(alpha_formula)) alpha_formula <- stats::as.formula(paste("~ Z +", rhsX))
  if (is.null(fz_formula))    fz_formula    <- stats::as.formula(paste("~", rhsX))
  if (is.null(fy_formula))    fy_formula    <- stats::as.formula(paste("~", rhsX))
  if (any(all.vars(fy_formula) == "Z")) {
    stop("`fy_formula` must not involve `Z`: p(y | u; eta) is the full-data ",
         "outcome density given the covariates only (Sun, Miao and ",
         "Wickramarachchi 2025, Section 2.3).", call. = FALSE)
  }
  if (is.null(support)) support <- sort(unique(Y0[R == 1]))
  family <- if (outcome_type == "binary") "binomial" else "gaussian"

  Dal <- siv_model_matrix(alpha_formula, df)
  Dz  <- siv_model_matrix(fz_formula, df)
  Dy  <- siv_model_matrix(fy_formula, df)
  bs  <- siv_basis_fit(basis, Y0[R == 1], X)
  Bobs <- siv_basis_eval(bs, Y0, X)
  q <- ncol(Bobs)

  n_al <- ncol(Dal); n_ze <- ncol(Dz); n_et <- siv_ydens_npar(family, Dy, support)
  blocks <- list(alpha = seq_len(n_al)); k <- n_al
  blocks$zeta  <- k + seq_len(n_ze); k <- k + n_ze
  blocks$eta   <- k + seq_len(n_et); k <- k + n_et
  blocks$gamma <- k + seq_len(q);    k <- k + q
  blocks$mu    <- k + 1L; k <- k + 1L        # mu = E(Y^{(0)})
  blocks$c0    <- k + 1L; k <- k + 1L        # c0 = E[(1-A) Y]
  blocks$pA    <- k + 1L; k <- k + 1L        # pA = Pr(A = 1)
  blocks$mu1   <- k + 1L; k <- k + 1L        # mu1 = E(Y | A = 1)
  npar <- k
  pnames <- c(paste0("alpha.", colnames(Dal)),
              paste0("zeta.", colnames(Dz)),
              paste0("eta", seq_len(n_et)),
              paste0("gamma.", colnames(Bobs)),
              "mu", "c0", "pA", "mu1")

  ## eq. (2.2)
  pi_at <- function(par, yv) {
    expit(as.vector(Dal %*% par[blocks$alpha]) +
            as.vector(siv_basis_eval(bs, yv, X) %*% par[blocks$gamma]))
  }

  estfun <- function(par) {
    pihat <- expit(as.vector(Dal %*% par[blocks$alpha]) +
                     as.vector(Bobs %*% par[blocks$gamma]))
    rw <- R / pihat
    pz <- expit(as.vector(Dz %*% par[blocks$zeta]))

    ## eq. (2.10), first line
    m1 <- (rw - 1) * Dal
    ## eq. (2.10), second line
    m2 <- Dz * (Z - pz)
    ## eq. (2.10), third line
    m3 <- siv_ydens_score(family, par[blocks$eta], Dy, Y0, support, rw)

    ## E{b(Y, X) | U; eta} under the full-data outcome density p(y | u; eta)
    EB <- siv_tilt_moments(family, par[blocks$eta], Dy, support, bs, X,
                           rep(0, q),
                           function(yv) siv_basis_eval(bs, yv, X),
                           sgn = 1, n_nodes = n_nodes)
    EB <- EB$num / EB$den
    ## d(W) - dbar(W) = {b(Y,X) - E(b|U)} {Z - p(Z=1|U)}   (Examples 2-4, S4)
    dW <- (Bobs - EB) * (Z - pz)
    qY <- Y0 - par[blocks$mu]

    if (!augmented) {
      ## Section 3.4: g~(O; nu, c = 0, d) = R q(W; mu, 0, d)/pi
      m4 <- rw * dW
      m5 <- rw * qY
    } else {
      ## Section 3.2 / Result 1:
      ##   g = R q/pi + {1 - R/pi} E(q | R = 0, X; nu),  with eq. (3.11)
      pif <- function(yv) pi_at(par, yv)
      m0B <- siv_sun2025_m0(family, par[blocks$eta], Dy, support, X, pif,
                            function(yv) siv_basis_eval(bs, yv, X), n_nodes)
      m0Y <- siv_sun2025_m0(family, par[blocks$eta], Dy, support, X, pif,
                            function(yv) matrix(yv, ncol = 1L), n_nodes)
      aug <- 1 - rw
      m4 <- rw * dW + aug * ((m0B - EB) * (Z - pz))
      m5 <- rw * qY + aug * (as.vector(m0Y) - par[blocks$mu])
    }
    ## Moments needed to map E(Y^{(0)}) to the ATT
    m6 <- R * Y0 - par[blocks$c0]
    m7 <- A - par[blocks$pA]
    m8 <- A * (Y - par[blocks$mu1])
    cbind(m1, m2, m3, m4, cbind(m5), cbind(m6), cbind(m7), cbind(m8))
  }

  par0 <- numeric(npar)
  mr <- suppressWarnings(stats::glm.fit(Dal, R, family = stats::binomial()))
  par0[blocks$alpha] <- ifelse(is.finite(mr$coefficients), mr$coefficients, 0)
  mz <- suppressWarnings(stats::glm.fit(Dz, Z, family = stats::binomial()))
  par0[blocks$zeta] <- ifelse(is.finite(mz$coefficients), mz$coefficients, 0)
  par0[blocks$eta] <- siv_ydens_start(family, Dy, Y0, which(R == 1), support)
  par0[blocks$gamma] <- 0
  par0[blocks$mu] <- mean(Y0[R == 1])
  par0[blocks$c0] <- mean(R * Y0)
  par0[blocks$pA] <- mean(A)
  par0[blocks$mu1] <- mean(Y[A == 1])

  list(estfun = estfun, par0 = par0, blocks = blocks, names = pnames,
       free = rep(TRUE, npar), n = n, outcome_type = outcome_type,
       basis = bs, support = support, augmented = augmented,
       Dal = Dal, Dz = Dz, Dy = Dy, R = R, Y0 = Y0, A = A, Y = Y)
}

#' Sun, Miao and Wickramarachchi (2025) doubly robust separable-IV estimator
#'
#' Implements the doubly robust estimator of Sun, Miao and Wickramarachchi
#' (2025), adapted to the causal ATT setting through the correspondence
#' \eqn{R = 1 - A}, \eqn{Y = Y^{(0)}}, in which the treatment-free potential
#' outcome is "missing" exactly for the treated.  Inference uses the sandwich
#' covariance of Propositions 1-2 of that paper (Supplement Sections S2-S3).
#'
#' @section Parameterisation:
#' The distinguishing feature of this method, relative to Liu et al. (2021) and
#' Sun et al. (2021), is the *alternative factorisation* of Section 2.3:
#' \deqn{p(r, w;\cdot) = \pi(w;\alpha,\gamma)^r\{1-\pi(w;\alpha,\gamma)\}^{1-r}
#'       p(y \mid u;\eta)\,p(z \mid u;\zeta)\,p(u),}
#' in which the *full-data* outcome density \eqn{p(y \mid u;\eta)} -- not the
#' complete-case density \eqn{p(y \mid R = 1, x)} -- is modelled, so that
#' \eqn{\eta} and \eqn{\zeta} are variation independent under the exclusion
#' restriction.  Consequently `fy_formula` must not involve `Z`.
#'
#' The estimator is doubly robust in the union model of Section 3: it is
#' consistent when \eqn{\alpha(x;\alpha)} is correctly specified and *either*
#' \eqn{p(z\mid u;\zeta)} *or* \eqn{p(y\mid u;\eta)} is correctly specified.
#'
#' @section Estimating equations:
#' `augmented = FALSE` (the default, and the estimator used in the paper's own
#' simulation study and application, Section 3.4 and Example 4) solves
#' eq. (2.10) jointly with
#' \deqn{P_n\left[\frac{R}{\pi(W;\alpha,\gamma)}
#'   \{b(Y,X) - E_\eta(b\mid U)\}\{Z - p(Z=1\mid U;\zeta)\}\right] = 0}
#' and takes the Hajek estimator of \eqn{E(Y^{(0)})} from
#' \eqn{P_n[R(Y-\mu)/\pi] = 0}.  With the default linear basis and no
#' covariate interaction this is exactly Example 4 with \eqn{d(w) = yz} and
#' \eqn{d_{yz}(U) = 1}; with a richer basis the outcome moments \eqn{v_1(Y)} of
#' Example 3 are used, and \eqn{d_{yz}(U)} is taken to be the covariate part of
#' the basis so that the system is exactly identified (a documented
#' implementation choice, see `vignette`-free note in
#' `Method_Documentation/SIV_Sun2025.md`).
#'
#' `augmented = TRUE` adds the non-respondent augmentation of Result 1 /
#' Section 3.2, using the conditional density of eq. (3.11).
#'
#' @param Y Numeric outcome vector.
#' @param A Binary treatment indicator (1 = treated); \eqn{R = 1 - A}.
#' @param Z Binary instrument.
#' @param X Matrix or data frame of covariates \eqn{U}.
#' @param outcome_type `"binary"` or `"continuous"`; selects the working model
#'   \eqn{p(y\mid u;\eta)} (logistic, or the homoscedastic normal of
#'   Supplement S5.2).
#' @param basis A [siv_basis_spec()] for the selection-bias function
#'   \eqn{h(y,x;\gamma)}.  The default reproduces \eqn{h(y,x;\gamma)=\gamma y},
#'   the specification used in Section 4 of the paper.
#' @param alpha_formula One-sided formula for \eqn{\alpha(x;\alpha)}
#'   (default `~ Z + <covariates>`, matching
#'   \eqn{\pi(w) = \mathrm{expit}\{(h_1(u), z)\alpha + \gamma y\}}).
#' @param fz_formula One-sided formula for \eqn{p(z\mid u;\zeta)}.
#' @param fy_formula One-sided formula for \eqn{p(y\mid u;\eta)}; must not
#'   involve `Z`.
#' @param augmented Logical; use the influence-function based estimator of
#'   Section 3.2 (`TRUE`) or the computationally simpler estimator of
#'   Section 3.4 (`FALSE`, the default and the one used in the paper).
#' @param support Outcome support for the binomial working model.
#' @param n_nodes Number of Gauss-Hermite nodes for a continuous outcome.
#' @param start Optional named list of starting values, named after the
#'   parameter blocks (`alpha`, `zeta`, `eta`, `gamma`, `mu`, `c0`, `pA`,
#'   `mu1`).  Supplying it disables the automatic warm start.
#' @param warm_start Logical.  When the tilt basis has more than one column,
#'   first fit the scalar model `h(y, x; gamma) = gamma * y` of Section 4 and
#'   project its solution onto the richer basis to initialise the full system.
#' @param level Confidence level.
#' @param tol,maxit Control parameters of the solver.  The default `tol = 1e-7`
#'   matches the convergence criterion stated in the note to Table S3 of the
#'   supplement.
#'
#' @return An object of class `"SIV_Sun"` with `estimate` containing `mu`
#'   (\eqn{E(Y^{(0)})}), `theta` (\eqn{E(Y^{(0)}\mid A=1)}), `ATT`, the tilting
#'   parameters `gamma.*` and the remaining working-model parameters; the
#'   sandwich `vcov`; `se`; `converged`; and fitted nuisance quantities.
#'
#' @references
#' Sun, B., Miao, W. and Wickramarachchi, D. S.  On doubly robust estimation
#' with nonignorable missing data using instrumental variables.
#' *Statistica Sinica*, SS-2023-0383, and its Supplementary Material.
#'
#' @examples
#' set.seed(8)
#' n <- 600; X <- cbind(X1 = rnorm(n)); Z <- rbinom(n, 1, plogis(0.2 * X[, 1]))
#' Y0 <- rbinom(n, 1, 0.5); A <- rbinom(n, 1, plogis(-2 + 4 * Z + 0.3 * Y0))
#' Y <- ifelse(A, 1, Y0)
#' fit <- SIV_Sun(Y, A, Z, X, outcome_type = "binary")
#' coef(fit)[c("mu", "theta", "ATT")]
#' @export
SIV_Sun <- function(Y, A, Z, X = NULL,
                    outcome_type = c("binary", "continuous"),
                    basis = siv_basis_spec("linear"),
                    alpha_formula = NULL, fz_formula = NULL, fy_formula = NULL,
                    augmented = FALSE,
                    support = NULL, n_nodes = 20L, start = NULL,
                    warm_start = TRUE,
                    level = 0.95, tol = 1e-7, maxit = 200L) {
  cl <- match.call()
  ## Two-stage warm start; see the note in `SIV_Liu()`.  The scalar model
  ## h(y, x; gamma) = gamma * y is exactly the specification of Section 4 of the
  ## paper and is well conditioned; its fit initialises the richer model.
  if (isTRUE(warm_start) && is.null(start)) {
    A0 <- as.numeric(A)
    Xm <- if (is.null(X)) matrix(numeric(0), length(Y), 0) else as.matrix(X)
    bs_full <- siv_basis_fit(basis, Y[A0 == 0], Xm)
    if (siv_basis_dim(bs_full) > 1L) {
      f0 <- try(SIV_Sun(Y, A, Z, X, outcome_type = outcome_type,
                        basis = siv_basis_scalar(bs_full),
                        alpha_formula = alpha_formula, fz_formula = fz_formula,
                        fy_formula = fy_formula, augmented = augmented,
                        support = support, n_nodes = n_nodes, start = NULL,
                        warm_start = FALSE, level = level, tol = tol,
                        maxit = maxit), silent = TRUE)
      if (!inherits(f0, "try-error") && isTRUE(f0$converged)) {
        b0 <- f0$blocks; p0 <- f0$par_full
        start <- list(
          alpha = p0[b0$alpha], zeta = p0[b0$zeta], eta = p0[b0$eta],
          gamma = siv_project_psi(bs_full, Y[A0 == 0],
                                  Xm[A0 == 0, , drop = FALSE],
                                  p0[b0$gamma], bs_full$y_ref),
          mu = p0[b0$mu], c0 = p0[b0$c0], pA = p0[b0$pA], mu1 = p0[b0$mu1])
      }
    }
  }
  su <- siv_sun2025_setup(Y, A, Z, X, outcome_type = outcome_type,
                          basis = basis, alpha_formula = alpha_formula,
                          fz_formula = fz_formula, fy_formula = fy_formula,
                          augmented = augmented, support = support,
                          n_nodes = n_nodes)
  estfun <- su$estfun; par0 <- su$par0; pnames <- su$names; bl <- su$blocks
  if (!is.null(start)) {
    for (nmb in intersect(names(start), names(bl))) {
      v <- start[[nmb]]
      if (!is.null(v) && length(v) == length(bl[[nmb]])) par0[bl[[nmb]]] <- v
    }
  }
  sol <- siv_solve_ee(function(p) colMeans(estfun(p)), par0, tol = tol,
                      maxit = maxit)
  par <- sol$par
  sw <- siv_sandwich(estfun, par, names = pnames)

  ## Delta method for  theta = (mu - c0)/pA  and  ATT = mu1 - theta.
  mu <- par[bl$mu]; c0 <- par[bl$c0]; pA <- par[bl$pA]; mu1 <- par[bl$mu1]
  theta <- (mu - c0) / pA
  att <- mu1 - theta
  p <- length(par)
  G <- matrix(0, 2L, p)
  G[1L, bl$mu] <-  1 / pA
  G[1L, bl$c0] <- -1 / pA
  G[1L, bl$pA] <- -(mu - c0) / pA^2
  G[2L, ] <- -G[1L, ]
  G[2L, bl$mu1] <- G[2L, bl$mu1] + 1
  Vd <- G %*% sw$vcov %*% t(G)
  Cross <- sw$vcov %*% t(G)

  V <- rbind(cbind(sw$vcov, Cross), cbind(t(Cross), Vd))
  V <- (V + t(V)) / 2
  est <- c(par, theta = theta, ATT = att)
  nm <- c(pnames, "theta", "ATT")
  names(est) <- nm; dimnames(V) <- list(nm, nm)
  se <- sqrt(pmax(diag(V), 0))
  q95 <- stats::qnorm(1 - (1 - level) / 2)

  out <- list(
    estimate = est, vcov = V, se = se,
    ci = cbind(lower = est - q95 * se, upper = est + q95 * se),
    level = level, converged = sol$converged, criterion = sol$criterion,
    checks = c(sw$checks, stats::setNames(siv_check_vcov(V, length(est)),
                                          paste0("final_", names(siv_check_vcov(V, length(est)))))),
    bread = sw$bread, meat = sw$meat, jacobian = sw$jacobian,
    estfun = sw$estfun, estfun_fn = estfun, blocks = bl, par_full = par,
    nuisance = list(
      pi = expit(as.vector(su$Dal %*% par[bl$alpha]) +
                   as.vector(siv_basis_eval(su$basis, su$Y0, as.matrix(
                     if (is.null(X)) matrix(numeric(0), su$n, 0) else X)) %*%
                     par[bl$gamma])),
      pZ1 = expit(as.vector(su$Dz %*% par[bl$zeta]))),
    n = su$n, outcome_type = su$outcome_type, basis = su$basis,
    augmented = augmented,
    estimator = if (augmented) "IF-based (Section 3.2)" else "simplified DR (Section 3.4)",
    se_type = "sandwich (stacked estimating equations; Propositions 1-2 / Supplement S2-S3)",
    method = "Sun, Miao & Wickramarachchi (2025) doubly robust separable IV",
    target = c(mu = "E(Y0)", theta = "E(Y0 | A = 1)", ATT = "E(Y1 - Y0 | A = 1)"),
    call = cl
  )
  class(out) <- c("SIV_Sun", "siv_fit")
  out
}
