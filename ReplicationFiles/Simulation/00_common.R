## Shared, non-method-specific utilities for the version-2 simulations.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x
}

if (!requireNamespace("SIV", quietly = TRUE)) {
  stop(
    "Install SIV and all required dependencies first:\n  Rscript Code/Install.R",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(SIV))

siv_env_int <- function(name, default, lower = 1L) {
  ans <- suppressWarnings(as.integer(Sys.getenv(name) %||% as.character(default)))
  if (length(ans) != 1L || is.na(ans) || ans < lower) {
    stop(sprintf("%s must be an integer >= %d.", name, lower), call. = FALSE)
  }
  ans
}

## Comma-separated list of positive integers, e.g. SIM_N="500,1000,2000".
siv_env_int_vector <- function(name, default, lower = 1L) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) value <- paste(default, collapse = ",")
  ans <- suppressWarnings(as.integer(trimws(
    strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(ans) || anyNA(ans) || any(ans < lower)) {
    stop(sprintf("%s must be a comma-separated list of integers >= %d.",
                 name, lower), call. = FALSE)
  }
  ans
}

siv_env_types <- function(name = "SIM_TYPE", default = "both") {
  value <- tolower(Sys.getenv(name) %||% default)
  if (identical(value, "both")) return(c("continuous", "binary"))
  value <- strsplit(value, ",", fixed = TRUE)[[1L]]
  value <- trimws(value)
  if (!length(value) || !all(value %in% c("continuous", "binary"))) {
    stop(sprintf("%s must be continuous, binary, or both.", name), call. = FALSE)
  }
  unique(value)
}

## ---------------------------------------------------------------------------
## Simulation DGPs. These live here, not in SIV_Package: they are properties
## of the simulation designs rather than of any estimator implementation.
## ---------------------------------------------------------------------------

## EDIT THE DGP HERE. All observed-data generation, potential-treatment draws,
## oracle nuisances, and Monte Carlo truths below use this specification.
siv_dgp_forall <- function(type = c("continuous", "binary"), p = 3L,
                           violation = c("siv", "mild", "severe")) {
  type <- match.arg(type)
  violation <- match.arg(violation)
  p <- as.integer(p)
  if (length(p) != 1L || is.na(p) || p < 1L) {
    stop("p must be a positive integer.", call. = FALSE)
  }
  
  draw_x <- function(n) {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- sprintf("X%02d", seq_len(p))
    X
  }
  x_index <- function(X) rowSums(as.matrix(X))
  p_z1 <- function(W) plogis(0.1 * W)
  
  
  log_gamma <- function(z, W){ 
    -2.75 + 2.75 * z + 0.1 * W
  } 
  
  alpha_shape <- function(y) {
    if (type == "continuous") {
      c0 <- 1
      c1 <- 0.5
      c1*((1-1/(y+1+c0)^2)*as.numeric(y< -2-c0) +
            ((y+1+c0)^2-1)*as.numeric(-2-c0<=y & y< -1-c0) +
            (-1)* as.numeric(-1-c0<=y & y< -1) +
            ((y+1)^2-1)*as.numeric(-1<=y & y<=0) +
            (1-1/(y+1)^2)*as.numeric(0<y))
    } else {
      y
    }
  }
  
  log_delta <- function(y0, W, X = NULL) {
    if (type == "continuous") {
      (1+2.5*tanh(2.25*W))*alpha_shape(y0)
    } else {
      -(1+3*tanh(2*X[,1])) * alpha_shape(y0)
    }
  }
  
  if (type == "continuous") {
    mean_y0 <- function(W) {
      +0.3-0.2 * W
    } 

    draw_y0 <- function(W, X = NULL) mean_y0(W) + rnorm(length(W))
    draw_y1 <- function(W, X = NULL){ 
      1 + 1.0 * W + 0.25 * rnorm(length(W))
    } 
    density_y0 <- function(grid, W, X = NULL) {
      vapply(grid, function(y) dnorm(y, mean_y0(W), 1),
             numeric(length(W)))
    }
    
    
  } else { ## Binary
    prob_y0 <- function(W, X = NULL){ 
      plogis(-W)
    } 
    prob_y1 <- function(W, X = NULL){ 
      plogis(2+W)
    } 
    
    draw_y0 <- function(W, X = NULL) {
      rbinom(length(W), 1, prob_y0(W, X))
    }
    draw_y1 <- function(W, X = NULL) {
      rbinom(length(W), 1, prob_y1(W, X))
    }
    
    density_y0 <- function(grid, W, X = NULL) {
      py1 <- prob_y0(W, X)
      vapply(grid, function(y) if (y == 1) py1 else 1 - py1,
             numeric(length(W)))
    }
    
  }
  
  
  
  ## -------------------------------------------------------------------------
  ## Departure from separability, used by Falsification/01_falsification_test.R.
  
  
  violation_parameters <- switch(
    violation,
    siv    = list(kappa = 0),
    mild   = list(kappa = 3),
    severe = list(kappa = 4.5)
  )
  
  b_shape <- function(y0) {
    4 * y0 * (1 - y0) * exp(-(y0 - 0.5)^2)
  }
  
  log_alpha <- function(y0, z, W, X = NULL) {
    log_delta(y0, W, X) -
      (2 * z - 1) * violation_parameters$kappa * b_shape(y0)
  }
  alpha <- function(y0, z, W, X = NULL) exp(log_alpha(y0, z, W, X))
  
  list(
    type = type, p = p, draw_x = draw_x, x_index = x_index,
    p_z1 = p_z1, draw_y0 = draw_y0, draw_y1 = draw_y1,
    density_y0 = density_y0, alpha_shape = alpha_shape,
    log_delta = log_delta, log_gamma = log_gamma,
    violation = violation, violation_parameters = violation_parameters,
    baseline_alpha_shape = alpha_shape, b_shape = b_shape,
    log_alpha = log_alpha, alpha = alpha
  )
}

## Compatibility helper used by a few diagnostics.
siv_alpha_shape <- function(y, type = c("continuous", "binary")) {
  type <- match.arg(type)
  siv_dgp_forall(type)$alpha_shape(y)
}

## Generic sampler. Do not repeat DGP formulas here.
siv_dgp <- function(n, type = c("continuous", "binary"), p = 3L,
                    violation = c("siv", "mild", "severe")) {
  type <- match.arg(type)
  violation <- match.arg(violation)
  spec <- siv_dgp_forall(type, p, violation)
  X <- spec$draw_x(n)
  W <- spec$x_index(X)
  pZ1 <- spec$p_z1(W)
  Z <- rbinom(n, 1, pZ1)
  Y0 <- spec$draw_y0(W, X)
  Y1 <- spec$draw_y1(W, X)
  ## With `violation = "siv"` (the default) log_alpha == log_delta for both z,
  ## so this is exactly the separable design used by the estimation scripts.
  logalpha_z0 <- spec$log_alpha(Y0, 0, W, X)
  logalpha_z1 <- spec$log_alpha(Y0, 1, W, X)
  logdelta <- spec$log_delta(Y0, W, X)
  loggamma_z0 <- spec$log_gamma(0, W)
  loggamma_z1 <- spec$log_gamma(1, W)
  pA0 <- plogis(logalpha_z0 + loggamma_z0)
  pA1 <- plogis(logalpha_z1 + loggamma_z1)
  treatment_rank <- runif(n)
  A0 <- as.integer(treatment_rank <= pA0)
  A1 <- as.integer(treatment_rank <= pA1)
  A <- ifelse(Z == 1, A1, A0)
  pA <- ifelse(Z == 1, pA1, pA0)
  loggamma <- ifelse(Z == 1, loggamma_z1, loggamma_z0)
  Y <- (1 - A) * Y0 + A * Y1
  list(Y = Y, A = A, Z = Z, X = X, Y0 = Y0, Y1 = Y1, W = W,
       A0 = A0, A1 = A1, pA = pA, pA0 = pA0, pA1 = pA1, pZ1 = pZ1,
       logdelta = logdelta, loggamma = loggamma, type = type,
       violation = violation,
       alpha = exp(ifelse(Z == 1, logalpha_z1, logalpha_z0)),
       alpha_z0 = exp(logalpha_z0), alpha_z1 = exp(logalpha_z1),
       logalpha = ifelse(Z == 1, logalpha_z1, logalpha_z0),
       logalpha_z0 = logalpha_z0, logalpha_z1 = logalpha_z1)
}

## Generic oracle calculation from the same central DGP specification.
siv_oracle_nuisance <- function(X, grid, weights,
                                type = c("continuous", "binary"),
                                violation = c("siv", "mild", "severe")) {
  type <- match.arg(type)
  violation <- match.arg(violation)
  X <- as.matrix(X)
  n <- nrow(X); M <- length(grid)
  spec <- siv_dgp_forall(type, ncol(X), violation)
  W <- spec$x_index(X)
  pZ1 <- spec$p_z1(W)
  fy <- matrix(spec$density_y0(grid, W, X), nrow = n, ncol = M)
  Wm <- matrix(weights, n, M, byrow = TRUE)
  logalpha_z0 <- matrix(vapply(
    grid, function(y) spec$log_alpha(rep(y, n), 0, W, X), numeric(n)),
    nrow = n, ncol = M)
  logalpha_z1 <- matrix(vapply(
    grid, function(y) spec$log_alpha(rep(y, n), 1, W, X), numeric(n)),
    nrow = n, ncol = M)
  ## The separable part only; equal to alpha for both z when violation = "siv".
  delta <- exp(matrix(vapply(
    grid, function(y) spec$log_delta(rep(y, n), W, X), numeric(n)),
    nrow = n, ncol = M))
  loggamma_z0 <- spec$log_gamma(0, W)
  loggamma_z1 <- spec$log_gamma(1, W)
  gz0 <- exp(loggamma_z0)
  gz1 <- exp(loggamma_z1)
  ## Normalize away the finite-grid tail error, then build the joint density the
  ## falsification moments ask for directly:
  ##   f(Y0=y, A=0 | z, X) = f(Y0=y | X) Pr(A=0 | Y0=y, z, X).
  ## For a separable alpha this reproduces the closed forms fY/(1 + delta gamma)
  ## and the corresponding Pr(A=1|z,X) exactly.
  fY0 <- fy / rowSums(fy * Wm)
  pA_grid_z0 <- plogis(sweep(logalpha_z0, 1L, loggamma_z0, "+"))
  pA_grid_z1 <- plogis(sweep(logalpha_z1, 1L, loggamma_z1, "+"))
  fYA0_z0 <- fY0 * (1 - pA_grid_z0)
  fYA0_z1 <- fY0 * (1 - pA_grid_z1)
  qA0_z0 <- rowSums(fYA0_z0 * Wm)
  qA0_z1 <- rowSums(fYA0_z1 * Wm)
  list(pZ1 = pZ1, pA1_z0 = 1 - qA0_z0, pA1_z1 = 1 - qA0_z1,
       fY_z0 = fYA0_z0 / qA0_z0, fY_z1 = fYA0_z1 / qA0_z1,
       fYA0_z0 = fYA0_z0, fYA0_z1 = fYA0_z1,
       delta = delta,
       alpha_z0 = exp(logalpha_z0), alpha_z1 = exp(logalpha_z1),
       logalpha_z0 = logalpha_z0, logalpha_z1 = logalpha_z1,
       gamma_z0 = gz0, gamma_z1 = gz1, fY0 = fY0)
}

## Validation DGPs from the competing papers, used by 02_check_competitors.R.
siv_dgp_liu2021 <- function(n) {
  C1 <- rbinom(n, 1, 0.4)
  C2 <- rbinom(n, 1, 0.6)
  Z <- rbinom(n, 1, plogis(0.2 + 0.4 * C1 - 0.5 * C2))
  Y0 <- rbinom(n, 1, plogis(0.6 + 0.8 * C1 - 2 * C2))
  Y1 <- rbinom(n, 1, plogis(0.7 - 0.3 * C1))
  A <- rbinom(n, 1,
              plogis(0.4 + 2 * Z + 0.8 * C1 - 0.6 * Y0 - 1.6 * C1 * Z))
  Y <- Y0 * (1 - A) + Y1 * A
  list(Y = Y, A = A, Z = Z, X = cbind(C1 = C1, C2 = C2),
       Y0 = Y0, Y1 = Y1, psi_true = -0.6)
}

siv_dgp_liu2021_cont <- function(n) {
  C1 <- rbinom(n, 1, 0.4)
  C2 <- rbinom(n, 1, 0.6)
  Z <- rbinom(n, 1, plogis(0.7 + 0.8 * C1 - C2))
  Y0 <- rnorm(n, 0.5 + C1 + 3 * C2, 1)
  Y1 <- rnorm(n, 1.1 - 1.3 * C1, 1)
  A <- rbinom(n, 1,
              plogis(-0.2 - 3 * Z - 3 * C1 + 0.3 * Y0 + 4 * C1 * Z))
  Y <- Y0 * (1 - A) + Y1 * A
  list(Y = Y, A = A, Z = Z, X = cbind(C1 = C1, C2 = C2),
       Y0 = Y0, Y1 = Y1, psi_true = 0.3)
}

siv_dgp_sun2021 <- function(n) {
  X1 <- rbinom(n, 1, 0.4)
  X2 <- rbinom(n, 1, 0.6)
  Z <- rbinom(n, 1,
              plogis(0.4 + 0.9 * X1 - 0.7 * X2 - 0.8 * X1 * X2))
  Y <- rbinom(n, 1, plogis(1 - 1.2 * X1 + 1.5 * X2))
  R <- rbinom(n, 1,
              plogis(-1.5 + 2.5 * Z + 0.8 * X1 - 1.2 * X2 + 1.8 * Y))
  list(Y = Y, Yobs = ifelse(R == 1, Y, NA_real_), R = R, Z = Z,
       X = cbind(X1 = X1, X2 = X2), gamma_true = 1.8, mu_true = 0.769)
}

siv_dgp_sun2025 <- function(n, type = c("binary", "continuous"),
                            gamma_true = 2) {
  type <- match.arg(type)
  if (type == "binary") {
    L <- chol(matrix(c(1, 0.2, 0.2, 1), 2, 2))
    U <- matrix(rnorm(2 * n), n, 2) %*% L
  } else {
    rtn <- function(k) qnorm(runif(
      k, pnorm(-1), pnorm(1)))
    U <- cbind(rtn(n), rtn(n))
  }
  U1 <- U[, 1L]; U2 <- U[, 2L]
  Z <- rbinom(n, 1,
              plogis(1 + 2 * U1 - U2 - 0.8 * U1 * U2))
  if (type == "binary") {
    Y <- rbinom(n, 1, plogis(0.5 - 2 * U1 + U2))
    R <- rbinom(n, 1,
                plogis(2 - 3 * Z + 0.8 * U1 + U2 + gamma_true * Y))
  } else {
    Y <- rnorm(n, 0.5 - 2 * U1 + U2, sqrt(0.25))
    R <- rbinom(n, 1,
                plogis(2 - 4 * Z + 0.8 * U1 + U2 + gamma_true * Y))
  }
  list(Y = Y, Yobs = ifelse(R == 1, Y, NA_real_), R = R, Z = Z,
       X = cbind(U1 = U1, U2 = U2), gamma_true = gamma_true, type = type)
}

## Automatic ATT truth from a large structural draw. pA integrates out the
## treatment Bernoulli draw, so this is less noisy than mean(effect[A == 1]).
siv_att_from_draw <- function(d) {
  effect <- d$Y1 - d$Y0
  sum(d$pA * effect) / sum(d$pA)
}

siv_att_truth <- function(type = c("continuous", "binary"),
                          n = 1000000L, seed = 20260902L, p = 3L) {
  type <- match.arg(type)
  set.seed(seed + match(type, c("continuous", "binary")))
  siv_att_from_draw(siv_dgp(n, type, p))
}

## Reference estimators. These deliberately remain outside SIV_Package.
siv_fit_2sls <- function(Y, A, Z, X, compute_se = TRUE) {
  X <- as.matrix(X)
  D <- cbind(`(Intercept)` = 1, A = A, X)
  Q <- cbind(`(Intercept)` = 1, Z = Z, X)
  qd <- crossprod(Q, D)
  b <- as.vector(solve(qd, crossprod(Q, Y)))
  if (!compute_se) return(c(estimate = b[2L], se = NA_real_))
  resid <- Y - as.vector(D %*% b)
  s2 <- sum(resid^2) / (length(Y) - ncol(D))
  iq <- solve(qd)
  V <- s2 * iq %*% crossprod(Q) %*% t(iq)
  c(estimate = b[2L], se = sqrt(V[2L, 2L]))
}

siv_fit_linear_model <- function(Y, A, Z, X) {
  dat <- data.frame(Y = Y, A = A, Z = Z, as.data.frame(X))
  rhs <- paste(c("A", "Z", colnames(X)), collapse = " + ")
  fit <- lm(as.formula(paste("Y ~", rhs)), data = dat)
  c(estimate = unname(coef(fit)["A"]),
    se = unname(sqrt(diag(vcov(fit)))["A"]))
}

## The ignorability benchmark: the cross-fitted EIF estimator of tau_Ign,
##
##   tau_Ign = sum_k sum_{i in I_k} {A_i - (1-A_i) Pr^(-k)(A_i=1|L_i)
##                                          / Pr^(-k)(A_i=0|L_i)}
##                                 {Y_i - E^(-k)(Y_i | A_i=1, L_i)}
##             / sum_i A_i ,
##
## with L = (Z, X) and both nuisances fitted by `glm()` on the other K-1 folds.
##
## It is the ATT one would report if A were assumed independent of Y(0) given L
## -- i.e. ignoring the selection on the treatment-free outcome that the
## separable-IV model allows.  Its bias against the true ATT is therefore the
## price of that assumption, which is exactly what the comparison is for.
##
## Note that the outcome regression cancels from the estimand: for any g(L),
## E[(1-A) {Pr(A=1|L)/Pr(A=0|L)} g(L)] = E[Pr(A=1|L) g(L)] = E[A g(L)],
## so the two g-terms cancel and only the treatment model has to be correct.
## The regression affects efficiency, not consistency; it is included as
## written above, fitted among the treated.
siv_fit_ignorability <- function(Y, A, Z, X, crossfit = 2L, trim = 1e-2) {
  X <- as.matrix(X)
  n <- length(Y)
  crossfit <- max(1L, as.integer(crossfit))
  df <- data.frame(Y = Y, A = as.integer(A), Z = Z, as.data.frame(X))
  rhs <- paste(c("Z", colnames(X)), collapse = " + ")
  f_a <- stats::as.formula(paste("A ~", rhs))
  f_y <- stats::as.formula(paste("Y ~", rhs))
  
  ## Folds stratified on (A, Z) so that every training set contains treated
  ## units (needed by the outcome regression) and both instrument arms.
  if (crossfit <= 1L) {
    folds <- list(seq_len(n))
  } else {
    cell <- interaction(df$A, df$Z, drop = TRUE)
    lab <- integer(n)
    for (lv in levels(cell)) {
      idx <- which(cell == lv)
      lab[idx] <- sample(rep_len(seq_len(crossfit), length(idx)))
    }
    folds <- split(seq_len(n), lab)
  }
  
  pA <- numeric(n)
  m1 <- numeric(n)
  for (te in folds) {
    tr <- if (length(folds) == 1L) te else setdiff(seq_len(n), te)
    fa <- stats::glm(f_a, data = df[tr, , drop = FALSE],
                     family = stats::binomial())
    pA[te] <- stats::predict(fa, newdata = df[te, , drop = FALSE],
                             type = "response")
    tr1 <- tr[df$A[tr] == 1L]
    fy <- stats::glm(f_y, data = df[tr1, , drop = FALSE],
                     family = stats::gaussian())
    m1[te] <- stats::predict(fy, newdata = df[te, , drop = FALSE],
                             type = "response")
  }
  ## Pr(A=1|L) -> 1 makes the odds diverge; bound it as SIV_LATE() bounds its
  ## own propensities.
  pA <- pmin(pmax(pA, trim), 1 - trim)
  
  num <- (A - (1 - A) * pA / (1 - pA)) * (Y - m1)
  estimate <- sum(num) / sum(A)
  eif <- (num - estimate * A) / mean(A)
  c(estimate = estimate, se = stats::sd(eif) / sqrt(n))
}

siv_compute_truths <- function(type, n = 300000L, seed = 20260902L) {
  set.seed(seed + match(type, c("continuous", "binary")))
  d <- siv_dgp(n, type)
  effect <- d$Y1 - d$Y0
  late_weight <- d$pA1 - d$pA0
  expected_y <- d$Y0 + d$pA * effect
  c(
    ATT = siv_att_from_draw(d),
    LATE = sum(late_weight * effect) / sum(late_weight),
    TSLS_projection = unname(siv_fit_2sls(
      expected_y, d$pA, d$Z, d$X, compute_se = FALSE)["estimate"])
  )
}

siv_summarise_results <- function(raw) {
  groups <- split(raw, interaction(raw$type, raw$method, drop = TRUE))
  out <- lapply(groups, function(g) {
    ## Preserve every finite estimate in the numerical summary. Diagnostics
    ## are flags, not missing-data instructions: an estimate that failed A3/A5
    ## must remain visible alongside n_ok/n_fail. Only genuine computation
    ## errors (non-finite estimate) are omitted from numerical summaries.
    computed <- is.finite(g$estimate)
    with_se <- computed & is.finite(g$se)
    ok <- g$ok & computed
    est <- g$estimate[computed]
    se <- g$se[with_se]
    truth <- g$truth[1L]
    data.frame(
      type = g$type[1L], method = g$method[1L], estimand = g$estimand[1L],
      truth = truth, n_computed = sum(computed),
      n_ok = sum(ok), n_fail = sum(!g$ok),
      mean = if (length(est)) mean(est) else NA_real_,
      bias = if (length(est) && is.finite(truth)) mean(est) - truth else NA_real_,
      empirical_sd = if (length(est) > 1L) sd(est) else NA_real_,
      rmse = if (length(est) && is.finite(truth))
        sqrt(mean((est - truth)^2)) else NA_real_,
      average_se = if (length(se)) mean(se) else NA_real_,
      coverage95 = if (any(with_se) && is.finite(truth))
        mean(abs(g$estimate[with_se] - truth) <=
               qnorm(0.975) * g$se[with_se]) else NA_real_,
      mean_seconds = mean(g$seconds),
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans[order(ans$type, ans$estimand, ans$method), ]
}
