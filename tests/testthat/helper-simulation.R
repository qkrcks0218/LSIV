## Test-only data generators. These are fixtures for estimator tests and are
## deliberately not part of the installed SIV namespace. The simulation source
## used for analysis lives in Simulation_Version2_based_on_SIV_Package/00_common.R.

test_alpha_shape <- function(y, type = c("continuous", "binary")) {
  type <- match.arg(type)
  if (type == "binary") return(2 * y)
  inside <- (y >= 0) & (y <= 2)
  out <- numeric(length(y))
  out[inside] <- y[inside] * (y[inside] - 2)
  out[!inside] <- 1 - (y[!inside] - 1)^(-2)
  out
}

test_dgp <- function(n, type = c("continuous", "binary"), p = 3L) {
  type <- match.arg(type)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- sprintf("X%02d", seq_len(p))
  W <- rowSums(X)
  Z <- stats::rbinom(n, 1, stats::plogis(0.1 * W))
  if (type == "continuous") {
    Y0 <- 0.1 * W + stats::rnorm(n)
    Y1 <- 1 + 0.1 * W + 0.25 * stats::rnorm(n)
  } else {
    Y0 <- stats::rbinom(n, 1, stats::plogis(0.2 * W - 0.25))
    Y1 <- stats::rbinom(n, 1, stats::plogis(0.2 * W + 0.25))
  }
  loggamma <- -2 + 4 * Z + 0.1 * W
  logdelta <- (0.4 + 0.05 * W) * test_alpha_shape(Y0, type)
  A <- stats::rbinom(n, 1, stats::plogis(logdelta + loggamma))
  Y <- (1 - A) * Y0 + A * Y1
  list(Y = Y, A = A, Z = Z, X = X, Y0 = Y0, Y1 = Y1, W = W,
       logdelta = logdelta, loggamma = loggamma, type = type)
}

test_oracle_nuisance <- function(X, grid, weights,
                                 type = c("continuous", "binary")) {
  type <- match.arg(type)
  X <- as.matrix(X)
  n <- nrow(X); M <- length(grid)
  W <- rowSums(X)
  pZ1 <- stats::plogis(0.1 * W)
  fy <- matrix(0, n, M)
  if (type == "continuous") {
    for (m in seq_len(M)) fy[, m] <- stats::dnorm(grid[m], 0.1 * W, 1)
  } else {
    py1 <- stats::plogis(0.2 * W - 0.25)
    for (m in seq_len(M)) fy[, m] <- ifelse(grid[m] == 1, py1, 1 - py1)
  }
  delta <- exp(outer(0.4 + 0.05 * W, test_alpha_shape(grid, type)))
  Wm <- matrix(weights, n, M, byrow = TRUE)
  gz0 <- exp(-2 + 0.1 * W)
  gz1 <- exp(2 + 0.1 * W)
  num0 <- fy / (1 + delta * gz0)
  num1 <- fy / (1 + delta * gz1)
  fY_z0 <- num0 / rowSums(num0 * Wm)
  fY_z1 <- num1 / rowSums(num1 * Wm)
  pA1_z0 <- rowSums(fy * (delta * gz0) / (1 + delta * gz0) * Wm) /
    rowSums(fy * Wm)
  pA1_z1 <- rowSums(fy * (delta * gz1) / (1 + delta * gz1) * Wm) /
    rowSums(fy * Wm)
  list(pZ1 = pZ1, pA1_z0 = pA1_z0, pA1_z1 = pA1_z1,
       fY_z0 = fY_z0, fY_z1 = fY_z1, delta = delta,
       gamma_z0 = gz0, gamma_z1 = gz1, fY0 = fy)
}

test_dgp_liu2021 <- function(n) {
  C1 <- stats::rbinom(n, 1, 0.4)
  C2 <- stats::rbinom(n, 1, 0.6)
  Z <- stats::rbinom(n, 1, stats::plogis(0.2 + 0.4 * C1 - 0.5 * C2))
  Y0 <- stats::rbinom(n, 1, stats::plogis(0.6 + 0.8 * C1 - 2 * C2))
  Y1 <- stats::rbinom(n, 1, stats::plogis(0.7 - 0.3 * C1))
  A <- stats::rbinom(n, 1,
                     stats::plogis(0.4 + 2 * Z + 0.8 * C1 - 0.6 * Y0 - 1.6 * C1 * Z))
  Y <- Y0 * (1 - A) + Y1 * A
  list(Y = Y, A = A, Z = Z, X = cbind(C1 = C1, C2 = C2),
       Y0 = Y0, Y1 = Y1, psi_true = -0.6)
}

test_dgp_liu2021_cont <- function(n) {
  C1 <- stats::rbinom(n, 1, 0.4)
  C2 <- stats::rbinom(n, 1, 0.6)
  Z <- stats::rbinom(n, 1, stats::plogis(0.7 + 0.8 * C1 - C2))
  Y0 <- stats::rnorm(n, 0.5 + C1 + 3 * C2, 1)
  Y1 <- stats::rnorm(n, 1.1 - 1.3 * C1, 1)
  A <- stats::rbinom(n, 1,
                     stats::plogis(-0.2 - 3 * Z - 3 * C1 + 0.3 * Y0 + 4 * C1 * Z))
  Y <- Y0 * (1 - A) + Y1 * A
  list(Y = Y, A = A, Z = Z, X = cbind(C1 = C1, C2 = C2),
       Y0 = Y0, Y1 = Y1, psi_true = 0.3)
}

test_dgp_sun2021 <- function(n) {
  X1 <- stats::rbinom(n, 1, 0.4)
  X2 <- stats::rbinom(n, 1, 0.6)
  Z <- stats::rbinom(n, 1,
                     stats::plogis(0.4 + 0.9 * X1 - 0.7 * X2 - 0.8 * X1 * X2))
  Y <- stats::rbinom(n, 1, stats::plogis(1 - 1.2 * X1 + 1.5 * X2))
  R <- stats::rbinom(n, 1,
                     stats::plogis(-1.5 + 2.5 * Z + 0.8 * X1 - 1.2 * X2 + 1.8 * Y))
  list(Y = Y, Yobs = ifelse(R == 1, Y, NA_real_), R = R, Z = Z,
       X = cbind(X1 = X1, X2 = X2), gamma_true = 1.8, mu_true = 0.769)
}

test_dgp_sun2025 <- function(n, type = c("binary", "continuous"),
                             gamma_true = 2) {
  type <- match.arg(type)
  if (type == "binary") {
    L <- chol(matrix(c(1, 0.2, 0.2, 1), 2, 2))
    U <- matrix(stats::rnorm(2 * n), n, 2) %*% L
  } else {
    rtn <- function(k) stats::qnorm(stats::runif(
      k, stats::pnorm(-1), stats::pnorm(1)))
    U <- cbind(rtn(n), rtn(n))
  }
  U1 <- U[, 1L]; U2 <- U[, 2L]
  Z <- stats::rbinom(n, 1,
                     stats::plogis(1 + 2 * U1 - U2 - 0.8 * U1 * U2))
  if (type == "binary") {
    Y <- stats::rbinom(n, 1, stats::plogis(0.5 - 2 * U1 + U2))
    R <- stats::rbinom(n, 1,
                       stats::plogis(2 - 3 * Z + 0.8 * U1 + U2 + gamma_true * Y))
  } else {
    Y <- stats::rnorm(n, 0.5 - 2 * U1 + U2, sqrt(0.25))
    R <- stats::rbinom(n, 1,
                       stats::plogis(2 - 4 * Z + 0.8 * U1 + U2 + gamma_true * Y))
  }
  list(Y = Y, Yobs = ifelse(R == 1, Y, NA_real_), R = R, Z = Z,
       X = cbind(U1 = U1, U2 = U2), gamma_true = gamma_true, type = type)
}
