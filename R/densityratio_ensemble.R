## KLIEP/uLSIF density-ratio ensemble for
##   R_Y(y, x) = f(y | A = 0, Z = 1, x) / f(y | A = 0, Z = 0, x).
##
## This follows Simulation_Version1/1.Sim_Conti_SepIV.R, lines 1040 through 1464.
## KLIEP and uLSIF use the same high-distance-quantile bandwidth library as the
## replication code, are refitted within each outer training fold for inner
## validation, and are ensembled after conversion to candidate Z = 1
## conditional densities.  The original code normalises at the density level;
## doing that is essential because KLIEP/uLSIF estimate a joint (Y, X) ratio,
## whose X-only factor must disappear from each conditional density row.

## Internal: squared Euclidean distances without relying on an unexported
## densityratio helper.
siv_kliep_ulsif_squared_distance <- function(x, centers) {
  x <- as.matrix(x)
  centers <- as.matrix(centers)
  ans <- outer(rowSums(x^2), rowSums(centers^2), "+") -
    2 * tcrossprod(x, centers)
  pmax(ans, 0)
}

## Internal: construct a common-bandwidth candidate range from distance
## quantiles.
siv_kliep_ulsif_sigma_candidates <- function(
    x1, x0, kernel_num = 100L,
    probs = c(0.5, 0.75, 0.9)) {
  x1 <- as.matrix(x1)
  x0 <- as.matrix(x0)
  kernel_num <- min(as.integer(kernel_num), nrow(x1))
  centers <- if (kernel_num >= nrow(x1)) {
    x1
  } else {
    x1[sample.int(nrow(x1), kernel_num), , drop = FALSE]
  }
  sigma_factor <- stats::IQR(x1[, 1L]) / stats::IQR(x0[, 1L])
  if (!is.finite(sigma_factor) || sigma_factor <= 0) sigma_factor <- 1
  sigma <- sqrt(stats::quantile(
    siv_kliep_ulsif_squared_distance(x1, centers), probs = probs,
    names = FALSE, na.rm = TRUE
  )) * sigma_factor
  sigma <- unique(as.numeric(sigma[is.finite(sigma) & sigma > 1e-8]))
  if (!length(sigma)) sigma <- 1
  list(sigma = sigma, factor = sigma_factor)
}

## Internal: common transformation for separate Y and X bandwidths.  X is
## centred and scaled using the pooled outer-training sample, giving every
## nondegenerate X column unit standard deviation.  Y is only centred because
## its scale is controlled directly by the Y bandwidth.
siv_kliep_ulsif_transform_info <- function(x1, x0) {
  all <- rbind(as.data.frame(x1), as.data.frame(x0))
  y_center <- mean(all[[1L]])
  if (ncol(all) == 1L) {
    return(list(y_center = y_center, x_center = numeric(0),
                x_scale = numeric(0)))
  }
  X <- as.matrix(all[-1L])
  x_center <- colMeans(X)
  x_scale <- apply(X, 2L, stats::sd)
  x_scale[!is.finite(x_scale) | x_scale <= 0] <- 1
  list(y_center = y_center, x_center = x_center, x_scale = x_scale)
}

siv_kliep_ulsif_transform <- function(data, transform, h_y, h_x) {
  data <- as.data.frame(data)
  out <- data.frame(Y = (data[[1L]] - transform$y_center) / h_y)
  if (ncol(data) > 1L) {
    X <- sweep(as.matrix(data[-1L]), 2L, transform$x_center, "-")
    X <- sweep(X, 2L, transform$x_scale * h_x, "/")
    out <- cbind(out, as.data.frame(X, check.names = FALSE))
    names(out) <- names(data)
  }
  out
}

## Internal: pooled outer-training standardisation for a common bandwidth.
## Every Y/X coordinate is put on a unit-standard-deviation scale before the
## scalar bandwidth is selected, fitted, and used for prediction.
siv_kliep_ulsif_common_transform_info <- function(x1, x0) {
  all <- rbind(as.data.frame(x1), as.data.frame(x0))
  center <- colMeans(all)
  scale <- apply(all, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  list(center = center, scale = scale)
}

siv_kliep_ulsif_common_transform <- function(data, transform) {
  data <- as.data.frame(data)
  out <- sweep(as.matrix(data), 2L, transform$center, "-")
  out <- sweep(out, 2L, transform$scale, "/")
  out <- as.data.frame(out, check.names = FALSE)
  names(out) <- names(data)
  out
}

## Internal: separate bandwidth ranges for Y and standardised X.  Distance
## quantiles for each axis define the two-dimensional CV grid.
siv_kliep_ulsif_bandwidth_candidates <- function(
    x1, x0, kernel_num = 100L,
    probs = c(0.5, 0.75, 0.9)) {
  x1 <- as.data.frame(x1)
  x0 <- as.data.frame(x0)
  all <- rbind(x1, x0)
  transform <- siv_kliep_ulsif_transform_info(x1, x0)
  kernel_num <- min(as.integer(kernel_num), nrow(x1), nrow(all))
  center_index <- if (kernel_num >= nrow(all)) {
    seq_len(nrow(all))
  } else {
    sample.int(nrow(all), kernel_num)
  }
  centers <- all[center_index, , drop = FALSE]
  dy2 <- outer(x1[[1L]], centers[[1L]], "-")^2
  y_factor <- stats::IQR(x1[[1L]]) / stats::IQR(x0[[1L]])
  if (!is.finite(y_factor) || y_factor <= 0) y_factor <- 1
  h_y <- sqrt(stats::quantile(
    dy2, probs = probs, names = FALSE, na.rm = TRUE
  )) * y_factor
  if (ncol(all) > 1L) {
    X1 <- sweep(as.matrix(x1[-1L]), 2L, transform$x_center, "-")
    X1 <- sweep(X1, 2L, transform$x_scale, "/")
    Xc <- sweep(as.matrix(centers[-1L]), 2L, transform$x_center, "-")
    Xc <- sweep(Xc, 2L, transform$x_scale, "/")
    dx2 <- siv_kliep_ulsif_squared_distance(X1, Xc)
    h_x <- sqrt(stats::quantile(
      dx2, probs = probs, names = FALSE, na.rm = TRUE
    ))
  } else {
    h_x <- 1
  }
  h_y <- unique(as.numeric(h_y[is.finite(h_y) & h_y > 1e-8]))
  h_x <- unique(as.numeric(h_x[is.finite(h_x) & h_x > 1e-8]))
  if (!length(h_y)) h_y <- 1
  if (!length(h_x)) h_x <- 1
  list(
    grid = expand.grid(h_y = h_y, h_x = h_x),
    h_y = h_y, h_x = h_x, probs = probs, y_factor = y_factor,
    transform = transform, centers = centers
  )
}

## Internal: fit one anisotropic density-ratio candidate using only the public
## densityratio interface.  Dividing Y by h_y and each standardised X column by
## h_x is exactly the separate-bandwidth Gaussian kernel.
siv_kliep_ulsif_anisotropic_fit <- function(
    method, x1, x0, bandwidth, transform, kernel_num = 100L,
    fold = 5L, centers = NULL, cv = FALSE) {
  h_y <- unname(bandwidth[["h_y"]])
  h_x <- unname(bandwidth[["h_x"]])
  tx1 <- siv_kliep_ulsif_transform(x1, transform, h_y, h_x)
  tx0 <- siv_kliep_ulsif_transform(x0, transform, h_y, h_x)
  tc <- if (is.null(centers)) NULL else {
    siv_kliep_ulsif_transform(centers, transform, h_y, h_x)
  }
  kernel_num <- min(as.integer(kernel_num), nrow(tx1))
  fit <- if (identical(method, "kliep")) {
    tryCatch(
      densityratio::kliep(
        tx1, tx0, scale = NULL, sigma = 1, ncenters = kernel_num,
        centers = tc, cv = cv,
        nfold = max(2L, min(as.integer(fold), nrow(tx1))),
        progressbar = FALSE
      ),
      error = function(e) NULL
    )
  } else {
    tryCatch(
      densityratio::ulsif(
        tx1, tx0, scale = NULL, sigma = 1, ncenters = kernel_num,
        centers = tc, progressbar = FALSE
      ),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) return(NULL)
  structure(
    list(fit = fit, bandwidth = c(h_y = h_y, h_x = h_x),
         transform = transform),
    class = "siv_kliep_ulsif_anisotropic"
  )
}

## Internal: fit the two direct density-ratio candidates with separate Y and X
## bandwidth selection.
siv_kliep_ulsif_fit <- function(x1, x0, kernel_num = 100L, fold = 5L,
                                sigma = NULL,
                                cv_grid = c(0.5, 0.75, 0.9),
                                methods = c("kliep", "ulsif")) {
  x1 <- as.data.frame(x1)
  x0 <- as.data.frame(x0)
  if (!identical(names(x1), names(x0))) names(x0) <- names(x1)
  methods <- unique(match.arg(methods, c("kliep", "ulsif"),
                              several.ok = TRUE))
  kernel_num <- min(as.integer(kernel_num), nrow(x1))
  fold <- max(2L, min(as.integer(fold), nrow(x1)))
  if (!is.null(sigma)) {
    stop("Fixed scalar `sigma` is unavailable for the separate-bandwidth fit.",
         call. = FALSE)
  }
  bandwidth_info <- siv_kliep_ulsif_bandwidth_candidates(
    x1, x0, kernel_num, probs = cv_grid
  )
  candidate_grid <- bandwidth_info$grid
  cv_seed <- sample.int(.Machine$integer.max, 1L)
  fits <- vector("list", length(methods))
  names(fits) <- methods
  selected_bandwidth <- matrix(
    NA_real_, length(methods), 2L,
    dimnames = list(methods, c("h_y", "h_x"))
  )
  cv_score <- stats::setNames(
    replicate(length(methods), rep(NA_real_, nrow(candidate_grid)),
              simplify = FALSE), methods
  )
  for (method in names(fits)) {
    candidate_fits <- vector("list", nrow(candidate_grid))
    for (j in seq_len(nrow(candidate_grid))) {
      if (identical(method, "kliep")) set.seed(cv_seed)
      candidate_fits[[j]] <- siv_kliep_ulsif_anisotropic_fit(
        method, x1, x0, candidate_grid[j, ], bandwidth_info$transform,
        kernel_num = kernel_num, fold = fold,
        centers = bandwidth_info$centers,
        cv = identical(method, "kliep")
      )
      if (is.null(candidate_fits[[j]])) next
      object <- candidate_fits[[j]]$fit
      cv_score[[method]][j] <- if (identical(method, "kliep")) {
        max(object$cv_score)
      } else {
        -min(object$cv_score)
      }
    }
    if (all(!is.finite(cv_score[[method]]))) return(NULL)
    best <- which.max(cv_score[[method]])
    fits[[method]] <- candidate_fits[[best]]
    selected_bandwidth[method, ] <- as.numeric(candidate_grid[best, ])
  }
  list(
    fits = fits,
    sigma_candidates = candidate_grid,
    sigma_factor = bandwidth_info$y_factor,
    sigma_selected = selected_bandwidth[, "h_x"],
    selected_bandwidth = selected_bandwidth,
    bandwidth_range = list(h_y = bandwidth_info$h_y,
                           h_x = bandwidth_info$h_x,
                           probs = bandwidth_info$probs),
    cv_score = cv_score,
    transform = bandwidth_info$transform
  )
}

## Internal: scalar-bandwidth CV fit.  Y and every X coordinate are first
## standardised using the pooled outer-training sample, after which the same
## sigma is selected and applied to every coordinate.
siv_kliep_ulsif_YXcommonbw_fit <- function(
    x1, x0, kernel_num = 100L, fold = 5L,
    cv_grid = c(0.5, 0.75, 0.9),
    methods = c("kliep", "ulsif")) {
  x1 <- as.data.frame(x1)
  x0 <- as.data.frame(x0)
  if (!identical(names(x1), names(x0))) names(x0) <- names(x1)
  methods <- unique(match.arg(methods, c("kliep", "ulsif"),
                              several.ok = TRUE))
  kernel_num <- min(as.integer(kernel_num), nrow(x1))
  fold <- max(2L, min(as.integer(fold), nrow(x1)))
  transform <- siv_kliep_ulsif_common_transform_info(x1, x0)
  tx1 <- siv_kliep_ulsif_common_transform(x1, transform)
  tx0 <- siv_kliep_ulsif_common_transform(x0, transform)
  sigma_info <- siv_kliep_ulsif_sigma_candidates(
    tx1, tx0, kernel_num, probs = cv_grid
  )

  fits <- lapply(methods, function(method) {
    fit <- if (identical(method, "kliep")) {
      tryCatch(
        densityratio::kliep(
          tx1, tx0, scale = NULL, ncenters = kernel_num,
          sigma = sigma_info$sigma, cv = TRUE, nfold = fold,
          progressbar = FALSE
        ),
        error = function(e) NULL
      )
    } else {
      tryCatch(
        densityratio::ulsif(
          tx1, tx0, scale = NULL, ncenters = kernel_num,
          sigma = sigma_info$sigma, progressbar = FALSE
        ),
        error = function(e) NULL
      )
    }
    if (is.null(fit)) return(NULL)
    structure(
      list(fit = fit, transform = transform),
      class = "siv_kliep_ulsif_common"
    )
  })
  names(fits) <- methods
  if (any(vapply(fits, is.null, logical(1L)))) return(NULL)
  sigma_selected <- vapply(methods, function(method) {
    unname(fits[[method]]$fit$sigma_opt)
  }, numeric(1L))
  list(
    fits = fits,
    sigma_candidates = sigma_info$sigma,
    sigma_factor = sigma_info$factor,
    sigma_selected = sigma_selected,
    selected_bandwidth = NULL,
    bandwidth_range = NULL,
    transform = transform
  )
}

## Internal: fast KLIEP/uLSIF library with the median heuristic.  Both
## candidates standardise Y and X on the pooled outer-training sample before
## using one common median bandwidth.  uLSIF retains its inexpensive lambda
## selection.
siv_kliep_ulsif_median_fit <- function(
    x1, x0, kernel_num = 100L,
    methods = c("kliep", "ulsif")) {
  x1 <- as.data.frame(x1)
  x0 <- as.data.frame(x0)
  if (!identical(names(x1), names(x0))) names(x0) <- names(x1)
  methods <- unique(match.arg(methods, c("kliep", "ulsif"),
                              several.ok = TRUE))
  kernel_num <- min(as.integer(kernel_num), nrow(x1))
  transform <- siv_kliep_ulsif_common_transform_info(x1, x0)
  tx1 <- siv_kliep_ulsif_common_transform(x1, transform)
  tx0 <- siv_kliep_ulsif_common_transform(x0, transform)

  fits <- lapply(methods, function(method) {
    fit <- if (identical(method, "kliep")) {
      tryCatch(
        densityratio::kliep(
          tx1, tx0, scale = NULL, ncenters = kernel_num, nsigma = 1L,
          cv = FALSE, progressbar = FALSE
        ),
        error = function(e) NULL
      )
    } else {
      tryCatch(
        densityratio::ulsif(
          tx1, tx0, scale = NULL, ncenters = kernel_num, nsigma = 1L,
          progressbar = FALSE
        ),
        error = function(e) NULL
      )
    }
    if (is.null(fit)) return(NULL)
    structure(
      list(fit = fit, transform = transform),
      class = "siv_kliep_ulsif_common"
    )
  })
  names(fits) <- methods
  if (any(vapply(fits, is.null, logical(1L)))) return(NULL)
  sigma_selected <- vapply(methods, function(method) {
    fit <- fits[[method]]$fit
    if (identical(method, "kliep")) {
      unname(fit$sigma[1L])
    } else {
      unname(fit$sigma_opt)
    }
  }, numeric(1L))
  list(
    fits = fits,
    sigma_candidates = sigma_selected,
    sigma_factor = NA_real_,
    sigma_selected = sigma_selected,
    transform = transform
  )
}

## Internal: median-heuristic fit with distinct bandwidths for Y and the
## standardised X coordinates.
siv_kliep_ulsif_median_separate_fit <- function(
    x1, x0, kernel_num = 100L,
    methods = c("kliep", "ulsif")) {
  x1 <- as.data.frame(x1)
  x0 <- as.data.frame(x0)
  if (!identical(names(x1), names(x0))) names(x0) <- names(x1)
  methods <- unique(match.arg(methods, c("kliep", "ulsif"),
                              several.ok = TRUE))
  info <- siv_kliep_ulsif_bandwidth_candidates(
    x1, x0, kernel_num = kernel_num, probs = 0.5
  )
  bandwidth <- info$grid[1L, , drop = FALSE]
  fits <- lapply(methods, function(method) {
    siv_kliep_ulsif_anisotropic_fit(
      method, x1, x0, bandwidth = bandwidth,
      transform = info$transform, kernel_num = kernel_num,
      centers = info$centers, cv = FALSE
    )
  })
  names(fits) <- methods
  if (any(vapply(fits, is.null, logical(1L)))) return(NULL)
  selected_bandwidth <- matrix(
    rep(as.numeric(bandwidth), each = length(methods)),
    nrow = length(methods), byrow = FALSE,
    dimnames = list(methods, c("h_y", "h_x"))
  )
  list(
    fits = fits,
    sigma_candidates = bandwidth,
    sigma_factor = info$y_factor,
    sigma_selected = selected_bandwidth[, "h_x"],
    selected_bandwidth = selected_bandwidth,
    bandwidth_range = list(h_y = info$h_y, h_x = info$h_x,
                           probs = info$probs),
    transform = info$transform
  )
}

## Internal: refit one Version1 candidate with the bandwidth selected on the
## complete outer training fold.  This is the role of lines 1096 through 1138
## and 1196 through 1236 in the original script.
siv_kliep_ulsif_refit <- function(method, x1, x0, sigma, transform,
                                  kernel_num, fold) {
  x1 <- as.data.frame(x1)
  x0 <- as.data.frame(x0)
  if (!identical(names(x1), names(x0))) names(x0) <- names(x1)
  tx1 <- siv_kliep_ulsif_common_transform(x1, transform)
  tx0 <- siv_kliep_ulsif_common_transform(x0, transform)
  kernel_num <- min(as.integer(kernel_num), nrow(x1))
  fit <- if (identical(method, "kliep")) {
    tryCatch(
      densityratio::kliep(
        tx1, tx0, scale = NULL, ncenters = kernel_num, sigma = sigma,
        cv = FALSE,
        nfold = max(2L, min(as.integer(fold), nrow(x1))),
        progressbar = FALSE
      ),
      error = function(e) NULL
    )
  } else {
    tryCatch(
      densityratio::ulsif(
        tx1, tx0, scale = NULL, ncenters = kernel_num, sigma = sigma,
        progressbar = FALSE
      ),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) return(NULL)
  structure(
    list(fit = fit, transform = transform),
    class = "siv_kliep_ulsif_common"
  )
}

## Internal: refit one separate-bandwidth candidate for ensemble-weight CV.
siv_kliep_ulsif_anisotropic_refit <- function(
    method, x1, x0, bandwidth, transform, kernel_num, fold) {
  siv_kliep_ulsif_anisotropic_fit(
    method, x1, x0, bandwidth = bandwidth, transform = transform,
    kernel_num = kernel_num, fold = fold, centers = NULL, cv = FALSE
  )
}

siv_kliep_ulsif_predict <- function(object, newdata) {
  if (inherits(object, "siv_kliep_ulsif_common")) {
    newdata <- siv_kliep_ulsif_common_transform(newdata, object$transform)
    object <- object$fit
  } else if (inherits(object, "siv_kliep_ulsif_anisotropic")) {
    newdata <- siv_kliep_ulsif_transform(
      newdata, object$transform,
      h_y = object$bandwidth[["h_y"]],
      h_x = object$bandwidth[["h_x"]]
    )
    object <- object$fit
  }
  ans <- tryCatch(
    as.numeric(stats::predict(object, newdata = as.data.frame(newdata))),
    error = function(e) rep(NA_real_, nrow(newdata))
  )
  ans
}

## Internal: fit one direct density-ratio method on every outer training fold.
siv_single_density_ratio_grid <- function(
    Y, A, Z, X, grid, weights, folds, baseline_density,
    method = c("kliep", "ulsif"),
    bandwidth = c("cv_separate", "median_separate", "cv_common",
                  "median_common"),
    cv_grid = c(0.5, 0.75, 0.9),
    kernel_num = 100L, fold = 5L, dr_floor = 1e-4,
    density_floor = 1e-4, min_arm = 10L) {
  method <- match.arg(method)
  bandwidth <- match.arg(bandwidth)
  if (!requireNamespace("densityratio", quietly = TRUE)) {
    stop("`ratio_method = \"", method,
         "\"` requires the `densityratio` package.", call. = FALSE)
  }
  n <- length(Y)
  M <- length(grid)
  X <- as.matrix(X)
  xnames <- colnames(X)
  if (is.null(xnames)) xnames <- sprintf("X%02d", seq_len(ncol(X)))
  df <- data.frame(Y = Y, X, check.names = FALSE)
  names(df) <- c("Y", xnames)
  baseline_density <- siv_density_normalize(
    baseline_density, weights, density_floor
  )
  ratio <- matrix(NA_real_, n, M)
  arm_sizes <- matrix(NA_integer_, length(folds), 2L,
                      dimnames = list(NULL, c("A0Z1", "A0Z0")))
  details <- vector("list", length(folds))

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (!length(tr)) tr <- te
    i1 <- intersect(tr, which(A == 0 & Z == 1))
    i0 <- intersect(tr, which(A == 0 & Z == 0))
    arm_sizes[k, ] <- c(length(i1), length(i0))
    if (length(i1) < min_arm || length(i0) < min_arm) {
      return(NULL)
    }
    fit <- switch(
      bandwidth,
      cv_separate = siv_kliep_ulsif_fit(
        df[i1, , drop = FALSE], df[i0, , drop = FALSE],
        kernel_num = kernel_num, fold = fold, cv_grid = cv_grid,
        methods = method
      ),
      median_separate = siv_kliep_ulsif_median_separate_fit(
        df[i1, , drop = FALSE], df[i0, , drop = FALSE],
        kernel_num = kernel_num, methods = method
      ),
      cv_common = siv_kliep_ulsif_YXcommonbw_fit(
        df[i1, , drop = FALSE], df[i0, , drop = FALSE],
        kernel_num = kernel_num, fold = fold, cv_grid = cv_grid,
        methods = method
      ),
      median_common = siv_kliep_ulsif_median_fit(
        df[i1, , drop = FALSE], df[i0, , drop = FALSE],
        kernel_num = kernel_num, methods = method
      )
    )
    if (is.null(fit)) return(NULL)
    nte <- length(te)
    big <- data.frame(
      Y = rep(grid, each = nte),
      X[rep(te, times = M), , drop = FALSE],
      check.names = FALSE
    )
    names(big) <- c("Y", xnames)
    prediction <- siv_kliep_ulsif_predict(fit$fits[[method]], big)
    if (any(!is.finite(prediction))) return(NULL)
    f1 <- siv_ratio_candidate_density(
      matrix(prediction, nte, M),
      baseline_density[te, , drop = FALSE], weights,
      density_floor, dr_floor
    )
    ratio[te, ] <- f1 / baseline_density[te, , drop = FALSE]
    details[[k]] <- fit
  }
  if (anyNA(ratio) || any(!is.finite(ratio))) return(NULL)
  list(
    ratio = ratio,
    weights = matrix(1, length(folds), 1L,
                     dimnames = list(NULL, method)),
    mse_weights = matrix(1, length(folds), 1L,
                         dimnames = list(NULL, method)),
    kl_weights = matrix(1, length(folds), 1L,
                        dimnames = list(NULL, method)),
    sigma_candidates = lapply(details, `[[`, "sigma_candidates"),
    sigma_factor = vapply(details, function(x) x$sigma_factor, numeric(1L)),
    selected_sigma = lapply(details, `[[`, "sigma_selected"),
    selected_bandwidth = lapply(details, `[[`, "selected_bandwidth"),
    bandwidth_range = lapply(details, `[[`, "bandwidth_range"),
    arm_sizes = arm_sizes,
    candidates = method,
    ratio_method = method,
    ratio_bandwidth = bandwidth,
    cv_grid = if (grepl("^cv_", bandwidth)) cv_grid else NULL,
    method = paste("cross-fitted", method, "density ratio")
  )
}

## Internal: stratified inner folds on the entire outer training fold.  Thus
## both density-ratio arms occur in every training/validation split, as in the
## cell-balanced Version1 splitting scheme.
siv_ratio_inner_folds <- function(index, A, Z, nfold) {
  nfold <- max(2L, as.integer(nfold))
  local <- siv_make_folds(length(index), A[index], Z[index], nfold)
  lapply(local, function(ii) index[ii])
}

## Internal: turn a joint density-ratio candidate into a conditional Z = 1
## density.  Row normalisation removes the unknown X-only factor.
siv_ratio_candidate_density <- function(ratio, f0, weights,
                                        density_floor, dr_floor) {
  ratio <- as.matrix(ratio)
  ratio[!is.finite(ratio)] <- dr_floor
  candidate <- ratio * f0
  candidate <- siv_density_normalize(candidate, weights,
                                     dr_floor * density_floor)
  candidate
}

## Internal: fast cross-fitted median-heuristic KLIEP/uLSIF ensemble.  Both
## direct ratio candidates are always used.  Their bandwidths come from the
## median heuristic, while their ensemble weights minimize inner-validation
## MSPE (also known as ISE) at the conditional-density level.
siv_kliep_ulsif_median_ratio_grid <- function(
    Y, A, Z, X, grid, weights, folds, baseline_density,
    kernel_num = 100L, fold = 5L, dr_floor = 1e-4,
    density_floor = 1e-4, min_arm = 10L, require_complete = TRUE) {
  if (!requireNamespace("densityratio", quietly = TRUE)) {
    stop("`ratio_bandwidth = \"median_common\"` requires the ",
         "`densityratio` package.", call. = FALSE)
  }
  n <- length(Y)
  M <- length(grid)
  X <- as.matrix(X)
  xnames <- colnames(X)
  if (is.null(xnames)) xnames <- sprintf("X%02d", seq_len(ncol(X)))
  df <- data.frame(Y = Y, X, check.names = FALSE)
  names(df) <- c("Y", xnames)
  baseline_density <- siv_density_normalize(
    baseline_density, weights, density_floor
  )

  RY <- matrix(NA_real_, n, M)
  ensemble_weights <- mspe_weights <- matrix(
    NA_real_, length(folds), 2L,
    dimnames = list(NULL, c("kliep", "ulsif"))
  )
  selected_sigma <- matrix(
    NA_real_, length(folds), 2L,
    dimnames = list(NULL, c("kliep", "ulsif"))
  )
  sigma_candidates <- vector("list", length(folds))
  arm_sizes <- matrix(
    NA_integer_, length(folds), 2L,
    dimnames = list(NULL, c("A0Z1", "A0Z0"))
  )

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (!length(tr)) tr <- te
    i1 <- intersect(tr, which(A == 0 & Z == 1))
    i0 <- intersect(tr, which(A == 0 & Z == 0))
    arm_sizes[k, ] <- c(length(i1), length(i0))
    if (length(i1) < min_arm || length(i0) < min_arm) {
      return(NULL)
    }

    full <- siv_kliep_ulsif_median_fit(
      df[i1, , drop = FALSE], df[i0, , drop = FALSE],
      kernel_num = kernel_num
    )
    if (is.null(full)) return(NULL)
    sigma_candidates[[k]] <- full$sigma_candidates
    selected_sigma[k, ] <- full$sigma_selected

    ## Estimate only the library weights by inner validation.  In contrast to
    ## the slow path, no bandwidth candidate search is repeated here.
    inner <- siv_ratio_inner_folds(tr, A, Z, fold)
    cv_ratio <- list(
      kliep = matrix(NA_real_, length(tr), M,
                     dimnames = list(as.character(tr), NULL)),
      ulsif = matrix(NA_real_, length(tr), M,
                     dimnames = list(as.character(tr), NULL))
    )
    for (v in seq_along(inner)) {
      va <- inner[[v]]
      inn_tr <- setdiff(tr, va)
      vi1 <- intersect(inn_tr, which(A == 0 & Z == 1))
      vi0 <- intersect(inn_tr, which(A == 0 & Z == 0))
      if (length(vi1) < min_arm || length(vi0) < min_arm) return(NULL)
      nva <- length(va)
      big <- data.frame(
        Y = rep(grid, each = nva),
        X[rep(va, times = M), , drop = FALSE],
        check.names = FALSE
      )
      names(big) <- c("Y", xnames)
      for (method in c("kliep", "ulsif")) {
        inner_fit <- siv_kliep_ulsif_refit(
          method, df[vi1, , drop = FALSE], df[vi0, , drop = FALSE],
          sigma = full$sigma_selected[[method]],
          transform = full$transform,
          kernel_num = kernel_num, fold = fold
        )
        if (is.null(inner_fit)) return(NULL)
        pred <- siv_kliep_ulsif_predict(inner_fit, big)
        if (any(!is.finite(pred))) return(NULL)
        cv_ratio[[method]][as.character(va), ] <- matrix(pred, nva, M)
      }
    }

    f0_cv <- baseline_density[tr, , drop = FALSE]
    density_cv <- lapply(cv_ratio, function(r) {
      siv_ratio_candidate_density(
        r[as.character(tr), , drop = FALSE], f0_cv, weights,
        density_floor, dr_floor
      )
    })
    target_pos <- which(A[tr] == 0 & Z[tr] == 1)
    if (!length(target_pos)) return(NULL)
    W <- matrix(weights, length(target_pos), M, byrow = TRUE)
    H <- matrix(0, 2L, 2L)
    for (j in 1:2) {
      for (l in j:2) {
        H[j, l] <- H[l, j] <- mean(rowSums(
          density_cv[[j]][target_pos, , drop = FALSE] *
            density_cv[[l]][target_pos, , drop = FALSE] * W
        ))
      }
    }
    observed <- cbind(
      kliep = siv_density_at_observed(
        density_cv$kliep[target_pos, , drop = FALSE], grid, Y[tr][target_pos]
      ),
      ulsif = siv_density_at_observed(
        density_cv$ulsif[target_pos, , drop = FALSE], grid, Y[tr][target_pos]
      )
    )
    ew <- siv_ratio_ensemble_weights(
      H, colMeans(observed)
    )
    ensemble_weights[k, ] <- ew$weight
    mspe_weights[k, ] <- ew$mspe

    nte <- length(te)
    big <- data.frame(
      Y = rep(grid, each = nte),
      X[rep(te, times = M), , drop = FALSE],
      check.names = FALSE
    )
    names(big) <- c("Y", xnames)
    candidate_test <- lapply(full$fits, function(object) {
      pred <- siv_kliep_ulsif_predict(object, big)
      if (any(!is.finite(pred))) return(NULL)
      siv_ratio_candidate_density(
        matrix(pred, nte, M), baseline_density[te, , drop = FALSE],
        weights, density_floor, dr_floor
      )
    })
    if (any(vapply(candidate_test, is.null, logical(1L)))) return(NULL)
    f1 <- ew$weight[1L] * candidate_test$kliep +
      ew$weight[2L] * candidate_test$ulsif
    f1 <- siv_density_normalize(f1, weights, density_floor)
    RY[te, ] <- f1 / baseline_density[te, , drop = FALSE]
  }

  checked <- if (isTRUE(require_complete)) {
    seq_len(n)
  } else {
    sort(unique(unlist(folds, use.names = FALSE)))
  }
  if (anyNA(RY[checked, , drop = FALSE]) ||
      any(!is.finite(RY[checked, , drop = FALSE]))) return(NULL)
  list(
    ratio = RY,
    weights = ensemble_weights,
    mspe_weights = mspe_weights,
    mse_weights = mspe_weights,
    kl_weights = NULL,
    sigma_candidates = sigma_candidates,
    sigma_factor = rep(NA_real_, length(folds)),
    selected_sigma = selected_sigma,
    arm_sizes = arm_sizes,
    candidates = c("kliep", "ulsif"),
    ratio_method = "kliep_ulsif",
    ratio_bandwidth = "median_common",
    ensemble_loss = "mspe",
    method = paste("cross-fitted CV-weighted median-heuristic",
                   "KLIEP/uLSIF density ensemble")
  )
}

## Internal: the full Version1 KLIEP/uLSIF density-ratio ensemble on the SIV
## outcome grid.  Returns NULL when a fold cannot support both candidates, so
## SIV() can use its documented logistic fallback.
siv_kliep_ulsif_ratio_grid <- function(
    Y, A, Z, X, grid, weights, folds, baseline_density,
    kernel_num = 100L, fold = 5L, dr_floor = 1e-4,
    density_floor = 1e-4, min_arm = 10L, require_complete = TRUE,
    bandwidth_mode = c("cv_separate", "median_separate", "cv_common"),
    cv_grid = c(0.5, 0.75, 0.9)) {
  bandwidth_mode <- match.arg(bandwidth_mode)
  if (!requireNamespace("densityratio", quietly = TRUE)) {
    stop("`ratio_method = \"kliep_ulsif\"` requires the `densityratio` package.",
         call. = FALSE)
  }
  n <- length(Y)
  M <- length(grid)
  X <- as.matrix(X)
  xnames <- colnames(X)
  if (is.null(xnames)) xnames <- sprintf("X%02d", seq_len(ncol(X)))
  df <- data.frame(Y = Y, X, check.names = FALSE)
  names(df) <- c("Y", xnames)
  baseline_density <- siv_density_normalize(
    baseline_density, weights, density_floor
  )

  RY <- matrix(NA_real_, n, M)
  ensemble_weights <- mspe_weights <- matrix(
    NA_real_, length(folds), 2L,
    dimnames = list(NULL, c("kliep", "ulsif"))
  )
  selected_sigma <- matrix(
    NA_real_, length(folds), 2L,
    dimnames = list(NULL, c("kliep", "ulsif"))
  )
  sigma_candidates <- vector("list", length(folds))
  selected_bandwidth <- bandwidth_range <- vector("list", length(folds))
  sigma_factor <- numeric(length(folds))
  arm_sizes <- matrix(NA_integer_, length(folds), 2L,
                      dimnames = list(NULL, c("A0Z1", "A0Z0")))

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)
    if (!length(tr)) tr <- te
    i1 <- intersect(tr, which(A == 0 & Z == 1))
    i0 <- intersect(tr, which(A == 0 & Z == 0))
    arm_sizes[k, ] <- c(length(i1), length(i0))
    if (length(i1) < min_arm || length(i0) < min_arm) {
      return(NULL)
    }

    full <- if (identical(bandwidth_mode, "cv_separate")) {
      siv_kliep_ulsif_fit(
        df[i1, , drop = FALSE], df[i0, , drop = FALSE],
        kernel_num = kernel_num, fold = fold, cv_grid = cv_grid
      )
    } else if (identical(bandwidth_mode, "median_separate")) {
      siv_kliep_ulsif_median_separate_fit(
        df[i1, , drop = FALSE], df[i0, , drop = FALSE],
        kernel_num = kernel_num
      )
    } else {
      siv_kliep_ulsif_YXcommonbw_fit(
        df[i1, , drop = FALSE], df[i0, , drop = FALSE],
        kernel_num = kernel_num, fold = fold, cv_grid = cv_grid
      )
    }
    if (is.null(full)) return(NULL)
    sigma_candidates[[k]] <- full$sigma_candidates
    sigma_factor[k] <- full$sigma_factor
    selected_sigma[k, ] <- full$sigma_selected
    selected_bandwidth[k] <- list(full$selected_bandwidth)
    bandwidth_range[k] <- list(full$bandwidth_range)

    inner <- siv_ratio_inner_folds(tr, A, Z, fold)
    cv_ratio <- list(
      kliep = matrix(NA_real_, length(tr), M,
                     dimnames = list(as.character(tr), NULL)),
      ulsif = matrix(NA_real_, length(tr), M,
                     dimnames = list(as.character(tr), NULL))
    )
    for (v in seq_along(inner)) {
      va <- inner[[v]]
      inn_tr <- setdiff(tr, va)
      vi1 <- intersect(inn_tr, which(A == 0 & Z == 1))
      vi0 <- intersect(inn_tr, which(A == 0 & Z == 0))
      if (length(vi1) < min_arm || length(vi0) < min_arm) return(NULL)
      nva <- length(va)
      big <- data.frame(
        Y = rep(grid, each = nva),
        X[rep(va, times = M), , drop = FALSE],
        check.names = FALSE
      )
      names(big) <- c("Y", xnames)
      for (method in c("kliep", "ulsif")) {
        inner_fit <- if (!identical(bandwidth_mode, "cv_common")) {
          siv_kliep_ulsif_anisotropic_refit(
            method, df[vi1, , drop = FALSE], df[vi0, , drop = FALSE],
            bandwidth = full$selected_bandwidth[method, ],
            transform = full$transform,
            kernel_num = kernel_num, fold = fold
          )
        } else {
          siv_kliep_ulsif_refit(
            method, df[vi1, , drop = FALSE], df[vi0, , drop = FALSE],
            sigma = full$sigma_selected[[method]],
            transform = full$transform,
            kernel_num = kernel_num, fold = fold
          )
        }
        if (is.null(inner_fit)) return(NULL)
        pred <- siv_kliep_ulsif_predict(inner_fit, big)
        if (any(!is.finite(pred))) return(NULL)
        cv_ratio[[method]][as.character(va), ] <- matrix(pred, nva, M)
      }
    }

    f0_cv <- baseline_density[tr, , drop = FALSE]
    density_cv <- lapply(cv_ratio, function(r) {
      siv_ratio_candidate_density(
        r[as.character(tr), , drop = FALSE], f0_cv, weights,
        density_floor, dr_floor
      )
    })
    target_pos <- which(A[tr] == 0 & Z[tr] == 1)
    if (!length(target_pos)) return(NULL)
    W <- matrix(weights, length(target_pos), M, byrow = TRUE)
    H <- matrix(0, 2L, 2L)
    for (j in 1:2) {
      for (l in j:2) {
        H[j, l] <- H[l, j] <- mean(rowSums(
          density_cv[[j]][target_pos, , drop = FALSE] *
            density_cv[[l]][target_pos, , drop = FALSE] * W
        ))
      }
    }
    observed <- cbind(
      kliep = siv_density_at_observed(
        density_cv$kliep[target_pos, , drop = FALSE], grid, Y[tr][target_pos]
      ),
      ulsif = siv_density_at_observed(
        density_cv$ulsif[target_pos, , drop = FALSE], grid, Y[tr][target_pos]
      )
    )
    ew <- siv_ratio_ensemble_weights(
      H, colMeans(observed)
    )
    ensemble_weights[k, ] <- ew$weight
    mspe_weights[k, ] <- ew$mspe

    nte <- length(te)
    big <- data.frame(
      Y = rep(grid, each = nte),
      X[rep(te, times = M), , drop = FALSE],
      check.names = FALSE
    )
    names(big) <- c("Y", xnames)
    candidate_test <- lapply(full$fits, function(object) {
      pred <- siv_kliep_ulsif_predict(object, big)
      if (any(!is.finite(pred))) return(NULL)
      siv_ratio_candidate_density(
        matrix(pred, nte, M), baseline_density[te, , drop = FALSE],
        weights, density_floor, dr_floor
      )
    })
    if (any(vapply(candidate_test, is.null, logical(1L)))) return(NULL)
    f1 <- ew$weight[1L] * candidate_test$kliep +
      ew$weight[2L] * candidate_test$ulsif
    f1 <- siv_density_normalize(f1, weights, density_floor)
    ratio <- f1 / baseline_density[te, , drop = FALSE]
    ratio <- pmax(ratio, dr_floor)
    f1 <- siv_density_normalize(
      ratio * baseline_density[te, , drop = FALSE], weights, density_floor
    )
    RY[te, ] <- f1 / baseline_density[te, , drop = FALSE]
  }

  checked <- if (isTRUE(require_complete)) {
    seq_len(n)
  } else {
    sort(unique(unlist(folds, use.names = FALSE)))
  }
  if (anyNA(RY[checked, , drop = FALSE]) ||
      any(!is.finite(RY[checked, , drop = FALSE]))) return(NULL)
  list(
    ratio = RY,
    weights = ensemble_weights,
    mspe_weights = mspe_weights,
    mse_weights = mspe_weights,
    kl_weights = NULL,
    sigma_candidates = sigma_candidates,
    sigma_factor = sigma_factor,
    selected_sigma = selected_sigma,
    selected_bandwidth = selected_bandwidth,
    bandwidth_range = bandwidth_range,
    arm_sizes = arm_sizes,
    candidates = c("kliep", "ulsif"),
    ratio_method = "kliep_ulsif",
    ratio_bandwidth = bandwidth_mode,
    ensemble_loss = "mspe",
    cv_grid = if (grepl("^cv_", bandwidth_mode)) cv_grid else NULL,
    method = if (identical(bandwidth_mode, "cv_separate")) {
      paste("cross-fitted separate-bandwidth KLIEP/uLSIF",
            "density ensemble")
    } else if (identical(bandwidth_mode, "median_separate")) {
      paste("cross-fitted median separate-bandwidth KLIEP/uLSIF",
            "density ensemble")
    } else {
      paste("cross-fitted common Y/X-bandwidth KLIEP/uLSIF",
            "density ensemble")
    }
  )
}

siv_kliep_ulsif_YXcommonbw_ratio_grid <- function(...) {
  siv_kliep_ulsif_ratio_grid(..., bandwidth_mode = "cv_common")
}
