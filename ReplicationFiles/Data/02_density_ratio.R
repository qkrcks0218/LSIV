files <- sort(list.files(path="./Ratio",
                         pattern = "JC_SIV_Density_B[0-9]+\\.rds$",
                         full.names = TRUE))
S <- length(files)
ref <- readRDS(files[1L])
grid <- ref$grid
n <- nrow(ref$fY_a0z0)
M <- ncol(ref$fY_a0z0)
rm(ref)
message(sprintf("%d cross-fitting repetition(s); N = %d, %d outcome values.",
                S, n, M))

row_median <- if (requireNamespace("matrixStats", quietly = TRUE)) {
  function(m) matrixStats::rowMedians(m)
} else {
  function(m) apply(m, 1L, stats::median)
}

## One density at a time keeps the memory use at N * M * S doubles.
median_density <- function(name) {
  stack <- matrix(NA_real_, n * M, S)
  for (s in seq_len(S)) {
    obj <- readRDS(files[s])
    if (!identical(dim(obj[[name]]), c(n, M)) ||
        !isTRUE(all.equal(obj$grid, grid))) {
      stop(files[s], " does not match the grid or dimension of ", files[1L],
           call. = FALSE)
    }
    stack[, s] <- as.vector(obj[[name]])
  }
  matrix(row_median(stack), n, M)
}

cells <- c("(0,0)" = "fY_a0z0", "(0,1)" = "fY_a0z1",
           "(1,0)" = "fY_a1z0", "(1,1)" = "fY_a1z1")
f_med <- list()
f_med[["(0,0)"]] <- median_density(cells[["(0,0)"]])
for (az in names(cells)[-1L]) f_med[[az]] <- median_density(cells[[az]])
ratio <- lapply(f_med, function(f) f / f_med[["(0,0)"]])

summary_table <- t(vapply(ratio, function(r) {
  stats::quantile(r, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
}, numeric(5L)))
colnames(summary_table) <- c("Min", "Q1", "Q2", "Q3", "Max")

out <- data.frame(az = rownames(summary_table), summary_table,
                  num_crossfit = S, row.names = NULL, check.names = FALSE)
write.csv(out, "JC_SIV_DensityRatio_Summary.csv", row.names = FALSE)
saveRDS(list(grid = grid, files = basename(files), f_med = f_med,
             ratio = ratio),
        "JC_SIV_DensityRatio.rds")

cat(sprintf("\nSummary statistics of r_i(y,a,z) across y in {%g,...,%g} and ",
            min(grid), max(grid)),
    sprintf("i in {1,...,%d}; median over %d cross-fitting repetition(s).\n\n",
            n, S), sep = "")
print(format(as.data.frame(summary_table), nsmall = 4L, digits = 4L),
      quote = FALSE)
