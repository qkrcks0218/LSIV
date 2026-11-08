## `in_file` / `out_dir` may already be set by a driver script; fall back to the
## merged file in the working directory so this script also runs on its own.
if (!exists("in_file")) in_file <- "Merge_Sim_1snp_F.csv"
if (!exists("out_dir")) out_dir <- "."

raw <- read.csv(in_file)
need <- c("method", "truth", "estimate", "se", "type", "n")
miss <- setdiff(need, names(raw))
if (length(miss)) {
  stop("`", basename(in_file), "` is missing column(s): ",
       paste(miss, collapse = ", "), ".", call. = FALSE)
}
## `ok` / `computed` are optional; treat them as all-TRUE when absent so that a
## file produced by an older version of 01_compare_methods.R still summarises.
if (is.null(raw$ok)) raw$ok <- TRUE
if (is.null(raw$computed)) raw$computed <- TRUE
if (is.null(raw$seconds)) raw$seconds <- NA_real_
raw$ok <- as.logical(raw$ok) & !is.na(raw$ok)

## Bias is formed row by row: each batch carries its own Monte Carlo truth.
raw$bias <- raw$estimate - raw$truth
raw$covered <- is.finite(raw$estimate) & is.finite(raw$se) &
  abs(raw$bias) <= stats::qnorm(0.975) * raw$se

## Order the methods the way the comparison is read, with anything unexpected
## appended rather than dropped.
method_order <- c("LSIV","LSIV_one_sided","LSIV_LATE")
raw$method <- factor(raw$method,
                     levels = c(intersect(method_order, unique(raw$method)),
                                setdiff(unique(raw$method), method_order)))

## ---------------------------------------------------------------------------
## Summary table
## ---------------------------------------------------------------------------

lsiv_summarise_merged <- function(d, scale = 10) {
  keys <- unique(d[, c("type", "n", "method")])
  keys <- keys[order(keys$type, keys$n, keys$method), , drop = FALSE]
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    g <- d[d$type == keys$type[i] & d$n == keys$n[i] &
             d$method == keys$method[i], , drop = FALSE]
    fin <- is.finite(g$estimate)
    ok <- fin & g$ok
    stat <- function(idx) {
      b <- g$bias[idx]; s <- g$se[idx]; e <- g$estimate[idx]
      c(reps = length(b),
        bias = if (length(b)) scale * mean(b) else NA_real_,
        med_bias = if (length(b)) scale * stats::median(b) else NA_real_,
        ese = if (length(e) > 1L) scale * stats::sd(e) else NA_real_,
        mad_ese = if (length(e) > 1L) scale * stats::mad(e) else NA_real_,
        ase = if (any(is.finite(s))) scale * mean(s[is.finite(s)]) else NA_real_,
        med_ase = if (any(is.finite(s)))
          scale * stats::median(s[is.finite(s)]) else NA_real_,
        rmse = if (length(b)) scale * sqrt(mean(b^2)) else NA_real_,
        coverage = if (length(b)) mean(g$covered[idx]) else NA_real_)
    }
    all_s <- stat(which(fin))
    ok_s <- stat(which(ok))
    data.frame(
      type = keys$type[i], n = keys$n[i], method = as.character(keys$method[i]),
      reps = nrow(g), n_computed = sum(fin),
      n_ok = sum(ok), n_fail = sum(!g$ok),
      bias = all_s[["bias"]], med_bias = all_s[["med_bias"]],
      ese = all_s[["ese"]], mad_ese = all_s[["mad_ese"]],
      ase = all_s[["ase"]], med_ase = all_s[["med_ase"]],
      rmse = all_s[["rmse"]], coverage = all_s[["coverage"]],
      bias_ok = ok_s[["bias"]], ese_ok = ok_s[["ese"]],
      ase_ok = ok_s[["ase"]], rmse_ok = ok_s[["rmse"]],
      coverage_ok = ok_s[["coverage"]],
      median_seconds = if (any(is.finite(g$seconds)))
        stats::median(g$seconds[is.finite(g$seconds)]) else NA_real_,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

SCALE <- 10
summary_tab <- lsiv_summarise_merged(raw, scale = SCALE)

cat(sprintf("\nSummary of %s  (%d rows, %d batches)\n",
            basename(in_file), nrow(raw), length(unique(raw$rep)) *
              length(unique(raw$n)) * length(unique(raw$type))))
cat(sprintf("Bias / ESE / ASE / RMSE are multiplied by %d, as in ",
            SCALE))
cat("Simulation_Version1/3.Summarize.R\n")
cat("truth varies by batch; bias is computed row by row as estimate - truth\n\n")

cat("--- all computed replications ---\n")
print(summary_tab[, c("type", "n", "method", "n_computed", "n_fail",
                      "bias", "med_bias", "ese", "mad_ese", "ase", "med_ase",
                      "rmse", "coverage", "median_seconds")],
      digits = 4, row.names = FALSE)

if (any(summary_tab$n_fail > 0L)) {
  cat("\n--- restricted to replications passing the A3/A5 diagnostics (ok) ---\n")
  print(summary_tab[, c("type", "n", "method", "n_ok",
                        "bias_ok", "ese_ok", "ase_ok", "rmse_ok",
                        "coverage_ok")],
        digits = 4, row.names = FALSE)
  cat("\nMethods with n_fail > 0 have replications whose fixed point degenerated;\n")
  cat("their all-replication mean is driven by those rows, so read med_bias and\n")
  cat("the ok-only block alongside it.\n")
}

summary_file <- file.path(out_dir, "compare_methods_1snp_summary.csv")
utils::write.csv(summary_tab, summary_file, row.names = FALSE)

## ---------------------------------------------------------------------------
## Figure: bias by method, ONE PANEL PER OUTCOME TYPE.  Each panel carries three
## groups of boxes, one per sample size, and inside a group the three compared
## methods appear side by side in `PLOT_METHODS` order, so a group is read as
## "all methods at this N" and the three groups show what happens as N grows.
## The statistics underneath use the same columns as the boxes, with a colour
## swatch identifying each estimator within every sample-size group.
## ---------------------------------------------------------------------------

## Methods to plot, in the order they appear inside a group, and their labels.
PLOT_METHODS <- c("LSIV","LSIV_one_sided","LSIV_LATE")
METHOD_LABEL <- c(LSIV = expression(hat(tau)), 
                  LSIV_one_sided = expression(hat(tau)["ATT,1snp"]),
                  LSIV_LATE = expression(hat(tau)["LATE"]))
## Set to FALSE for bare boxplots without the tabulated statistics underneath.
SHOW_STATS <- TRUE

methods <- intersect(PLOT_METHODS, as.character(unique(raw$method)))
TYPE_ORDER <- c("continuous", "binary")
types <- c(intersect(TYPE_ORDER, unique(raw$type)),
           setdiff(unique(raw$type), TYPE_ORDER))
ns <- sort(unique(raw$n))

## Grey, white, and light blue distinguish the estimators by lightness and hue.
## Dark borders keep the white boxes visible against the background.
COL <- c(LSIV = "grey45", LSIV_one_sided = "white",
         LSIV_LATE = "#8CC6EC")[methods]

## Geometry: `GAP` is the blank space, in box widths, left between two sample
## sizes.  Boxes inside a group sit at consecutive integers.
GAP <- 1
at_of <- function(ni, mi) (ni - 1L) * (length(methods) + GAP) + mi
AT <- as.vector(t(outer(seq_along(ns), seq_along(methods), at_of)))
GROUP_MID <- vapply(seq_along(ns),
                    function(ni) mean(at_of(ni, seq_along(methods))), numeric(1))
XL <- c(min(AT) - 0.9, max(AT) + 0.9)

draw_panel <- function(tp) {
  d <- raw[raw$type == tp & raw$method %in% methods & is.finite(raw$bias), ,
           drop = FALSE]
  bl <- vector("list", length(ns) * length(methods))
  cols <- character(length(bl))
  k <- 0L
  for (ni in seq_along(ns)) for (mi in seq_along(methods)) {
    k <- k + 1L
    bl[[k]] <- d$bias[d$n == ns[ni] & d$method == methods[mi]]
    cols[k] <- COL[methods[mi]]
  }

  ## One y-window for the whole panel, so the three sample sizes are directly
  ## comparable, with head-room for the method key.
  rng <- range(unlist(bl), na.rm = TRUE)
  sp <- max(diff(rng), 1e-6)
  if(tp=="continuous"){
    YL <- c(-1,1)*0.65
  } else {
    YL <- c(-1,1)*0.2
  }

  graphics::par(mar = c(0.5, 4.2, 4, 0.25))
  graphics::boxplot(lapply(bl, function(v) pmin(pmax(v, YL[1L]), YL[2L])),
                    at = AT, xlim = XL, ylim = YL, axes = FALSE, 
                    outline = TRUE, range = 1.5, col = cols,
                    boxwex = 0.6, medlwd = 1, pch=19, cex=0.25)
  graphics::abline(h = 0, col = 2, lty = 2)
  graphics::axis(2, las = 1, cex.axis = 0.9, at=round(seq(YL[1],YL[2],length=5),3))
  graphics::title(ylab = "Bias", line = 2.8)
  # graphics::box(bty = "l")

  ## Sample-size headers over each group, and the outcome type as the title.
  graphics::mtext(sprintf("N=%s", format(ns, big.mark = ",", trim = TRUE)),
                  side = 3, at = GROUP_MID, line = 0.25, cex = 0.9)
  graphics::mtext(sprintf("%s Y", tools::toTitleCase(tp)), side = 3,
                  line = 2.6, cex = 1.25, font = 1)

  ## Method key, drawn once inside the bottom-right corner of the binary panel
  ## as a 2 x 3 grid.  `legend()` fills column by column, so the methods are
  ## reordered to read row by row (first row = first three methods).
  ## `x.intersp` is the gap between a colour box and its label.
  if (tp == "binary") {
    key_cex <- 1.25
    key_ncol <- 3L
    key_ord <- as.vector(matrix(seq_along(methods), ncol = key_ncol,
                                byrow = TRUE))
    graphics::legend("bottomright", inset = c(0.01, -0.02),
                     legend = METHOD_LABEL[methods][key_ord],
                     fill = COL[methods][key_ord], border = "grey30",
                     ncol = key_ncol, x.intersp = 0.4,
                     bg = "white", box.col = "white", cex = key_cex)
  }
}

draw_stats <- function(tp) {
  graphics::par(mar = c(0.2, 4.2, 0.6, 0.25))
  graphics::plot.new()
  graphics::plot.window(xlim = XL, ylim = c(0, 1))

  labs <- c(sprintf("Bias (x%d)", SCALE), sprintf("ESE (x%d)", SCALE),
            sprintf("ASE (x%d)", SCALE), "Coverage")
  ## Row labels are right-aligned into the left margin so they never run into
  ## the N = 1,000 column.
  yp <- seq(0.72, 0.08, length.out = length(labs))
  graphics::text(XL[1L] - 0.15, yp, labs, adj = c(1, 0.5),
                 xpd = NA, cex = 0.85)
  graphics::points(AT, rep(0.93, length(AT)), pch = 22, cex = 2.2, lwd = 1,
                   bg = rep(COL[methods], times = length(ns)), col = "grey20")
  st <- summary_tab[summary_tab$type == tp, , drop = FALSE]
  for (ni in seq_along(ns)) for (mi in seq_along(methods)) {
    r <- st[st$n == ns[ni] & st$method == methods[mi], , drop = FALSE]
    if (!nrow(r)) next
    vals <- c(sprintf("%.3f", r$bias), sprintf("%.3f", r$ese),
              sprintf("%.3f", r$ase), sprintf("%.3f", r$coverage))
    graphics::text(at_of(ni, mi), yp, vals, cex = 0.9, col = "black")
  }
}

draw_figure <- function() {
  if (!SHOW_STATS) {
    graphics::layout(matrix(seq_along(types), nrow = 1L))
    for (tp in types) draw_panel(tp)
    return(invisible(NULL))
  }

  panel_ids <- 2L * seq_along(types) - 1L
  graphics::layout(rbind(panel_ids, rep(2L * length(types) + 1L, length(types)),
                          panel_ids + 1L), heights = c(3, 0.5, 1.2))
  for (tp in types) {
    draw_panel(tp)
    draw_stats(tp)
  }

  ## Heading spanning all outcome panels; the swatches are above the values.
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  graphics::text(0.5, 0.5, expression("Summary Statistics of Estimates"),
                 cex = 1.8, xpd = NA)
}

png_file <- file.path(out_dir, "Main_Sim_1snp.png")
grDevices::png(png_file, width = 11, height = 4.5, units = "in", res = 500)
draw_figure()
dev.off()
