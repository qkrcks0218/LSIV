DGP_VALUES <- c("SIV", "mild", "severe")
DGP_LABELS <- c(SIV = "SIV holds", mild = "Mild violation",
                severe = "Severe violation")

raw <- read.csv("Merge_Fal_F.csv")

## Accept files produced before the simulation column was renamed from
## `scenario` to `dgp`.
if (!"dgp" %in% names(raw) && "scenario" %in% names(raw)) {
  raw$dgp <- raw$scenario
}
if (!"dgp" %in% names(raw)) {
  stop("Input must contain a `dgp` column.", call. = FALSE)
}
raw$dgp <- ifelse(tolower(raw$dgp) == "siv", "SIV", tolower(raw$dgp))

N_VALUES <- sort(unique(raw$n))
DGP_VALUES <- DGP_VALUES[DGP_VALUES %in% raw$dgp]

## Older per-batch files predate the estimated-nuisance columns, so a missing
## column summarises as NA rather than aborting the merge.
col <- function(g, name) if (is.null(g[[name]])) rep(NA_real_, nrow(g)) else g[[name]]
finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
finite_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

summarize_group <- function(g) {
  ok <- !is.na(g$reject)
  k <- sum(g$reject[ok]); m <- sum(ok)
  p <- if (m) k / m else NA_real_
  mcse <- if (m) sqrt(p * (1 - p) / m) else NA_real_
  p_value <- col(g, "p_value")
  statistic <- col(g, "statistic")
  critical_value <- col(g, "critical_value")
  marginal_first_stage <- col(g, "marginal_first_stage")
  reversal_probability <- col(g, "reversal_probability")
  converged <- col(g, "converged")
  ## A Wilson interval remains informative when there are zero rejections.
  z <- stats::qnorm(0.975)
  center <- if (m) (p + z^2 / (2 * m)) / (1 + z^2 / m) else NA_real_
  half <- if (m) z * sqrt(p * (1 - p) / m + z^2 / (4 * m^2)) /
    (1 + z^2 / m) else NA_real_
  data.frame(
    dgp = g$dgp[1L],
    dgp_label = unname(DGP_LABELS[g$dgp[1L]]),
    n = g$n[1L], reps = nrow(g), n_ok = m, n_failed = sum(!ok),
    rejections = k, rejection_probability = p, mcse = mcse,
    wilson_lower = max(0, center - half), wilson_upper = min(1, center + half),
    mean_p_value = finite_mean(p_value[ok]),
    mean_statistic = finite_mean(statistic[ok]),
    mean_critical_value = finite_mean(critical_value[ok]),
    mean_marginal_first_stage = finite_mean(marginal_first_stage[ok]),
    mean_reversal_probability = finite_mean(reversal_probability[ok]),
    ## With Super Learner nuisances the fit can fail to converge or hit the
    ## Assumption (A3)/(A5) diagnostics, so a rejection probability is only
    ## interpretable alongside these.
    n_converged = if (all(is.na(converged[ok]))) NA_integer_ else
      sum(converged[ok] == 1, na.rm = TRUE),
    mean_att = finite_mean(col(g, "att")[ok]),
    mean_g_floor_frac = finite_mean(col(g, "g_floor_frac")[ok]),
    min_relevance = finite_min(col(g, "min_relevance")[ok]),
    mean_seconds = finite_mean(col(g, "seconds")), stringsAsFactors = FALSE)
}

groups <- split(raw, interaction(raw$dgp, raw$n, drop = TRUE))
summary <- do.call(rbind, lapply(groups, summarize_group))
rownames(summary) <- NULL
summary$dgp <- factor(summary$dgp, levels = DGP_VALUES)
summary <- summary[order(summary$dgp, summary$n), ]
summary$dgp <- as.character(summary$dgp)

colors <- c(SIV = "#1B6CA8", mild = "#E69F00", severe = "#C0392B")
shapes <- c(SIV = 16, mild = 17, severe = 15)

png("Falsification.png",height=5,width=10,unit="in",res=500)
par(mar=c(3.5,3.5,2,0.5))
draw_power_curve <- function() {
  graphics::plot(range(N_VALUES), c(-0.05, 1.05), type = "n",
                 xlab = "", ylab = "",
                 yaxs = "i",
                 main = "")
  
  title(main="Falsification Test")
  title(xlab="Sample size N", ylab="Rejection probability",line=2.5)
  graphics::abline(h = 0.05, lty = 2, col = "gray40")
  for (dgp in DGP_VALUES) {
    z <- summary[summary$dgp == dgp, ]
    # graphics::segments(z$n, z$wilson_lower, z$n, z$wilson_upper,
    #                    col = colors[dgp])
    graphics::lines(z$n, z$rejection_probability, type = "b", lwd = 2,
                    pch = shapes[dgp], col = colors[dgp])
  }
  graphics::legend("topleft", legend = unname(DGP_LABELS[DGP_VALUES]),
                   col = colors[DGP_VALUES], pch = shapes[DGP_VALUES],
                   lty = 1, lwd = 2, bty = "n")
  text(4000,0.08,"Nominal 0.05",pos=2)
}
draw_power_curve()
dev.off()
