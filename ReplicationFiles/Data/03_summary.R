LEVEL <- 0.95

raw <- read.csv("Merge_JC_Result.csv")
zcrit <- stats::qnorm(1 - (1 - LEVEL) / 2)

combine <- function(g) {
  keep <- is.finite(g$estimate)
  est <- g$estimate[keep]
  se <- g$se[keep]
  med <- if (length(est)) stats::median(est) else NA_real_
  if (grepl("^QTT_", g$estimand[1L])) {
    se_med <- NA_real_
    lower <- stats::median(g$ci95_lower[keep])
    upper <- stats::median(g$ci95_upper[keep])
  } else {
    se_med <- if (length(est) && all(is.finite(se))) {
      sqrt(stats::median(se^2 + (est - med)^2))
    } else NA_real_
    lower <- med - zcrit * se_med
    upper <- med + zcrit * se_med
  }
  data.frame(method = g$method[1L], estimand = g$estimand[1L],
             num_crossfit = length(est), 
             estimate = med, se = se_med,
             ci95_lower = lower, ci95_upper = upper,
             stringsAsFactors = FALSE)
}

groups <- split(raw, list(raw$method, raw$estimand), drop = TRUE)
summary <- do.call(rbind, lapply(groups, combine))
method_order <- c("SIV", "SIV_Liu", "SIV_Sun", "SIV_LATE", "2SLS",
                  "2SLS_no_X", "Ign", "linear model")
summary <- summary[order(match(summary$method, method_order),
                         summary$estimand), ]
rownames(summary) <- NULL
write.csv(summary, "JC_Summary.csv", row.names = FALSE)

tests <- read.csv("Merge_JC_Falsification.csv")
p_values <- tests$falsification_p_value[is.finite(tests$falsification_p_value)]
median_p <- if (length(p_values)) stats::median(p_values) else NA_real_
falsification <- data.frame(
  num_crossfit = length(p_values),
  median_p_value = median_p,
  reject = median_p <= 1 - LEVEL,
  median_statistic = stats::median(tests$falsification_statistic, na.rm = TRUE),
  median_critical_value = stats::median(tests$falsification_critical_value,
                                        na.rm = TRUE)
)
write.csv(falsification, "JC_Falsification_Summary.csv", row.names = FALSE)

print(transform(summary,
                estimate = round(estimate, 2), se = round(se, 2),
                ci95_lower = round(ci95_lower, 2),
                ci95_upper = round(ci95_upper, 2)))
cat("\nFalsification test\n")
print(falsification)
