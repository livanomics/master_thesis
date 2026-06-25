# == 01. Information ==========================================================
#' Functions for the feature analysis of gene lists: normality and Wilcoxon
#' tests, Cliff's delta effect size and group comparisons


# == 02. Effect size ==========================================================
# Cliff's delta for two independent samples
cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA)
  r  <- rank(c(x, y))
  rx <- sum(r[seq_along(x)])
  U  <- rx - length(x) * (length(x) + 1) / 2
  (2 * U) / (length(x) * length(y)) - 1
}

# Qualitative magnitude of Cliff's delta
delta_label <- function(d) {
  if (is.na(d)) return(NA)
  ad <- abs(d)
  if (ad < 0.147) "negligible"
  else if (ad < 0.33)  "small"
  else if (ad < 0.474) "medium"
  else "large"
}


# == 03. Statistical tests ====================================================
# Shapiro-Wilk p-value (needs >=3 finite values, capped at 5000 for the test)
shapiro_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA)
  if (length(x) > 5000) x <- sample(x, 5000)
  shapiro.test(x)$p.value
}

# Wilcoxon test between two named groups for one variable
test_groups <- function(data, group_a, group_b, v) {
  x <- data[data$group == group_a, v]; x <- x[is.finite(x)]
  y <- data[data$group == group_b, v]; y <- y[is.finite(y)]
  w <- wilcox.test(x, y)
  d <- cliffs_delta(x, y)
  data.frame(
    group_a  = group_a,
    group_b  = group_b,
    variable = v,
    median_a = round(median(x), 2),
    median_b = round(median(y), 2),
    p_value  = signif(w$p.value, 3),
    cliff_delta = round(d, 3),
    effect      = delta_label(d))
}

# Run Shapiro-Wilk per group and variable
run_normality <- function(data, vars) {
  normality <- data.frame()
  for (g in unique(data$group)) {
    for (v in names(vars)) {
      x <- data[data$group == g, v]
      p <- shapiro_safe(x)
      normality <- rbind(normality, data.frame(
        group     = g,
        variable  = v,
        n         = sum(is.finite(x)),
        shapiro_p = signif(p, 3),
        normal    = ifelse(is.na(p), NA, p > 0.05)))
    }
  }
  normality
}

# Run Wilcoxon + Cliff's delta (UP vs DOWN per condition, and each list vs siC)
# with BH correction across this family of tests
run_wilcoxon <- function(data, vars) {
  results <- data.frame()
  diff_groups <- setdiff(unique(data$group), "siC")
  for (v in names(vars)) {
    # UP vs DOWN within each condition
    for (smp in c("siUAP56", "siBRG1")) {
      a <- paste0(smp, "_UP"); b <- paste0(smp, "_DOWN")
      if (a %in% data$group && b %in% data$group)
        results <- rbind(results, test_groups(data, a, b, v))
    }
    # each list vs siC
    if ("siC" %in% data$group)
      for (g in diff_groups)
        results <- rbind(results, test_groups(data, g, "siC", v))
  }
  # FDR correction within this family of tests
  results$p_adj       <- signif(p.adjust(results$p_value, method = "BH"), 3)
  results$significant <- results$p_adj < 0.05
  results
}

