# =============================================================================
# AUC Comparison Statistics: Method × Host Effect Size and CI
# Purpose : Statistically compare AUC distributions across four classification
#           feature types (Total, MF, RF, Core) for GC and CRC separately.
#           Uses Kruskal-Wallis + Dunn post-hoc for independent comparisons,
#           and seed-matched paired Wilcoxon for within-host method contrasts.
#           All effect sizes are rank-biserial r (directional) with 95% bootstrap CI.
#
# Design  : host (GC, CRC) × method (Total, MF, RF, Core), n = 50 seeds/cell
# Effect size interpretation:
#   |r| < 0.10       negligible
#   |r| 0.10 – 0.29  small
#   |r| 0.30 – 0.49  moderate
#   |r| >= 0.50      large
#
# INPUT  : mouth_mdeep_rf.csv — long-format AUC table (cols: host, method, auc, seed)
# OUTPUT : result_kruskal_effectsize.csv
#          result_delta_auc.csv        (seed-matched Hodges-Lehmann + r)
#          result_dunn_pairwise.csv    (Dunn-Bonferroni + Hodges-Lehmann + r)
#          result_host_comparison.csv  (GC vs CRC per method)
#          auc_method_comparison.pdf
# =============================================================================

# ============================================================
#  AUC Comparison Statistics: Effect Size and CI by Host × Method
#  Design: host(GC, CRC) × method(Total, MF, RF, Core), n=50/cell
#  Rank-biserial r computed directly from W statistic (directional)
#     independent:  r = 1 - 2W / (n1*n2)
#     paired:       r = 1 - 2W / (n*(n+1)/2)
#     → r > 0 when hodges_hl > 0 (group1 > group2)
# ============================================================

library(tidyverse)
library(rstatix)    # kruskal_effsize, dunn_test
library(ggplot2)
library(ggpubr)

# Helper: rank-biserial r (independent samples, directional)
r_rb_indep <- function(x1, x2, nboot = 1000, conf.level = 0.95) {
  wt  <- wilcox.test(x1, x2, conf.int = TRUE, conf.level = conf.level,
                     exact = FALSE)
  n1  <- length(x1); n2 <- length(x2)
  r   <- 1 - (2 * as.numeric(wt$statistic)) / (n1 * n2)

  set.seed(42)
  boot_r <- replicate(nboot, {
    b1 <- sample(x1, replace = TRUE)
    b2 <- sample(x2, replace = TRUE)
    wt_b <- wilcox.test(b1, b2, exact = FALSE)
    1 - (2 * as.numeric(wt_b$statistic)) / (n1 * n2)
  })
  list(
    hodges_hl   = as.numeric(wt$estimate),
    hl_ci_lower = wt$conf.int[1],
    hl_ci_upper = wt$conf.int[2],
    r           = r,
    r_ci_lower  = quantile(boot_r, (1 - conf.level) / 2),
    r_ci_upper  = quantile(boot_r, 1 - (1 - conf.level) / 2)
  )
}

# Helper: rank-biserial r (paired samples, directional)
r_rb_paired <- function(x1, x2, nboot = 1000, conf.level = 0.95) {
  wsr <- wilcox.test(x1, x2, paired = TRUE,
                     conf.int = TRUE, conf.level = conf.level, exact = FALSE)
  n   <- length(x1)
  r   <- 1 - (2 * as.numeric(wsr$statistic)) / (n * (n + 1) / 2)

  set.seed(42)
  delta  <- x1 - x2
  boot_r <- replicate(nboot, {
    d_b  <- sample(delta, replace = TRUE)
    wt_b <- wilcox.test(d_b, mu = 0, exact = FALSE)
    1 - (2 * as.numeric(wt_b$statistic)) / (n * (n + 1) / 2)
  })
  list(
    hodges_hl   = as.numeric(wsr$estimate),
    hl_ci_lower = wsr$conf.int[1],
    hl_ci_upper = wsr$conf.int[2],
    r           = r,
    r_ci_lower  = quantile(boot_r, (1 - conf.level) / 2),
    r_ci_upper  = quantile(boot_r, 1 - (1 - conf.level) / 2),
    p_value     = wsr$p.value
  )
}


# ── 0. Load data ────────────────────────────────────────────
df <- read.csv("mouth_mdeep_rf.csv") %>%
  mutate(
    host   = factor(host,   levels = c("GC", "CRC")),
    method = factor(method, levels = c("Total", "MF", "RF", "Core"))
  )

cat("=== Data summary ===\n")
df %>%
  group_by(host, method) %>%
  summarise(n = n(), mean_auc = mean(auc), sd_auc = sd(auc), .groups = "drop") %>%
  print()


# ============================================================
# 1. Kruskal-Wallis: global test of method differences (per host)
#    + rank epsilon squared effect size + 95% bootstrap CI
# ============================================================

cat("\n\n=== [1] Kruskal-Wallis: method comparison per host ===\n")

kw_results <- df %>%
  group_by(host) %>%
  group_modify(~ {
    kt   <- kruskal.test(auc ~ method, data = .x)
    eps2 <- kruskal_effsize(.x, auc ~ method,
                            ci = TRUE, conf.level = 0.95, nboot = 1000)
    tibble(
      statistic = kt$statistic,
      df        = kt$parameter,
      p.value   = kt$p.value,
      epsilon2  = eps2$effsize,
      magnitude = eps2$magnitude,
      ci_lower  = eps2$conf.low,
      ci_upper  = eps2$conf.high
    )
  }) %>%
  ungroup()

kw_results %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  print()


# ============================================================
# 2. Seed-matched pairwise ΔAUC (paired)
#    Hodges-Lehmann estimator + 95% CI + rank-biserial r (directional)
# ============================================================

cat("\n\n=== [2] Seed-matched delta AUC + Hodges-Lehmann + r (per host) ===\n")

method_pairs <- combn(levels(df$method), 2, simplify = FALSE)

delta_results <- df %>%
  group_by(host) %>%
  group_modify(~ {
    host_data <- .x
    map_dfr(method_pairs, function(pair) {
      m1   <- pair[1]; m2 <- pair[2]
      auc1 <- host_data$auc[host_data$method == m1]
      auc2 <- host_data$auc[host_data$method == m2]

      res   <- r_rb_paired(auc1, auc2)
      p_adj <- p.adjust(res$p_value, method = "bonferroni",
                        n = length(method_pairs))

      tibble(
        comparison   = paste0(m1, " - ", m2),
        group1       = m1,
        group2       = m2,
        n_seeds      = length(auc1),
        hl_estimate  = res$hodges_hl,
        hl_ci_lower  = res$hl_ci_lower,
        hl_ci_upper  = res$hl_ci_upper,
        r            = res$r,
        r_ci_lower   = res$r_ci_lower,
        r_ci_upper   = res$r_ci_upper,
        p_value      = res$p_value,
        p_bonferroni = p_adj
      )
    })
  }) %>%
  ungroup()

delta_results %>%
  select(host, comparison, hl_estimate, hl_ci_lower, hl_ci_upper,
         r, r_ci_lower, r_ci_upper, p_bonferroni) %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  print(n = Inf)


# ============================================================
# 3. Dunn post-hoc: pairwise method comparisons (independent)
#    Hodges-Lehmann + 95% CI + rank-biserial r (directional)
# ============================================================

cat("\n\n=== [3] Dunn post-hoc + Hodges-Lehmann + r (per host) ===\n")

dunn_results <- df %>%
  group_by(host) %>%
  group_modify(~ {
    dunn  <- dunn_test(.x, auc ~ method, p.adjust.method = "bonferroni")
    pairs <- dunn %>% select(group1, group2)

    rb_list <- map_dfr(seq_len(nrow(pairs)), function(i) {
      g1  <- pairs$group1[i]; g2 <- pairs$group2[i]
      x1  <- .x$auc[.x$method == g1]
      x2  <- .x$auc[.x$method == g2]
      res <- r_rb_indep(x1, x2)
      tibble(
        group1      = g1,
        group2      = g2,
        hodges_hl   = res$hodges_hl,
        hl_ci_lower = res$hl_ci_lower,
        hl_ci_upper = res$hl_ci_upper,
        r           = res$r,
        r_ci_lower  = res$r_ci_lower,
        r_ci_upper  = res$r_ci_upper
      )
    })

    dunn %>%
      select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
      left_join(rb_list, by = c("group1", "group2"))
  }) %>%
  ungroup()

dunn_results %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  print(n = Inf)


# ============================================================
# 4. GC vs CRC comparison (per method, independent)
#    Hodges-Lehmann + 95% CI + rank-biserial r (directional)
# ============================================================

cat("\n\n=== [4] GC vs CRC comparison (per method) ===\n")

host_compare <- df %>%
  group_by(method) %>%
  group_modify(~ {
    gc  <- .x$auc[.x$host == "GC"]
    crc <- .x$auc[.x$host == "CRC"]
    res <- r_rb_indep(gc, crc)
    wt  <- wilcox.test(gc, crc, exact = FALSE)

    tibble(
      n_GC        = length(gc),
      n_CRC       = length(crc),
      median_GC   = median(gc),
      median_CRC  = median(crc),
      hodges_hl   = res$hodges_hl,
      hl_ci_lower = res$hl_ci_lower,
      hl_ci_upper = res$hl_ci_upper,
      r           = res$r,
      r_ci_lower  = res$r_ci_lower,
      r_ci_upper  = res$r_ci_upper,
      p.value     = wt$p.value,
      p.adj       = p.adjust(wt$p.value, method = "bonferroni", n = 4)
    )
  }) %>%
  ungroup()

host_compare %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  print()


# ============================================================
# 5. Save results
# ============================================================

write.csv(kw_results,    "result_kruskal_effectsize.csv", row.names = FALSE)
write.csv(delta_results, "result_delta_auc.csv",          row.names = FALSE)
write.csv(dunn_results,  "result_dunn_pairwise.csv",      row.names = FALSE)
write.csv(host_compare,  "result_host_comparison.csv",    row.names = FALSE)

cat("\nResults saved:\n")
cat("  - result_kruskal_effectsize.csv\n")
cat("  - result_delta_auc.csv\n")
cat("  - result_dunn_pairwise.csv\n")
cat("  - result_host_comparison.csv\n")


# ============================================================
# 6. Visualization
# ============================================================

p <- ggplot(df, aes(x = method, y = auc, fill = method)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = 16, outlier.size = 1.5) +
  stat_compare_means(method = "kruskal.test", label = "p.format",
                     label.y = max(df$auc) * 1.05) +
  facet_wrap(~ host, ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  labs(title    = "AUC distribution by method and host",
       subtitle = "Kruskal-Wallis p-value shown; post-hoc: Dunn-Bonferroni",
       x = "Method", y = "AUC") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("auc_method_comparison.pdf", p, width = 10, height = 5)
cat("  - auc_method_comparison.pdf\n")
