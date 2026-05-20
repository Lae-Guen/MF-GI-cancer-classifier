# =============================================================================
# Multivariate Analysis: MF Index ~ Group + Alcohol + Exercise + BMI
# Purpose : Confirm that MF index elevation in GC/CRC remains significant
#           after adjustment for lifestyle confounders (alcohol, exercise, BMI).
#           Addresses reviewer concern regarding potential confounding.
#
# Model    : Linear regression with log1p-transformed MF index as outcome.
#            log1p transformation applied because the MF index is right-skewed
#            with structural zeros (Box-Cox optimal lambda ≈ 0).
# Approach : Complete-case analysis; HC as reference group.
#            HC3 robust standard errors reported if heteroscedasticity detected
#            (Breusch-Pagan test p < 0.05).
#
# Variable coding (from lifestyle questionnaire):
#   ALCO_A_MOD : 1 = Non-drinker (< once/month)  2 = Drinker
#   EXER_A_MOD : 1 = No regular exercise          2 = Regular exercise
#   SMOK_A_MOD : 1 = Never  2 = Current  3 = Former smoker
#   Group      : HC (reference), MD, GC, CRC
#
# INPUT  : data/mf_index_meta_15000.csv — person-level MF index with metadata
# OUTPUT : MF_multivariate_coef_table.csv — coefficient table (Supplementary)
#          MF_model_diagnostics.pdf       — residual diagnostic plots
#          MF_index_adjusted_EMM.pdf      — estimated marginal means plot
# =============================================================================

library(tidyverse)
library(emmeans)   # estimated marginal means
library(car)       # Anova() Type III SS, vif()
library(lmtest)    # bptest() heteroscedasticity test
library(sandwich)  # vcovHC() robust standard errors
library(broom)     # tidy() coefficient table
library(ggpubr)

# -----------------------------------------------------------------------------
# 1. Load data and recode covariates
# -----------------------------------------------------------------------------
df <- read.csv("data/mf_index_meta_15000.csv", stringsAsFactors = FALSE)

df <- df %>%
  mutate(
    # Recode to binary (0/1) for cleaner model interpretation
    Alcohol  = ifelse(ALCO_A_MOD == 2, 1L, 0L),  # 1 = drinker
    Exercise = ifelse(EXER_A_MOD == 2, 1L, 0L),  # 1 = regular exercise

    # Set HC as the reference group
    Group = factor(Group, levels = c("HC", "MD", "GC", "CRC")),

    # log1p transformation to address right skew with structural zeros
    logMF = log1p(MTG.index)
  )

# Complete-case analysis: exclude rows with missing covariates
df_cc <- df %>% filter(!is.na(BMI), !is.na(Alcohol), !is.na(Exercise))

cat("=== Complete-case N per group ===\n")
print(table(df_cc$Group))

# -----------------------------------------------------------------------------
# 2. Unadjusted model (baseline)
# -----------------------------------------------------------------------------
m0 <- lm(logMF ~ Group, data = df_cc)
cat("\n=== Model 0: Unadjusted ===\n")
print(summary(m0))

# -----------------------------------------------------------------------------
# 3. Primary adjusted model: Group + Alcohol + Exercise + BMI
# -----------------------------------------------------------------------------
m1 <- lm(logMF ~ Group + Alcohol + Exercise + BMI, data = df_cc)

cat("\n=== Model 1: Adjusted (Group + Alcohol + Exercise + BMI) ===\n")
print(summary(m1))

# Type III ANOVA: appropriate for unbalanced designs with covariates
cat("\n--- Type III ANOVA ---\n")
print(Anova(m1, type = "III"))

# -----------------------------------------------------------------------------
# 4. Model diagnostics
# -----------------------------------------------------------------------------
cat("\n=== Variance Inflation Factors (multicollinearity; all should be < 5) ===\n")
print(vif(m1))

bp <- bptest(m1)
cat(sprintf("\nBreusch-Pagan test: BP = %.3f, p = %.4f\n", bp$statistic, bp$p.value))

# Report HC3 robust SEs if heteroscedasticity is detected
if (bp$p.value < 0.05) {
  cat("Heteroscedasticity detected; reporting HC3 robust standard errors.\n")
  cat("\n--- Coefficient table with HC3 robust SEs ---\n")
  print(coeftest(m1, vcov = vcovHC(m1, type = "HC3")))
}

# Residual diagnostic plots
pdf("MF_model_diagnostics.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
plot(m1, main = "Model 1: Residual Diagnostics")
dev.off()
cat("Diagnostic plots saved: MF_model_diagnostics.pdf\n")

# -----------------------------------------------------------------------------
# 5. Estimated marginal means (EMM): adjusted group comparisons
#    Covariates (Alcohol, Exercise, BMI) held at their observed means
# -----------------------------------------------------------------------------
emm <- emmeans(m1, specs = "Group")

cat("\n=== Estimated Marginal Means (log1p scale) ===\n")
print(emm)

# Pairwise contrasts vs HC (Bonferroni correction)
cat("\n=== Pairwise contrasts vs HC (Bonferroni-corrected) ===\n")
contr_vs_HC <- contrast(emm, method = "trt.vs.ctrl", ref = "HC",
                        adjust = "bonferroni")
print(summary(contr_vs_HC, infer = TRUE))

# All pairwise contrasts
cat("\n=== All pairwise contrasts (Bonferroni-corrected) ===\n")
print(summary(pairs(emm, adjust = "bonferroni"), infer = TRUE))

# Back-transform EMMs to original MF index scale
emm_df <- as.data.frame(summary(emm)) %>%
  mutate(
    MF_adj       = expm1(emmean),
    MF_adj_lower = expm1(lower.CL),
    MF_adj_upper = expm1(upper.CL)
  )

cat("\n=== Back-transformed EMMs (original MF index scale) ===\n")
print(emm_df[, c("Group", "MF_adj", "MF_adj_lower", "MF_adj_upper")])

# -----------------------------------------------------------------------------
# 6. Sensitivity analysis: HC vs GC vs CRC only (excluding MD group)
# -----------------------------------------------------------------------------
df_cancer <- df_cc %>%
  filter(Group %in% c("HC", "GC", "CRC")) %>%
  mutate(Group = droplevels(Group))

m2 <- lm(logMF ~ Group + Alcohol + Exercise + BMI, data = df_cancer)

cat("\n=== Sensitivity model (HC + GC + CRC only, MD excluded) ===\n")
print(summary(m2))
print(Anova(m2, type = "III"))

emm2 <- emmeans(m2, specs = "Group")
cat("\nContrasts vs HC (sensitivity model, Bonferroni-corrected):\n")
print(summary(contrast(emm2, method = "trt.vs.ctrl", ref = "HC",
                       adjust = "bonferroni"), infer = TRUE))

# -----------------------------------------------------------------------------
# 7. Coefficient table for supplementary material
# -----------------------------------------------------------------------------
coef_table <- tidy(m1, conf.int = TRUE) %>%
  mutate(
    term = case_match(term,
      "(Intercept)" ~ "Intercept",
      "GroupMD"     ~ "Group: MD vs HC",
      "GroupGC"     ~ "Group: GC vs HC",
      "GroupCRC"    ~ "Group: CRC vs HC",
      "Alcohol"     ~ "Alcohol (drinker vs non-drinker)",
      "Exercise"    ~ "Exercise (regular vs none)",
      "BMI"         ~ "BMI",
      .default = term
    ),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  select(term, estimate, std.error, conf.low, conf.high, statistic, p.value, significance) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

cat("\n=== Coefficient Table (Supplementary Table S2) ===\n")
print(coef_table)

write.csv(coef_table, "MF_multivariate_coef_table.csv", row.names = FALSE)
cat("Coefficient table saved: MF_multivariate_coef_table.csv\n")

# -----------------------------------------------------------------------------
# 8. Visualization: adjusted group means ± 95% CI (back-transformed)
# -----------------------------------------------------------------------------
emm_df <- emm_df %>%
  mutate(Group = factor(Group, levels = c("HC", "MD", "GC", "CRC")))

# Build significance bracket data from contrasts vs HC
sig_data <- summary(contr_vs_HC, infer = TRUE) %>%
  as.data.frame() %>%
  filter(p.value < 0.05) %>%
  mutate(
    group1   = "HC",
    group2   = sub(" - HC", "", contrast),
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  )

p <- ggplot(emm_df, aes(x = Group, y = MF_adj, color = Group)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = MF_adj_lower, ymax = MF_adj_upper),
                width = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("HC"  = "#75ca9d",
                                 "MD"  = "#ffeb33",
                                 "GC"  = "#ca75dd",
                                 "CRC" = "#eaa4a4")) +
  labs(
    title    = "MF Index: Adjusted Group Comparison",
    subtitle = paste("Estimated marginal means ± 95% CI",
                     "(adjusted for Alcohol, Exercise, BMI;",
                     "back-transformed from log1p scale)", sep = "\n"),
    x = NULL,
    y = "Adjusted MF Index"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey45", size = 10)
  )

# Add significance brackets for comparisons with p < 0.05
if (nrow(sig_data) > 0) {
  y_max   <- max(emm_df$MF_adj_upper, na.rm = TRUE)
  y_steps <- y_max * seq(1.15, by = 0.15, length.out = nrow(sig_data))
  sig_data <- sig_data %>% mutate(y.position = y_steps)

  p <- p + stat_pvalue_manual(
    sig_data, label = "p.signif",
    xmin = "group1", xmax = "group2",
    y.position = "y.position", tip.length = 0.01
  )
}

ggsave("MF_index_adjusted_EMM.pdf", p, width = 5.5, height = 5)
cat("Figure saved: MF_index_adjusted_EMM.pdf\n")

# -----------------------------------------------------------------------------
# 9. Summary for methods/results reporting
# -----------------------------------------------------------------------------
s   <- summary(m1)
cfs <- coef(s)
f_p <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)

cat("\n")
cat("================================================================\n")
cat("  SUMMARY FOR RESULTS REPORTING\n")
cat("================================================================\n")
cat(sprintf("  Complete-case N = %d (HC=%d, MD=%d, GC=%d, CRC=%d)\n",
    nrow(df_cc),
    sum(df_cc$Group == "HC"), sum(df_cc$Group == "MD"),
    sum(df_cc$Group == "GC"), sum(df_cc$Group == "CRC")))
cat(sprintf("  Model: R² = %.3f, Adj.R² = %.3f, F p-value = %.2e\n\n",
    s$r.squared, s$adj.r.squared, f_p))

cat("  Disease group effects vs HC (adjusted):\n")
for (g in c("GroupMD", "GroupGC", "GroupCRC")) {
  b  <- cfs[g, "Estimate"]
  se <- cfs[g, "Std. Error"]
  pv <- cfs[g, "Pr(>|t|)"]
  pl <- ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.4f", pv))
  cat(sprintf("    %-8s  beta = %+.3f (SE %.3f), %s\n",
              sub("Group", "", g), b, se, pl))
}

cat("\n  Lifestyle covariate effects:\n")
for (cv in c("Alcohol", "Exercise", "BMI")) {
  b  <- cfs[cv, "Estimate"]
  pv <- cfs[cv, "Pr(>|t|)"]
  pl <- ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.4f", pv))
  cat(sprintf("    %-10s  beta = %+.3f, %s\n", cv, b, pl))
}
cat("================================================================\n")
