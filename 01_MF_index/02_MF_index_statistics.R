# =============================================================================
# MF Index: Group Comparison Statistics and Visualization
# Purpose : Compute person-level MF index from shared ASV counts,
#           perform pre-specified pairwise Wilcoxon rank-sum tests, report
#           rank-biserial correlation (r) with 95% bootstrap CIs, and
#           visualize group differences and lifestyle associations.
#
# Statistical approach:
#   - No multiple testing correction applied to pre-specified hypothesis-driven
#     comparisons (HC vs MD, HC vs GC, HC vs CRC, GC vs CRC), consistent with
#     the paper's Statistics section.
#
# Effect size interpretation (rank-biserial r):
#   |r| < 0.10       negligible
#   |r| 0.10 – 0.29  small
#   |r| 0.30 – 0.49  moderate
#   |r| >= 0.50      large
#
# INPUT  : data/yonsei_r15000.RDS       — rarefied Yonsei phyloseq object
#          data/yonsei.mtg.csv          — shared ASV table from script 01
# OUTPUT : data/mf_index_meta_15000.csv — person-level MF index with metadata
#          mf_index_wilcoxon_effectsize.csv
#          mf_index_descriptive_stats.csv
#          mf_index.png                 — violin plot of group comparisons
# =============================================================================

library(phyloseq)
library(microbiome)
library(tidyverse)
library(rstatix)   # wilcox_test(), wilcox_effsize()
library(ggpubr)
library(Hmisc)     # rcorr() for Spearman correlation

# -----------------------------------------------------------------------------
# 1. Compute person-level MF index
# -----------------------------------------------------------------------------
# MF index (%) = (sum of shared fecal ASV counts / rarefaction depth) × 100
# Samples without any shared ASVs (paired but no MF bacteria detected) are
# assigned an MF index of 0. Unpaired samples are excluded.

d   <- readRDS("data/yonsei_r15000.RDS")
mtg <- read.csv("data/yonsei.mtg.csv")

# Aggregate shared fecal counts per sample and compute MF index
mtg_index <- aggregate(Fecal_count ~ Feces_sample, data = mtg, sum) %>%
  mutate(MTG.index = (Fecal_count / 15000) * 100)

# Extract fecal sample metadata
meta <- microbiome::meta(d) %>%
  filter(Type == "Feces") %>%
  mutate(Feces_sample = rownames(.)) %>%
  select(Feces_sample, Group, BMI, SMOK_A_MOD, ALCO_A_MOD, EXER_A_MOD, MED9)

# Left join: samples with no shared ASVs receive NA → recode to 0
mtg_index <- left_join(meta, mtg_index, by = "Feces_sample") %>%
  mutate(MTG.index = replace_na(MTG.index, 0))

write.csv(mtg_index, "data/mf_index_meta_15000.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Load data and set group factor order
# -----------------------------------------------------------------------------
df <- read.csv("data/mf_index_meta_15000.csv")
df$Group <- factor(df$Group, levels = c("HC", "MD", "GC", "CRC"))

# -----------------------------------------------------------------------------
# 3. Pre-specified pairwise Wilcoxon rank-sum tests (no correction)
# -----------------------------------------------------------------------------
comparisons <- list(
  c("HC", "MD"),
  c("HC", "GC"),
  c("HC", "CRC"),
  c("GC", "CRC")
)

wilcox_results <- rstatix::wilcox_test(
  data            = df,
  formula         = MTG.index ~ Group,
  comparisons     = comparisons,
  p.adjust.method = "none"
)

# -----------------------------------------------------------------------------
# 4. Effect size: rank-biserial correlation r with 95% bootstrap CI
# -----------------------------------------------------------------------------
effsize_results <- rstatix::wilcox_effsize(
  data        = df,
  formula     = MTG.index ~ Group,
  comparisons = comparisons,
  ci          = TRUE,
  conf.level  = 0.95,
  nboot       = 5000
)

# -----------------------------------------------------------------------------
# 5. Merge Wilcoxon results and effect sizes into a single table
# -----------------------------------------------------------------------------
result_table <- wilcox_results %>%
  select(group1, group2, statistic, p) %>%
  left_join(
    effsize_results %>% select(group1, group2, effsize, magnitude, conf.low, conf.high),
    by = c("group1", "group2")
  ) %>%
  rename(
    W         = statistic,
    p_value   = p,
    r         = effsize,
    r_CI_low  = conf.low,
    r_CI_high = conf.high
  ) %>%
  mutate(
    p_value   = format(round(p_value, 4), nsmall = 4),
    r         = round(r, 3),
    r_CI_low  = round(r_CI_low, 3),
    r_CI_high = round(r_CI_high, 3)
  )

print(result_table)

# -----------------------------------------------------------------------------
# 6. Descriptive statistics: mean ± SD and median [IQR] per group
# -----------------------------------------------------------------------------
desc_stats <- df %>%
  group_by(Group) %>%
  summarise(
    n      = n(),
    mean   = round(mean(MTG.index), 3),
    sd     = round(sd(MTG.index), 3),
    median = round(median(MTG.index), 3),
    Q1     = round(quantile(MTG.index, 0.25), 3),
    Q3     = round(quantile(MTG.index, 0.75), 3),
    .groups = "drop"
  )

print(desc_stats)

# -----------------------------------------------------------------------------
# 7. Export results
# -----------------------------------------------------------------------------
write.csv(result_table, "mf_index_wilcoxon_effectsize.csv",  row.names = FALSE)
write.csv(desc_stats,   "mf_index_descriptive_stats.csv",    row.names = FALSE)
cat("Results saved: mf_index_wilcoxon_effectsize.csv\n")

# -----------------------------------------------------------------------------
# 8. Visualization: violin plot of MF index by group
# -----------------------------------------------------------------------------
my_comparisons <- list(c("HC", "MD"), c("HC", "GC"), c("HC", "CRC"))

p <- ggviolin(mtg_index, x = "Group", y = "MTG.index",
              fill = "Group",
              palette = c("#75ca9d", "#ffeb33", "#ca75dd", "#eaa4a4"),
              add = "boxplot",
              add.params = list(fill = "white", alpha = 0.5)) +
  stat_compare_means(comparisons = my_comparisons, size = 5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

ggsave("mf_index.png", p, width = 10, height = 12, dpi = 300, units = "cm")
cat("Figure saved: mf_index.png\n")

# -----------------------------------------------------------------------------
# 9. Lifestyle association: MF index vs alcohol consumption (CRC subgroup)
# -----------------------------------------------------------------------------
filtered_alco <- mtg_index %>% filter(!is.na(ALCO_A_MOD), Group == "CRC")

ggboxplot(filtered_alco, x = "ALCO_A_MOD", y = "MTG.index",
          fill = "#D3D3D3") +
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.7) +
  scale_color_manual(values = c("#eaa4a4")) +
  stat_compare_means(size = 5)

# -----------------------------------------------------------------------------
# 10. Lifestyle association: MF index vs antidiabetic medication (MD subgroup)
# -----------------------------------------------------------------------------
# MED9: antidiabetic medication use (binary: 0 = no, 1 = yes)
filtered_med <- mtg_index %>% filter(!is.na(MED9), Group == "MD")

ggboxplot(filtered_med, x = "MED9", y = "MTG.index",
          fill = "#D3D3D3") +
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.7) +
  scale_color_manual(values = c("#ffeb33")) +
  stat_compare_means(size = 5)

# -----------------------------------------------------------------------------
# 11. Spearman correlation: MF-shared bacteria abundance vs clinical indicators
# -----------------------------------------------------------------------------
data_corr <- read.csv("data/yonsei.mtg.csv") %>%
  mutate(
    Rel.abund.mouth = Oral_count  / 15000,
    Rel.abund.feces = Fecal_count / 15000,
    orallog  = log10(Rel.abund.mouth),
    fecallog = log10(Rel.abund.feces)
  ) %>%
  select(Feces_sample, Species, orallog, fecallog)

meta_full <- microbiome::meta(d) %>%
  filter(Type == "Feces") %>%
  mutate(Feces_sample = rownames(.))

data_corr <- left_join(data_corr, meta_full, by = "Feces_sample")

# Clinical indicators to test
clinical_vars <- c("TG", "HDL", "LDL", "CHO", "BMI", "HbA1c", "FBS",
                   "SBP", "DBP", "WBC", "PLT", "globulin", "GGT", "LD",
                   "hsCRP", "FreeT4", "TSH", "BUN.Creatinineratio",
                   "Albumin.Globulin", "AFP", "CA19.9", "VitD")

# Keep only species with >= 50 non-missing fecal log-abundance observations
valid_species <- data_corr %>%
  group_by(Species) %>%
  filter(sum(!is.na(fecallog)) >= 50) %>%
  pull(Species) %>%
  unique()

# Run Spearman correlation for each valid species
significant_list <- list()

for (sp in valid_species) {
  sp_data <- data_corr %>% filter(Species == sp)
  if (nrow(sp_data) < 10) next

  cor_mat <- Hmisc::rcorr(
    as.matrix(sp_data[, c("fecallog", clinical_vars)]),
    type = "spearman"
  )

  rho_vals <- cor_mat$r[1, -1]
  p_vals   <- cor_mat$P[1, -1]
  sig_idx  <- which(p_vals < 0.05)

  if (length(sig_idx) > 0) {
    significant_list[[length(significant_list) + 1]] <- data.frame(
      Species        = sp,
      Clinical_Metric = names(rho_vals)[sig_idx],
      Spearman_rho   = rho_vals[sig_idx],
      p_value        = p_vals[sig_idx]
    )
  }
}

sig_corr_df <- do.call(rbind, significant_list)
write.csv(sig_corr_df, "fecallog_spearman_clinical.csv", row.names = FALSE)
cat("Spearman correlation results saved: fecallog_spearman_clinical.csv\n")

# Scatter plot: example visualization (DBP vs oral log-abundance)
sig_species <- unique(sig_corr_df$Species)
plot_data <- data_corr %>% filter(Species %in% sig_species)

ggscatter(plot_data, x = "DBP", y = "orallog",
          color = "Group",
          palette = c("#75ca9d", "#ffeb33", "#ca75dd", "#eaa4a4"),
          size = 3, shape = 1,
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.coeff.args = list(method = "spearman", label.sep = "\n")) +
  facet_wrap(~ Species, scales = "free", nrow = 1)
