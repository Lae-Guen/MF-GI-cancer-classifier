# =============================================================================
# Stage-Stratified Analysis of MF-Shared Bacteria (Figure 3A/3B)
# Purpose : Test whether gut abundance of MF-transmitted species differs
#           across cancer stages (HC vs Early vs Late) for GC and CRC,
#           and visualize the stage-dependent enrichment patterns.
#
# Stage grouping:
#   HC    — healthy control (Stage = "HC" in metadata)
#   Early — pathological Stage 0-2
#   Late  — pathological Stage 3-4
#
# Statistics:
#   Pairwise Wilcoxon rank-sum tests (HC vs Early, HC vs Late, Early vs Late)
#   FDR correction (Benjamini-Hochberg) applied globally across all species
#   within each cancer type (not per comparison pair).
#   Rank-biserial correlation (r) with 95% bootstrap CI (nboot = 5,000).
#
# INPUT  : data/yonsei.mtg.csv — shared ASV table with Stage and Group columns
# OUTPUT : fig3A_GC_stage_wilcoxon_effectsize.csv
#          fig3B_CRC_stage_wilcoxon_effectsize.csv
#          Stage-stratified box plots (Figure 3A and 3B)
# =============================================================================

library(tidyverse)
library(rstatix)
library(ggpubr)
library(purrr)

data <- read.csv("data/yonsei.mtg.csv")

# Color palettes for stage groupings
gc_palette  <- c("HC" = "#75ca9d", "Early" = "#FF8080", "Late" = "#800000")
crc_palette <- c("HC" = "#75ca9d", "Early" = "#BD8EF5", "Late" = "#521B88")

# -----------------------------------------------------------------------------
# Helper: assign cancer stage labels
# -----------------------------------------------------------------------------
assign_stage <- function(df, cancer_group) {
  df %>%
    filter(Group %in% c("HC", cancer_group)) %>%
    mutate(
      Rel.abund.feces = (Fecal_count / 15000) * 100,
      Cancer.stage = case_when(
        Stage == "HC"                                ~ "HC",
        suppressWarnings(as.numeric(Stage)) <= 2    ~ "Early",
        suppressWarnings(as.numeric(Stage)) >= 3    ~ "Late",
        TRUE                                         ~ NA_character_
      ),
      Cancer.stage = factor(Cancer.stage, levels = c("HC", "Early", "Late"))
    ) %>%
    filter(!is.na(Cancer.stage))
}

df_gc  <- assign_stage(data, "GC")
df_crc <- assign_stage(data, "CRC")

# -----------------------------------------------------------------------------
# Core function: pairwise Wilcoxon + rank-biserial r per species
# -----------------------------------------------------------------------------
run_stage_stats <- function(df, cancer_label,
                             comparisons = list(c("HC","Early"),
                                                c("HC","Late"),
                                                c("Early","Late")),
                             boot_n = 5000) {
  cat("\n=====", cancer_label, "=====\n")

  # Retain only species with >= 2 observations in at least 2 stage groups
  ok_species <- df %>%
    group_by(Species, Cancer.stage) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Species) %>%
    summarise(
      n_groups      = sum(n > 0),
      min_n_nonzero = min(n[n > 0]),
      .groups = "drop"
    ) %>%
    filter(n_groups >= 2, min_n_nonzero >= 2) %>%
    pull(Species)

  df_ok <- df %>% filter(Species %in% ok_species)

  # Summary statistics per species per stage
  sum_stats <- df_ok %>%
    group_by(Species, Cancer.stage) %>%
    summarise(
      n      = n(),
      median = median(Rel.abund.feces, na.rm = TRUE),
      q1     = quantile(Rel.abund.feces, 0.25, na.rm = TRUE),
      q3     = quantile(Rel.abund.feces, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  # Pairwise Wilcoxon tests with global BH correction across all species
  pw_wilcox <- df_ok %>%
    group_by(Species) %>%
    pairwise_wilcox_test(Rel.abund.feces ~ Cancer.stage,
                         p.adjust.method = "none") %>%
    ungroup() %>%
    mutate(q_value = p.adjust(p, method = "BH")) %>%
    transmute(Species, group1, group2, p_value = p, q_value)

  # Rank-biserial r with 95% bootstrap CI for each comparison pair × species
  # Uses tryCatch to skip species with single-level factor (edge case)
  effsize_list <- map_dfr(comparisons, function(comp) {
    ok_sp <- df_ok %>%
      filter(Cancer.stage %in% comp) %>%
      group_by(Species, Cancer.stage) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(Species) %>%
      summarise(
        has_g1 = any(Cancer.stage == comp[1] & n >= 2),
        has_g2 = any(Cancer.stage == comp[2] & n >= 2),
        .groups = "drop"
      ) %>%
      filter(has_g1, has_g2) %>%
      pull(Species)

    if (length(ok_sp) == 0) return(NULL)

    map_dfr(ok_sp, function(sp) {
      d_sp <- df_ok %>%
        filter(Species == sp, Cancer.stage %in% comp) %>%
        mutate(Cancer.stage = factor(as.character(Cancer.stage), levels = comp)) %>%
        droplevels()

      if (nlevels(d_sp$Cancer.stage) < 2) {
        message("  [skip: single level] ", sp, " (", comp[1], " vs ", comp[2], ")")
        return(NULL)
      }

      tryCatch({
        d_sp %>%
          wilcox_effsize(
            Rel.abund.feces ~ Cancer.stage,
            ci         = TRUE,
            conf.level = 0.95,
            nboot      = boot_n
          ) %>%
          transmute(
            Species   = sp,
            group1    = comp[1],
            group2    = comp[2],
            r         = effsize,
            r_CI_low  = conf.low,
            r_CI_high = conf.high,
            magnitude
          )
      }, error = function(e) {
        message("  [skip: error] ", sp, " (", comp[1], " vs ", comp[2], "): ",
                conditionMessage(e))
        NULL
      })
    })
  })

  # Merge Wilcoxon results, effect sizes, and summary statistics
  results <- pw_wilcox %>%
    left_join(effsize_list, by = c("Species", "group1", "group2")) %>%
    left_join(sum_stats %>% rename(group1 = Cancer.stage,
                                   n_g1 = n, median_g1 = median,
                                   q1_g1 = q1, q3_g1 = q3),
              by = c("Species", "group1")) %>%
    left_join(sum_stats %>% rename(group2 = Cancer.stage,
                                   n_g2 = n, median_g2 = median,
                                   q1_g2 = q1, q3_g2 = q3),
              by = c("Species", "group2")) %>%
    mutate(
      delta_median  = round(median_g2 - median_g1, 6),
      log2FC_median = log2((median_g2 + 1e-6) / (median_g1 + 1e-6)),
      r         = round(r, 3),
      r_CI_low  = round(r_CI_low, 3),
      r_CI_high = round(r_CI_high, 3),
      p_value   = signif(p_value, 3),
      q_value   = signif(q_value, 3)
    ) %>%
    arrange(group1, group2, q_value)

  results
}

# -----------------------------------------------------------------------------
# Run statistics: GC (Figure 3A) and CRC (Figure 3B)
# -----------------------------------------------------------------------------
results_gc  <- run_stage_stats(df_gc,  cancer_label = "GC  (Figure 3A)")
results_crc <- run_stage_stats(df_crc, cancer_label = "CRC (Figure 3B)")

# Print top results for reporting
cat("\n--- GC: HC vs Early (top 10 by q-value) ---\n")
results_gc %>%
  filter(group1 == "HC", group2 == "Early") %>%
  select(Species, n_g1, median_g1, n_g2, median_g2,
         delta_median, p_value, q_value, r, r_CI_low, r_CI_high, magnitude) %>%
  slice_min(q_value, n = 10) %>%
  print(width = 150)

cat("\n--- GC: HC vs Late (top 10 by q-value) ---\n")
results_gc %>%
  filter(group1 == "HC", group2 == "Late") %>%
  select(Species, n_g1, median_g1, n_g2, median_g2,
         delta_median, p_value, q_value, r, r_CI_low, r_CI_high, magnitude) %>%
  slice_min(q_value, n = 10) %>%
  print(width = 150)

write.csv(results_gc,  "fig3A_GC_stage_wilcoxon_effectsize.csv",  row.names = FALSE)
write.csv(results_crc, "fig3B_CRC_stage_wilcoxon_effectsize.csv", row.names = FALSE)
cat("Saved: fig3A_GC_stage_wilcoxon_effectsize.csv\n")
cat("Saved: fig3B_CRC_stage_wilcoxon_effectsize.csv\n")

# -----------------------------------------------------------------------------
# Stage-stratified visualization (Figure 3A and 3B box plots)
# -----------------------------------------------------------------------------
data_plot <- data %>%
  mutate(
    Rel.abund.feces = (Fecal_count / 15000) * 100,
    fecallog = log10(Rel.abund.feces),
    Cancer.stage = case_when(
      Stage == "HC" ~ "HC",
      as.numeric(Stage) <= 2 ~ "Early",
      as.numeric(Stage) >= 3 ~ "Late"
    )
  )

# Species to display for GC (Figure 3A)
gc_species <- c("Streptococcus parasanguinis",
                "Veillonella_A atypica",
                "Streptococcus vestibularis",
                "Streptococcus anginosus")

# Species to display for CRC (Figure 3B)
crc_species <- c("Peptostreptococcus stomatis",
                 "Streptococcus anginosus",
                 "Streptococcus parasanguinis",
                 "Granulicatella adiacens")

my_comparisons <- list(c("HC", "Early"), c("HC", "Late"))

# Figure 3A: GC
df_gc_plot <- data_plot %>%
  filter(Group %in% c("HC", "GC"), Species %in% gc_species,
         !is.na(Cancer.stage)) %>%
  mutate(
    Species      = factor(Species, levels = gc_species),
    Cancer.stage = factor(Cancer.stage, levels = c("HC", "Early", "Late"))
  )

p_gc <- ggboxplot(df_gc_plot, x = "Cancer.stage", y = "fecallog",
                  fill = "Cancer.stage",
                  palette = unname(gc_palette)) +
  stat_compare_means(
    comparisons = my_comparisons, size = 3,
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05),
                       symbols   = c("****", "***", "**", "*"))
  ) +
  facet_wrap(~ Species, scales = "free", nrow = 1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave("fig3A_GC_stage_boxplot.pdf", p_gc, width = 9.5, height = 5)

# Figure 3B: CRC
df_crc_plot <- data_plot %>%
  filter(Group %in% c("HC", "CRC"), Species %in% crc_species,
         !is.na(Cancer.stage)) %>%
  mutate(
    Species      = factor(Species, levels = crc_species),
    Cancer.stage = factor(Cancer.stage, levels = c("HC", "Early", "Late"))
  )

p_crc <- ggboxplot(df_crc_plot, x = "Cancer.stage", y = "fecallog",
                   fill = "Cancer.stage",
                   palette = unname(crc_palette)) +
  stat_compare_means(
    comparisons = my_comparisons, size = 3,
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05),
                       symbols   = c("****", "***", "**", "*"))
  ) +
  facet_wrap(~ Species, scales = "free", nrow = 1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave("fig3B_CRC_stage_boxplot.pdf", p_crc, width = 9.5, height = 5)

cat("Figures saved: fig3A_GC_stage_boxplot.pdf, fig3B_CRC_stage_boxplot.pdf\n")
