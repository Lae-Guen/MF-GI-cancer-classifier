# =============================================================================
# MF Species Profiling
# Purpose : Compute species-level MF index and oral relative abundance for
#           each disease group and generate the scatter plot shown in Figure 2
#           (oral mean abundance vs. MF index, point size = prevalence).
#
# Key steps:
#   1. Aggregate fecal ASV abundance to species level (denominator for MF index)
#   2. Compute species-level MF index per group (shared fecal counts / mean fecal
#      species abundance, normalized to rarefaction depth)
#   3. Compute oral prevalence of each MF-transmitted species per group
#   4. Compute mean oral relative abundance of MF-transmitted species per group
#   5. Produce Figure 2C-F scatter plots
#
# INPUT  : data/yonsei_r15000.RDS  — rarefied Yonsei phyloseq object
#          data/yonsei.mtg.csv     — shared ASV table from script 01
# OUTPUT : feces.species_means.csv  — mean fecal abundance per species
#          <Group>.mtg.index.species.csv — species-level MF index per group
#          oral<Group>.mtg.rel.csv  — mean oral relative abundance per group
#          mtg.oral<Group>.preval.csv — oral prevalence per group
#          <Group>.rel.mtg.prev.species.csv — combined data for Figure 2
# =============================================================================

library(phyloseq)
library(microViz)   # ps_filter(), otu_tibble(), tax_tibble(), aggregate_taxa()
library(microbiome)
library(tidyverse)
library(ggpubr)

d <- readRDS("data/yonsei_r15000.RDS")
mtg <- read.csv("data/yonsei.mtg.csv")

# Group-specific color palette (consistent across all figures)
group_palette <- c("HC" = "#75ca9d", "MD" = "#ffeb33",
                   "GC" = "#ca75dd", "CRC" = "#eaa4a4")

# -----------------------------------------------------------------------------
# Step 1. Compute mean fecal species abundance (denominator for MF index)
# -----------------------------------------------------------------------------
feces_ps <- ps_filter(d, Type == "Feces")
feces_spe <- aggregate_taxa(feces_ps, "Species")

otu_data <- otu_tibble(feces_spe)
write.csv(otu_data, "feces.species.csv", row.names = FALSE)

feces_mat <- read.csv("feces.species.csv", header = TRUE, row.names = 1)
species_means <- rowMeans(feces_mat)
feces_means_df <- data.frame(Species = rownames(feces_mat), Avg.mean = species_means)
write.csv(feces_means_df, "feces.species_means.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Helper function: process one disease group
# Computes MF index, oral prevalence, oral relative abundance,
# and merges into a single data frame for plotting.
# -----------------------------------------------------------------------------
process_group <- function(group_label, n_samples_oral) {

  # -- Species-level MF index ------------------------------------------------
  gr_mtg <- mtg %>% filter(Group == group_label)
  mtg_idx <- aggregate(Fecal_count ~ Species, data = gr_mtg, sum) %>%
    left_join(feces_means_df, by = "Species") %>%
    mutate(MTG       = Fecal_count / Avg.mean,
           MTG.index = (MTG / 15000) * 100) %>%
    select(Species, MTG.index)
  write.csv(mtg_idx, paste0(group_label, ".mtg.index.species.csv"), row.names = FALSE)

  # -- Oral prevalence of MF-transmitted species ------------------------------
  oral_ps <- ps_filter(d, Type == "Mouth", Group == group_label)
  oral_otu <- otu_tibble(oral_ps)
  tax_info <- tax_tibble(oral_ps) %>% select(FeatureID, Genus)

  # Get ASV IDs for this group's MF-transmitted taxa
  mtg_asvs <- gr_mtg %>%
    distinct(FeatureID, .keep_all = TRUE)

  # Maximum oral prevalence per species (across all ASVs assigned to that species)
  sp_asv_map <- mtg_asvs %>%
    group_by(Genus) %>%
    summarise(FeatureIDs = list(FeatureID), .groups = "drop") %>%
    filter(!is.na(Genus))

  asv_prev <- oral_otu %>%
    rowwise() %>%
    mutate(Prevalence = sum(c_across(-FeatureID) > 0)) %>%
    select(FeatureID, Prevalence)

  sp_prevalence <- sp_asv_map %>%
    rowwise() %>%
    mutate(Max_Prevalence = max(
      asv_prev %>% filter(FeatureID %in% FeatureIDs) %>% pull(Prevalence),
      na.rm = TRUE
    )) %>%
    mutate(Preval = (Max_Prevalence / n_samples_oral) * 100) %>%
    select(Genus, Preval)
  write.csv(sp_prevalence,
            paste0("mtg.oral", group_label, ".preval.csv"), row.names = FALSE)

  # -- Mean oral relative abundance of MF-transmitted species ----------------
  mtg_asv_ids <- gr_mtg$FeatureID %>% unique()
  oral_mtg_otu <- oral_otu %>% filter(FeatureID %in% mtg_asv_ids)
  oral_mtg_otu <- left_join(oral_mtg_otu, tax_info, by = "FeatureID") %>%
    group_by(Genus) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
    mutate(across(where(is.numeric), ~ . / 15000))

  oral_rel <- oral_mtg_otu %>%
    pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance") %>%
    group_by(Genus) %>%
    summarise(Oralmean = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    rename_with(~ paste0("Oralmean_", group_label), "Oralmean")
  write.csv(oral_rel, paste0("oral", group_label, ".mtg.rel.csv"), row.names = FALSE)

  # -- Merge and return combined data frame -----------------------------------
  mf <- mtg_idx %>%
    left_join(oral_rel, by = "Species") %>%
    left_join(sp_prevalence, by = c("Species" = "Genus")) %>%
    left_join(feces_means_df, by = "Species") %>%
    mutate(orallog = log10(.data[[paste0("Oralmean_", group_label)]]))
  write.csv(mf, paste0(group_label, ".rel.mtg.prev.species.csv"), row.names = FALSE)

  return(mf)
}

# -----------------------------------------------------------------------------
# Step 2. Process each group (sample counts from the Yonsei cohort)
# -----------------------------------------------------------------------------
mf_hc  <- process_group("HC",  n_samples_oral = 129)
mf_md  <- process_group("MD",  n_samples_oral = 213)
mf_gc  <- process_group("GC",  n_samples_oral = 77)
mf_crc <- process_group("CRC", n_samples_oral = 86)

# -----------------------------------------------------------------------------
# Step 3. Figure 2C-F: oral abundance vs MF index scatter plots
# Inclusion criteria: orallog >= -4 (mean oral relative abundance >= 0.01%)
#                     MTG.index >= 2
# -----------------------------------------------------------------------------
plot_mf_scatter <- function(mf_df, group_label, label_col = "Species") {
  filtered <- mf_df %>% filter(orallog >= -4, MTG.index >= 2)

  ggscatter(mf_df, x = "orallog", y = "MTG.index",
            color = "grey", size = "Preval", label = NULL) +
    geom_text(data = filtered, aes(label = .data[[label_col]]),
              vjust = 0.5, hjust = 0.5, size = 3) +
    geom_vline(xintercept = -4, linetype = "dotted", color = "#2F2F2F") +
    geom_hline(yintercept = 2,  linetype = "dotted", color = "#2F2F2F") +
    labs(
      title = group_label,
      x = "Mean oral relative abundance (log10)",
      y = "MF index (species level)"
    )
}

p_hc  <- plot_mf_scatter(mf_hc,  "HC")
p_md  <- plot_mf_scatter(mf_md,  "MD")
p_gc  <- plot_mf_scatter(mf_gc,  "GC")
p_crc <- plot_mf_scatter(mf_crc, "CRC")

# Save individual panels (to be assembled in Illustrator as Figure 2C-F)
ggsave("fig2C_HC_mf_scatter.pdf",  p_hc,  width = 5, height = 5)
ggsave("fig2D_MD_mf_scatter.pdf",  p_md,  width = 5, height = 5)
ggsave("fig2E_GC_mf_scatter.pdf",  p_gc,  width = 5, height = 5)
ggsave("fig2F_CRC_mf_scatter.pdf", p_crc, width = 5, height = 5)

cat("Figure 2 panels saved.\n")
