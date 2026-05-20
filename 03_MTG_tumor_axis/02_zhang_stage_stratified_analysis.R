# =============================================================================
# Stage-Stratified Transmission Index Analysis (Zhang Cohort)
# Purpose : Test whether the five transmission indices (MTT, MTP, TTG, PTG,
#           MTG) and inter-pathway ASV overlap ratios differ between early-
#           and late-stage cancer in the Zhang cohort (GC and CRC separately).
#
# zhang.xlsx sheet structure:
#   Sheet 1 — MTG (Mouth-to-Feces) shared ASV feature table
#   Sheet 2 — MTT (Mouth-to-Tumor) shared ASV feature table
#   Sheet 3 — MTP (Mouth-to-Paratumor) shared ASV feature table
#   Sheet 4 — TTG (Tumor-to-Feces) shared ASV feature table
#   Sheet 5 — PTG (Paratumor-to-Feces) shared ASV feature table
# Each sheet was produced by applying the logic in 01_compute_MF_index.R
# to the respective sample-type pair in the Zhang dataset.
#
# Stage grouping:
#   Early — TNM stage 0-2
#   Late  — TNM stage 3-4
#
# INPUT  : data/zhang.xlsx        — transmission ASV feature tables (5 sheets)
#          data/zhang.mtg.meta.csv — Zhang sample metadata with TNM_stage
#          data/*.mtgplot.zhang.csv — species-level index tables per pathway
# OUTPUT : data/gc.trans.index.stage.csv  — stage-annotated transmission indices (GC)
#          data/crc.trans.index.stage.csv — stage-annotated transmission indices (CRC)
#          data/gc.pathway.ratio.csv / data/crc.pathway.ratio.csv
#          Stage-stratified violin plots and pathway ratio box plots
# =============================================================================

library(readxl)
library(tidyverse)
library(ggpubr)

# Color palettes (consistent with Figure 3)
gc_palette  <- c("Early" = "#FF8080", "Late" = "#800000")
crc_palette <- c("Early" = "#BD8EF5", "Late" = "#521B88")

# =============================================================================
# Part 1. Stage-stratified transmission index (MTT, MTP, TTG, PTG)
# =============================================================================

# ── Helper: load one transmission sheet, compute index, annotate stage ──────
load_trans_index <- function(xlsx_path, sheet_num, sample_type,
                              group_label, index_label) {
  dat <- read_excel(xlsx_path, sheet = sheet_num) %>%
    filter(Group == group_label)

  index_df <- aggregate(Fecal_count ~ SampleID2, data = dat, sum) %>%
    mutate(MTG.index = (Fecal_count / 15000) * 100)

  # Retrieve TNM_stage from metadata (Tumor-type rows only)
  meta <- read.csv("data/zhang.mtg.meta.csv") %>%
    filter(Group == group_label, Type == sample_type) %>%
    select(SampleID2, TNM_stage)

  index_stage <- left_join(index_df, meta, by = "SampleID2") %>%
    filter(!is.na(TNM_stage)) %>%
    mutate(
      Cancer.stage = case_when(
        TNM_stage <= 2 ~ "Early",
        TNM_stage >= 3 ~ "Late"
      ),
      Cancer.stage = factor(Cancer.stage, levels = c("Early", "Late")),
      Model = index_label
    )
  return(index_stage)
}

xlsx_path <- "data/zhang.xlsx"

# -- GC: combine all four indices into one long data frame -------------------
gc_mtt <- load_trans_index(xlsx_path, 2, "Tumor", "GC",  "MTT")
gc_mtp <- load_trans_index(xlsx_path, 3, "Para",  "GC",  "MTP")
gc_ttg <- load_trans_index(xlsx_path, 4, "Tumor", "GC",  "TTG")
gc_ptg <- load_trans_index(xlsx_path, 5, "Para",  "GC",  "PTG")

gc_combined <- bind_rows(gc_mtt, gc_mtp, gc_ttg, gc_ptg) %>%
  mutate(Model = factor(Model, levels = c("MTT", "MTP", "TTG", "PTG")))

write.csv(gc_combined, "data/gc.trans.index.stage.csv", row.names = FALSE)

# -- CRC: combine all four indices -------------------------------------------
crc_mtt <- load_trans_index(xlsx_path, 2, "Tumor", "CRC", "MTT")
crc_mtp <- load_trans_index(xlsx_path, 3, "Para",  "CRC", "MTP")
crc_ttg <- load_trans_index(xlsx_path, 4, "Tumor", "CRC", "TTG")
crc_ptg <- load_trans_index(xlsx_path, 5, "Para",  "CRC", "PTG")

crc_combined <- bind_rows(crc_mtt, crc_mtp, crc_ttg, crc_ptg) %>%
  mutate(Model = factor(Model, levels = c("MTT", "MTP", "TTG", "PTG")))

write.csv(crc_combined, "data/crc.trans.index.stage.csv", row.names = FALSE)
cat("Saved: data/gc.trans.index.stage.csv, data/crc.trans.index.stage.csv\n")

# -- Violin plots: stage comparison of each transmission index ----------------
plot_stage_violin <- function(dat, palette, title_label) {
  ggviolin(dat, x = "Cancer.stage", y = "MTG.index",
           fill = "Cancer.stage", palette = palette,
           add = "boxplot",
           add.params = list(fill = "white", alpha = 0.5)) +
    stat_compare_means(size = 3) +
    facet_wrap(~ Model, scales = "free", nrow = 1) +
    labs(title = title_label) +
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.position = "none")
}

p_gc_stage  <- plot_stage_violin(gc_combined,  gc_palette,  "GC — Stage comparison")
p_crc_stage <- plot_stage_violin(crc_combined, crc_palette, "CRC — Stage comparison")

ggsave("fig4C_GC_stage_violin.pdf",  p_gc_stage,  width = 10, height = 5)
ggsave("fig4D_CRC_stage_violin.pdf", p_crc_stage, width = 10, height = 5)

# =============================================================================
# Part 2. Inter-pathway ASV overlap ratios
# Purpose : Quantify the fraction of ASVs in each pathway (MTT, MTP, TTG, PTG)
#           that are shared with another pathway (e.g., MTT ∩ MTG / all MTT).
#           This reveals how much of the tissue-associated oral microbiome
#           overlaps with the direct mouth-to-feces transmission pool.
# =============================================================================

# ── Helper: compute pairwise ASV overlap ratio ──────────────────────────────
# Returns: for each sample, what fraction of pathway_b's Fecal_count
#          corresponds to ASVs that are ALSO present in pathway_a
calculate_asv_ratio <- function(pathway_a, pathway_b, sample_ids) {
  sapply(sample_ids, function(sid) {
    a_asvs <- pathway_a$ASV[pathway_a$SampleID2 == sid]
    b_data  <- pathway_b[pathway_b$SampleID2 == sid, ]

    common_sum <- b_data %>%
      filter(ASV %in% a_asvs) %>%
      pull(Fecal_count) %>%
      sum()

    total_sum <- sum(b_data$Fecal_count)

    if (total_sum > 0) (common_sum / total_sum) * 100 else NA_real_
  })
}

# ── Process one group (GC or CRC) ───────────────────────────────────────────
compute_pathway_ratios <- function(group_label, meta_type = "Tumor",
                                    out_file, palette) {
  meta <- read.csv("data/zhang.mtg.meta.csv") %>%
    filter(Group == group_label, Type == meta_type)
  sample_ids <- unique(meta$SampleID2)

  mtg <- read_excel(xlsx_path, sheet = 1) %>% filter(SampleID2 %in% sample_ids)
  mtt <- read_excel(xlsx_path, sheet = 2) %>% filter(SampleID2 %in% sample_ids)
  mtp <- read_excel(xlsx_path, sheet = 3) %>% filter(SampleID2 %in% sample_ids)
  ttg <- read_excel(xlsx_path, sheet = 4) %>% filter(SampleID2 %in% sample_ids)
  ptg <- read_excel(xlsx_path, sheet = 5) %>% filter(SampleID2 %in% sample_ids)

  # Pathway overlap ratios (6 combinations)
  ratios_df <- data.frame(
    SampleID2    = sample_ids,
    mtt_mtg      = calculate_asv_ratio(mtt, mtg, sample_ids),  # MTT ∩ MTG / MTG
    mtp_mtg      = calculate_asv_ratio(mtp, mtg, sample_ids),  # MTP ∩ MTG / MTG
    mtg_mtt      = calculate_asv_ratio(mtg, mtt, sample_ids),  # MTG ∩ MTT / MTT
    mtg_mtp      = calculate_asv_ratio(mtg, mtp, sample_ids),  # MTG ∩ MTP / MTP
    mtg_ttg      = calculate_asv_ratio(mtg, ttg, sample_ids),  # MTG ∩ TTG / TTG
    mtg_ptg      = calculate_asv_ratio(mtg, ptg, sample_ids)   # MTG ∩ PTG / PTG
  )

  write.csv(ratios_df, out_file, row.names = FALSE)
  cat("Saved:", out_file, "\n")

  # Annotate cancer stage
  stage_info <- read.csv("data/zhang.mtg.meta.csv") %>%
    filter(Group == group_label) %>%
    filter(!is.na(TNM_stage)) %>%
    mutate(Cancer.stage = case_when(TNM_stage <= 2 ~ "Early",
                                    TNM_stage >= 3 ~ "Late")) %>%
    distinct(SampleID2, Cancer.stage)

  ratio_stage <- left_join(ratios_df, stage_info, by = "SampleID2") %>%
    filter(!is.na(Cancer.stage)) %>%
    pivot_longer(cols = -c(SampleID2, Cancer.stage),
                 names_to = "pathway", values_to = "percentage") %>%
    mutate(
      Cancer.stage = factor(Cancer.stage, levels = c("Early", "Late")),
      pathway = factor(pathway, levels = c("mtt_mtg", "mtp_mtg",
                                           "mtg_mtt", "mtg_mtp",
                                           "mtg_ttg", "mtg_ptg"))
    )

  # Box plots: stage comparison of pathway ratios
  p_box <- ggboxplot(ratio_stage, x = "Cancer.stage", y = "percentage",
                     color = "Cancer.stage", palette = palette,
                     add = "jitter", shape = "Cancer.stage") +
    stat_compare_means(size = 3.5) +
    facet_wrap(~ pathway, scales = "free", nrow = 2) +
    labs(title = paste(group_label, "— Pathway overlap ratios by stage")) +
    theme(axis.title.x = element_blank())

  return(p_box)
}

p_gc_ratio  <- compute_pathway_ratios("GC",  out_file = "data/gc.pathway.ratio.csv",
                                       palette = gc_palette)
p_crc_ratio <- compute_pathway_ratios("CRC", out_file = "data/crc.pathway.ratio.csv",
                                       palette = crc_palette)

ggsave("figS_GC_pathway_ratio_stage.pdf",  p_gc_ratio,  width = 8, height = 6)
ggsave("figS_CRC_pathway_ratio_stage.pdf", p_crc_ratio, width = 8, height = 6)

# =============================================================================
# Part 3. Species-level MTG/MTT/MTP/TTG/PTG association scatter plots
# Purpose : Identify species whose oral-fecal and oral-tissue transmission
#           indices are correlated, validating the tissue-mediated route.
#           Uses pre-computed species-level index tables (from 04_MF_species_profiling.R
#           applied to the Zhang dataset).
# =============================================================================

load_species_index <- function(prefix, cancer_label) {
  mtg  <- read.csv(paste0("data/", cancer_label, ".mtgplot.zhang.csv"))
  mtt  <- read.csv(paste0("data/", cancer_label, ".mttplot.zhang.csv"))
  mtp  <- read.csv(paste0("data/", cancer_label, ".mtpplot.zhang.csv"))
  ttg  <- read.csv(paste0("data/", cancer_label, ".ttgplot.zhang.csv"))
  ptg  <- read.csv(paste0("data/", cancer_label, ".ptgplot.zhang.csv"))
  prev <- read.csv("data/mtg.oralCRC.preval.csv")

  # Retain only species present across all five transmission tables
  common_species <- Reduce(intersect, list(mtg$Species, mtt$Species,
                                            mtp$Species, ttg$Species, ptg$Species))
  mtg  <- filter(mtg,  Species %in% common_species)
  mtt  <- filter(mtt,  Species %in% common_species)
  prev <- filter(prev, Species %in% common_species)

  merged <- left_join(mtg, mtt, by = "Species") %>%
    left_join(prev, by = "Species")

  out_file <- paste0("data/", tolower(cancer_label), ".mtg.mtt.csv")
  write.csv(merged, out_file, row.names = FALSE)
  cat("Saved:", out_file, "\n")

  return(merged)
}

# -- CRC: MTG vs MTT/MTT index scatter plot -----------------------------------
crc_assoc <- load_species_index(prefix = "crc", cancer_label = "crc")

# Labeled threshold: MTG.index >= 0.2 AND MTT.index >= 0.1
labeled <- crc_assoc %>% filter(MTG.index >= 0.2, MTT.index >= 0.1)

# Bubble plot: point size = oral prevalence (Preval)
p_scatter_bubble <- ggscatter(crc_assoc, x = "MTG.index", y = "MTT.index",
                               color = "grey", size = "Preval", label = NULL) +
  geom_text(data = labeled, aes(label = Species),
            vjust = 0.5, hjust = 0.5, size = 3) +
  geom_vline(xintercept = 0.2, linetype = "dotted", color = "#2F2F2F") +
  geom_hline(yintercept = 0.1, linetype = "dotted", color = "#2F2F2F") +
  labs(x = "MTG index (species level)",
       y = "MTT index (species level)")

# Spearman correlation scatter plot
p_scatter_corr <- ggscatter(crc_assoc, x = "MTG.index", y = "MTT.index",
                              color = "black", size = 3, shape = 1,
                              add = "reg.line",
                              add.params = list(color = "blue", fill = "lightgray"),
                              conf.int = TRUE,
                              cor.coef = TRUE,
                              cor.coeff.args = list(method = "spearman",
                                                    label.sep = "\n")) +
  labs(x = "MTG index (species level)",
       y = "MTT index (species level)")

ggsave("fig4_CRC_MTG_MTT_scatter.pdf",      p_scatter_bubble, width = 6, height = 5)
ggsave("fig4_CRC_MTG_MTT_correlation.pdf",  p_scatter_corr,   width = 5, height = 5)

cat("All stage-stratified and pathway analysis figures saved.\n")
