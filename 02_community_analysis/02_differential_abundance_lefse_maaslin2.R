# =============================================================================
# Differential Abundance Analysis: LEfSe and MaAsLin2
# Purpose : Identify microbial taxa that are significantly enriched in one
#           disease group relative to another, using two complementary methods:
#           (1) LEfSe (non-parametric + LDA effect size, via microeco)
#           (2) MaAsLin2 (multivariable linear model, FDR-corrected)
#           Taxa consistently identified by both methods at the same direction
#           are reported as high-confidence signals (Figure S2, S3).
#
# Comparisons: HC vs GC, HC vs CRC, GC vs CRC (oral and fecal separately)
# Taxonomy:    Family level (primary); Species level (secondary/exploratory)
# FDR threshold: p.adj < 0.05 (LEfSe), q < 0.10 (MaAsLin2)
#
# INPUT  : data/yonsei_r15000.RDS — rarefied Yonsei phyloseq object
#          Table_S2.Differentially_enriched_taxa_related_to_Figure1.fam.xlsx
#            (pre-computed summary; Sheet 1 = LEfSe, Sheet 2 = MaAsLin2)
# OUTPUT : lefse.<site>.<comparison>.<level>.csv  — LEfSe results per comparison
#          MaAsLin2_<site><comparison>_<level>/   — MaAsLin2 output directory
#          MaAsLin2_<site><comparison>_<level>.csv — MaAsLin2 significant results
#          Table_S2_consensus_Family_Species.xlsx  — consensus taxa (both methods)
# =============================================================================

library(phyloseq)
library(microbiome)
library(microeco)
library(file2meco)
library(Maaslin2)
library(tidyverse)
library(readxl)
library(writexl)
library(ggpubr)
library(magrittr)

d <- readRDS("data/yonsei_r15000.RDS")

# -----------------------------------------------------------------------------
# 1. LEfSe: pairwise group comparisons at Family level (Figure S2)
# -----------------------------------------------------------------------------
# All six pairwise comparisons (oral + fecal × HC-GC, HC-CRC, GC-CRC)
run_lefse <- function(ps_subset, site_label, comparison_label, taxa_level = "Family") {
  meco <- phyloseq2meco(ps_subset)
  meco$cal_abund()
  t1 <- trans_diff$new(dataset = meco, method = "lefse",
                       group = "Group", alpha = 0.05, taxa_level = taxa_level)
  res <- t1$res_diff
  out_file <- paste0("lefse.", site_label, ".", comparison_label, ".",
                     tolower(taxa_level), ".csv")
  write.csv(res, out_file, row.names = FALSE)
  cat("LEfSe saved:", out_file, "\n")
  return(res)
}

# Oral samples
oral_ps <- ps_filter(d, Type == "Mouth")

run_lefse(ps_filter(oral_ps, Group %in% c("HC", "GC")),
          site_label = "oral", comparison_label = "HC-GC")
run_lefse(ps_filter(oral_ps, Group %in% c("HC", "CRC")),
          site_label = "oral", comparison_label = "HC-CRC")
run_lefse(ps_filter(oral_ps, Group %in% c("GC", "CRC")),
          site_label = "oral", comparison_label = "GC-CRC")

# Fecal samples
feces_ps <- ps_filter(d, Type == "Feces")

run_lefse(ps_filter(feces_ps, Group %in% c("HC", "GC")),
          site_label = "feces", comparison_label = "HC-GC")
run_lefse(ps_filter(feces_ps, Group %in% c("HC", "CRC")),
          site_label = "feces", comparison_label = "HC-CRC")
run_lefse(ps_filter(feces_ps, Group %in% c("GC", "CRC")),
          site_label = "feces", comparison_label = "GC-CRC")

# -----------------------------------------------------------------------------
# 2. MaAsLin2: multivariable linear model at Family level (Figure S3)
# -----------------------------------------------------------------------------
run_maaslin2 <- function(ps_subset, site_label, comparison_label,
                         ref_group, taxa_level = "Family",
                         q_threshold = 0.10) {

  # Agglomerate to specified taxonomic level
  ps_agg <- tax_glom(ps_subset, taxrank = taxa_level, NArm = TRUE)

  # Build clean taxonomy labels, resolving NAs and duplicates
  tax_df <- as.data.frame(tax_table(ps_agg))
  tax_df$Taxa_clean <- tax_df[[taxa_level]]
  tax_df$Taxa_clean[is.na(tax_df$Taxa_clean) | tax_df$Taxa_clean == ""] <- "Unclassified"
  tax_df$Taxa_clean <- make.unique(tax_df$Taxa_clean)

  # OTU matrix: rows = samples, columns = taxa
  otu_mat <- as(otu_table(ps_agg), "matrix")
  if (taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)
  colnames(otu_mat) <- tax_df$Taxa_clean

  meta_df <- as(sample_data(ps_agg), "data.frame")
  otu_mat  <- otu_mat[rownames(meta_df), , drop = FALSE]

  # Total sum scaling (TSS) normalization
  otu_tss <- sweep(otu_mat, 1, rowSums(otu_mat), "/")
  otu_tss[!is.finite(otu_tss)] <- 0

  input_meta <- meta_df %>%
    mutate(Group = relevel(factor(Group), ref = ref_group)) %>%
    as.data.frame()

  input_data <- as.data.frame(otu_tss, check.names = FALSE)
  stopifnot(nrow(input_data) == nrow(input_meta))

  out_dir <- paste0("MaAsLin2_", site_label, comparison_label, "_",
                    tolower(taxa_level))

  fit <- Maaslin2(
    input_data      = input_data,
    input_metadata  = input_meta,
    output          = out_dir,
    fixed_effects   = "Group",
    normalization   = "NONE",
    transform       = "LOG",
    analysis_method = "LM",
    correction      = "BH",
    max_significance = q_threshold
  )

  res_sig <- read_tsv(file.path(out_dir, "significant_results.tsv"),
                      show_col_types = FALSE)
  out_csv <- paste0("MaAsLin2_", site_label, comparison_label, "_",
                    tolower(taxa_level), ".csv")
  write.csv(res_sig, out_csv, row.names = FALSE)
  cat("MaAsLin2 saved:", out_csv, "\n")
  return(res_sig)
}

# Fecal HC vs GC (reference = HC)
run_maaslin2(ps_filter(feces_ps, Group %in% c("HC", "GC")),
             site_label = "feces", comparison_label = "HC-GC", ref_group = "HC")

# Fecal HC vs CRC (reference = HC)
run_maaslin2(ps_filter(feces_ps, Group %in% c("HC", "CRC")),
             site_label = "feces", comparison_label = "HC-CRC", ref_group = "HC")

# Fecal GC vs CRC (reference = GC)
run_maaslin2(ps_filter(feces_ps, Group %in% c("GC", "CRC")),
             site_label = "feces", comparison_label = "GC-CRC", ref_group = "GC")

# Oral HC vs GC (reference = HC)
run_maaslin2(ps_filter(oral_ps, Group %in% c("HC", "GC")),
             site_label = "oral", comparison_label = "HC-GC", ref_group = "HC")

# Oral HC vs CRC (reference = HC)
run_maaslin2(ps_filter(oral_ps, Group %in% c("HC", "CRC")),
             site_label = "oral", comparison_label = "HC-CRC", ref_group = "HC")

# Oral GC vs CRC (reference = GC)
run_maaslin2(ps_filter(oral_ps, Group %in% c("GC", "CRC")),
             site_label = "oral", comparison_label = "GC-CRC", ref_group = "GC")

# -----------------------------------------------------------------------------
# 3. Consensus: taxa significant in BOTH LEfSe AND MaAsLin2 (same direction)
# -----------------------------------------------------------------------------
# Reads from pre-compiled Table S2 (Excel); can also be built directly from
# the per-comparison CSVs above.

path <- "Table_S2.Differentially_enriched_taxa_related_to_Figure1.fam.xlsx"
lefse_all   <- read_excel(path, sheet = 1)
maaslin_all <- read_excel(path, sheet = 2)

alpha_thresh <- 0.05
p_col        <- "P.adj"

# --- Family level consensus ---
lefse_fam_sig <- lefse_all %>%
  filter(Rank == "Family") %>%
  mutate(Taxa_clean = str_extract(Taxa, "f__[^|]+") %>% str_remove("^f__")) %>%
  filter(!is.na(Taxa_clean), .data[[p_col]] < alpha_thresh) %>%
  transmute(Comparison, Site, Rank, Taxa_clean,
            LEfSe_Group = Group, LEfSe_LDA = LDA, LEfSe_p = .data[[p_col]])

maaslin_fam_sig <- maaslin_all %>%
  filter(Rank == "Family") %>%
  mutate(Taxa_clean = Taxa) %>%
  filter(!is.na(Taxa_clean), .data[[p_col]] < alpha_thresh) %>%
  transmute(Comparison, Site, Rank, Taxa_clean,
            MaAsLin2_Group = Group, MaAsLin2_Coef = Coef,
            MaAsLin2_Stderr = Stderr, MaAsLin2_p = .data[[p_col]])

family_consensus <- lefse_fam_sig %>%
  inner_join(maaslin_fam_sig, by = c("Comparison", "Site", "Rank", "Taxa_clean")) %>%
  filter(LEfSe_Group == MaAsLin2_Group) %>%                  # directional concordance
  arrange(Comparison, Site, desc(abs(MaAsLin2_Coef))) %>%
  mutate(Stratum = paste(Comparison, Site, Rank, sep = " | ")) %>%
  select(Stratum, Comparison, Site, Rank, Taxa_clean,
         LEfSe_Group, LEfSe_LDA, LEfSe_p,
         MaAsLin2_Group, MaAsLin2_Coef, MaAsLin2_Stderr, MaAsLin2_p)

# --- Species level consensus ---
lefse_spe_sig <- lefse_all %>%
  filter(Rank == "Species") %>%
  mutate(
    Taxa_raw   = str_extract(Taxa, "s__[^|]+") %>% str_remove("^s__") %>% str_squish(),
    Taxa_clean = str_replace(Taxa_raw, " ", ".")   # match MaAsLin2 "Genus.species" format
  ) %>%
  filter(!is.na(Taxa_clean), .data[[p_col]] < alpha_thresh) %>%
  transmute(Comparison, Site, Rank, Taxa_clean,
            LEfSe_Group = Group, LEfSe_LDA = LDA, LEfSe_p = .data[[p_col]])

maaslin_spe_sig <- maaslin_all %>%
  filter(Rank == "Species") %>%
  mutate(Taxa_clean = str_squish(Taxa)) %>%
  filter(!is.na(Taxa_clean), .data[[p_col]] < alpha_thresh) %>%
  transmute(Comparison, Site, Rank, Taxa_clean,
            MaAsLin2_Group = Group, MaAsLin2_Coef = Coef,
            MaAsLin2_Stderr = Stderr, MaAsLin2_p = .data[[p_col]])

species_consensus <- lefse_spe_sig %>%
  inner_join(maaslin_spe_sig, by = c("Comparison", "Site", "Rank", "Taxa_clean")) %>%
  filter(LEfSe_Group == MaAsLin2_Group) %>%
  arrange(Comparison, Site, desc(abs(MaAsLin2_Coef))) %>%
  mutate(Stratum = paste(Comparison, Site, Rank, sep = " | ")) %>%
  select(Stratum, Comparison, Site, Rank, Taxa_clean,
         LEfSe_Group, LEfSe_LDA, LEfSe_p,
         MaAsLin2_Group, MaAsLin2_Coef, MaAsLin2_Stderr, MaAsLin2_p)

write_xlsx(
  list(Family_consensus  = family_consensus,
       Species_consensus = species_consensus),
  path = "Table_S2_consensus_Family_Species.xlsx"
)
cat("Consensus table saved: Table_S2_consensus_Family_Species.xlsx\n")

# -----------------------------------------------------------------------------
# 4. Visualization: MaAsLin2 coefficient bar plot (Figure S3)
# -----------------------------------------------------------------------------
plot_maaslin_coef <- function(df,
                               top_n       = 20,
                               title       = NULL,
                               group_col   = "MaAsLin2_Group",
                               coef_col    = "MaAsLin2_Coef",
                               group_colors = c("HC" = "#75ca9d", "CRC" = "#eaa4a4")) {
  df_plot <- df %>%
    mutate(
      Taxa_label = str_replace_all(Taxa_clean, "_", " "),
      Taxa_label = str_wrap(Taxa_label, width = 28),
      "{coef_col}" := as.numeric(.data[[coef_col]]),
      "{group_col}" := factor(.data[[group_col]])
    ) %>%
    arrange(desc(.data[[coef_col]])) %>%
    slice_head(n = top_n) %>%
    mutate(Taxa_label = factor(Taxa_label, levels = rev(Taxa_label)))

  ggplot(df_plot, aes(x = Taxa_label, y = .data[[coef_col]],
                      fill = .data[[group_col]])) +
    geom_col(width = 0.75, color = "black", linewidth = 0.3) +
    coord_flip() +
    labs(x = NULL, y = "MaAsLin2 coefficient", title = title, fill = NULL) +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5))
}

# Example: fecal HC vs CRC (species level)
maaslin_crc_feces_spe <- maaslin_spe_sig %>%
  filter(Comparison == "CRC - HC", Site == "Feces",
         MaAsLin2_Coef < 1) %>%
  mutate(MaAsLin2_Group = factor(MaAsLin2_Group, levels = c("HC", "CRC")))

p_coef <- plot_maaslin_coef(maaslin_crc_feces_spe,
                             title = "Fecal differential species: HC vs CRC")
ggsave("figS3_maaslin2_feces_HC_CRC_species.pdf", p_coef, width = 6, height = 4)

cat("Figure S3 saved: figS3_maaslin2_feces_HC_CRC_species.pdf\n")
