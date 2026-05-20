# =============================================================================
# Sequence Similarity-Based MF Feature Evaluation and Extraction
# Purpose : (1) Evaluate VSEARCH hit counts across cohorts and identity
#               thresholds using heatmaps (Figure S7).
#           (2) Extract MF-transmitted ASVs matched in external phyloseq
#               objects and build abundance / binary feature tables for
#               external validation.
#
# INPUT  : v3v4_identity/identity_hit_summary_v3v4.tsv — VSEARCH hit summary
#          identity_hit_summary.tsv                     — all-cohort summary
#          data/crc3_Wang_r15000.RDS                    — Wang cohort phyloseq
#          wang_ASV_1.00.txt                            — matched ASV IDs (100%)
#          query_to_wang_1.00.tsv                       — query-target mapping
# OUTPUT : Heatmap PDFs (matched ASVs, log10, normalized ratio, connectivity)
#          identity_hit_summary_*_with_metrics.csv
#          fcrc_wang_1.00_abundance.csv / binary.csv
#          russo_1.00_abundance.csv / binary.csv
# =============================================================================

library(phyloseq)
library(microbiome)
library(microViz)
library(tidyverse)
library(magrittr)

# =============================================================================
# Part 1. VSEARCH Hit Evaluation — V3-V4 Cohorts (Russo, Flemer, Wang)
# =============================================================================

# ── 1A. Raw and log10 hit count heatmaps ────────────────────────────────────
df_v3v4 <- read.delim("v3v4_identity/identity_hit_summary_v3v4.tsv",
                       sep = "\t", header = TRUE) %>%
  mutate(
    Cohort   = factor(Cohort, levels = c("flemer", "russo", "wang")),
    Identity = factor(Identity,
                      levels = c(1.00, 0.99, 0.98, 0.97, 0.96, 0.95),
                      labels = c("100%", "99%", "98%", "97%", "96%", "95%"))
  )

# Raw count heatmap
p_raw_v3v4 <- ggplot(df_v3v4,
                      aes(x = Identity, y = Cohort,
                          fill = Matched_external_ASVs)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Matched_external_ASVs), size = 4) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Matched external ASVs by cohort and identity threshold",
       x = "Sequence identity cutoff", y = "Cohort",
       fill = "Matched\nASVs") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("heatmap_matched_external_ASVs_raw.pdf", p_raw_v3v4, width = 7, height = 4.5)

# Log10 heatmap
df_v3v4 <- df_v3v4 %>%
  mutate(log10_Matched = log10(Matched_external_ASVs + 1))

p_log_v3v4 <- ggplot(df_v3v4,
                      aes(x = Identity, y = Cohort, fill = log10_Matched)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Matched_external_ASVs), size = 4) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Matched external ASVs (log10 scale)",
       x = "Sequence identity cutoff", y = "Cohort",
       fill = "log10(count+1)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("heatmap_matched_external_ASVs_log10.pdf", p_log_v3v4, width = 7, height = 4.5)
cat("V3-V4 heatmaps saved.\n")

# =============================================================================
# Part 2. All-Cohort Hit Evaluation with Normalized Metrics
# Cohorts: zhang, baxter, zackular, zeller (fecal-only)
#          russo, flemer, wang, zhang     (oral V3-V4)
# =============================================================================

df_all <- read.delim("identity_hit_summary.tsv", sep = "\t", header = TRUE)

# Total ASV counts per cohort (from QIIME2 feature table inspection)
cohort_totals <- tibble(
  Cohort = c("zhang", "baxter", "zackular", "zeller"),
  Total_external_ASVs = c(43810, 5229, 2914, 5420)
)

n_query_asvs <- 327   # number of transmitted query ASVs

df_all2 <- df_all %>%
  left_join(cohort_totals, by = "Cohort") %>%
  mutate(
    Identity = factor(
      Identity,
      levels = c(1.00, 0.99, 0.98, 0.96, 0.94, 0.90, 0.85, 0.80, 0.70),
      labels = c("100%", "99%", "98%", "96%", "94%", "90%", "85%", "80%", "70%")
    ),
    Cohort = factor(Cohort, levels = c("zhang", "baxter", "zackular", "zeller")),
    # Fraction of external ASVs that were matched
    Matched_ASV_ratio       = Matched_external_ASVs / Total_external_ASVs,
    # Average number of query ASVs hitting each matched external ASV
    Pairs_per_matched_ASV   = Query_to_target_pairs / Matched_external_ASVs,
    # Average matched external ASVs per query transmitted ASV
    Matched_ASVs_per_query  = Matched_external_ASVs / n_query_asvs
  )

write.csv(df_all2, "identity_hit_summary_with_normalized_metrics.csv", row.names = FALSE)

# ── Heatmap helper ───────────────────────────────────────────────────────────
make_heatmap <- function(df, fill_col, label_expr, title_str, fill_label,
                          fill_color = "red", out_pdf, w = 9, h = 5) {
  p <- ggplot(df, aes(x = Identity, y = Cohort, fill = .data[[fill_col]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = !!rlang::parse_expr(label_expr)), size = 4) +
    scale_fill_gradient(low = "white", high = fill_color) +
    labs(title = title_str, x = "Sequence identity cutoff",
         y = "Cohort", fill = fill_label) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_pdf, p, width = w, height = h)
  cat("Saved:", out_pdf, "\n")
}

df_all3 <- df_all2 %>%
  mutate(log10_Matched = log10(Matched_external_ASVs + 1),
         Identity = Identity %<>%
           factor(., levels = c("70%","80%","85%","90%","94%","96%","98%","99%","100%")))

make_heatmap(df_all2, "Matched_external_ASVs",
             "Matched_external_ASVs",
             "Matched external ASVs across all cohorts",
             "Matched\nASVs", "red",
             "allcohort_heatmap_matched_ASVs_raw.pdf")

make_heatmap(df_all3, "log10_Matched",
             "Matched_external_ASVs",
             "Matched external ASVs (log10 scale)",
             "log10(count+1)", "red",
             "allcohort_heatmap_matched_ASVs_log10.pdf")

make_heatmap(df_all2, "Matched_ASV_ratio",
             "sprintf('%.3f', Matched_ASV_ratio)",
             "Normalized matched-ASV ratio",
             "Matched ASV\nratio", "darkorange",
             "allcohort_heatmap_matched_ASV_ratio.pdf")

make_heatmap(df_all2, "Pairs_per_matched_ASV",
             "sprintf('%.2f', Pairs_per_matched_ASV)",
             "Connectivity per matched external ASV",
             "Pairs per\nmatched ASV", "purple",
             "allcohort_heatmap_pairs_per_matched_ASV.pdf")

make_heatmap(df_all2, "Matched_ASVs_per_query",
             "sprintf('%.2f', Matched_ASVs_per_query)",
             paste0("Matched external ASVs per transmitted query ASV (n = ", n_query_asvs, ")"),
             "Matched ASVs\nper query", "forestgreen",
             "allcohort_heatmap_matched_ASVs_per_query.pdf")

# =============================================================================
# Part 3. Extract MF-Matched ASVs from External Phyloseq and Build Feature Tables
# =============================================================================

# ── Helper: orient OTU table to sample × feature ────────────────────────────
get_otu_df <- function(ps_obj) {
  otu <- as.data.frame(otu_table(ps_obj))
  if (taxa_are_rows(ps_obj)) otu <- t(otu)
  as.data.frame(otu)
}

# ── Example: Wang CRC cohort at 100% sequence identity ──────────────────────
target_asvs <- readLines("wang_ASV_1.00.txt") %>%
  trimws() %>%
  .[nchar(.) > 0]

cat("Target ASV count:", length(target_asvs), "\n")

d      <- readRDS("data/crc3_Wang_r15000.RDS")
ps     <- d %>%
  ps_filter(Type == "Feces") %>%
  ps_filter(Group %in% c("HC", "CRC"))

matched_asvs <- intersect(taxa_names(ps), target_asvs)
cat("Matched ASVs in Wang phyloseq:", length(matched_asvs), "/",
    length(target_asvs), "\n")

# Relative abundance normalization, then prune to matched ASVs
ps_rel      <- transform_sample_counts(ps, function(x) x / sum(x))
ps_filtered <- prune_taxa(matched_asvs, ps_rel)

# ASV-level abundance and binary tables
asv_abund  <- get_otu_df(ps_filtered)
asv_binary <- asv_abund; asv_binary[asv_binary > 0] <- 1L

# Append sample metadata
meta_df <- microbiome::meta(ps_filtered) %>%
  mutate(SampleID = rownames(.)) %>%
  select(SampleID, Group)

asv_abund$SampleID  <- rownames(asv_abund)
asv_abund  <- left_join(asv_abund,  meta_df, by = "SampleID")
asv_binary$SampleID <- rownames(asv_binary)
asv_binary <- left_join(asv_binary, meta_df, by = "SampleID")

write.csv(asv_abund,  "fcrc_wang_1.00_abundance.csv", row.names = FALSE)
write.csv(asv_binary, "fcrc_wang_1.00_binary.csv",    row.names = FALSE)
cat("Saved: fcrc_wang_1.00_abundance.csv, fcrc_wang_1.00_binary.csv\n")

# ── Species and Genus level aggregation (optional, uncomment if needed) ──────
# ps_species <- tax_glom(ps_filtered, taxrank = "Species", NArm = FALSE)
# ps_genus   <- tax_glom(ps_filtered, taxrank = "Genus",   NArm = FALSE)
# (Add taxonomy-label renaming and export as needed)
