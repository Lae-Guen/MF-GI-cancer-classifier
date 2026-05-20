# =============================================================================
# Alpha and Beta Diversity Analysis
# Purpose : Characterize oral and fecal microbial community structure across
#           disease groups (HC, MD, GC, CRC) and assess within-individual
#           oral-fecal diversity correlations (Figure 1B, S1B-E).
#
# Analyses:
#   Alpha diversity : Shannon index (microbiome::alpha)
#   Beta diversity  : Bray-Curtis PCoA (microeco::trans_beta),
#                     PERMANOVA (manova_all = TRUE)
#
# INPUT  : data/yonsei_r15000.RDS — rarefied Yonsei phyloseq object
# OUTPUT : yonsei.alpha.csv — per-sample Shannon diversity with metadata
#          PCoA plots and group distance box plots
# =============================================================================

library(phyloseq)
library(microbiome)
library(microeco)
library(file2meco)
library(tidyverse)
library(ggpubr)

d <- readRDS("data/yonsei_r15000.RDS")

group_palette <- c("HC"  = "#75ca9d",
                   "MD"  = "#ffeb33",
                   "GC"  = "#ca75dd",
                   "CRC" = "#eaa4a4")

# -----------------------------------------------------------------------------
# 1. Alpha diversity: Shannon index
# -----------------------------------------------------------------------------
meta  <- microbiome::meta(d)
alpha <- microbiome::alpha(d, index = "all")

shannon       <- alpha %>% select(diversity_shannon)
shannon$id    <- rownames(shannon)
meta$id       <- rownames(meta)
meta          <- left_join(meta, shannon, by = "id")

write.csv(meta, "yonsei.alpha.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Within-individual oral-fecal diversity correlation (Figure 1B, S1D)
# -----------------------------------------------------------------------------
# Reshape to wide format: one row per individual with oral and fecal Shannon values
meta_o <- meta %>%
  filter(Type == "Mouth") %>%
  rename(oshannon = diversity_shannon)

meta_f <- meta %>%
  filter(Type == "Feces") %>%
  rename(fshannon = diversity_shannon)

shannon_paired <- left_join(meta_o, meta_f, by = "SubjectID")

# Spearman correlation scatter plot, faceted by group (Figure 1B)
p_corr <- ggscatter(shannon_paired,
                    x = "fshannon", y = "oshannon",
                    color = "black", shape = 21, size = 2,
                    add = "reg.line",
                    add.params = list(color = "blue", fill = "lightgray"),
                    conf.int = TRUE) +
  stat_cor(method = "spearman", size = 6) +
  facet_wrap(~ factor(Group, levels = c("HC", "MD", "GC", "CRC"))) +
  theme(strip.text.x = element_text(size = 20))

ggsave("fig1B_oral_fecal_shannon_correlation.pdf", p_corr,
       width = 10, height = 8)

# -----------------------------------------------------------------------------
# 3. Alpha diversity comparison across groups (Figure S1B/C)
# -----------------------------------------------------------------------------
my_comparisons <- list(c("HC", "MD"), c("HC", "GC"), c("HC", "CRC"))

p_alpha <- ggboxplot(meta, x = "Group", y = "diversity_shannon",
                     fill = "Group",
                     palette = group_palette,
                     order = c("HC", "MD", "GC", "CRC")) +
  stat_compare_means(comparisons = my_comparisons, size = 4) +
  facet_wrap(~ factor(Type, levels = c("Mouth", "Feces"))) +
  theme(
    legend.position  = "none",
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 14),
    axis.title.y     = element_text(size = 14)
  ) +
  labs(y = "Shannon diversity index")

ggsave("figS1BC_shannon_by_group.pdf", p_alpha, width = 8, height = 5)

# -----------------------------------------------------------------------------
# 4. Beta diversity: Bray-Curtis PCoA via microeco (Figure 1D/E, S1E)
# -----------------------------------------------------------------------------
# Analysis is run separately for oral and fecal samples.
# Family-level taxonomy is used for ordination (consistent with paper methods).

run_pcoa <- function(ps, sample_type) {
  cat("\n--- Beta diversity:", sample_type, "---\n")

  # Convert to microtable and aggregate to Family level
  meco_obj <- phyloseq2meco(ps)
  meco_obj$cal_abund()
  meco_fam <- meco_obj$merge_taxa(taxa = "Family")
  meco_fam$tax_table <- meco_fam$tax_table[meco_fam$tax_table$Family != "f__", ]
  meco_fam$tidy_dataset()
  rownames(meco_fam$otu_table)  <- meco_fam$tax_table[rownames(meco_fam$otu_table), "Family"]
  rownames(meco_fam$tax_table)  <- meco_fam$tax_table[, "Family"]

  # Compute Bray-Curtis distance and run PCoA
  meco_fam$cal_betadiv(unifrac = FALSE)
  t1 <- trans_beta$new(dataset = meco_fam, group = "Group", measure = "bray")

  set.seed(111)
  t1$cal_ordination(ordination = "PCoA")

  set.seed(111)
  t1$cal_manova(manova_all = TRUE)
  cat("PERMANOVA result:\n")
  print(t1$res_manova)

  # PCoA scatter plot with group ellipses
  p_pcoa <- t1$plot_ordination(
    plot_color       = "Group",
    plot_type        = c("point", "ellipse"),
    plot_group_order = c("HC", "MD", "GC", "CRC"),
    ellipse_chull_alpha = 0.05,
    point_alpha      = 0.7
  ) +
    theme_bw() +
    scale_color_manual(values = group_palette) +
    theme(
      axis.title.x  = element_text(size = 20),
      axis.title.y  = element_text(size = 20),
      legend.title  = element_text(size = 18),
      legend.text   = element_text(size = 16)
    )

  ggsave(paste0("fig1_PCoA_", sample_type, "_BrayCurtis.pdf"),
         p_pcoa, width = 8, height = 6)

  # Within-group distance box plot and Wilcoxon test
  t1$cal_group_distance(within_group = TRUE)
  t1$cal_group_distance_diff(method = "wilcox")

  # PCoA axis scores box plot for group comparison
  pcoa_scores <- t1$res_ordination$scores
  my_comp_pcoa <- list(c("HC", "GC"), c("GC", "CRC"))

  p_pco1 <- ggboxplot(pcoa_scores, x = "Group", y = "PCo1",
                      fill = "Group",
                      palette = group_palette,
                      order = c("HC", "MD", "GC", "CRC")) +
    stat_compare_means(
      comparisons = my_comp_pcoa,
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05),
                         symbols   = c("****", "***", "**", "*")),
      size = 4
    ) +
    theme(
      legend.position = "none",
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(size = 18),
      axis.title.y    = element_text(size = 18)
    ) +
    labs(y = "Bray-Curtis PCo1")

  ggsave(paste0("fig1_PCo1_", sample_type, "_boxplot.pdf"),
         p_pco1, width = 5, height = 5)

  return(t1)
}

# Run for oral samples (Figure 1D) and fecal samples (Figure 1E)
oral_ps  <- ps_filter(d, Type == "Mouth")
feces_ps <- ps_filter(d, Type == "Feces")

pcoa_oral  <- run_pcoa(oral_ps,  sample_type = "oral")
pcoa_feces <- run_pcoa(feces_ps, sample_type = "fecal")

cat("Beta diversity analysis complete.\n")
