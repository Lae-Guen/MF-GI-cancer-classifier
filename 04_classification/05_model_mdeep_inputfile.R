# =============================================================================
# MDeep Input File Preparation
# Purpose : Build the three input arrays required by MDeep for each
#           train/test iteration (50 seeds, seeds 100–149):
#             X_train.npy / X_eval.npy — ASV abundance matrix (samples × ASVs)
#             Y_train.npy / Y_eval.npy — binary outcome labels (0 = HC, 1 = disease)
#             c.npy                    — phylogenetic correlation matrix C = exp(-2ρD)
#
#   ρ (rho) = 2 was selected based on prior tuning experiments.
#   The phylogenetic distance matrix D is derived from the cophenetic
#   distance of the phylogenetic tree stored in the phyloseq object.
#   Train/test split: 70/30 stratified by outcome (createDataPartition).
#
# INPUT  : data/yonsei_r15000.RDS — rarefied Yonsei phyloseq object
#          data/yonsei.mtg.csv    — MF-transmitted ASV table
# OUTPUT : MDeep/data/<experiment>/iter_<seed>/{X_train, X_eval, Y_train,
#                                               Y_eval, c}.npy
# =============================================================================

library(phyloseq)
library(microbiome)
library(microViz)
library(tidyverse)
library(caret)    # createDataPartition
library(RcppCNPy) # npySave

d <- readRDS("data/yonsei_r15000.RDS")

# -----------------------------------------------------------------------------
# 1. Subset to MF-transmitted ASVs and target disease group
# -----------------------------------------------------------------------------
# GC example (change ps_filter Group for CRC or other comparisons)
asv_ids <- read.csv("data/yonsei.mtg.csv")$ASV %>% unique()
ps_mtg  <- prune_taxa(asv_ids, d)

ps_filtered <- ps_mtg %>%
  ps_filter(Type  == "Mouth") %>%
  ps_filter(Group %in% c("HC", "GC"))

# -----------------------------------------------------------------------------
# 2. Prevalence filter and relative abundance transformation
# -----------------------------------------------------------------------------
# Retain ASVs present in >= 10% of samples to reduce noise
ps_filtered <- ps_filtered %>%
  tax_filter(min_prevalence = 0.1, verbose = TRUE) %>%
  microbiome::transform("compositional")

# -----------------------------------------------------------------------------
# 3. Build feature matrix X (samples × ASVs)
# -----------------------------------------------------------------------------
X_data <- as(otu_table(ps_filtered), "matrix")
if (taxa_are_rows(ps_filtered)) X_data <- t(X_data)

# -----------------------------------------------------------------------------
# 4. Build phylogenetic correlation matrix C = exp(-2ρD)
# -----------------------------------------------------------------------------
# D: pairwise cophenetic (patristic) distance between ASVs
# ρ: evolutionary rate parameter (rho = 2 selected by prior cross-validation)
rho    <- 2
D_mat  <- cophenetic(phy_tree(ps_filtered))
D_mat  <- D_mat[colnames(X_data), colnames(X_data)]  # align to X columns
C_mat  <- exp(-2 * rho * D_mat)

# -----------------------------------------------------------------------------
# 5. Build outcome vector Y (0 = HC, 1 = GC)
# -----------------------------------------------------------------------------
ps_filtered <- ps_filtered %>%
  ps_mutate(
    Y_label = case_when(
      Group == "HC" ~ 0L,
      Group == "GC" ~ 1L,
      TRUE          ~ NA_integer_
    )
  )

y_data <- as.integer(sample_data(ps_filtered)$Y_label)

# -----------------------------------------------------------------------------
# 6. Save per-seed train/test splits (seeds 100–149, n = 50)
# -----------------------------------------------------------------------------
base_dir <- "MDeep/data/GC_mouth_mtg"

for (seed in 100:149) {
  set.seed(seed)
  train_idx <- createDataPartition(y_data, p = 0.7, list = FALSE)
  dir_path  <- file.path(base_dir, paste0("iter_", seed))
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  npySave(file.path(dir_path, "X_train.npy"), X_data[train_idx, ])
  npySave(file.path(dir_path, "Y_train.npy"), y_data[train_idx])
  npySave(file.path(dir_path, "X_eval.npy"),  X_data[-train_idx, ])
  npySave(file.path(dir_path, "Y_eval.npy"),  y_data[-train_idx])

  # C matrix is the same across all seeds (fixed ρ and tree)
  npySave(file.path(dir_path, "c.npy"), C_mat)
}

cat("MDeep input files saved to:", base_dir, "\n")
cat("Seeds:", 100, "–", 149, "| n_features:", ncol(X_data),
    "| n_samples:", nrow(X_data), "\n")
