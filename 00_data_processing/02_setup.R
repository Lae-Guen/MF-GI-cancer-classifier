# =============================================================================
# Setup: Load Required R Packages
# Purpose : Install (if necessary) and load all packages used across the
#           analysis pipeline. Run this script before any other script.
# R version: 4.1.0
# =============================================================================

# Set working directory to the repository root.
# Update this path to match your local environment.
setwd("path/to/MF-GI-cancer-classifier")

# -----------------------------------------------------------------------------
# Core packages
# -----------------------------------------------------------------------------
library(phyloseq)    # v1.42.0 — phylogenetic microbiome data structures
library(microbiome)  # v1.23.1 — alpha/beta diversity, microbiome utilities
library(microViz)    # — phyloseq visualization helpers (ps_filter, otu_tibble)
library(qiime2R)     # — import QIIME2 .qza artifacts into phyloseq
library(microeco)    # v1.6.0 — beta diversity, LEfSe, ordination via trans_* classes
library(file2meco)   # — phyloseq-to-microtable conversion (phyloseq2meco)
library(tidyverse)   # v2.0.0 — data wrangling and visualization (dplyr, ggplot2, etc.)
library(magrittr)    # — pipe operator (%>%, %<>%)
library(ggpubr)      # v0.6.2 — publication-ready ggplot2 extensions
library(Maaslin2)    # v1.12.0 — multivariable association testing (MaAsLin2)
library(mikropml)    # v1.7.0 — machine learning pipeline for microbiome classification
