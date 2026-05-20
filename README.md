# MF-GI-Cancer-Classifier

**Oral microbiome signatures of gut transmission enable robust, non-invasive classification of gastrointestinal cancers**

Jang LG, Huh JW, et al. — *Cell Host & Microbe* (under review)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20302778.svg)](https://doi.org/10.5281/zenodo.20302778)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R version](https://img.shields.io/badge/R-4.1.0-blue.svg)](https://www.r-project.org/)

---

## Overview

This repository contains all analysis code for the manuscript. We propose using the **MF (Mouth-to-Feces) index**—a measure of oral-to-gut microbial transmission—as a non-invasive feature set for classifying gastric cancer (GC) and colorectal cancer (CRC). The pipeline covers:

1. **MF index computation** — ASV-level co-occurrence between matched oral and fecal samples
2. **Community analysis** — alpha/beta diversity, LEfSe, MaAsLin2
3. **Mouth–Tumor–Gut (MTG) axis** — transmission pathway decomposition in the Zhang cohort
4. **Classification** — internal (Yonsei) and external (7 public cohorts) validation using Random Forest, glmnet, and MDeep

---

## Repository Structure

```
MF-GI-cancer-classifier/
│
├── README.md
├── LICENSE
├── sessionInfo.txt                    # R session info for reproducibility
│
├── 00_data_processing/
│   ├── 01_qiime2_dada2_pipeline.sh.md   # QIIME2 DADA2 processing pipeline
│   ├── 02_setup.R                        # Package loading (run first)
│   └── 03_phyloseq_import_rarefaction.R  # Import .qza → phyloseq, rarefy to 15,000
│
├── 01_MF_index/
│   ├── 01_compute_MF_index.R             # Compute MF index per individual
│   ├── 02_MF_index_statistics.R          # Wilcoxon tests, effect sizes, visualization
│   ├── 03_MF_index_multivariate_analysis.R  # Multivariate regression with lifestyle covariates
│   ├── 04_MF_species_profiling.R         # Species-level MF index and oral prevalence
│   └── 05_MF_stage_stratified_analysis.R # Stage-stratified (Early/Late) analysis
│
├── 02_community_analysis/
│   ├── 01_diversity_ordination.R         # Alpha/beta diversity, PCoA, PERMANOVA
│   └── 02_differential_abundance_lefse_maaslin2.R  # LEfSe + MaAsLin2 consensus
│
├── 03_MTG_tumor_axis/
│   ├── 01_zhang_transmission_index.R     # MTT/MTP/TTG/PTG/MTG indices (Zhang cohort)
│   └── 02_zhang_stage_stratified_analysis.R  # Stage comparison, pathway ratios
│
└── 04_classification/
    ├── 01_feature_construction_prevalence.R      # Build MF / Total / T-MF feature tables
    ├── 02_MF_feature_extraction_using_VSEARCH.md # VSEARCH pipeline (bash)
    ├── 03_sequence_similarity_based_mf_feature.R # VSEARCH hit evaluation + ASV extraction
    ├── 04_rename_collapse_asvs.R                 # Rename external ASVs to Yonsei query IDs
    ├── 05_model_mdeep_inputfile.R                # Build MDeep input arrays (.npy)
    ├── 06_run_model_Mdeep.md                     # MDeep training and evaluation pipeline
    ├── 07_merging_external_dataset.R             # Merge Yonsei + external feature tables
    ├── 08_functions_external_test.R              # External validation functions (multi-algorithm)
    ├── 09_run_scenario_comparison.R              # Scenario × algorithm benchmarking
    ├── 10_run_internal_rf_optimized.R            # Internal validation (RF, 100 seeds)
    ├── 11_functions_external_rf_optimized.R      # External RF validation functions
    ├── 12_run_external_rf_optimized.R            # External RF validation runner
    ├── 13_auc_stats_effectsize_ci.R              # AUC comparison statistics
    └── 14_bootstrap_for_auc_calculation_of_FOBT.R  # FOBT benchmark
```

---

## Data Availability

Raw sequencing data are deposited at NCBI Sequence Read Archive (SRA) and DDBJ.

### Internal (Yonsei) cohort

| Cohort | Sample type | Disease groups | n (subjects) | Accession |
|--------|------------|----------------|--------------|-----------|
| Yonsei | Oral + Fecal (paired) | HC, MD, GC, CRC | 507 | PRJNA1248208 |

### External validation cohorts

| Cohort | Sample type | Disease | Region | Accession |
|--------|------------|---------|--------|-----------|
| Zhang  | Oral + Fecal + Tumor + Para (paired) | GC, CRC | V4 | PRJNA778008 |
| Uchino | Oral + Fecal (paired) | CRC | V1-V2 | PRJDB11845 |
| Russo  | Oral + Fecal (paired) | CRC | V3-V4 | ENA (Russo et al. 2020) |
| Wang   | Oral + Fecal (unpaired) | CRC | V3-V4 | PRJNA698923 |
| Flemer | Oral (unpaired) | CRC | V3-V4 | PRJEB50080 |
| Zackular | Fecal (unpaired) | CRC | V3-V4 | PRJNA389927 |
| Zeller | Fecal (unpaired) | CRC | V3-V4 | PRJEB6070 |

Pre-processed QIIME2 artifacts (.qza) for all cohorts are provided via Zenodo:
**[https://doi.org/10.5281/zenodo.20302778](https://doi.org/10.5281/zenodo.20302778)**

---

## Requirements

### R (version 4.1.0)

Install all required packages by running:

```r
source("00_data_processing/02_setup.R")
```

| Package | Version | Purpose |
|---------|---------|---------|
| phyloseq | 1.42.0 | Microbiome data structure |
| microbiome | 1.23.1 | Alpha/beta diversity utilities |
| microViz | — | phyloseq helpers (ps_filter, otu_tibble) |
| qiime2R | — | Import QIIME2 .qza into R |
| microeco | 1.6.0 | LEfSe, PCoA, PERMANOVA |
| file2meco | — | phyloseq → microtable conversion |
| tidyverse | 2.0.0 | Data wrangling and visualization |
| ggpubr | 0.6.2 | Publication-ready plots |
| Maaslin2 | 1.12.0 | Multivariable association testing |
| mikropml | 1.7.0 | Machine learning pipeline |
| pROC | — | ROC and AUC computation |
| future / future.apply | — | Parallel execution |
| rstatix | 0.7.3 | Statistical testing helpers |
| emmeans | — | Estimated marginal means |
| car | — | Type III ANOVA, VIF |
| RcppCNPy | — | Save .npy arrays for MDeep |

> **Important:** `mikropml` v1.7.0 is required. Earlier versions (≤ 1.6.1) have a
> different argument signature for `get_partition_indices()`. See `sessionInfo.txt`
> for the full version list used in this study.

### External tools

| Tool | Version | Purpose |
|------|---------|---------|
| QIIME2 | 2021.02 | Raw read processing (Step 00) |
| VSEARCH | 2.x | ASV sequence identity matching (Step 04) |

### Python (MDeep benchmarking only)

```bash
mamba create -n mdeep -c anaconda -c conda-forge \
  python=3.6 tensorflow=1.12.0 numpy=1.16 scipy=1.2 \
  pandas=0.24 scikit-learn=0.20 matplotlib=3.0 h5py=2.9 -y
conda activate mdeep
git clone https://github.com/alfredyewang/MDeep.git
```

---

## Usage

Scripts are numbered and should be run in order within each module.

### Step 0 — Raw read processing (QIIME2 + DADA2)

```bash
# See 00_data_processing/01_qiime2_dada2_pipeline.sh.md for full commands
# Truncation lengths: forward 260 bp, reverse 230 bp (V3-V4 region)
```

```r
source("00_data_processing/02_setup.R")
source("00_data_processing/03_phyloseq_import_rarefaction.R")
# Rarefied to 15,000 reads/sample (minimum depth across all samples)
# Output: data/*_r15000.RDS for each cohort
```

### Step 1 — MF index computation

```r
source("01_MF_index/01_compute_MF_index.R")
# Identifies ASVs with 100% 16S identity in matched oral-fecal pairs
# Output: data/yonsei.mtg.csv (shared ASV table)

source("01_MF_index/02_MF_index_statistics.R")
source("01_MF_index/03_MF_index_multivariate_analysis.R")
source("01_MF_index/04_MF_species_profiling.R")
source("01_MF_index/05_MF_stage_stratified_analysis.R")
```

### Step 2 — Community analysis

```r
source("02_community_analysis/01_diversity_ordination.R")
source("02_community_analysis/02_differential_abundance_lefse_maaslin2.R")
```

### Step 3 — MTG tumor axis (Zhang cohort)

```r
source("03_MTG_tumor_axis/01_zhang_transmission_index.R")
source("03_MTG_tumor_axis/02_zhang_stage_stratified_analysis.R")
```

### Step 4 — Classification

```r
# 4-1. Build feature tables
source("04_classification/01_feature_construction_prevalence.R")

# 4-2. VSEARCH — see 02_MF_feature_extraction_using_VSEARCH.md (bash)

# 4-3. Process VSEARCH results and build external feature tables
source("04_classification/03_sequence_similarity_based_mf_feature.R")
source("04_classification/04_rename_collapse_asvs.R")

# 4-4. Merge Yonsei + external datasets
source("04_classification/07_merging_external_dataset.R")

# 4-5. Internal validation (Random Forest, 100 seeds)
source("04_classification/10_run_internal_rf_optimized.R")

# 4-6. External validation (RF, fixed group-partitioned split)
source("04_classification/12_run_external_rf_optimized.R")

# 4-7. FOBT benchmark
source("04_classification/14_bootstrap_for_auc_calculation_of_FOBT.R")

# 4-8. AUC comparison statistics
source("04_classification/13_auc_stats_effectsize_ci.R")
```

---

## Key Methods

### MF index
Shared ASVs are identified by 100% 16S rRNA sequence identity between matched oral and fecal samples of the same individual. The MF index (%) is:

```
MF index = (sum of shared fecal ASV counts / rarefaction depth) × 100
```

Samples with no shared ASVs detected (paired but MF bacteria absent) receive an MF index of 0. Unpaired samples are excluded.

### Random Forest classifier
- Package: `mikropml` v1.7.0 (`method = "rf"`, uses `randomForest` package)
- `ntree` = 500; `mtry` tuned dynamically based on feature count (√p, √p × 1.5, √p × 2)
- Internal validation: 100-seed repeated 70/30 stratified split
- External validation: fixed group-partitioned split (`groups` / `group_partitions` in `run_ml()`)
- Patient-level AUC: mean predicted probability across 100 seeds
- Classification threshold: Youden's J index on patient-level ROC

### External validation feature construction
MF-transmitted ASV sequences from Yonsei are matched to external cohort representative sequences using VSEARCH (`--usearch_global`, `--iddef 2`, `--strand both`). External ASV IDs matching at 100% identity are extracted, renamed to Yonsei query ASV IDs, and abundance values summed for multi-hit cases. Cohort-wise z-score normalization is applied to abundance features before classification.

### VSEARCH parameters
- Identity range tested: 70%–100%
- External cohorts (V3-V4): Russo, Flemer, Wang
- External cohorts (V4 core region): Zhang, Baxter, Zackular, Zeller

---

## Reproducibility

A fixed random seed is applied at each stage:
- Rarefaction: `rngseed = 1`
- Classification: seeds 100–199 (internal); seeds 100–149 (external/MDeep)
- Effect size bootstrap: `set.seed(42)`

For the exact R package versions used in this study, see [`sessionInfo.txt`](sessionInfo.txt).

```r
# To regenerate sessionInfo.txt:
source("00_data_processing/02_setup.R")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

Server environment: Linux (Ubuntu 20.04), parallel execution via `future::multicore`
with ~20 workers. Set `OMP_NUM_THREADS=1` to prevent nested parallelism conflicts.

---

## Citation

If you use this code or data, please cite:

> Jang LG, Huh JW, et al. Oral microbiome signatures of gut transmission enable
> robust, non-invasive classification of gastrointestinal cancers.
> *Cell Host & Microbe* (under review).

Code and data archived at Zenodo:
> [https://doi.org/10.5281/zenodo.20302778](https://doi.org/10.5281/zenodo.20302778)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

- **Jihyun F. Kim** (corresponding author) — jfk1@yonsei.ac.kr
- Department of Systems Biology, Yonsei University, Seoul, Republic of Korea
