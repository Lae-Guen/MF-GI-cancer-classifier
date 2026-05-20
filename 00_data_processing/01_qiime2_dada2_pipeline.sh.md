# QIIME2 / DADA2 Pipeline
# Purpose : Raw 16S rRNA amplicon read processing for all cohorts
#           (Yonsei internal + external validation cohorts)
# Region  : V3-V4 (primers 337F / 806R)
# Platform: Illumina MiSeq, 2 × 300 bp paired-end
# Software: QIIME2 v2021.02, DADA2 v1.18.0, Greengenes2 v2022.10

# ─────────────────────────────────────────────
# Step 1. Quality Control and Denoising (DADA2)
# ─────────────────────────────────────────────
# Truncation lengths were chosen based on quality score profiles
# (forward: retain positions 1-260; reverse: retain positions 1-230)
# Output files are named with the truncation lengths appended (e.g., 260230)

$ qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 1_demux.qza \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 230 \
  --o-table yonsei.gg2-filt-table260230.qza \
  --o-representative-sequences rep-seq.qza \
  --o-denoising-stats stats.qza

# ─────────────────────────────────────────────
# Step 2. Taxonomic Assignment (Greengenes2)
# ─────────────────────────────────────────────

# 2-1. Extract simulated amplicon reads from the Greengenes2 reference database
#      matching the V3-V4 primer set used in this study
$ qiime feature-classifier extract-reads \
  --i-sequences 2022.10.backbone.full-length.fna.qza \
  --p-f-primer CCTACGGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --o-reads gg2-V3V4-seq.qza \
  --p-n-jobs 4

# 2-2. Train a naive Bayes classifier on the extracted reference reads
$ qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads gg2-V3V4-seq.qza \
  --i-reference-taxonomy 2022.10.backbone.tax.qza \
  --o-classifier gg2_classifier_V3V4.qza

# 2-3. Classify representative sequences using the trained classifier
$ qiime feature-classifier classify-sklearn \
  --i-reads rep-seq.qza \
  --i-classifier gg2_classifier_V3V4.qza \
  --o-classification yonsei.taxonomy-gg2-260230.qza

# ─────────────────────────────────────────────
# Notes on external cohorts
# ─────────────────────────────────────────────
# Pre-processed QIIME2 artifacts (.qza) for external cohorts
# (Zhang, Uchino, Russo, Wang, Flemer, Zackular, Zeller) are provided
# in the data/ directory. These were processed using the same DADA2
# parameters and Greengenes2 classifier described above.
# Cohort-specific truncation lengths are encoded in the filename suffix
# (e.g., zhang.gg2-filt-table370.qza = trunc-len-f 370).
