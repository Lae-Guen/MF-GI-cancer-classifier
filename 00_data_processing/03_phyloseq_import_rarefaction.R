# =============================================================================
# Phyloseq Import and Rarefaction
# Purpose : Import QIIME2 .qza artifacts for all cohorts into phyloseq objects
#           and rarefy each to a uniform sequencing depth of 15,000 reads.
#
# INPUT  : QIIME2 feature table (.qza), taxonomy (.qza), and metadata (.txt)
#          for each cohort (located in the data/ directory)
# OUTPUT : Rarefied phyloseq objects (one per cohort), saved as .RDS files
#
# Rarefaction depth : 15,000 reads/sample (minimum depth across all samples)
# Replacement       : TRUE (standard practice for rarefaction)
# Random seed       : 1 (fixed for reproducibility)
# =============================================================================

source("00_data_processing/02_setup.R")

# -----------------------------------------------------------------------------
# Helper function: import a QIIME2 artifact set and rarefy
# -----------------------------------------------------------------------------
import_and_rarefy <- function(features, taxonomy, metadata,
                               depth = 15000, seed = 1) {
  ps <- qiime2R::qza_to_phyloseq(
    features = features,
    taxonomy = taxonomy,
    metadata = metadata
  )
  ps_rare <- ps %>%
    rarefy_even_depth(sample.size = depth, rngseed = seed, replace = TRUE)
  return(ps_rare)
}

# -----------------------------------------------------------------------------
# Yonsei cohort (internal; paired oral-fecal, n = 507)
# -----------------------------------------------------------------------------
yonsei <- import_and_rarefy(
  features = "data/yonsei.gg2-filt-table260230.qza",
  taxonomy = "data/yonsei.taxonomy-gg2-260230.qza",
  metadata = "data/yonsei.meta.txt"
)
saveRDS(yonsei, "data/yonsei_r15000.RDS")

# -----------------------------------------------------------------------------
# Zhang cohort (external; paired oral-fecal, GC + CRC)
# -----------------------------------------------------------------------------
zhang <- import_and_rarefy(
  features = "data/zhang.gg2-filt-table370.qza",
  taxonomy = "data/zhang.gg2-taxonomy370.qza",
  metadata = "data/zhang.gc-crc-meta.txt"
)
saveRDS(zhang, "data/zhang_r15000.RDS")

# -----------------------------------------------------------------------------
# Uchino cohort (external; paired oral-fecal, CRC; V1-V2 region)
# NOTE: ASV-based VSEARCH matching is not performed for this cohort
#       due to incompatible sequencing region (V1-V2 vs V3-V4 in Yonsei)
# -----------------------------------------------------------------------------
uchino <- import_and_rarefy(
  features = "data/uchino.dada-filt-table230210.qza",
  taxonomy = "data/uchino.taxonomy-gg230210.qza",
  metadata = "data/uchino.CRCmetadata.txt"
)
saveRDS(uchino, "data/uchino_r15000.RDS")

# -----------------------------------------------------------------------------
# Russo cohort (external; paired oral-fecal, CRC; saliva)
# -----------------------------------------------------------------------------
russo <- import_and_rarefy(
  features = "data/Russo-gg2-table270230.qza",
  taxonomy = "data/Russo-gg2-taxonomy270230.qza",
  metadata = "data/Russo-gg2-meta.txt"
)
saveRDS(russo, "data/russo_r15000.RDS")

# -----------------------------------------------------------------------------
# Wang cohort (external; oral and fecal unpaired, CRC; saliva)
# -----------------------------------------------------------------------------
wang <- import_and_rarefy(
  features = "data/Wang-gg2-table270250.qza",
  taxonomy = "data/Wang-gg2-taxonomy270250.qza",
  metadata = "data/Wang-gg2-meta.txt"
)
saveRDS(wang, "data/wang_r15000.RDS")

# -----------------------------------------------------------------------------
# Flemer cohort (external; oral unpaired, CRC; buccal swab)
# -----------------------------------------------------------------------------
flemer <- import_and_rarefy(
  features = "data/Flemer-gg2-table.qza",
  taxonomy = "data/Flemer-gg2-taxonomy.qza",
  metadata = "data/Flemer.meta.txt"
)
saveRDS(flemer, "data/flemer_r15000.RDS")

# -----------------------------------------------------------------------------
# Zackular cohort (external; fecal unpaired, CRC)
# -----------------------------------------------------------------------------
zackular <- import_and_rarefy(
  features = "data/Zackular-gg2-table.220220.qza",
  taxonomy = "data/Zackular-gg2-taxonomy220220.qza",
  metadata = "data/Zackular-gg2-meta.txt"
)
saveRDS(zackular, "data/zackular_r15000.RDS")

# -----------------------------------------------------------------------------
# Zeller cohort (external; fecal unpaired, CRC)
# -----------------------------------------------------------------------------
zeller <- import_and_rarefy(
  features = "data/Zeller-gg2-table.220210.qza",
  taxonomy = "data/Zeller-gg2-taxonomy220210.qza",
  metadata = "data/Zeller.meta.txt"
)
saveRDS(zeller, "data/zeller_r15000.RDS")
