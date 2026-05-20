# =============================================================================
# Rename and Collapse External ASVs to Yonsei Query ASV IDs
# Purpose : The VSEARCH pipeline produces a query-to-target mapping file
#           (query_to_<cohort>_<pid>.tsv) that links Yonsei MF-transmitted
#           query ASV IDs to matched external cohort ASV IDs.
#           This script renames external ASV columns to their corresponding
#           query IDs and collapses (sums) multiple external ASVs that map
#           to the same query ASV, enabling cross-cohort feature matrix merging.
#
# Processing logic:
#   - External ASVs with a mapping entry  → renamed to query ASV ID, then
#     summed if multiple external ASVs share the same query ID
#   - External ASVs without a mapping entry → retained under their original ID
#
# INPUT  : fcrc_wang_1.00_abundance.csv  — abundance feature table (external)
#          fcrc_wang_1.00_binary.csv     — binary feature table (external)
#          query_to_wang_1.00.tsv        — query-to-target mapping (col1=query, col2=external)
# OUTPUT : fcrc_wang_1.00_abundance_renamed.csv
#          fcrc_wang_1.00_binary_renamed.csv
#          merge_summary.csv             — number of external ASVs collapsed per query ASV
# =============================================================================

library(dplyr)
library(tidyr)

# =============================================================================
# Helper: rename and collapse one feature table (abundance or binary)
# =============================================================================
rename_and_collapse <- function(abundance_file, mapping_file,
                                 meta_cols  = c("SampleID", "Group"),
                                 is_binary  = FALSE,
                                 out_file   = NULL) {

  df      <- read.csv(abundance_file, check.names = FALSE)
  mapping <- read.table(mapping_file, sep = "\t", header = FALSE,
                         col.names = c("query_asv", "external_asv"),
                         stringsAsFactors = FALSE)

  # Separate metadata columns from ASV feature columns
  asv_cols  <- setdiff(colnames(df), meta_cols)
  meta_df   <- df[, meta_cols,  drop = FALSE]
  asv_df    <- df[, asv_cols,   drop = FALSE]

  # Identify external ASVs that have a mapping entry
  mapping_present  <- mapping[mapping$external_asv %in% asv_cols, ]
  external_in_map  <- unique(mapping_present$external_asv)
  external_not_map <- setdiff(asv_cols, external_in_map)

  cat("External ASVs in abundance table AND mapping:", length(external_in_map), "\n")
  cat("External ASVs NOT in mapping (kept as-is):  ", length(external_not_map), "\n")

  # ── Rename mapped external ASVs to query ASV IDs ────────────────────────
  asv_mapped   <- asv_df[, external_in_map, drop = FALSE]
  rename_vec   <- setNames(mapping_present$query_asv, mapping_present$external_asv)
  colnames(asv_mapped) <- rename_vec[colnames(asv_mapped)]

  # ── Collapse (sum/OR) columns sharing the same query ASV ID ─────────────
  query_names    <- colnames(asv_mapped)
  unique_queries <- unique(query_names)

  cat("Unique query ASVs after rename:", length(unique_queries), "\n")

  collapsed_list <- lapply(unique_queries, function(q) {
    cols <- which(query_names == q)
    if (length(cols) == 1) {
      if (is_binary) {
        df_out <- as.data.frame(as.integer(asv_mapped[, cols] >= 1))
      } else {
        df_out <- asv_mapped[, cols, drop = FALSE]
      }
      colnames(df_out) <- q
      return(df_out)
    } else {
      # Multiple external ASVs → sum (abundance) or binary OR
      mat    <- asv_mapped[, cols, drop = FALSE]
      summed <- rowSums(mat, na.rm = TRUE)
      if (is_binary) summed <- as.integer(summed >= 1)
      df_out <- data.frame(summed)
      colnames(df_out) <- q
      return(df_out)
    }
  })

  collapsed_df <- do.call(cbind, collapsed_list)
  rownames(collapsed_df) <- NULL

  cat("Collapsed feature matrix:", nrow(collapsed_df), "samples ×",
      ncol(collapsed_df), "query ASVs\n")

  # ── Unmapped external ASVs (retain original IDs) ─────────────────────────
  asv_unmapped <- asv_df[, external_not_map, drop = FALSE]

  # ── Final output: meta + renamed/collapsed + unmapped ────────────────────
  final_df <- cbind(meta_df, collapsed_df, asv_unmapped)

  if (!is.null(out_file)) {
    write.csv(final_df, out_file, row.names = FALSE, quote = TRUE)
    cat("Saved:", out_file, "\n")
  }

  return(final_df)
}

# =============================================================================
# Wang CRC cohort — 100% identity
# =============================================================================

# Abundance
rename_and_collapse(
  abundance_file = "fcrc_wang_1.00_abundance.csv",
  mapping_file   = "query_to_wang_1.00.tsv",
  meta_cols      = c("SampleID", "Group"),
  is_binary      = FALSE,
  out_file       = "fcrc_wang_1.00_abundance_renamed.csv"
)

# Binary
rename_and_collapse(
  abundance_file = "fcrc_wang_1.00_binary.csv",
  mapping_file   = "query_to_wang_1.00.tsv",
  meta_cols      = c("SampleID", "Group"),
  is_binary      = TRUE,
  out_file       = "fcrc_wang_1.00_binary_renamed.csv"
)

# =============================================================================
# Collapse diagnostic: how many external ASVs merged per query ASV
# =============================================================================
mapping_wang <- read.table("query_to_wang_1.00.tsv", sep = "\t", header = FALSE,
                            col.names = c("query_asv", "external_asv"))
abundance_wang <- read.csv("fcrc_wang_1.00_abundance.csv", check.names = FALSE)
asv_cols_wang  <- setdiff(colnames(abundance_wang), c("SampleID", "Group"))

mapping_present <- mapping_wang[mapping_wang$external_asv %in% asv_cols_wang, ]

merge_summary <- mapping_present %>%
  group_by(query_asv) %>%
  summarise(
    n_external_asvs_collapsed = n(),
    external_asvs             = paste(external_asv, collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_external_asvs_collapsed))

write.csv(merge_summary, "merge_summary.csv", row.names = FALSE)
cat("Saved: merge_summary.csv\n")
cat("\nTop 10 query ASVs by number of collapsed external ASVs:\n")
print(head(merge_summary[, c("query_asv", "n_external_asvs_collapsed")], 10))
