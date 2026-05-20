# =============================================================================
# Feature Table Construction for Classification
# Purpose : Build classification input tables from the MF shared ASV table
#           (yonsei.mtg.csv) and from the total microbiome (phyloseq object).
#           Three feature set types are produced:
#             (1) MF-shared feature table  — abundance and binary formats
#             (2) Total feature table      — oral or fecal community at genus level
#             (3) T-MF feature table       — total features with MF ASVs removed
#
#   Feature tables are output in "Group-first" wide format:
#     rows = samples, columns = Group + taxon abundance/binary
#
# INPUT  : data/yonsei.mtg.csv   — shared ASV table (from 01_compute_MF_index.R)
#          data/yonsei_r15000.RDS — rarefied Yonsei phyloseq object
# OUTPUT : <prefix>.<site>.<rank>.mtg.csv      — MF-shared abundance tables
#          <prefix>.<site>.<rank>.mtg.binary.csv — MF-shared binary tables
#          ogc.yonsei.gen.core5.csv             — total oral genus feature table
#          ogc.core-mtg.rf.csv                  — T-MF feature table
# =============================================================================

library(phyloseq)
library(microbiome)
library(microViz)
library(tidyverse)

# =============================================================================
# Part 1. MF-shared feature tables
# =============================================================================

# ── Helper: build one abundance or binary RF input table ────────────────────
# Aggregates MF-shared ASV counts to a specified taxonomic rank,
# converts to relative abundance (%), and pivots to wide format.
make_rf_table <- function(df, tax_rank, count_col, sample_col,
                           depth  = 15000,
                           prefix = "crc.yonsei",
                           outdir = ".") {
  message("Building feature table: ", tax_rank, " / ", count_col)

  sub_df <- df %>%
    select(SampleID2, Group,
           all_of(sample_col), all_of(tax_rank), all_of(count_col)) %>%
    filter(!is.na(.data[[tax_rank]]), .data[[tax_rank]] != "")

  # Aggregate counts within each sample × taxon combination
  agg_df <- sub_df %>%
    group_by(SampleID2, Group, .data[[sample_col]], .data[[tax_rank]]) %>%
    summarise(Total = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(RelAbund = (Total / depth) * 100)

  # Pivot to sample × taxon wide format
  wide_df <- agg_df %>%
    select(SampleID2, Group, Taxon = all_of(tax_rank), RelAbund) %>%
    pivot_wider(names_from  = Taxon,
                values_from = RelAbund,
                values_fill = 0) %>%
    select(-SampleID2)  # remove SampleID2; Group is the first column

  site_label <- ifelse(count_col == "Oral_count", "oral", "feces")
  rank_label <- tolower(substr(tax_rank, 1, 3))  # spe/gen/fam
  out_file   <- file.path(outdir,
                           paste0(prefix, ".", site_label, ".", rank_label, ".mtg.csv"))
  write.csv(wide_df, out_file, row.names = FALSE)
  message("Saved: ", out_file)
  return(wide_df)
}

# ── Part 1A: HC vs CRC MF-shared feature tables ─────────────────────────────
mtg <- read.csv("data/yonsei.mtg.csv", check.names = FALSE) %>%
  filter(Group %in% c("HC", "CRC"))

tax_ranks  <- c("Genus")
count_info <- list(
  oral  = list(count_col = "Oral_count",  sample_col = "Oral_sample"),
  feces = list(count_col = "Fecal_count", sample_col = "Feces_sample")
)

rf_tables <- list()
for (tax_rank in tax_ranks) {
  for (site in names(count_info)) {
    obj_name <- paste(site, tax_rank, sep = "_")
    rf_tables[[obj_name]] <- make_rf_table(
      df        = mtg,
      tax_rank  = tax_rank,
      count_col = count_info[[site]]$count_col,
      sample_col = count_info[[site]]$sample_col,
      depth     = 15000,
      prefix    = "crc.yonsei",
      outdir    = "."
    )
  }
}

# ── Part 1B: GC stage-stratified feature tables (HC vs Early vs Late) ────────
# Cancer.stage coding:
#   HC    — healthy control
#   Early — Stage 0-2
#   Late  — Stage 3-4
mtg_gc <- read.csv("data/yonsei.mtg.csv", stringsAsFactors = FALSE) %>%
  filter(Group %in% c("HC", "GC")) %>%
  mutate(
    Stage_num    = suppressWarnings(as.numeric(Stage)),
    Cancer.stage = case_when(
      Group == "HC"                                        ~ "HC",
      Group == "GC" & !is.na(Stage_num) & Stage_num <= 2  ~ "Early",
      Group == "GC" & !is.na(Stage_num) & Stage_num >= 3  ~ "Late",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Cancer.stage))

cat("Stage distribution:\n")
print(table(mtg_gc$Cancer.stage, useNA = "ifany"))

# Helper: stage-stratified input table
make_stage_input_table <- function(data, rank_col,
                                    site         = c("mouth", "feces"),
                                    stage_keep   = c("HC", "Early", "Late"),
                                    depth        = 15000) {
  site      <- match.arg(site)
  count_col <- ifelse(site == "mouth", "Oral_count",  "Fecal_count")
  abund_col <- ifelse(site == "mouth", "Rel.abund.mouth", "Rel.abund.feces")

  df <- data %>%
    filter(Cancer.stage %in% stage_keep,
           !is.na(.data[[rank_col]]), .data[[rank_col]] != "")

  agg <- df %>%
    group_by(SampleID2, Cancer.stage, .data[[rank_col]]) %>%
    summarise(Total_count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(!!abund_col := (Total_count / depth) * 100)

  wide <- agg %>%
    select(SampleID2, Cancer.stage, !!sym(rank_col), !!sym(abund_col)) %>%
    pivot_wider(names_from  = !!sym(rank_col),
                values_from = !!sym(abund_col),
                values_fill = 0)
  return(wide)
}

# Generate all rank × site combinations for GC stage classification
tax_ranks_stage <- c("Species", "Genus", "Family")
sites_stage     <- c("mouth", "feces")

for (rank in tax_ranks_stage) {
  for (site in sites_stage) {
    out_df   <- make_stage_input_table(mtg_gc, rank_col = rank,
                                       site = site,
                                       stage_keep = c("HC", "Late"),
                                       depth = 15000)
    obj_name <- paste0("gc.hc-late.", tolower(substr(rank, 1, 3)), ".mf.", site)
    out_file <- paste0(obj_name, ".csv")
    write.csv(out_df, out_file, row.names = FALSE)
    cat("Saved:", out_file, "\n")
  }
}

# ── Part 1C: Binary conversion ───────────────────────────────────────────────
# Convert abundance feature tables to binary (presence/absence) format.
# All files matching "fit_late.csv" in the current directory are processed.
convert_to_binary <- function(file) {
  message("Converting to binary: ", file)
  df        <- read.csv(file, check.names = FALSE)
  group_col <- df$Group
  binary_df <- df %>%
    select(-Group) %>%
    mutate(across(everything(), ~ ifelse(. > 0, 1L, 0L)))
  binary_df <- cbind(Group = group_col, binary_df)
  out_file  <- gsub(".fit_late.csv", ".fit_late.binary.csv", file)
  write.csv(binary_df, out_file, row.names = FALSE)
}

files <- list.files(pattern = "fit_late\\.csv$")
lapply(files, convert_to_binary)

# =============================================================================
# Part 2. Total feature table (oral genus-level, core microbiome)
# =============================================================================
# Core filter: prevalence >= 5% across all oral samples in the HC + GC subset.
# Abundance: compositional (relative abundance, TSS-normalized).

d    <- readRDS("data/yonsei_r15000.RDS")
gc_d <- d %>%
  ps_filter(Group %in% c("HC", "GC"))

oral_d  <- ps_filter(gc_d, Type == "Mouth")
# fecal_d <- ps_filter(gc_d, Type == "Feces")  # uncomment for fecal total feature

# Filter to core taxa (>= 5% prevalence) and aggregate to genus level
oral_core <- core(oral_d, detection = 2, prevalence = 5 / 100) %>%
  aggregate_taxa("Genus") %>%
  microbiome::transform("compositional")

# Pivot to sample × genus wide format and append Group metadata
otu_wide <- otu_tibble(oral_core) %>%
  pivot_longer(cols = -FeatureID, names_to = "SampleID", values_to = "Value") %>%
  pivot_wider(names_from = FeatureID, values_from = Value)

meta_oral <- microbiome::meta(oral_core) %>%
  mutate(SampleID = rownames(.)) %>%
  select(SampleID, Group)

total_df <- left_join(otu_wide, meta_oral, by = "SampleID")
write.csv(total_df, "ogc.yonsei.gen.core5.csv", row.names = FALSE)
cat("Saved: ogc.yonsei.gen.core5.csv\n")

# =============================================================================
# Part 3. T-MF feature table
# Purpose : Total microbiome features with MF-transmitted ASVs removed,
#           to isolate the classification signal of non-MF bacteria.
# =============================================================================
# Load the list of MF-transmitted ASVs (at >= 25% prevalence threshold)
mtg_25 <- read.csv("data/yonsei.mtg.25.csv")
mtg_asvs <- unique(mtg_25$ASV)

# Remove MF ASV columns from the total feature table
total_tmf <- total_df[, !(names(total_df) %in% mtg_asvs)]
write.csv(total_tmf, "ogc.core-mtg.rf.csv", row.names = FALSE)
cat("Saved: ogc.core-mtg.rf.csv (T-MF feature table, MF ASVs removed)\n")
