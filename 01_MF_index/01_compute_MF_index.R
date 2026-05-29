# =============================================================================
# Compute MF (Mouth-to-Feces) Transmission Index
# Purpose : Identify ASVs shared between matched oral and fecal samples
#           (100% 16S rRNA sequence identity within each individual) and
#           compute the MF index as the proportion of shared ASV abundance
#           relative to total fecal ASV abundance.
#
# MF index formula:
#   MF index (%) = (sum of shared ASV counts in feces / total fecal reads) × 100
#
# INPUT  : data/yonsei_r15000.RDS  — rarefied Yonsei phyloseq object
#          data/yonsei.meta.txt    — sample metadata with SampleID, Type, Group
# OUTPUT : yonsei.meta_checked.txt — metadata with SubjectID and pairing status
#          yonsei_mtg15000_<SubjectID>_results.txt — per-subject shared ASV tables
#          yonsei.mtg.csv          — combined shared ASV table across all subjects
# =============================================================================

library(phyloseq)
library(dplyr)
library(readr)   # parse_number()

# Load rarefied phyloseq object
d <- readRDS("data/yonsei_r15000.RDS")

# -----------------------------------------------------------------------------
# Step 1. Assign SubjectIDs and detect paired samples
# -----------------------------------------------------------------------------
# SampleID format: YYMMDD_<individual_code>[S|F]
#   S = oral (Mouth), F = fecal (Feces)
# PersonKey: SampleID with the trailing S/F removed → identifies one individual

meta <- read.table("./yonsei.meta.txt",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

meta <- meta %>%
  mutate(
    PersonKey = sub("[SF]$", "", SampleID),          # strip trailing S/F
    DateKey   = substr(PersonKey, 1, 6),             # 6-digit date prefix (YYMMDD)
    WithinKey = sub("^\\d{6}_", "", PersonKey),      # individual code after date
    WithinNum = parse_number(WithinKey)              # numeric part for sorting
  )

# Build a unique SubjectID (Yonsei1, Yonsei2, ...) per individual
id_map <- meta %>%
  distinct(PersonKey, DateKey, WithinNum, WithinKey) %>%
  arrange(as.integer(DateKey), WithinNum, WithinKey, PersonKey) %>%
  mutate(SubjectID = paste0("Yonsei", row_number())) %>%
  select(PersonKey, SubjectID)

meta <- meta %>%
  left_join(id_map, by = "PersonKey")

# -----------------------------------------------------------------------------
# Step 2. Classify each subject as paired / mouth_only / feces_only
# -----------------------------------------------------------------------------
pair_status <- meta %>%
  group_by(SubjectID, PersonKey) %>%
  summarise(
    n_samples = n(),
    has_mouth = any(Type == "Mouth"),
    has_feces = any(Type == "Feces"),
    pairing   = case_when(
      has_mouth & has_feces  ~ "paired",
      has_mouth & !has_feces ~ "mouth_only",
      !has_mouth & has_feces ~ "feces_only",
      TRUE                   ~ "unknown"
    ),
    .groups = "drop"
  )

# Retain only paired subjects for MF index computation
meta_paired <- meta %>%
  left_join(pair_status %>% select(SubjectID, pairing), by = "SubjectID") %>%
  filter(pairing == "paired")

# Summary of pairing status across subjects
table(pair_status$pairing)

write.table(meta, file = "data/yonsei.meta_checked.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# -----------------------------------------------------------------------------
# Step 3. Compute shared ASV table per subject
# -----------------------------------------------------------------------------
# An ASV is considered "shared" if it has ≥ 1 read in BOTH the oral
# and the matched fecal sample of the same individual.
# The MF index is then the ratio of shared fecal counts to total fecal reads.

metadata <- read.table("data/yonsei.meta_checked.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(metadata) <- metadata$SampleID
sample_data(d) <- sample_data(metadata)

unique_subject_ids <- unique(metadata$SubjectID)
results_list <- list()

for (subject_id in unique_subject_ids) {
  tryCatch({

    # Retrieve metadata rows for this subject
    subject_rows <- metadata[metadata$SubjectID == subject_id, ]
    if (nrow(subject_rows) == 0) {
      cat("No metadata found for SubjectID:", subject_id, "\n")
      next
    }

    subject_meta   <- subject_rows[1, ]
    sex_info       <- subject_meta$Sex
    dob_info       <- subject_meta$DOB
    group_info     <- subject_meta$Group

    # Subset phyloseq object to this subject's samples
    subset_ps <- subset_samples(d, SubjectID == subject_id)

    # Identify ASVs present in both oral and fecal samples (shared ASVs)
    asv_table     <- as.data.frame(otu_table(subset_ps))
    shared_asvs   <- asv_table[rowSums(asv_table >= 1) == 2, ]  # present in both samples
    taxa_info     <- as.data.frame(tax_table(subset_ps))

    # Identify the oral and fecal sample IDs for this subject
    oral_id  <- as.character(
      metadata[metadata$SubjectID == subject_id & metadata$Type == "Mouth", "SampleID"])
    fecal_id <- as.character(
      metadata[metadata$SubjectID == subject_id & metadata$Type == "Feces", "SampleID"])

    if (length(oral_id) == 0 || length(fecal_id) == 0) {
      cat("Oral or fecal sample missing for SubjectID:", subject_id, "\n")
      next
    }

    # Compile the shared ASV table with taxonomy and metadata
    if (nrow(shared_asvs) >= 1) {
      result_df <- cbind(
        ASV         = rownames(shared_asvs),
        shared_asvs[, c(oral_id, fecal_id)],
        taxa_info[rownames(shared_asvs), ]
      )
      result_df$SubjectID    <- subject_id
      result_df$Oral_sample  <- oral_id
      result_df$Feces_sample <- fecal_id
      result_df$Sex          <- sex_info
      result_df$DOB          <- dob_info
      result_df$Group        <- group_info

      # Rename abundance columns for clarity
      colnames(result_df)[2] <- "Oral_count"
      colnames(result_df)[3] <- "Fecal_count"

      results_list[[subject_id]] <- result_df

      # Save per-subject result
      out_file <- paste0("data/yonsei_mtg15000_", subject_id, "_results.txt")
      write.table(result_df, file = out_file,
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }

  }, error = function(e) {
    cat("Error for SubjectID:", subject_id, "-", e$message, "\n")
  })
}

# Combine all per-subject results into a single data frame and save
all_mtg <- do.call(rbind, results_list)
write.csv(all_mtg, "data/yonsei.mtg.csv", row.names = FALSE)

cat("MF shared ASV table saved: data/yonsei.mtg.csv\n")
cat("Total subjects processed:", length(results_list), "\n")
