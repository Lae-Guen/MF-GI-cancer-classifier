# =============================================================================
# Merge Yonsei (Training) and External (Test) Feature Tables
# Purpose : Align external cohort feature tables to the Yonsei training set
#           column structure (adding missing features as 0), then concatenate.
#           For abundance features, cohort-wise z-score normalization is applied
#           to remove batch effects before classification.
#           Binary features do not require normalization.
#
# Merging strategy:
#   - Features present in train but absent in test → filled with 0
#   - Features present in test but absent in train → discarded
#   - Column order in the merged file follows the training set exactly
#
# INPUT  : fyonsei_1.00_binary.csv              — Yonsei training set
#          fcrc.zhang.seq1.00_binary_fobt.csv   — external test set
#          fcrc_Zhang-Yonsei_gen_abund_late.csv — pre-merged abundance file
#                                                 (for z-score normalization)
# OUTPUT : fcrc_zhang-yonsei_seq1.00_binary_fobt.csv — merged binary dataset
#          fcrc_Zhang-Yonsei_gen_abund_late_zscore.csv — z-score normalized file
# =============================================================================

library(tidyverse)

# =============================================================================
# Helper: align test set columns to match training set
# =============================================================================
# Features missing from test → set to 0
# Features in test not in train → dropped
align_test_to_train <- function(train, test,
                                 meta_cols = c("Group")) {
  train_features   <- setdiff(colnames(train), meta_cols)
  test_features    <- setdiff(colnames(test),  meta_cols)

  # Add missing feature columns to test (fill with 0)
  missing_features <- setdiff(train_features, test_features)
  for (col in missing_features) test[[col]] <- 0

  # Add missing meta columns as NA (e.g., if test has no "set" column yet)
  missing_meta <- setdiff(meta_cols, colnames(test))
  for (col in missing_meta) test[[col]] <- NA

  # Reorder test to match training column order exactly
  test <- test %>% select(all_of(colnames(train)))
  return(test)
}

# =============================================================================
# Part 1. Merge binary feature tables
# =============================================================================
train <- read.csv("fyonsei_1.00_binary.csv",               check.names = FALSE)
test  <- read.csv("fcrc.zhang.seq1.00_binary_fobt.csv",    check.names = FALSE)

test_aligned <- align_test_to_train(train, test, meta_cols = c("Group"))

# Label cohort membership
train$set        <- "t1"   # Yonsei = training
test_aligned$set <- "v1"   # external = test

merged <- bind_rows(train, test_aligned)
write.csv(merged, "fcrc_zhang-yonsei_seq1.00_binary_fobt.csv", row.names = FALSE)
cat("Saved: fcrc_zhang-yonsei_seq1.00_binary_fobt.csv\n")
cat("  Train (t1):", sum(merged$set == "t1"),
    "  Test  (v1):", sum(merged$set == "v1"), "\n")

# =============================================================================
# Part 2. Cohort-wise z-score normalization for abundance feature tables
# Purpose : Remove systematic abundance differences between cohorts (batch effect)
#           without altering within-cohort relative relationships.
#           Applied independently per cohort and per feature:
#             z = (x - mean_cohort) / sd_cohort
#           Features with sd = 0 within a cohort are assigned sd = 1e-10
#           to avoid division by zero.
# =============================================================================
df <- read.csv("fcrc_Zhang-Yonsei_gen_abund_late.csv", check.names = FALSE)

feature_cols <- setdiff(colnames(df), c("Group", "set"))

cat("=== Input data ===\n")
cat("Rows:", nrow(df), "| Features:", length(feature_cols), "\n")
cat("Cohort (set) distribution:",
    paste(names(table(df$set)), table(df$set), sep = "=", collapse = ", "), "\n")

df_z <- df

for (cohort in unique(df$set)) {
  mask <- df$set == cohort
  for (feat in feature_cols) {
    vals <- df[[feat]][mask]
    m    <- mean(vals, na.rm = TRUE)
    s    <- sd(vals,   na.rm = TRUE)
    if (is.na(s) || s == 0) s <- 1e-10
    df_z[[feat]][mask] <- (vals - m) / s
  }
}

# Verify: each cohort × feature should have mean ≈ 0 and SD ≈ 1
cat("\n=== Z-score verification (per cohort) ===\n")
for (cohort in unique(df$set)) {
  means <- colMeans(df_z[df_z$set == cohort, feature_cols])
  sds   <- apply(df_z[df_z$set == cohort, feature_cols], 2, sd)
  cat(sprintf("Cohort [%s] → max|mean|: %.2e  |  mean(SD): %.4f\n",
              cohort, max(abs(means)), mean(sds, na.rm = TRUE)))
}

out_file <- "fcrc_Zhang-Yonsei_gen_abund_late_zscore.csv"
write.csv(df_z, out_file, row.names = FALSE)
cat("\nSaved:", out_file, "\n")
cat("  Columns:", ncol(df_z), "| Rows:", nrow(df_z), "\n")
