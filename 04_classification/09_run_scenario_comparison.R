# =============================================================================
# Scenario × Algorithm Comparison — External Validation (Zackular Cohort)
# Purpose : Systematically evaluate all combinations of VSEARCH identity
#           thresholds (70%–100%) × feature types (abundance/binary) ×
#           classification algorithms (glmnet, RF, svmRadial, xgbTree).
#           Identifies the optimal feature construction strategy for CRC
#           classification using MF-transmitted bacteria.
#
# Source  : 08_functions_external_test.R (provides run_external_analysis())
# INPUT   : fcrc_zackular_seq<pid>_<type>.csv — one file per scenario
#           (type = abundance_zscore or binary)
# OUTPUT  : zackular_seqs_scenario_method_comparison.csv — ranked comparison table
#           Per-scenario: .<method>.seed_auc_*.csv, .patient_level_*.csv, .ROC.png
# =============================================================================

source("08_functions_external_test.R")

# =============================================================================
# Scenario definitions: all identity thresholds × feature types
# =============================================================================
# Identity thresholds tested: 100%, 99%, 98%, 96%, 94%, 90%, 85%, 80%, 70%
# Feature types: abundance (cohort-wise z-scored) and binary (presence/absence)

identity_thresholds <- c("0.70", "0.80", "0.85", "0.90", "0.94",
                          "0.96", "0.98", "0.99", "1.00")

scenarios <- c(
  lapply(identity_thresholds, function(pid) {
    list(file = paste0("fcrc_zackular_seq", pid, "_abundance_zscore.csv"),
         name = paste0("fcrc_zackular_seq", pid, "_abundance"))
  }),
  lapply(identity_thresholds, function(pid) {
    list(file = paste0("fcrc_zackular_seq", pid, "_binary.csv"),
         name = paste0("fcrc_zackular_seq", pid, "_binary"))
  })
)

# =============================================================================
# Algorithm set (comment out any method to exclude from comparison)
# NOTE: xgbTree is excluded from the final paper analyses due to CPU load.
#       It is included here for benchmarking purposes only.
# =============================================================================
methods <- c("glmnet", "rf", "svmRadial", "xgbTree")

# =============================================================================
# Run all scenario × method combinations
# =============================================================================
all_results <- list()

for (scen in scenarios) {
  for (m in methods) {

    prefix <- paste0(scen$name, ".", m)
    cat("\n##################################################\n")
    cat("## ", prefix, "\n")
    cat("##################################################\n")

    res <- tryCatch({
      run_external_analysis(
        input_file         = paste0("./", scen$file),
        output_prefix      = prefix,
        method             = m,
        use_inverse_weight = FALSE
      )
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(res)) {
      auc_row <- res$patient_level_benchmark %>% filter(Metric == "AUC")
      all_results[[prefix]] <- tibble(
        Scenario        = scen$name,
        Method          = m,
        Patient_AUC     = auc_row$Estimate,
        Patient_AUC_CI  = auc_row$Estimate_CI,
        Seed_AUC_mean   = res$auc_seed_summary$mean_AUC,
        Seed_AUC_sd     = res$auc_seed_summary$sd_AUC
      )
    }
  }
}

# =============================================================================
# Final ranked comparison table
# =============================================================================
comparison <- bind_rows(all_results) %>%
  arrange(desc(Patient_AUC))

cat("\n\n========================================\n")
cat("=== FINAL SCENARIO × METHOD COMPARISON ===\n")
cat("========================================\n\n")
print(comparison, n = 50)

write.csv(comparison, "zackular_seqs_scenario_method_comparison.csv", row.names = FALSE)
cat("Saved: zackular_seqs_scenario_method_comparison.csv\n")
