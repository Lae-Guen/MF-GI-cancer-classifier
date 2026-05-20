# =============================================================================
# External Validation — Random Forest (Optimized)
# Purpose : Apply the optimized RF pipeline from 11_functions_external_rf_optimized.R
#           to external cohort scenarios using a fixed group-partitioned
#           train (Yonsei) / test (external) split.
#           Add scenarios to the list below and run.
#
# Source  : 11_functions_external_rf_optimized.R
# INPUT   : Merged train+test CSV files (produced by 07_merging_external_dataset.R)
# OUTPUT  : external_optimized_comparison.csv — ranked patient-level AUC table
#           Per-scenario: *.rf_best_hp.csv, *.patient_level_*.csv, *.patient_level_ROC.pdf
# =============================================================================

## ============================================================
## External Validation – RF Only (Optimized)
## run_external_rf_optimized.R
## ============================================================

source("functions_external_rf_optimized.R")

# ----------------------------------------------------------
# Scenario configuration (RF only; add file paths below)
# ----------------------------------------------------------
scenarios <- list(
  list(file = "",
       name = "")
)

# ----------------------------------------------------------
# Run all scenarios and collect results
# ----------------------------------------------------------
all_results <- list()

for (scen in scenarios) {

  prefix <- paste0(scen$name, ".rf")

  cat("\n##################################################\n")
  cat("##", prefix, "\n")
  cat("##################################################\n")

  res <- tryCatch({
    run_external_rf(
      input_file         = paste0("./", scen$file),
      output_prefix      = prefix,
      use_inverse_weight = FALSE
    )
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(res)) {
    auc_row  <- res$patient_level_benchmark %>% filter(Metric == "AUC")
    sens_row <- res$patient_level_benchmark %>% filter(Metric == "Sensitivity")
    spec_row <- res$patient_level_benchmark %>% filter(Metric == "Specificity")
    ppv_row  <- res$patient_level_benchmark %>% filter(Metric == "PPV")
    npv_row  <- res$patient_level_benchmark %>% filter(Metric == "NPV")

    top_hp <- res$hp_frequency[1, ]

    all_results[[prefix]] <- tibble(
      Scenario        = scen$name,
      Method          = "rf",
      Cutoff          = res$optimal_cutoff,
      Best_mtry       = top_hp$mtry,
      HP_n_selected   = top_hp$n_selected,
      Seed_AUC_mean   = res$auc_seed_summary$mean_AUC,
      Seed_AUC_sd     = res$auc_seed_summary$sd_AUC,
      Patient_AUC     = auc_row$Estimate,
      Patient_AUC_CI  = auc_row$Estimate_CI,
      Patient_Sens    = sens_row$Estimate,
      Patient_Spec    = spec_row$Estimate,
      Patient_PPV     = ppv_row$Estimate,
      Patient_NPV     = npv_row$Estimate
    )
  }
}

# ----------------------------------------------------------
# Final ranked comparison table
# ----------------------------------------------------------
comparison <- bind_rows(all_results) %>% arrange(desc(Patient_AUC))

cat("\n\n========================================\n")
cat("=== EXTERNAL VALIDATION – RF COMPARISON ===\n")
cat("========================================\n\n")
print(comparison, n = 50)

write.csv(comparison, "external_optimized_comparison.csv", row.names = FALSE)
cat("\nSaved: external_optimized_comparison.csv\n")
