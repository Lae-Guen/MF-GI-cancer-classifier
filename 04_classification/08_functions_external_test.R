# =============================================================================
# External Validation Functions — Algorithm Benchmarking
# Purpose : Provide core functions for evaluating multiple classification
#           algorithms (RF, glmnet, svmRadial, xgbTree) and identity
#           thresholds in the external validation setting.
#           Designed to be sourced by 09_run_scenario_comparison.R.
#
# Key design:
#   - Fixed train/test split via `groups` and `group_partitions` in mikropml
#   - 50 seeds (100–149) run in parallel via future.apply
#   - Patient-level probability = mean of 50-seed predicted probabilities
#   - Optimal classification threshold = Youden's J index on patient-level ROC
#   - pROC::coords() list return handled via [[1]][1]
# =============================================================================

library(mikropml)
library(tidyverse)
library(future)
library(future.apply)
library(doFuture)
library(pROC)

# ── Parallel backend setup ───────────────────────────────────────────────────
# OMP_NUM_THREADS = 1 prevents nested parallelism conflicts on shared servers
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
doFuture::registerDoFuture()
future::plan(future::multicore, workers = 20)

# =============================================================================
# Helper functions
# =============================================================================

# Binomial exact 95% CI (Clopper-Pearson)
# Returns c(estimate, lower, upper); returns NaN when n = 0
calc_binom_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(estimate = NaN, lower = NaN, upper = NaN))
  bt <- binom.test(x, n, conf.level = conf.level)
  c(estimate = unname(bt$estimate),
    lower    = bt$conf.int[1],
    upper    = bt$conf.int[2])
}

# Format "estimate (lower–upper)" string
fmt_ci <- function(est, low, up, digits = 3) {
  paste0(round(est, digits), " (", round(low, digits), "\u2013", round(up, digits), ")")
}

# Patient-level performance metrics (sensitivity, specificity, PPV, NPV, AUC)
# with 95% Clopper-Pearson CIs for binary metrics and 5000-boot CI for AUC
evaluate_patient_level <- function(df_pred,
                                    truth_col      = "Actual_Group",
                                    prob_col       = "Mean_Prob",
                                    pred_col       = "Pred_Class",
                                    positive_class = "CRC",
                                    negative_class = "HC") {
  truth <- factor(df_pred[[truth_col]], levels = c(negative_class, positive_class))
  pred  <- factor(df_pred[[pred_col]],  levels = c(negative_class, positive_class))
  prob  <- df_pred[[prob_col]]

  TP <- sum(truth == positive_class & pred == positive_class, na.rm = TRUE)
  FN <- sum(truth == positive_class & pred == negative_class, na.rm = TRUE)
  TN <- sum(truth == negative_class & pred == negative_class, na.rm = TRUE)
  FP <- sum(truth == negative_class & pred == positive_class, na.rm = TRUE)

  sens_ci <- calc_binom_ci(TP, TP + FN)
  spec_ci <- calc_binom_ci(TN, TN + FP)
  ppv_ci  <- calc_binom_ci(TP, TP + FP)
  npv_ci  <- calc_binom_ci(TN, TN + FN)

  roc_obj <- pROC::roc(response  = truth,
                        predictor = prob,
                        levels    = c(negative_class, positive_class),
                        direction = "<", quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci  <- pROC::ci.auc(roc_obj, method = "bootstrap",
                            boot.n = 1000, strata = truth)

  out <- tibble(
    Metric   = c("Sensitivity", "Specificity", "PPV", "NPV", "AUC"),
    Estimate = c(sens_ci["estimate"], spec_ci["estimate"],
                 ppv_ci["estimate"],  npv_ci["estimate"],  auc_val),
    Lower95  = c(sens_ci["lower"],    spec_ci["lower"],
                 ppv_ci["lower"],     npv_ci["lower"],     auc_ci[1]),
    Upper95  = c(sens_ci["upper"],    spec_ci["upper"],
                 ppv_ci["upper"],     npv_ci["upper"],     auc_ci[3])
  ) %>%
    mutate(Estimate_CI = fmt_ci(Estimate, Lower95, Upper95))

  list(summary = out, roc_obj = roc_obj)
}

# =============================================================================
# Main external validation function
# =============================================================================
# Runs 50-seed parallel classification using a fixed group-based train/test split.
# Supports RF, glmnet, svmRadial, and xgbTree via the `method` argument.
run_external_analysis <- function(input_file,
                                   output_prefix,
                                   method          = "rf",
                                   train_label     = "t1",
                                   test_label      = "v1",
                                   positive_class  = "CRC",
                                   negative_class  = "HC",
                                   seeds           = 100:149,
                                   use_inverse_weight = FALSE) {
  cat("========================================\n")
  cat("Input:  ", input_file, "\n")
  cat("Method: ", method, "\n")
  cat("Output: ", output_prefix, "\n")
  cat("========================================\n")

  df <- read.csv(input_file, check.names = FALSE)
  df[is.na(df)] <- 0
  df$Group     <- factor(df$Group, levels = c(negative_class, positive_class))
  df$Sample_ID <- paste0("Sample_", seq_len(nrow(df)))
  rownames(df) <- df$Sample_ID

  stopifnot(all(c("Group", "set") %in% colnames(df)))

  feature_cols <- setdiff(colnames(df), c("Sample_ID", "Group", "set"))
  model_df     <- df %>% select(Sample_ID, Group, all_of(feature_cols))
  grps         <- df$set

  # Auto-adjust kfold to avoid errors with very few features
  n_features <- length(feature_cols)
  kfold_use  <- min(5L, n_features)
  cat("Features:", n_features, " | kfold:", kfold_use, "\n")
  cat("Train (", train_label, "):", sum(df$set == train_label), "\n")
  cat("Test  (", test_label,  "):", sum(df$set == test_label),  "\n")

  # Class weights for imbalanced data (optional)
  if (use_inverse_weight) {
    train_df   <- model_df[df$set == train_label, ]
    class_cnt  <- table(train_df$Group)
    class_w    <- setNames(as.numeric(max(class_cnt) / class_cnt),
                           names(class_cnt))
    safe_weights <- class_w[as.character(train_df$Group)]
    names(safe_weights) <- NULL
  } else {
    safe_weights <- NULL
  }

  # ── 50-seed parallel runs ─────────────────────────────────────────────────
  results <- future_lapply(seeds, function(seed) {

    ml_args <- list(
      dataset         = model_df %>% select(-Sample_ID),
      method          = method,
      outcome_colname = "Group",
      seed            = seed,
      kfold           = kfold_use,
      groups          = grps,
      group_partitions = list(train = train_label, test = test_label)
    )
    if (!is.null(safe_weights)) ml_args$weights <- safe_weights
    ml_result <- do.call(run_ml, ml_args)

    test_dat  <- ml_result$test_data
    prob_test <- predict(ml_result$trained_model, newdata = test_dat, type = "prob")

    roc_seed <- pROC::roc(
      response  = factor(test_dat$Group, levels = c(negative_class, positive_class)),
      predictor = prob_test[[positive_class]],
      levels    = c(negative_class, positive_class),
      direction = "<", quiet = TRUE
    )
    auc_seed <- as.numeric(pROC::auc(roc_seed))

    pred_df <- tibble(
      Sample_ID    = rownames(test_dat),
      Seed         = seed,
      Actual_Group = test_dat$Group,
      Prob_Target  = prob_test[[positive_class]]
    )

    list(
      patient_preds         = pred_df,
      auc_seed              = tibble(Seed = seed, AUC = auc_seed),
      trained_model_results = as.data.frame(ml_result$trained_model$results)
    )
  }, future.seed = TRUE)

  # ── Collect seed-level results ────────────────────────────────────────────
  auc_seed_tbl  <- bind_rows(lapply(results, `[[`, "auc_seed"))
  pred_seed_tbl <- bind_rows(lapply(results, `[[`, "patient_preds"))

  auc_seed_summary <- auc_seed_tbl %>%
    summarise(n_seed = n(), mean_AUC = mean(AUC), sd_AUC = sd(AUC),
              median_AUC = median(AUC), iqr_AUC = IQR(AUC),
              min_AUC = min(AUC), max_AUC = max(AUC))

  write.csv(auc_seed_tbl,     paste0(output_prefix, ".seed_auc_distribution.csv"), row.names = FALSE)
  write.csv(auc_seed_summary, paste0(output_prefix, ".seed_auc_summary.csv"),      row.names = FALSE)
  write.csv(pred_seed_tbl,    paste0(output_prefix, ".seed_patient_predictions.csv"), row.names = FALSE)

  # ── Patient-level aggregation + Youden's J optimal cutoff ─────────────────
  patient_mean <- pred_seed_tbl %>%
    group_by(Sample_ID, Actual_Group) %>%
    summarise(Mean_Prob = mean(Prob_Target, na.rm = TRUE), .groups = "drop")

  roc_cutoff <- pROC::roc(
    response  = factor(patient_mean$Actual_Group,
                        levels = c(negative_class, positive_class)),
    predictor = patient_mean$Mean_Prob,
    levels    = c(negative_class, positive_class),
    direction = "<", quiet = TRUE
  )
  # [[1]][1] handles pROC::coords() returning a list when ties exist
  youden_coords  <- coords(roc_cutoff, x = "best",
                            best.method = "youden", ret = "threshold")
  optimal_cutoff <- as.numeric(youden_coords[[1]][1])
  cat("Optimal cutoff (Youden's J):", round(optimal_cutoff, 4), "\n")

  patient_final <- patient_mean %>%
    mutate(
      Global_Cutoff = optimal_cutoff,
      Pred_Class    = factor(
        ifelse(Mean_Prob >= optimal_cutoff, positive_class, negative_class),
        levels = c(negative_class, positive_class)
      )
    )

  write.csv(patient_final, paste0(output_prefix, ".patient_level_predictions.csv"),
            row.names = FALSE)

  # ── Patient-level benchmark ───────────────────────────────────────────────
  final_eval  <- evaluate_patient_level(patient_final,
                                         positive_class = positive_class,
                                         negative_class = negative_class)
  final_table <- final_eval$summary
  write.csv(final_table, paste0(output_prefix, ".patient_level_benchmark.csv"),
            row.names = FALSE)

  # ── ROC plot (PNG) ────────────────────────────────────────────────────────
  auc_row <- final_table %>% filter(Metric == "AUC")
  png(paste0(output_prefix, ".patient_level_ROC.png"),
      width = 1800, height = 1600, res = 220)
  plot(final_eval$roc_obj,
       main = paste0(output_prefix, "\nPatient-level ROC"),
       xlab = "1 - Specificity", ylab = "Sensitivity", lwd = 3)
  text(x = 0.60, y = 0.20,
       labels = paste0("AUC = ", round(auc_row$Estimate, 3), "\n",
                        "95% CI = ", round(auc_row$Lower95, 3), " \u2013 ",
                        round(auc_row$Upper95, 3), "\n",
                        "Cutoff = ", round(optimal_cutoff, 3)),
       cex = 1.1)
  dev.off()

  # ── Summary output ────────────────────────────────────────────────────────
  cat("\n--- Results:", output_prefix, "---\n")
  cat("Method:", method, "\n")
  cat("Optimal cutoff:", round(optimal_cutoff, 4), "\n")
  cat("Seed-level AUC: mean =", round(auc_seed_summary$mean_AUC, 3),
      "± SD", round(auc_seed_summary$sd_AUC, 3), "\n")
  print(final_table %>% select(Metric, Estimate_CI))

  # ── Hyperparameter tuning results ─────────────────────────────────────────
  hyper_tbl <- bind_rows(lapply(results, `[[`, "trained_model_results"))
  write.csv(hyper_tbl, paste0(output_prefix, ".hyper_results.csv"), row.names = FALSE)

  return(list(
    auc_seed_distribution   = auc_seed_tbl,
    auc_seed_summary        = auc_seed_summary,
    patient_final           = patient_final,
    patient_level_benchmark = final_table
  ))
}
