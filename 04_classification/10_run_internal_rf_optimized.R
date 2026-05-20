# =============================================================================
# Internal Validation — Random Forest (Optimized)
# Purpose : Train and internally validate a Random Forest classifier across
#           100 seeds (70/30 stratified split) using mikropml.
#           Hyperparameter (mtry) grid is dynamically built based on feature
#           count. Patient-level AUC is derived from mean predicted probabilities
#           across 100 seeds, with Youden's J optimal classification threshold.
#           ROC curves include an SD ribbon (Mean ± SD across seeds).
#
# Key parameters:
#   - method    : "rf" (randomForest package via mikropml; mtry tuning only)
#   - ntree     : 500
#   - seeds     : 100:199 (100 seeds)
#   - kfold     : min(5, n_features)  — auto-adjusted to avoid errors
#   - workers   : 20 (future::multicore)
#   - OMP_NUM_THREADS = 1 (prevents nested parallelism on shared server)
#
# INPUT  : <input_file>.csv — feature table with Group column
# OUTPUT : <prefix>.rf_best_hp.csv / .rf_hp_frequency.csv
#          <prefix>.auc_dist.csv / .auc_summary.csv
#          <prefix>.pred.csv / .patient_predictions.csv
#          <prefix>.patient_benchmark.csv
#          <prefix>.feature_imp_summary.csv
#          <prefix>.patient_level_ROC.pdf
# =============================================================================

## ============================================================
## Internal Validation – Random Forest Only
## 70/30 random split × 100 seeds
## RF hyperparameter optimization (mtry + min.node.size)
## Youden's J optimal cutoff
## get_feature_importance() integration (aggregated across 100 models)
## ============================================================

library(mikropml)
library(tidyverse)
library(future)
library(future.apply)
library(doFuture)
library(pROC)
library(caret)        # createGrid, confusionMatrix
library(randomForest) # feature importance fallback

# ====== Limit internal multithreading (prevents conflicts with future parallelism) ======
options(rf.cores = 1)
Sys.setenv(OMP_NUM_THREADS      = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS      = 1)

doFuture::registerDoFuture()
future::plan(future::multicore, workers = 20)

# ----------------------------------------------------------
# Helper functions
# ----------------------------------------------------------
calc_binom_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(estimate = NaN, lower = NaN, upper = NaN))
  bt <- binom.test(x, n, conf.level = conf.level)
  c(estimate = unname(bt$estimate),
    lower    = bt$conf.int[1],
    upper    = bt$conf.int[2])
}

fmt_ci <- function(est, low, up, digits = 3) {
  paste0(round(est, digits), " (", round(low, digits), "–", round(up, digits), ")")
}

# ----------------------------------------------------------
# RF hyperparameter grid — built dynamically based on feature count
#
# NOTE: mikropml method = "rf" uses the randomForest package (not ranger)
#   → only mtry is tunable (no min.node.size support)
#   → pass hyperparameters as list(mtry = c(...))
#
# mtry strategy by feature count range:
#   20–60   : sqrt, sqrt*1.5, sqrt*2
#   61–120  : sqrt, sqrt*1.5, sqrt*2, log2+1
#   121–200 : sqrt, sqrt*1.5, sqrt*2, log2+1, p/10
# ----------------------------------------------------------
build_rf_hparams <- function(n_features) {

  sqrt_mtry <- max(1L, round(sqrt(n_features)))

  if (n_features <= 60) {
    mtry_vals <- unique(round(c(sqrt_mtry,
                                sqrt_mtry * 1.5,
                                sqrt_mtry * 2)))
  } else if (n_features <= 120) {
    mtry_vals <- unique(round(c(sqrt_mtry,
                                sqrt_mtry * 1.5,
                                sqrt_mtry * 2,
                                floor(log2(n_features)) + 1L)))
  } else {
    mtry_vals <- unique(round(c(sqrt_mtry,
                                sqrt_mtry * 1.5,
                                sqrt_mtry * 2,
                                floor(log2(n_features)) + 1L,
                                floor(n_features / 10))))
  }

  # Clip mtry to [1, n_features] and convert to integer
  mtry_vals <- sort(unique(as.integer(pmax(1L, pmin(mtry_vals, n_features)))))

  # Return as named list (mikropml hyperparameters format)
  list(mtry = mtry_vals)
}

# ----------------------------------------------------------
# Internal validation function (RF only)
# ----------------------------------------------------------
run_internal_rf <- function(input_file,
                            output_prefix,
                            positive_class = "CRC",
                            negative_class = "HC",
                            seeds          = 100:199,   # 100 seeds
                            ntree          = 500) {

  cat("========================================\n")
  cat("Input       :", input_file, "\n")
  cat("Method      : rf (optimized)\n")
  cat("Output prefix:", output_prefix, "\n")
  cat("ntree       :", ntree, "\n")
  cat("========================================\n")

  # ---- Load data ----
  df <- read.csv(input_file, check.names = FALSE)
  df[is.na(df)] <- 0
  df$Group     <- factor(df$Group, levels = c(negative_class, positive_class))
  df$Sample_ID <- paste0("Sample_", seq_len(nrow(df)))
  rownames(df) <- df$Sample_ID

  feature_cols <- setdiff(colnames(df), c("Sample_ID", "Group", "set"))
  n_features   <- length(feature_cols)
  kfold_use    <- min(5L, n_features)

  # ---- RF hyperparameters ----
  rf_hparams <- build_rf_hparams(n_features)

  cat("Features    :", n_features, "\n")
  cat("kfold       :", kfold_use, "(auto-adjusted)\n")
  cat("Samples     :", nrow(df), "\n")
  cat("mtry candidates:", paste(rf_hparams$mtry, collapse = ", "), "\n\n")

  model_df <- df

  # ---- Parallel execution across seeds ----
  results <- future_lapply(seeds, function(seed) {

    set.seed(seed)
    train_idx <- get_partition_indices(model_df$Group, training_frac = 0.70)

    # mikropml run_ml: randomForest-based RF (method = 'rf')
    # Pass hyperparameter grid directly for tuning
    ml_result <- run_ml(
      dataset         = model_df %>% select(-Sample_ID),
      method          = "rf",
      outcome_colname = "Group",
      seed            = seed,
      kfold           = kfold_use,
      training_frac   = train_idx,
      hyperparameters = rf_hparams,   # list(mtry = c(...))
      ntree           = ntree         # passed as extra argument (...) to randomForest
    )

    # ---- Test set prediction ----
    test_dat  <- ml_result$test_data
    prob_test <- predict(ml_result$trained_model,
                         newdata = test_dat, type = "prob")

    roc_seed <- pROC::roc(
      response  = factor(test_dat$Group,
                         levels = c(negative_class, positive_class)),
      predictor = prob_test[[positive_class]],
      levels    = c(negative_class, positive_class),
      direction = "<",
      quiet     = TRUE
    )

    auc_seed <- as.numeric(pROC::auc(roc_seed))

    # ---- Record best hyperparameter for this seed ----
    best_hp <- ml_result$trained_model$bestTune
    best_hp$Seed <- seed
    best_hp$AUC  <- auc_seed

    # ---- Prediction table ----
    pred_df <- tibble(
      Sample_ID    = rownames(test_dat),
      Seed         = seed,
      Actual_Group = test_dat$Group,
      Prob_CRC     = prob_test[[positive_class]]
    )

    # ---- Feature importance (permutation-based via mikropml) ----
    # mikropml::get_feature_importance() returns permutation-based importance
    # Passes outcome_colname, seed, and other required arguments
    imp_result <- tryCatch({
      get_feature_importance(
        trained_model   = ml_result$trained_model,
        train_data      = ml_result$trained_model$trainingData,
        test_data       = test_dat,
        outcome_colname = "Group",
        perf_metric_function = caret::multiClassSummary,
        perf_metric_name     = "AUC",
        class_probs          = TRUE,
        seed                 = seed
      )
    }, error = function(e) {
      # Fallback: use caret::varImp if get_feature_importance fails
      vi <- caret::varImp(ml_result$trained_model, scale = FALSE)$importance
      vi_df <- rownames_to_column(vi, var = "Feature")
      colnames(vi_df)[2] <- "importance"
      vi_df$Seed <- seed
      vi_df
    })

    # Normalize imp_result — standardize column names
    if ("feat" %in% colnames(imp_result)) {
      imp_df <- imp_result %>%
        select(Feature = feat,
               importance = contains("perf_metric"),
               everything()) %>%
        mutate(Seed = seed)
    } else {
      imp_df <- imp_result %>% mutate(Seed = seed)
    }

    list(
      best_hp      = as.data.frame(best_hp),
      auc_seed     = tibble(Seed = seed, AUC = auc_seed),
      patient_preds = pred_df,
      feature_imp  = imp_df
    )
  }, future.seed = TRUE)

  # ---- Collect results ----
  best_hp_tbl <- bind_rows(lapply(results, `[[`, "best_hp"))
  auc_tbl     <- bind_rows(lapply(results, `[[`, "auc_seed"))
  pred_tbl    <- bind_rows(lapply(results, `[[`, "patient_preds"))
  imp_tbl_raw <- bind_rows(lapply(results, `[[`, "feature_imp"))

  # ---- Aggregate feature importance across seeds ----
  # Handle flexible column names from get_feature_importance output
  imp_col <- intersect(c("importance", "Overall", "CRC", positive_class),
                       colnames(imp_tbl_raw))[1]

  feature_col <- intersect(c("Feature", "feat"), colnames(imp_tbl_raw))[1]

  imp_summary <- imp_tbl_raw %>%
    rename(Feature    = !!sym(feature_col),
           Importance = !!sym(imp_col)) %>%
    group_by(Feature) %>%
    summarise(
      mean_importance   = mean(Importance,  na.rm = TRUE),
      sd_importance     = sd(Importance,    na.rm = TRUE),
      median_importance = median(Importance, na.rm = TRUE),
      n_seeds_present   = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_importance))

  # ---- Hyperparameter selection frequency (mtry only for randomForest) ----
  hp_freq <- best_hp_tbl %>%
    count(mtry, name = "n_selected") %>%
    arrange(desc(n_selected))

  # ---- Save CSV outputs ----
  write.csv(best_hp_tbl,  paste0(output_prefix, ".rf_best_hp.csv"),       row.names = FALSE)
  write.csv(hp_freq,      paste0(output_prefix, ".rf_hp_frequency.csv"),   row.names = FALSE)
  write.csv(auc_tbl,      paste0(output_prefix, ".auc_dist.csv"),          row.names = FALSE)
  write.csv(pred_tbl,     paste0(output_prefix, ".pred.csv"),              row.names = FALSE)
  write.csv(imp_tbl_raw,  paste0(output_prefix, ".feature_imp_all.csv"),   row.names = FALSE)
  write.csv(imp_summary,  paste0(output_prefix, ".feature_imp_summary.csv"), row.names = FALSE)

  # ---- AUC summary statistics ----
  auc_summary <- auc_tbl %>%
    summarise(
      n_seed     = n(),
      mean_AUC   = mean(AUC,   na.rm = TRUE),
      sd_AUC     = sd(AUC,     na.rm = TRUE),
      median_AUC = median(AUC, na.rm = TRUE),
      iqr_AUC    = IQR(AUC,    na.rm = TRUE),
      min_AUC    = min(AUC,    na.rm = TRUE),
      max_AUC    = max(AUC,    na.rm = TRUE)
    )

  write.csv(auc_summary, paste0(output_prefix, ".auc_summary.csv"), row.names = FALSE)

  # ---- Patient-level aggregation + Youden's J optimal cutoff ----
  # SD_Prob: per-patient SD across 100 seeds → used for SD ribbon on ROC plot
  patient_mean <- pred_tbl %>%
    group_by(Sample_ID, Actual_Group) %>%
    summarise(
      Mean_Prob = mean(Prob_CRC, na.rm = TRUE),
      SD_Prob   = sd(Prob_CRC,   na.rm = TRUE),
      N_seeds   = n(),
      .groups   = "drop"
    ) %>%
    mutate(SD_Prob = replace_na(SD_Prob, 0))

  roc_for_cutoff <- pROC::roc(
    response  = factor(patient_mean$Actual_Group,
                       levels = c(negative_class, positive_class)),
    predictor = patient_mean$Mean_Prob,
    levels    = c(negative_class, positive_class),
    direction = "<",
    quiet     = TRUE
  )
  youden_coords  <- coords(roc_for_cutoff, x = "best",
                            best.method = "youden", ret = "threshold")
  optimal_cutoff <- as.numeric(youden_coords[[1]][1])

  patient_final <- patient_mean %>%
    mutate(
      Cutoff     = optimal_cutoff,
      Pred_Class = factor(
        ifelse(Mean_Prob >= optimal_cutoff, positive_class, negative_class),
        levels = c(negative_class, positive_class)
      )
    )

  # ---- Patient-level performance metrics ----
  truth <- factor(patient_final$Actual_Group, levels = c(negative_class, positive_class))
  pred  <- factor(patient_final$Pred_Class,   levels = c(negative_class, positive_class))

  TP <- sum(truth == positive_class & pred == positive_class, na.rm = TRUE)
  FN <- sum(truth == positive_class & pred == negative_class, na.rm = TRUE)
  TN <- sum(truth == negative_class & pred == negative_class, na.rm = TRUE)
  FP <- sum(truth == negative_class & pred == positive_class, na.rm = TRUE)

  sens_ci <- calc_binom_ci(TP, TP + FN)
  spec_ci <- calc_binom_ci(TN, TN + FP)
  ppv_ci  <- calc_binom_ci(TP, TP + FP)
  npv_ci  <- calc_binom_ci(TN, TN + FN)

  auc_val <- as.numeric(pROC::auc(roc_for_cutoff))
  auc_ci  <- pROC::ci.auc(roc_for_cutoff,
                            method = "bootstrap", boot.n = 5000, strata = truth)

  patient_benchmark <- tibble(
    Metric   = c("Sensitivity", "Specificity", "PPV", "NPV", "AUC"),
    Estimate = c(sens_ci["estimate"], spec_ci["estimate"],
                 ppv_ci["estimate"],  npv_ci["estimate"], auc_val),
    Lower95  = c(sens_ci["lower"],    spec_ci["lower"],
                 ppv_ci["lower"],     npv_ci["lower"],    auc_ci[1]),
    Upper95  = c(sens_ci["upper"],    spec_ci["upper"],
                 ppv_ci["upper"],     npv_ci["upper"],    auc_ci[3])
  ) %>%
    mutate(Estimate_CI = fmt_ci(Estimate, Lower95, Upper95))

  write.csv(patient_final,     paste0(output_prefix, ".patient_predictions.csv"), row.names = FALSE)
  write.csv(patient_benchmark, paste0(output_prefix, ".patient_benchmark.csv"),   row.names = FALSE)

  # ---- Compute upper/lower ROC bounds for SD ribbon visualization ----
  patient_upper <- patient_mean %>%
    mutate(Prob_upper = pmin(Mean_Prob + SD_Prob, 1))
  patient_lower <- patient_mean %>%
    mutate(Prob_lower = pmax(Mean_Prob - SD_Prob, 0))

  roc_upper <- pROC::roc(
    response  = factor(patient_upper$Actual_Group,
                       levels = c(negative_class, positive_class)),
    predictor = patient_upper$Prob_upper,
    levels    = c(negative_class, positive_class),
    direction = "<", quiet = TRUE
  )
  roc_lower <- pROC::roc(
    response  = factor(patient_lower$Actual_Group,
                       levels = c(negative_class, positive_class)),
    predictor = patient_lower$Prob_lower,
    levels    = c(negative_class, positive_class),
    direction = "<", quiet = TRUE
  )

  extract_roc_xy <- function(roc_obj) {
    fpr <- 1 - rev(roc_obj$specificities)
    tpr <- rev(roc_obj$sensitivities)
    list(fpr = fpr, tpr = tpr)
  }

  build_ribbon_polygon <- function(roc_up, roc_lo) {
    fpr_grid <- seq(0, 1, length.out = 200)
    xy_up <- extract_roc_xy(roc_up)
    xy_lo <- extract_roc_xy(roc_lo)
    tpr_up <- approx(xy_up$fpr, xy_up$tpr, xout = fpr_grid,
                     method = "linear", rule = 2)$y
    tpr_lo <- approx(xy_lo$fpr, xy_lo$tpr, xout = fpr_grid,
                     method = "linear", rule = 2)$y
    list(
      x = c(fpr_grid, rev(fpr_grid)),
      y = c(tpr_up,   rev(tpr_lo))
    )
  }

  ribbon_poly <- build_ribbon_polygon(roc_upper, roc_lower)

  # ---- Save ROC as PDF (vector graphics; editable in Adobe Illustrator) ----
  auc_row    <- patient_benchmark %>% filter(Metric == "AUC")
  ribbon_col <- adjustcolor("#E07B39", alpha.f = 0.25)

  pdf(paste0(output_prefix, ".patient_level_ROC.pdf"),
      width = 4, height = 4)
  par(mar = c(4.5, 4.5, 3, 1.5))

  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  # 1) SD ribbon
  polygon(ribbon_poly$x, ribbon_poly$y,
          col = ribbon_col, border = NA)

  # 2) Diagonal reference line
  abline(a = 0, b = 1, lty = 2, col = "grey50", lwd = 1)

  # 3) Main ROC curve (mean predicted probability)
  roc_main <- pROC::roc(
    response  = factor(patient_final$Actual_Group,
                       levels = c(negative_class, positive_class)),
    predictor = patient_final$Mean_Prob,
    levels    = c(negative_class, positive_class),
    direction = "<", quiet = TRUE
  )
  xy_mean <- extract_roc_xy(roc_main)
  lines(xy_mean$fpr, xy_mean$tpr, col = "#E07B39", lwd = 2.5)

  # 4) Upper/lower boundary lines (thin dashed)
  xy_up <- extract_roc_xy(roc_upper)
  xy_lo <- extract_roc_xy(roc_lower)
  lines(xy_up$fpr, xy_up$tpr, col = "#E07B39", lwd = 1, lty = 3)
  lines(xy_lo$fpr, xy_lo$tpr, col = "#E07B39", lwd = 1, lty = 3)

  # 5) Axes and labels
  axis(1, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), cex.axis = 0.85)
  axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), cex.axis = 0.85, las = 1)
  title(
    main     = output_prefix,
    xlab     = "1 - Specificity (FPR)",
    ylab     = "Sensitivity (TPR)",
    cex.main = 0.75, cex.lab = 0.9
  )
  box()

  # 6) AUC annotation text
  text(
    x      = 0.62, y = 0.18,
    labels = paste0(
      "AUC = ",    round(auc_row$Estimate, 3), "\n",
      "95% CI: ",  round(auc_row$Lower95,  3), " – ",
                   round(auc_row$Upper95,  3), "\n",
      "Cutoff = ", round(optimal_cutoff,   3)
    ),
    cex = 0.75, adj = c(0, 0)
  )

  dev.off()
  cat("\n--- Results:", output_prefix, "---\n")
  cat("Optimal cutoff:", round(optimal_cutoff, 4), "\n")
  cat("Seed-level AUC: mean =", round(auc_summary$mean_AUC, 3),
      "± SD", round(auc_summary$sd_AUC, 3),
      "| median =", round(auc_summary$median_AUC, 3), "\n")
  cat("Most selected HP (top 3):\n")
  print(head(hp_freq, 3))
  cat("Top 10 features (mean importance):\n")
  print(head(imp_summary %>% select(Feature, mean_importance, sd_importance), 10))
  cat("Patient-level:\n")
  print(patient_benchmark %>% select(Metric, Estimate_CI))
  cat("\n")

  return(list(
    auc_distribution  = auc_tbl,
    auc_summary       = auc_summary,
    predictions       = pred_tbl,
    patient_benchmark = patient_benchmark,
    optimal_cutoff    = optimal_cutoff,
    hp_frequency      = hp_freq,
    feature_importance = imp_summary
  ))
}

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
    run_internal_rf(
      input_file    = paste0("./", scen$file),
      output_prefix = prefix
    )
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(res)) {
    auc_row  <- res$patient_benchmark %>% filter(Metric == "AUC")
    sens_row <- res$patient_benchmark %>% filter(Metric == "Sensitivity")
    spec_row <- res$patient_benchmark %>% filter(Metric == "Specificity")
    ppv_row  <- res$patient_benchmark %>% filter(Metric == "PPV")
    npv_row  <- res$patient_benchmark %>% filter(Metric == "NPV")

    # Record most frequently selected hyperparameter
    top_hp <- res$hp_frequency[1, ]

    all_results[[prefix]] <- tibble(
      Scenario        = scen$name,
      Method          = "rf",
      Cutoff          = res$optimal_cutoff,
      Best_mtry       = res$hp_frequency$mtry[1],
      HP_n_selected   = res$hp_frequency$n_selected[1],
      Seed_AUC_mean   = res$auc_summary$mean_AUC,
      Seed_AUC_sd     = res$auc_summary$sd_AUC,
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
cat("=== INTERNAL VALIDATION – RF COMPARISON ===\n")
cat("========================================\n\n")
print(comparison, n = 50)

write.csv(comparison, "internal_optimized_comparison.csv", row.names = FALSE)
cat("\nSaved: internal_optimized_comparison.csv\n")
