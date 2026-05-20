## ============================================================
## External Validation – RF Only (Optimized)
## functions_external_rf_optimized.R
##
## Applies the same hyperparameter strategy as the internal validation pipeline:
##   - Dynamic mtry grid based on feature count
##   - min.node.size tuning (not applicable for randomForest; retained for reference)
##   - get_feature_importance() integration (aggregated across 100 models)
##   - Fixed train/test split maintained via groups / group_partitions
## ============================================================

library(mikropml)
library(tidyverse)
library(future)
library(future.apply)
library(doFuture)
library(pROC)
library(caret)

# ====== Limit internal multithreading ======
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
# Patient-level evaluation function (with positive/negative class arguments)
# ----------------------------------------------------------
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

  roc_obj <- pROC::roc(
    response  = truth,
    predictor = prob,
    levels    = c(negative_class, positive_class),
    direction = "<",
    quiet     = TRUE
  )

  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci  <- pROC::ci.auc(roc_obj, method = "bootstrap",
                           boot.n = 5000, strata = truth)

  out <- tibble(
    Metric   = c("Sensitivity", "Specificity", "PPV", "NPV", "AUC"),
    Estimate = c(sens_ci["estimate"], spec_ci["estimate"],
                 ppv_ci["estimate"],  npv_ci["estimate"], auc_val),
    Lower95  = c(sens_ci["lower"],    spec_ci["lower"],
                 ppv_ci["lower"],     npv_ci["lower"],    auc_ci[1]),
    Upper95  = c(sens_ci["upper"],    spec_ci["upper"],
                 ppv_ci["upper"],     npv_ci["upper"],    auc_ci[3])
  ) %>%
    mutate(Estimate_CI = fmt_ci(Estimate, Lower95, Upper95))

  list(summary = out, roc_obj = roc_obj)
}

# ----------------------------------------------------------
# RF hyperparameter grid — built dynamically based on feature count
#
# NOTE: mikropml method = "rf" uses the randomForest package (not ranger)
#   → only mtry is tunable
#   → pass as list(mtry = c(...))
#
# mtry ranges by feature count:
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

  mtry_vals <- sort(unique(as.integer(pmax(1L, pmin(mtry_vals, n_features)))))
  list(mtry = mtry_vals)
}

# ----------------------------------------------------------
# External validation main function (RF only, optimized)
# ----------------------------------------------------------
run_external_rf <- function(input_file,
                            output_prefix,
                            train_label        = "t1",
                            test_label         = "v1",
                            positive_class     = "CRC",
                            negative_class     = "HC",
                            seeds              = 100:199,
                            use_inverse_weight = FALSE,
                            ntree              = 500) {

  cat("========================================\n")
  cat("Input        :", input_file, "\n")
  cat("Method       : rf (optimized)\n")
  cat("Output prefix:", output_prefix, "\n")
  cat("ntree        :", ntree, "\n")
  cat("========================================\n")

  # ---- Load data ----
  df <- read.csv(input_file, check.names = FALSE)
  df[is.na(df)] <- 0

  df$Group     <- factor(df$Group, levels = c(negative_class, positive_class))
  df$Sample_ID <- paste0("Sample_", seq_len(nrow(df)))
  rownames(df) <- df$Sample_ID

  stopifnot(all(c("Group", "set", "Sample_ID") %in% colnames(df)))

  feature_cols <- setdiff(colnames(df), c("Sample_ID", "Group", "set"))
  model_df     <- df %>% select(Sample_ID, Group, all_of(feature_cols))
  grps         <- df$set

  n_features <- length(feature_cols)
  kfold_use  <- min(5L, n_features)

  # ---- RF hyperparameters ----
  rf_hparams <- build_rf_hparams(n_features)

  cat("Features     :", n_features, "\n")
  cat("kfold        :", kfold_use, "(auto-adjusted)\n")
  cat("Train (", train_label, "):", sum(df$set == train_label), "\n")
  cat("Test  (", test_label,  "):", sum(df$set == test_label),  "\n")
  cat("mtry candidates:", paste(rf_hparams$mtry, collapse = ", "), "\n\n")

  # ---- Class weights (computed from training subset) ----
  train_rows <- which(df$set == train_label)
  train_df   <- model_df[train_rows, , drop = FALSE]

  if (use_inverse_weight) {
    train_class_counts <- table(train_df$Group)
    class_w <- setNames(
      as.numeric(max(train_class_counts) / train_class_counts),
      names(train_class_counts)
    )
    safe_weights <- class_w[as.character(train_df$Group)]
    names(safe_weights) <- NULL
  } else {
    safe_weights <- NULL
  }

  # ---- Parallel execution across seeds ----
  results <- future_lapply(seeds, function(seed) {

    ml_result <- if (is.null(safe_weights)) {
      run_ml(
        dataset         = model_df %>% select(-Sample_ID),
        method          = "rf",
        outcome_colname = "Group",
        seed            = seed,
        kfold           = kfold_use,
        groups          = grps,
        group_partitions = list(
          train = c(train_label),
          test  = c(test_label)
        ),
        hyperparameters = rf_hparams    # list(mtry = c(...))
      )
    } else {
      run_ml(
        dataset         = model_df %>% select(-Sample_ID),
        method          = "rf",
        outcome_colname = "Group",
        seed            = seed,
        kfold           = kfold_use,
        groups          = grps,
        group_partitions = list(
          train = c(train_label),
          test  = c(test_label)
        ),
        weights         = safe_weights,
        hyperparameters = rf_hparams
      )
    }

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
    best_hp      <- ml_result$trained_model$bestTune
    best_hp$Seed <- seed
    best_hp$AUC  <- auc_seed

    # ---- Prediction table ----
    pred_df <- tibble(
      Sample_ID    = rownames(test_dat),
      Seed         = seed,
      Actual_Group = test_dat$Group,
      Prob_Target  = prob_test[[positive_class]]
    )

    # ---- Feature importance ----
    imp_result <- tryCatch({
      get_feature_importance(
        trained_model        = ml_result$trained_model,
        train_data           = ml_result$trained_model$trainingData,
        test_data            = test_dat,
        outcome_colname      = "Group",
        perf_metric_function = caret::multiClassSummary,
        perf_metric_name     = "AUC",
        class_probs          = TRUE,
        seed                 = seed
      )
    }, error = function(e) {
      # Fallback: use caret varImp if get_feature_importance fails
      vi    <- caret::varImp(ml_result$trained_model, scale = FALSE)$importance
      vi_df <- rownames_to_column(vi, var = "Feature")
      colnames(vi_df)[2] <- "importance"
      vi_df$Seed <- seed
      vi_df
    })

    # Standardize column names
    if ("feat" %in% colnames(imp_result)) {
      imp_df <- imp_result %>%
        select(Feature    = feat,
               importance = contains("perf_metric"),
               everything()) %>%
        mutate(Seed = seed)
    } else {
      imp_df <- imp_result %>% mutate(Seed = seed)
    }

    list(
      best_hp       = as.data.frame(best_hp),
      auc_seed      = tibble(Seed = seed, AUC = auc_seed),
      patient_preds = pred_df,
      feature_imp   = imp_df
    )
  }, future.seed = TRUE)

  # ---- Collect results ----
  best_hp_tbl <- bind_rows(lapply(results, `[[`, "best_hp"))
  auc_seed_tbl <- bind_rows(lapply(results, `[[`, "auc_seed"))
  pred_seed_tbl <- bind_rows(lapply(results, `[[`, "patient_preds"))
  imp_tbl_raw   <- bind_rows(lapply(results, `[[`, "feature_imp"))

  # ---- Aggregate feature importance across seeds ----
  imp_col     <- intersect(c("importance", "Overall", "CRC", positive_class),
                           colnames(imp_tbl_raw))[1]
  feature_col <- intersect(c("Feature", "feat"), colnames(imp_tbl_raw))[1]

  imp_summary <- imp_tbl_raw %>%
    rename(Feature    = !!sym(feature_col),
           Importance = !!sym(imp_col)) %>%
    group_by(Feature) %>%
    summarise(
      mean_importance   = mean(Importance,   na.rm = TRUE),
      sd_importance     = sd(Importance,     na.rm = TRUE),
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
  write.csv(best_hp_tbl,   paste0(output_prefix, ".rf_best_hp.csv"),          row.names = FALSE)
  write.csv(hp_freq,       paste0(output_prefix, ".rf_hp_frequency.csv"),      row.names = FALSE)
  write.csv(auc_seed_tbl,  paste0(output_prefix, ".seed_auc_distribution.csv"), row.names = FALSE)
  write.csv(pred_seed_tbl, paste0(output_prefix, ".seed_patient_predictions.csv"), row.names = FALSE)
  write.csv(imp_tbl_raw,   paste0(output_prefix, ".feature_imp_all.csv"),      row.names = FALSE)
  write.csv(imp_summary,   paste0(output_prefix, ".feature_imp_summary.csv"),  row.names = FALSE)

  # ---- AUC summary statistics ----
  auc_seed_summary <- auc_seed_tbl %>%
    summarise(
      n_seed     = n(),
      mean_AUC   = mean(AUC,   na.rm = TRUE),
      sd_AUC     = sd(AUC,     na.rm = TRUE),
      median_AUC = median(AUC, na.rm = TRUE),
      iqr_AUC    = IQR(AUC,    na.rm = TRUE),
      min_AUC    = min(AUC,    na.rm = TRUE),
      max_AUC    = max(AUC,    na.rm = TRUE)
    )
  write.csv(auc_seed_summary,
            paste0(output_prefix, ".seed_auc_summary.csv"), row.names = FALSE)

  # ---- Patient-level aggregation + Youden's J optimal cutoff ----
  # SD_Prob: per-patient SD across 100 seeds → used for SD ribbon on ROC plot
  patient_mean <- pred_seed_tbl %>%
    group_by(Sample_ID, Actual_Group) %>%
    summarise(
      Mean_Prob = mean(Prob_Target, na.rm = TRUE),
      SD_Prob   = sd(Prob_Target,   na.rm = TRUE),
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

  cat("Optimal cutoff (Youden's J):", round(optimal_cutoff, 4), "\n")

  patient_final <- patient_mean %>%
    mutate(
      Global_Cutoff = optimal_cutoff,
      Pred_Class    = factor(
        ifelse(Mean_Prob >= optimal_cutoff, positive_class, negative_class),
        levels = c(negative_class, positive_class)
      )
    )

  write.csv(patient_final,
            paste0(output_prefix, ".patient_level_predictions.csv"),
            row.names = FALSE)

  # ---- Patient-level benchmark ----
  final_eval  <- evaluate_patient_level(
    patient_final,
    positive_class = positive_class,
    negative_class = negative_class
  )
  final_table <- final_eval$summary

  write.csv(final_table,
            paste0(output_prefix, ".patient_level_benchmark.csv"),
            row.names = FALSE)

  # ---- Compute upper/lower ROC bounds for SD ribbon visualization ----
  # Build separate ROC objects using Mean ± SD clipped probabilities for ribbon
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

  # Extract (FPR, TPR) from pROC object and interpolate to a common FPR grid
  extract_roc_xy <- function(roc_obj) {
    # pROC: specificities = 1-FPR (descending), sensitivities = TPR
    fpr <- 1 - rev(roc_obj$specificities)
    tpr <- rev(roc_obj$sensitivities)
    list(fpr = fpr, tpr = tpr)
  }

  build_ribbon_polygon <- function(roc_up, roc_lo) {
    # Interpolate TPR to a common FPR grid, then build ribbon polygon coordinates
    fpr_grid <- seq(0, 1, length.out = 200)

    xy_up <- extract_roc_xy(roc_up)
    xy_lo <- extract_roc_xy(roc_lo)

    tpr_up <- approx(xy_up$fpr, xy_up$tpr, xout = fpr_grid,
                     method = "linear", rule = 2)$y
    tpr_lo <- approx(xy_lo$fpr, xy_lo$tpr, xout = fpr_grid,
                     method = "linear", rule = 2)$y

    # Polygon: upper boundary left-to-right, lower boundary right-to-left
    list(
      x = c(fpr_grid, rev(fpr_grid)),
      y = c(tpr_up,   rev(tpr_lo))
    )
  }

  ribbon_poly <- build_ribbon_polygon(roc_upper, roc_lower)

  # ---- Save ROC as PDF (vector graphics; editable in Adobe Illustrator) ----
  auc_row <- final_table %>% filter(Metric == "AUC")
  ribbon_col <- adjustcolor("#4878CF", alpha.f = 0.25)   # transparent version of main curve color

  pdf(paste0(output_prefix, ".patient_level_ROC.pdf"),
      width = 4, height = 4)          # in inches; freely resizable in Illustrator
  par(mar = c(4.5, 4.5, 3, 1.5))

  # 1) SD ribbon (shaded polygon)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  polygon(ribbon_poly$x, ribbon_poly$y,
          col = ribbon_col, border = NA)

  # 2) Diagonal reference line
  abline(a = 0, b = 1, lty = 2, col = "grey50", lwd = 1)

  # 3) Main ROC curve (based on patient mean predicted probability)
  xy_mean <- extract_roc_xy(final_eval$roc_obj)
  lines(xy_mean$fpr, xy_mean$tpr, col = "#4878CF", lwd = 2.5)

  # 4) Upper/lower boundary lines (thin dashed)
  xy_up <- extract_roc_xy(roc_upper)
  xy_lo <- extract_roc_xy(roc_lower)
  lines(xy_up$fpr, xy_up$tpr, col = "#4878CF", lwd = 1, lty = 3)
  lines(xy_lo$fpr, xy_lo$tpr, col = "#4878CF", lwd = 1, lty = 3)

  # 5) Axes and labels
  axis(1, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), cex.axis = 0.85)
  axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), cex.axis = 0.85, las = 1)
  title(
    main = output_prefix,
    xlab = "1 - Specificity (FPR)",
    ylab = "Sensitivity (TPR)",
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

  # ---- Console summary ----
  cat("\n--- Results:", output_prefix, "---\n")
  cat("Method       : rf (optimized)\n")
  cat("Optimal cutoff:", round(optimal_cutoff, 4), "\n")
  cat("Seed-level AUC: mean =", round(auc_seed_summary$mean_AUC, 3),
      "± SD", round(auc_seed_summary$sd_AUC, 3),
      "| median =", round(auc_seed_summary$median_AUC, 3), "\n")
  cat("Most selected HP (top 3):\n")
  print(head(hp_freq, 3))
  cat("Top 10 features (mean importance):\n")
  print(head(imp_summary %>% select(Feature, mean_importance, sd_importance), 10))
  cat("Patient-level:\n")
  print(final_table %>% select(Metric, Estimate_CI))
  cat("\n")

  return(list(
    auc_seed_distribution   = auc_seed_tbl,
    auc_seed_summary        = auc_seed_summary,
    patient_final           = patient_final,
    patient_level_benchmark = final_table,
    optimal_cutoff          = optimal_cutoff,
    hp_frequency            = hp_freq,
    feature_importance      = imp_summary
  ))
}
