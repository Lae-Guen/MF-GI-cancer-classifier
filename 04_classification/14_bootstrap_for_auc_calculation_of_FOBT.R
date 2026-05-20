# =============================================================================
# FOBT Benchmark: AUC and Performance Metrics with Bootstrap CI
# Purpose : Compute AUC and 95% bootstrap CI for FOBT as a binary classifier
#           (CRC vs HC), providing a clinical benchmarking reference for the
#           MF microbiome classifiers reported in Figure 5.
#
# FOBT coding in fobt.meta.csv:
#   FOBT = 1 → Negative (no blood detected)
#   FOBT = 2 → Positive (blood detected)
#   → Recoded to FOBT_numeric: 0 = Negative, 1 = Positive
#
# AUC 95% CI: 5000 stratified bootstrap replicates (pROC::ci.auc)
# Sensitivity/Specificity CI: 5000 stratified bootstrap (pROC::ci.coords)
#
# INPUT  : data/fobt.meta.csv — columns: Group (HC/CRC), FOBT (1/2)
# OUTPUT : Console output with AUC, 95% CI, sensitivity, specificity,
#          PPV, NPV, and accuracy (each with 95% bootstrap CI)
# =============================================================================

### bootstrap for auc calculation of FOBT ####

# Load FOBT metadata

data_fobt <- read.csv("./fobt.meta.csv")
data_fobt %>% colnames()
data_fobt <- data_fobt %>% select(Group, FOBT)

# Set Group factor levels: HC (negative) first, CRC (positive) second
data_fobt$Group <- factor(data_fobt$Group, levels = c("HC", "CRC"))

# Recode FOBT result: 1 = Negative (No), 2 = Positive (Yes)
data_fobt$FOBT_label <- factor(data_fobt$FOBT,
                               levels = c(1, 2),
                               labels = c("No", "Yes"))

data_fobt$FOBT_numeric <- as.numeric(data_fobt$FOBT_label) - 1

# Validate data and remove rows with missing Group or FOBT values
if (any(is.na(data_fobt$Group))) {
  warning("Group column contains NA values. Removing rows with NA values.")
  data_fobt <- data_fobt[!is.na(data_fobt$Group), ]
}
if (any(is.na(data_fobt$FOBT_numeric))) { # check the recoded numeric column
  warning("FOBT_numeric column contains NA values. Removing rows with NA values.")
  data_fobt <- data_fobt[!is.na(data_fobt$FOBT_numeric), ]
}

# Inspect data structure
print(str(data_fobt))
print(summary(data_fobt))

# Build ROC object
# response: true class label (Group)
# predictor: FOBT result encoded as numeric (0 = negative, 1 = positive)
roc_obj_fobt <- roc(response = data_fobt$Group, predictor = data_fobt$FOBT_numeric, levels = c("HC", "CRC"))

# pROC direction 'controls < cases': FOBT_numeric = 0 (negative) for HC,
# FOBT_numeric = 1 (positive) for CRC — correct direction.

# Plot ROC curve
plot(roc_obj_fobt, main = "ROC Curve for FOBT Result (CRC vs. HC)",
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)",
     col = "#2ca02c",
     lwd = 2)

# Annotate AUC on plot
auc_val_fobt <- auc(roc_obj_fobt)
text(x = 0.5, y = 0.3, paste("AUC:", round(auc_val_fobt, 3)), cex = 1.2, col = "#2ca02c")

# Compute 95% CI for AUC via 5000 stratified bootstrap replicates
ci_auc_val_fobt <- ci.auc(roc_obj_fobt, method = "bootstrap", boot.n = 5000, strata = data_fobt$Group)
text(x = 0.5, y = 0.25, paste0("95% CI (boot): [", round(ci_auc_val_fobt[1], 3), ", ", round(ci_auc_val_fobt[3], 3), "]"), cex = 1, col = "#2ca02c")
print("AUC 95% Confidence Interval (5000 stratified bootstrap replicates) for FOBT:")
print(ci_auc_val_fobt)

# FOBT is inherently binary: 'Positive' (1) serves as the threshold.
# Using 0.5 as the cutoff is equivalent to classifying by FOBT result directly.
threshold_fobt <- 0.5

cat("\n--- Results for FOBT Classification (Yes/No) ---\n")

ci_fobt_coords <- ci.coords(roc_obj_fobt, x = threshold_fobt, input = "threshold",
                            ret = c("sensitivity", "specificity", "ppv", "npv", "accuracy"),
                            method = "bootstrap", boot.n = 5000,
                            strata = data_fobt$Group)

print("Sensitivity (95% CI):")
print(ci_fobt_coords$sensitivity)
print("Specificity (95% CI):")
print(ci_fobt_coords$specificity)
print("Positive Predictive Value (95% CI):")
print(ci_fobt_coords$ppv)
print("Negative Predictive Value (95% CI):")
print(ci_fobt_coords$npv)
print("Accuracy (95% CI):")
print(ci_fobt_coords$accuracy)