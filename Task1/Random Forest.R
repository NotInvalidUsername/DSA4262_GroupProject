library(dplyr)
library(pROC)
library(PRROC)
library(caret)
library(ranger)
library(doParallel)
library(randomForest)

## Data
site_features_full <- read.csv("site_features_full.csv")
info <- read.csv("info.csv")
df0_full <- read.csv("train_df_full.csv")
df1_full_no_labels <- read.csv("df1_full_no_labels.csv")
df2_full_no_labels <- read.csv("df2_full_no_labels.csv")

##Train test split
unique_genes <- unique(df0_full$gene_id)
set.seed(42)
test_genes <- sample(unique_genes, floor(length(unique_genes) * 0.2))

df0_train <- df0_full %>% filter(!gene_id %in% test_genes)
df0_test  <- df0_full %>% filter(gene_id %in% test_genes)

cat(sprintf("Train genes: %d | Test genes: %d\n",
            n_distinct(df0_train$gene_id), n_distinct(df0_test$gene_id)))

df0_train$label <- as.numeric(as.character(df0_train$label))
df0_test$label  <- as.numeric(as.character(df0_test$label))

numeric_cols <- df0_train %>% 
  select(-label, -transcript, -position, -gene_id) %>%
  names()

df0_train[, numeric_cols] <- lapply(df0_train[, numeric_cols], function(x) as.numeric(as.character(x)))
df0_test[, numeric_cols]  <- lapply(df0_test[, numeric_cols],  function(x) as.numeric(as.character(x)))

cat("Train/test split complete. ID columns preserved.\n\n")

set.seed(42)

##Random forest
## Version 4
k <- 5
n_search <- 30

ntree_vals <- sample(700:900, n_search, replace = TRUE)
mtry_vals <- sample(2:3, n_search, replace = TRUE)
nodesize_vals <- rep(1, n_search)

classwt_vals <- runif(n_search, min = 45, max = 70) 

majority_prop_vals <- runif(n_search, min = 0.15, max = 0.50)  # down to 0.15

# Console summary
cat("REFINED PR-AUC OPTIMIZATION - ROUND 2\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")
cat("Refinements based on previous results:\n")
cat("Narrowed ntree range around 820 (700–900)\n")
cat("Focused mtry on 2–3 (best was 2)\n")
cat("nodesize fixed at 1 (deep trees)\n")
cat("classwt range 45–70 (higher than 48.8)\n")
cat("majority_prop range 0.15–0.50 (stronger undersampling)\n\n")

cat("Hyperparameter ranges:\n")
cat(sprintf("ntree: [%d, %d]\n", min(ntree_vals), max(ntree_vals)))
cat(sprintf("mtry: [%d, %d]\n", min(mtry_vals), max(mtry_vals)))
cat(sprintf("nodesize: %d\n", unique(nodesize_vals)))
cat(sprintf("classwt (positive): [%.1f, %.1f]\n", min(classwt_vals), max(classwt_vals)))
cat(sprintf("majority_prop: [%.2f, %.2f]\n\n", min(majority_prop_vals), max(majority_prop_vals)))

# ============================================================================
# FEATURE ENGINEERING - ADD INTERACTIONS
# ============================================================================
cat("Feature engineering with interactions...\n")

# NOW it's safe to exclude these columns
feature_cols_base <- df0_train %>%
  select(-label, -transcript, -position, -gene_id) %>%
  select(where(is.numeric)) %>%
  names()

# Remove near-zero variance
x_temp <- df0_train %>% select(all_of(feature_cols_base))
nzv <- nearZeroVar(x_temp, saveMetrics = TRUE)
feature_cols_filtered <- feature_cols_base[!nzv$nzv]

cat(sprintf("  Starting features: %d\n", length(feature_cols_base)))
cat(sprintf("  After NZV removal: %d\n", length(feature_cols_filtered)))

# Quick importance check
set.seed(42)
sample_size <- min(3000, nrow(df0_train))
sample_idx <- sample(nrow(df0_train), sample_size)
x_sample <- df0_train[sample_idx, ] %>% select(all_of(feature_cols_filtered))
y_sample <- as.factor(df0_train$label[sample_idx])

quick_rf <- randomForest(
  x = x_sample, 
  y = y_sample, 
  ntree = 826, 
  mtry = 2,
  importance = TRUE
)

importance_df <- data.frame(
  feature = rownames(quick_rf$importance),
  importance = quick_rf$importance[, "MeanDecreaseGini"]
) %>% arrange(desc(importance))

top_features <- head(importance_df$feature, min(8, length(feature_cols_filtered)))

cat("  Top features by importance:\n")
for (i in 1:min(5, length(top_features))) {
  cat(sprintf("    %d. %s (%.2f)\n", i, top_features[i], importance_df$importance[i]))
}

# Add interactions for top features
interaction_count <- 0
if (length(top_features) >= 3) {
  cat("\n  Adding interaction features...\n")
  for (i in 1:min(4, length(top_features)-1)) {
    for (j in (i+1):min(5, length(top_features))) {
      feat1 <- top_features[i]
      feat2 <- top_features[j]
      new_feat <- paste0("interact_", i, "_", j)
      
      df0_train[[new_feat]] <- 
        df0_train[[feat1]] * df0_train[[feat2]]
      
      feature_cols_filtered <- c(feature_cols_filtered, new_feat)
      interaction_count <- interaction_count + 1
    }
  }
  cat(sprintf("    Added %d interaction features\n", interaction_count))
}

# Add squared terms for top 5
for (i in 1:min(5, length(top_features))) {
  feat <- top_features[i]
  new_feat <- paste0(feat, "_sq")
  df0_train[[new_feat]] <- df0_train[[feat]]^2
  feature_cols_filtered <- c(feature_cols_filtered, new_feat)
}

cat(sprintf("  Final feature count: %d\n\n", length(feature_cols_filtered)))
feature_cols <- feature_cols_filtered

# ============================================================================
# STRATIFIED GENE-BASED FOLDS (FIXED - NO GENE LEAKAGE)
# ============================================================================
cat("\nCreating STRATIFIED gene-based folds...\n")

set.seed(42)

gene_labels <- df0_train %>%
  group_by(gene_id) %>%
  summarise(
    majority_label = as.numeric(mean(label) > 0.5),
    n_samples = n(),
    n_positive = sum(label == 1),
    .groups = 'drop'
  )

positive_genes <- gene_labels %>% filter(majority_label == 1) %>% pull(gene_id)
negative_genes <- gene_labels %>% filter(majority_label == 0) %>% pull(gene_id)

cat(sprintf("  Total genes: %d\n", nrow(gene_labels)))
cat(sprintf("  Positive-labeled genes: %d (%.1f%%)\n", 
            length(positive_genes), 100 * length(positive_genes) / nrow(gene_labels)))
cat(sprintf("  Negative-labeled genes: %d (%.1f%%)\n\n", 
            length(negative_genes), 100 * length(negative_genes) / nrow(gene_labels)))

gene_fold_assignment <- data.frame(
  gene_id = gene_labels$gene_id,
  majority_label = gene_labels$majority_label,
  fold = NA,
  stringsAsFactors = FALSE
)

# stratification
if (length(positive_genes) > 0) {
  pos_idx <- which(gene_fold_assignment$gene_id %in% positive_genes)
  gene_fold_assignment$fold[pos_idx] <- sample(rep(1:k, length.out = length(pos_idx)))
}

if (length(negative_genes) > 0) {
  neg_idx <- which(gene_fold_assignment$gene_id %in% negative_genes)
  gene_fold_assignment$fold[neg_idx] <- sample(rep(1:k, length.out = length(neg_idx)))
}

folds_balanced_rf <- gene_fold_assignment$fold[match(
  df0_train$gene_id, 
  gene_fold_assignment$gene_id
)]

cat("Verifying NO gene leakage...\n")
gene_fold_check <- df0_train %>%
  select(gene_id) %>%
  mutate(fold = folds_balanced_rf) %>%
  group_by(gene_id) %>%
  summarise(
    n_folds = n_distinct(fold),
    .groups = 'drop'
  )

genes_in_multiple_folds <- gene_fold_check %>% filter(n_folds > 1)

if (nrow(genes_in_multiple_folds) > 0) {
  stop("ERROR: Gene leakage detected! Fix fold assignment.")
} else {
  cat("SUCCESS: Each gene in exactly ONE fold\n\n")
}

cat("Fold distribution:\n")
for (i in 1:k) {
  genes_in_fold <- unique(df0_train$gene_id[folds_balanced_rf == i])
  samples_in_fold <- sum(folds_balanced_rf == i)
  pos_samples <- sum(df0_train$label[folds_balanced_rf == i] == 1)
  pos_genes <- sum(genes_in_fold %in% positive_genes)
  
  cat(sprintf("  Fold %d: %d genes (%d pos), %d samples (%d pos, %.1f%%)\n", 
              i, length(genes_in_fold), pos_genes, samples_in_fold, 
              pos_samples, 100 * pos_samples / samples_in_fold))
}
cat("\n")


# ============================================================================
# PARALLEL CROSS-VALIDATION WITH RANGER
# ============================================================================

n_cores <- max(1, parallel::detectCores())
cat(sprintf("Detected %d cores — enabling parallel CV (ranger).\n", n_cores))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cv_preds_balanced_rf <- list()

cat(sprintf("Starting optimized Parallel Ranger CV at %s\n", 
            format(Sys.time(), "%H:%M:%S")))

for (j in 1:n_search) {
  start <- Sys.time()
  
  ntree <- ntree_vals[j]
  mtry <- mtry_vals[j]
  nodesize <- nodesize_vals[j]
  classwt <- classwt_vals[j]
  majority_prop <- majority_prop_vals[j]
  
  cat(sprintf("\n[%2d/%d] ntree=%4d mtry=%d classwt=%.1f maj_prop=%.2f\n", 
              j, n_search, ntree, mtry, classwt, majority_prop))
  
  fold_results <- foreach(
    i = 1:k,
    .packages = c("ranger", "dplyr"),
    .export = c("df0_train", "feature_cols", "folds_balanced_rf")
  ) %dopar% {
    set.seed(42 + j * 10 + i)  
    
    train_idx <- which(folds_balanced_rf != i)
    test_idx  <- which(folds_balanced_rf == i)
    
    train_fold <- df0_train[train_idx, ]
    test_fold  <- df0_train[test_idx, ]
    
    x_train <- train_fold %>% select(all_of(feature_cols))
    y_train <- as.factor(train_fold$label)
    x_test  <- test_fold %>% select(all_of(feature_cols))
    
    class_weights <- c("0" = 1, "1" = classwt)
    
    label_counts <- table(y_train)
    n_minority <- as.numeric(label_counts["1"])
    n_majority <- as.numeric(label_counts["0"])
    
    sample_fraction <- c(
      "0" = min(1, round(n_majority * majority_prop) / n_majority),
      "1" = 1
    )
    
    rf_model <- ranger(
      dependent.variable.name = "label",
      data = data.frame(label = y_train, x_train),
      num.trees = ntree,
      mtry = mtry,
      min.node.size = nodesize,
      class.weights = class_weights,
      sample.fraction = sample_fraction,
      replace = TRUE,
      probability = TRUE,
      num.threads = 1  
    )
    
    preds <- predict(rf_model, data = x_test)$predictions[, "1"]
    
    data.frame(
      fold = i,
      test_idx = test_idx,
      preds = preds
    )
  }
  
  fold_results_df <- bind_rows(fold_results)
  cv_scores_balanced_rf <- numeric(nrow(df0_train))
  cv_scores_balanced_rf[fold_results_df$test_idx] <- fold_results_df$preds
  
  cv_preds_balanced_rf[[j]] <- data.frame(
    transcript = df0_train$transcript,
    position = df0_train$position,
    gene_id = df0_train$gene_id,
    label = df0_train$label,
    score = cv_scores_balanced_rf,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    classwt = classwt,
    majority_prop = majority_prop
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  cat(sprintf("   Completed in %.1fs (%.1f min)\n", elapsed, elapsed / 60))
  gc()  
}

stopCluster(cl)
cat(sprintf("\nParallel Ranger CV completed at %s\n", format(Sys.time(), "%H:%M:%S")))

# ============================================================================
# EVALUATE RANGER RESULTS
# ============================================================================
cv_results_ranger <- bind_rows(cv_preds_balanced_rf)

results_summary_ranger <- cv_results_ranger %>%
  group_by(ntree, mtry, nodesize, classwt, majority_prop) %>%
  summarise(
    auc_roc = as.numeric(auc(roc(label, score, quiet = TRUE))),
    auc_pr = {
      pr_curve <- pr.curve(
        scores.class0 = score[label == 1],
        scores.class1 = score[label == 0],
        curve = TRUE
      )
      pr_curve$auc.integral
    },
    auc_combined = (auc_roc + auc_pr) / 2,
    n_samples = n(),
    n_genes = n_distinct(gene_id),
    .groups = 'drop'
  ) %>%
  arrange(desc(auc_pr))

cat("\n", paste0(rep("=", 70), collapse = ""), "\n")
cat("RANGER RESULTS - ALL HYPERPARAMETER COMBINATIONS\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("All results sorted by PR-AUC:\n\n")
print(results_summary_ranger %>% 
        select(ntree, mtry, nodesize, classwt, majority_prop, 
               auc_roc, auc_pr, auc_combined), n = Inf)

# ============================================================================
# BEST HYPERPARAMETERS
# ============================================================================
cat("\n", paste0(rep("=", 70), collapse = ""), "\n")
cat("BEST HYPERPARAMETERS\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

best_by_pr <- results_summary_ranger %>% slice_max(auc_pr, n = 1)
best_by_roc <- results_summary_ranger %>% slice_max(auc_roc, n = 1)
best_by_combined <- results_summary_ranger %>% slice_max(auc_combined, n = 1)

cat("Best by PR-AUC (Precision-Recall):\n")
cat(sprintf("  ntree:          %d\n", best_by_pr$ntree))
cat(sprintf("  mtry:           %d\n", best_by_pr$mtry))
cat(sprintf("  nodesize:       %d\n", best_by_pr$nodesize))
cat(sprintf("  classwt:        %.2f\n", best_by_pr$classwt))
cat(sprintf("  majority_prop:  %.2f\n", best_by_pr$majority_prop))
cat(sprintf("  ROC-AUC:        %.6f\n", best_by_pr$auc_roc))
cat(sprintf("  PR-AUC:         %.6f\n", best_by_pr$auc_pr))
cat(sprintf("  Combined:       %.6f\n\n", best_by_pr$auc_combined))

cat("Best by ROC-AUC:\n")
cat(sprintf("  ntree:          %d\n", best_by_roc$ntree))
cat(sprintf("  mtry:           %d\n", best_by_roc$mtry))
cat(sprintf("  nodesize:       %d\n", best_by_roc$nodesize))
cat(sprintf("  classwt:        %.2f\n", best_by_roc$classwt))
cat(sprintf("  majority_prop:  %.2f\n", best_by_roc$majority_prop))
cat(sprintf("  ROC-AUC:        %.6f\n", best_by_roc$auc_roc))
cat(sprintf("  PR-AUC:         %.6f\n", best_by_roc$auc_pr))
cat(sprintf("  Combined:       %.6f\n\n", best_by_roc$auc_combined))

cat("Best by Combined AUC:\n")
cat(sprintf("  ntree:          %d\n", best_by_combined$ntree))
cat(sprintf("  mtry:           %d\n", best_by_combined$mtry))
cat(sprintf("  nodesize:       %d\n", best_by_combined$nodesize))
cat(sprintf("  classwt:        %.2f\n", best_by_combined$classwt))
cat(sprintf("  majority_prop:  %.2f\n", best_by_combined$majority_prop))
cat(sprintf("  ROC-AUC:        %.6f\n", best_by_combined$auc_roc))
cat(sprintf("  PR-AUC:         %.6f\n", best_by_combined$auc_pr))
cat(sprintf("  Combined:       %.6f\n\n", best_by_combined$auc_combined))

# ============================================================================
# IMPROVEMENT SUMMARY
# ============================================================================
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("IMPROVEMENT vs BASELINE\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("Your baseline PR-AUC: 0.297750\n")
cat(sprintf("New best PR-AUC:      %.6f\n", best_by_pr$auc_pr))
cat(sprintf("Absolute change:      %+.6f\n", best_by_pr$auc_pr - 0.297750))
cat(sprintf("Relative change:      %+.1f%%\n\n", 
            100 * (best_by_pr$auc_pr - 0.297750) / 0.297750))


# ============================================================================
# RECOMMENDATION
# ============================================================================
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("RECOMMENDED HYPERPARAMETERS\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Choose based on which metric matters most
if (best_by_pr$auc_pr > 0.30) {
  cat("Recommendation: Use PR-AUC optimized parameters\n")
  cat("(PR-AUC > 0.30 indicates good minority class prediction)\n\n")
  recommended <- best_by_pr
} else {
  cat("Recommendation: Use Combined AUC optimized parameters\n")
  cat("(Balances both ROC-AUC and PR-AUC performance)\n\n")
  recommended <- best_by_combined
}

cat("Final parameters:\n")
cat(sprintf("  ntree:          %d\n", recommended$ntree))
cat(sprintf("  mtry:           %d\n", recommended$mtry))
cat(sprintf("  nodesize:       %d\n", recommended$nodesize))
cat(sprintf("  classwt:        %.2f\n", recommended$classwt))
cat(sprintf("  majority_prop:  %.2f\n", recommended$majority_prop))
cat(sprintf("  ROC-AUC:        %.6f\n", recommended$auc_roc))
cat(sprintf("  PR-AUC:         %.6f\n", recommended$auc_pr))

cat("\nHyperparameter optimization complete!\n")

# Function to create engineered features
create_engineered_features <- function(df, top_features) {
  for (i in 1:min(4, length(top_features)-1)) {
    for (j in (i+1):min(5, length(top_features))) {
      feat1 <- top_features[i]
      feat2 <- top_features[j]
      new_feat <- paste0("interact_", i, "_", j)
      
      if (feat1 %in% names(df) && feat2 %in% names(df)) {
        df[[new_feat]] <- df[[feat1]] * df[[feat2]]
      }
    }
  }
  
  for (i in 1:min(5, length(top_features))) {
    feat <- top_features[i]
    new_feat <- paste0(feat, "_sq")
    
    if (feat %in% names(df) && !endsWith(feat, "_sq")) {  # Don't square already squared features
      df[[new_feat]] <- df[[feat]]^2
    }
  }
  
  return(df)
}

# ============================================================
# FINAL MODEL TRAINING AND TEST EVALUATION
# ============================================================

cat("\nTraining final model with recommended parameters...\n")
set.seed(42)

# Creating model with recommended hyperparameters
final_model <- ranger(
  dependent.variable.name = "label",
  data = train_df %>% select(label, all_of(feature_cols)),
  num.trees = recommended$ntree,
  mtry = recommended$mtry,
  min.node.size = recommended$nodesize,
  class.weights = c("0" = 1, "1" = recommended$classwt),
  sample.fraction = c("0" = recommended$majority_prop, "1" = 1),
  probability = TRUE
)

cat("Generating predictions on test set...\n")

test_preds <- predict(
  final_model,
  data = test_df %>% select(all_of(feature_cols))
)$predictions[, "1"]

df0_full$label <- as.factor(df0_full$label)

feature_cols_base <- df0_full %>%
  select(-label, -transcript, -position, -gene_id) %>%
  select(where(is.numeric)) %>%
  names()

df0_full = create_engineered_features(df0_full, top_features)

feature_cols <- df0_full %>%
  select(-label, -transcript, -position, -gene_id) %>%
  names()

set.seed(42)
df0_full_test <- create_engineered_features(df0_test, top_features)

x_test <- df0_full_test %>% select(all_of(feature_cols))
y_test <- df0_full_test$label

test_preds <- predict(final_model, data = x_test)$predictions[, "1"]

# ============================================================
# METRIC EVALUATION ON TEST SET
# ============================================================

roc_auc <- pROC::auc(pROC::roc(y_test, test_preds))
pr_auc  <- PRROC::pr.curve(
  scores.class0 = test_preds[y_test == 1],
  scores.class1 = test_preds[y_test == 0]
)$auc.integral

cat(sprintf("\nTest ROC-AUC: %.6f\n", roc_auc))
cat(sprintf("Test PR-AUC:  %.6f\n", pr_auc))

# Saving final model and hyperparameters
final_model <- ranger(
  dependent.variable.name = "label",
  data = train_df %>% select(label, all_of(feature_cols)),
  num.trees = recommended$ntree,
  mtry = recommended$mtry,
  min.node.size = recommended$nodesize,
  class.weights = c("0" = 1, "1" = recommended$classwt),
  sample.fraction = c("0" = recommended$majority_prop, "1" = 1),
  probability = TRUE
)

model_bundle <- list(
  model = final_model,
  params = list(
    ntree = recommended$ntree,
    mtry = recommended$mtry,
    nodesize = recommended$nodesize,
    classwt = recommended$classwt,
    majority_prop = recommended$majority_prop
  ),
  features = feature_cols
)

# hard-coded model bundle 
'# Combine model and its metadata
model_bundle <- list(
  model = final_model,
  params = list(
    ntree = 771,
    mtry = 2,
    nodesize = 1,
    classwt = 50.91,
    majority_prop = 0.44
  ),
  features = feature_cols
)
'
saveRDS(model_bundle, "final_rf_model_bundle.rds")

cat("Model saved to 'final_rf_model_bundle.rds'\n")

# --- To reload later (recommended params) ---
bundle <- readRDS("final_rf_model_bundle.rds")
final_model <- bundle$model
params <- bundle$params
feature_cols <- bundle$features


# ============================================================
# MODEL APPLICATION OF DATASET2
# ============================================================
df3_full_no_labels = read.csv("/Users/ryann_/Downloads/df1_full_no_labels.csv")
names(df3_full_no_labels) <- make.names(names(df3_full_no_labels))

new_data <- create_engineered_features(df3_full_no_labels, top_features)

preds <- predict(final_model,
                 data = new_data %>% select(all_of(feature_cols)))$predictions[, "1"]

cat("\nModel Ran Successfully.n")

submission <- data.frame(
  transcript_id = new_data$transcript,
  transcript_position = new_data$position,
  score = preds
)

write.csv(submission, file = paste0("/Users/ryann_/Downloads/", "rf_df1_predictions", ".csv"), row.names = FALSE)

write.csv(submission, "rf_df0_predictions.csv", row.names = FALSE, quote = FALSE)

cat("\nFinal model trained and test predictions saved to rf_df0_predictions.csv\n")

# ============================================================
# MODEL APPLICATION OF DATASET2
# ============================================================
new_data <- create_engineered_features(df1_full_no_labels, top_features)

preds <- predict(final_model,
                 data = new_data %>% select(all_of(feature_cols)))$predictions[, "1"]

cat("\nModel Ran Successfully.n")

submission <- data.frame(
  transcript_id = new_data$transcript,
  transcript_position = new_data$position,
  probability_score = preds
)

write.csv(submission, "rf_df1_predictions.csv", row.names = FALSE, quote = FALSE)

cat("\nFinal model trained and test predictions saved to rf_df1_predictions.csv\n")

# ============================================================
# MODEL APPLICATION OF DATASET2
# ============================================================
new_data <- create_engineered_features(df2_full_no_labels, top_features)

preds <- predict(final_model,
                 data = new_data %>% select(all_of(feature_cols)))$predictions[, "1"]

cat("\nModel Ran Successfully.n")

submission <- data.frame(
  transcript_id = new_data$transcript,
  transcript_position = new_data$position,
  probability_score = preds
)

write.csv(submission, "rf_df2_predictions.csv", row.names = FALSE, quote = FALSE)

cat("\nFinal model trained and test predictions saved to rf_df2_predictions.csv\n")
