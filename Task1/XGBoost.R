
## XGBoost

### dataset0

# Load Libraries & Set Cores
# load libraries
library(dplyr)
library(xgboost)
library(pROC)
library(PRROC)
library(data.table)
library(parallel)

num_cores <- parallel::detectCores(logical = TRUE)
Sys.setenv("OMP_NUM_THREADS" = num_cores)
message("Using ", num_cores, " CPU threads")


# Load & Prepare Datasets
site_features_full <- read.csv("site_features_full.csv")
info <- read.csv("info.csv")
df0_full <- read.csv("df0_full.csv") # from joining site_features_full and info

df0_full <- df0_full %>%
  select(-transcript, -position)


#### Cross-validation to find best hyperparameters (using df0_train)
# Split into train & test sets by gene_id
# Split by gene_id
unique_genes <- unique(df0_full$gene_id)
set.seed(42)
test_genes <- sample(unique_genes, floor(length(unique_genes) * 0.2))  # 20%

df0_train <- df0_full %>% filter(!gene_id %in% test_genes)
df0_test  <- df0_full %>% filter(gene_id %in% test_genes)

cat(sprintf("Train genes: %d | Test genes: %d\n",
            n_distinct(df0_train$gene_id), n_distinct(df0_test$gene_id)))

# drop gene_id before training
df0_train <- df0_train %>% select(-gene_id)
df0_test  <- df0_test %>% select(-gene_id)

# convert label to numeric type
df0_train$label <- as.numeric(as.character(df0_train$label))
df0_test$label  <- as.numeric(as.character(df0_test$label))

# ensure all features are numeric
x_cols <- setdiff(names(df0_train), "label")
df0_train[, x_cols] <- lapply(df0_train[, x_cols], function(x) as.numeric(as.character(x)))
df0_test[, x_cols]  <- lapply(df0_test[, x_cols],  function(x) as.numeric(as.character(x)))


# Define CV parameters & deal with class imbalance
n <- nrow(df0_train)
k <- 5 # number of folds
n_search <- 50 # number of random hyperparameter searches

# compute imbalance ratio for scale_pos_weight
prop_1 <- mean(df0_train$label == 1)
prop_0 <- 1 - prop_1
scale_pos_weight <- prop_0 / prop_1
cat(sprintf("Positive class ratio: %.3f | scale_pos_weight = %.2f\n\n", prop_1, scale_pos_weight))

# create stratified folds to prepare label proportions
set.seed(42)
folds <- df0_train %>%
  group_by(label) %>%
  mutate(fold = sample(rep(1:k, length.out = n()))) %>%
  ungroup() %>%
  pull(fold)

cat(sprintf("Total rows: %d | %d stratified folds | %d random searches\n\n",
            n, k, n_search))


# Generate random hyperparameter search grid
set.seed(42)
max_depth_vals <- sample(3:10, n_search, replace = TRUE)
eta_vals <- runif(n_search, 0.01, 0.2)
subsample_vals <- runif(n_search, 0.6, 1.0)
colsample_vals <- runif(n_search, 0.6, 1.0)
lambda_vals <- runif(n_search, 0, 2)
alpha_vals <- runif(n_search, 0, 2)
nrounds_vals <- sample(seq(100, 400, 50), n_search, replace = TRUE)
seed_vals <- sample(1:1e6, n_search)

cv_preds <- list()
cat(sprintf("Running %d random searches × %d folds\n", n_search, k))
cat(sprintf("Started at %s\n", format(Sys.time(), "%H:%M:%S")))


# Cross-Validation loop pt 1: random search loop
for (j in 1:n_search) {
  start <- Sys.time()
  
  # ensure reproducibility for each random search iteration
  iteration_seed <- seed_vals[j]
  set.seed(iteration_seed)
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = max_depth_vals[j],
    eta = eta_vals[j],
    subsample = subsample_vals[j],
    colsample_bytree = colsample_vals[j],
    lambda = lambda_vals[j],
    alpha = alpha_vals[j],
    scale_pos_weight = scale_pos_weight,
    nthread = num_cores 
  )
  
  nrounds <- nrounds_vals[j]
  cv_scores <- numeric(n)


# Cross-validation loop pt 2: Inner 5-Fold CV
  for (i in 1:k) {
    train_idx <- which(folds != i)
    test_idx  <- which(folds == i)
    
    train_fold <- df0_train[train_idx, ]
    test_fold  <- df0_train[test_idx, ]
    
    # convert labels to numeric for xgboost
    train_fold$label <- as.numeric(train_fold$label)
    test_fold$label  <- as.numeric(test_fold$label)
    
    cv_train <- xgb.DMatrix(as.matrix(train_fold %>% select(-label)),
                          label = train_fold$label)
    cv_test <- xgb.DMatrix(as.matrix(test_fold %>% select(-label)))
    
    # train XGBoost model on fold
    fit <- xgb.train(
      params = params,
      data = cv_train,
      nrounds = nrounds,
      verbose = 0
    )
    
    # store predictions
    preds <- predict(fit, cv_test)
    cv_scores[test_idx] <- preds
  }


# Cross-validation loop pt 3
  # store results for each random iteration
  cv_preds[[j]] <- data.frame(
    label = as.numeric(df0_train$label),
    score = cv_scores,
    max_depth = max_depth_vals[j],
    eta = eta_vals[j],
    subsample = subsample_vals[j],
    colsample_bytree = colsample_vals[j],
    lambda = lambda_vals[j],
    alpha = alpha_vals[j],
    nrounds = nrounds,
    iteration_seed = iteration_seed
  )
  
  elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)
  cat(sprintf("[%d/%d] depth=%d eta=%.3f nrounds=%d seed=%d done (%.1fs)\n",
              j, n_search, max_depth_vals[j], eta_vals[j], nrounds, iteration_seed, elapsed))
}

# combine all CV results
cv_results <- bind_rows(cv_preds)


# Evaluate CV performance, select best model
# Pick best hyperparameters from CV
cv_summary <- cv_results %>%
  group_by(iteration_seed, max_depth, eta, subsample, colsample_bytree, lambda, alpha, nrounds) %>%
  summarise(
    auc_roc = as.numeric(auc(roc(label, score, quiet = TRUE))),
    auc_pr = {
      if(length(unique(label)) == 2) {
        pr_curve <- PRROC::pr.curve(
          scores.class0 = score[label == 1],
          scores.class1 = score[label == 0],
          curve = TRUE
        )
        pr_curve$auc.integral
      } else NA_real_
    },
    .groups = 'drop'
  ) %>%
  arrange(desc(auc_roc))

print(cv_summary, n = Inf)

# select best-performing hyperparameters by ROC-AUC
best_by_roc <- cv_summary %>% slice_max(auc_roc, n = 1)
cat("\nBest hyperparameters selected by ROC-AUC:\n")
print(best_by_roc)

cat("\nCross-validation complete\n")


# Save best params and seeds
best_params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = best_by_roc$max_depth,
  eta = best_by_roc$eta,
  subsample = best_by_roc$subsample,
  colsample_bytree = best_by_roc$colsample_bytree,
  lambda = best_by_roc$lambda,
  alpha = best_by_roc$alpha,
  scale_pos_weight = scale_pos_weight,
  nthread = num_cores
)
best_nrounds <- best_by_roc$nrounds
iteration_seed <- best_by_roc$iteration_seed

set.seed(iteration_seed)  # ensure reproducibility


### Fit model to full dataset0 with best hyperparameters and save model

dtrain_full <- xgb.DMatrix(as.matrix(df0_full[, x_cols]),
                           label = as.numeric(df0_full$label))

full_model <- xgb.train(
  params = best_params,
  data = dtrain_full,
  nrounds = best_nrounds,
  verbose = 1
)

cat("\nFull model training complete\n")

# Save trained model as .xgb
xgb.save(full_model, "xgboost_full_model_df0.xgb")
cat("Model saved to: xgboost_full_model_df0.xgb\n")

# Generate predictions on df0_full
df0_full$score <- predict(full_model, as.matrix(df0_full[, x_cols]))

df0_full_orig <- read.csv("df0_full.csv")
df0_with_preds <- df0_full_orig %>%
  mutate(score = df0_full$score)

output_preds <- df0_with_preds %>%
  select(transcript, position, score) %>%
  rename(transcript_id = transcript, transcript_position = position)

write.csv(output_preds, "xgboost_predictions_df0.csv", row.names = FALSE)
cat("Predictions saved to: xgboost_predictions_df0.csv\n")


### Generate predictions for dataset1
# Load dataset
df1_full <- read.csv("df1_full_no_labels.csv")

# Ensure features match training set and are numeric
df1_full[, x_cols] <- lapply(df1_full[, x_cols], function(x) as.numeric(as.character(x)))

# Convert to DMatrix
dmat1 <- xgb.DMatrix(as.matrix(df1_full[, x_cols]))

# Generate predictions
df1_full$score <- predict(full_model, dmat1)

# Include transcript info if present
if ("transcript" %in% colnames(df1_full) & "position" %in% colnames(df1_full)) {
  output_preds_df1 <- df1_full %>%
    select(transcript, position, score) %>%
    rename(transcript_id = transcript,
           transcript_position = position)
} else {
  output_preds_df1 <- df1_full %>%
    mutate(row_id = row_number()) %>%
    select(row_id, score)
}

# Save predictions
write.csv(output_preds_df1, "xgboost_predictions_df1.csv", row.names = FALSE)
cat("Predictions for df1_full_no_labels saved to xgboost_predictions_df1.csv\n")


### Generate predictions for dataset2
# Load dataset
df2_full <- read.csv("df2_full_no_labels.csv")

# Ensure features match training set and are numeric
df2_full[, x_cols] <- lapply(df2_full[, x_cols], function(x) as.numeric(as.character(x)))

# Convert to DMatrix
dmat2 <- xgb.DMatrix(as.matrix(df2_full[, x_cols]))

# Generate predictions
df2_full$score <- predict(full_model, dmat2)

# Include transcript info if present
if ("transcript" %in% colnames(df2_full) & "position" %in% colnames(df2_full)) {
  output_preds_df2 <- df2_full %>%
    select(transcript, position, score) %>%
    rename(transcript_id = transcript,
           transcript_position = position)
} else {
  output_preds_df2 <- df2_full %>%
    mutate(row_id = row_number()) %>%
    select(row_id, score)
}

# Save predictions
write.csv(output_preds_df2, "xgboost_predictions_df2.csv", row.names = FALSE)
cat("Predictions for df2_full_no_labels saved to xgboost_predictions_df2.csv\n")



