#Info data
info <- read.csv("data.info.labelled", stringsAsFactors = FALSE)
head(info)

##Dataset0
parse_m6anet_jsonl_full <- function(path) {
  lines <- readLines(path)  # read the full file
  out <- list()
  idx <- 0L
  
  for (ln in lines) {
    dat <- fromJSON(ln, simplifyVector = FALSE)
    
    for (tx in names(dat)) {
      pos_list <- dat[[tx]]
      for (pos in names(pos_list)) {
        kmer_list <- pos_list[[pos]]
        for (kmer in names(kmer_list)) {
          reads <- kmer_list[[kmer]]
          if (length(reads) == 0) next
          mat <- do.call(rbind, reads)
          df <- as.data.frame(mat, stringsAsFactors = FALSE)
          
          names(df) <- c(
            "dwelling_time_p-1","stdev_p-1","mean_current_p-1",
            "dwelling_time_p0","stdev_p0","mean_current_p0",
            "dwelling_time_p+1","stdev_p+1","mean_current_p+1"
          )
          
          df$transcript <- tx
          df$position   <- as.integer(pos)
          df$kmer       <- kmer
          
          idx <- idx + 1L
          out[[idx]] <- df
        }
      }
    }
  }
  
  bind_rows(out) %>%
    select(transcript, position, kmer, everything())
}

df_full <- parse_m6anet_jsonl_full(
  "dataset0.json.gz"
)

# Aggregate dataset0
site_features_full <- df_full %>% #calculating average for each transcript and position
  group_by(transcript, position) %>%
  mutate(
    across(
      c(starts_with("dwelling_time"), starts_with("stdev"), starts_with("mean_current")),
      as.numeric
    )
  ) %>%
  summarise(
    across(starts_with("dwelling_time"), \(x) mean(x, na.rm = TRUE)),
    across(starts_with("stdev"),        \(x) mean(x, na.rm = TRUE)),
    across(starts_with("mean_current"), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

# Left join with info to get the labels
train_df_full <- site_features_full %>%
  left_join(info,
            by = c("transcript" = "transcript_id",
                   "position"   = "transcript_position"))

write.csv(train_df_full,
          file = "train_df_full.csv",
          row.names = FALSE)

df0_full <- read.csv("train_df_full.csv") # from joining site_features_full and info

# Load Libraries & Model
library(dplyr)
library(pROC)
library(PRROC)
library(ranger)

bundle <- readRDS("final_rf_model_bundle.rds")

final_model <- bundle$model
params <- bundle$params
feature_cols <- bundle$features

n_cores <- parallel::detectCores()
cat("Detected", n_cores, "cores — using all for ranger.\n")
options(ranger.num.threads = n_cores)

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

near_zero_var <- function(df, freq_cut = 95/5, unique_cut = 10) {
  nzv <- logical(ncol(df))
  names(nzv) <- names(df)
  
  for (i in seq_along(df)) {
    x <- df[[i]]
    if (is.numeric(x) || is.factor(x)) {
      tab <- table(x, useNA = "no")
      if (length(tab) <= 1) {
        nzv[i] <- TRUE
      } else {
        freq_ratio <- max(tab) / sort(tab, decreasing = TRUE)[2]
        percent_unique <- (length(tab) / length(x)) * 100
        nzv[i] <- freq_ratio > freq_cut & percent_unique < unique_cut
      }
    } else {
      nzv[i] <- FALSE
    }
  }
  
  data.frame(
    nzv = nzv,
    row.names = names(df)
  )
}

library(randomForest)

feature_cols_base <- df0_train %>%
  select(-label, -transcript, -position, -gene_id) %>%
  select(where(is.numeric)) %>%
  names()

# Remove near-zero variance features
x_temp <- df0_train %>% select(all_of(feature_cols_base))
nzv <- near_zero_var(x_temp)
feature_cols_filtered <- feature_cols_base[!nzv$nzv]

cat(sprintf("  Starting features: %d\n", length(feature_cols_base)))
cat(sprintf("  After NZV removal: %d\n", length(feature_cols_filtered)), "\n")

if (any(nzv$nzv)) {
  removed <- names(nzv$nzv)[nzv$nzv]
  cat("  Removed near-zero variance features:\n    ", paste(removed, collapse = ", "), "\n")
}

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

# Function to create engineered features
create_engineered_features <- function(df, top_features) {
  # Add interaction features
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
  
  # Add squared terms (ONLY ONCE)
  for (i in 1:min(5, length(top_features))) {
    feat <- top_features[i]
    new_feat <- paste0(feat, "_sq")
    
    if (feat %in% names(df) && !endsWith(feat, "_sq")) {  # Don't square already squared features
      df[[new_feat]] <- df[[feat]]^2
    }
  }
  
  return(df)
}

# Select only valid predictors (exclude leakage columns)
feature_cols_base <- df0_full %>%
  select(-label, -transcript, -position, -gene_id) %>%
  select(where(is.numeric)) %>%
  names()

df0_full = create_engineered_features(df0_full, top_features)

feature_cols <- df0_full %>%
  select(-label, -transcript, -position, -gene_id) %>%
  names()

set.seed(42)
new_test_set <- create_engineered_features(test_set, top_features)

x_test <- new_test_set %>% select(all_of(feature_cols))
y_test <- new_test_set$label

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
