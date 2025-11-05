
bundle <- readRDS("../models/final_rf_model_bundle.rds")
final_model <- bundle$model
params <- bundle$params
feature_cols <- bundle$features


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

# ============================================================
# MODEL APPLICATION OF DATASET0
# ============================================================
new_data = read.csv("df0_full_no_labels.csv")

new_data <- create_engineered_features(new_data, top_features)

preds <- predict(final_model,
                 data = new_data %>% select(all_of(feature_cols)))$predictions[, "1"]

cat("\nModel Ran Successfully.n")

submission <- data.frame(
  transcript_id = new_data$transcript,
  transcript_position = new_data$position,
  score = preds
)

write.csv(submission, "rf_df0_predictions.csv", row.names = FALSE, quote = FALSE)

cat("\nFinal model trained and test predictions saved to rf_df0_predictions.csv\n")

# ============================================================
# MODEL APPLICATION OF DATASET1
# ============================================================
new_data = read.csv("df1_full_no_labels.csv")

new_data <- create_engineered_features(new_data, top_features)

preds <- predict(final_model,
                 data = new_data %>% select(all_of(feature_cols)))$predictions[, "1"]

cat("\nModel Ran Successfully.n")

submission <- data.frame(
  transcript_id = new_data$transcript,
  transcript_position = new_data$position,
  score = preds
)

write.csv(submission, "rf_df1_predictions.csv", row.names = FALSE, quote = FALSE)

cat("\nFinal model trained and test predictions saved to rf_df1_predictions.csv\n")

# ============================================================
# MODEL APPLICATION OF DATASET2
# ============================================================
new_data = read.csv("df2_full_no_labels.csv")

new_data <- create_engineered_features(new_data, top_features)

preds <- predict(final_model,
                 data = new_data %>% select(all_of(feature_cols)))$predictions[, "1"]

cat("\nModel Ran Successfully.n")

submission <- data.frame(
  transcript_id = new_data$transcript,
  transcript_position = new_data$position,
  score = preds
)

write.csv(submission, "rf_df2_predictions.csv", row.names = FALSE, quote = FALSE)

cat("\nFinal model trained and test predictions saved to rf_df2_predictions.csv\n")

# ============================================================
# MODEL APPLICATION OF DATASET3
# ============================================================
new_data = read.csv("df3_full_no_labels.csv")

new_data <- create_engineered_features(new_data, top_features)

preds <- predict(final_model,
                 data = new_data %>% select(all_of(feature_cols)))$predictions[, "1"]

cat("\nModel Ran Successfully.n")

submission <- data.frame(
  transcript_id = new_data$transcript,
  transcript_position = new_data$position,
  score = preds
)

write.csv(submission, "rf_df3_predictions.csv", row.names = FALSE, quote = FALSE)

cat("\nFinal model trained and test predictions saved to rf_df3_predictions.csv\n")
