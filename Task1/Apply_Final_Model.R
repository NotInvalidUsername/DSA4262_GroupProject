
bundle <- readRDS("final_rf_model_bundle.rds")
final_model <- bundle$model
params <- bundle$params
feature_cols <- bundle$features


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
