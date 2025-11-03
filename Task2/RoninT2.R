
# Install any missing packages with dependencies
# install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)

# ======== Load Dependencies ======== 
library(dplyr)
library(jsonlite)
library(xgboost)
library(ranger)
library(randomForest)
library(readr)
library(purrr)
library(stringr)
library(data.table)
library(ggplot2)
library(tidyr)
library(UpSetR)
library(pheatmap)
library(tidytext)

# ======== Parse json Function ======== 
parse_m6anet_jsonl <- function(path) {
  message("Starting parse for:", path)
  lines <- readLines(path)
  out <- list()
  idx <- 0L
  t0 <- Sys.time()
  
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
  
  t_end <- Sys.time()
  message(sprintf("Total time: %.2f min",
                  as.numeric(difftime(t_end, t0, units = "mins"))))
  
  bind_rows(out) %>%
    select(transcript, position, kmer, everything())
  
}

# ======== Data Loading & Pre-processing ======== 
set.seed(42)

## Load Data
A549r5r1 <- parse_m6anet_jsonl("task2/SGNex_A549_directRNA_replicate5_run1.json") 
A549r6r1 <- parse_m6anet_jsonl("task2/SGNex_A549_directRNA_replicate6_run1.json") 
Hct116r3r1 <- parse_m6anet_jsonl("task2/SGNex_Hct116_directRNA_replicate3_run1.json")
Hct116r3r4 <- parse_m6anet_jsonl("task2/SGNex_Hct116_directRNA_replicate3_run4.json") 
Hct116r4r3 <- parse_m6anet_jsonl("task2/SGNex_Hct116_directRNA_replicate4_run3.json") 
HepG2r5r2 <- parse_m6anet_jsonl("task2/SGNex_HepG2_directRNA_replicate5_run2.json") 
HepG2r6r1 <- parse_m6anet_jsonl("task2/SGNex_HepG2_directRNA_replicate6_run1.json") 
K562r4r1 <- parse_m6anet_jsonl("task2/SGNex_K562_directRNA_replicate4_run1.json") 
K562r5r1 <- parse_m6anet_jsonl("task2/SGNex_K562_directRNA_replicate5_run1.json") 
K562r6r1 <- parse_m6anet_jsonl("task2/SGNex_K562_directRNA_replicate6_run1.json") 
MCF7r3r1 <- parse_m6anet_jsonl("task2/SGNex_MCF7_directRNA_replicate3_run1.json") 
MCF7r4r1 <- parse_m6anet_jsonl("task2/SGNex_MCF7_directRNA_replicate4_run1.json")

# Aggregation: Data Pre-processing Function
process_df <- function(df) {
  df %>%
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
}

# list of variable names
varnames <- c(
  "A549r5r1", "A549r6r1",
  "Hct116r3r1", "Hct116r3r4", "Hct116r4r3",
  "HepG2r5r2", "HepG2r6r1",
  "K562r4r1", "K562r5r1", "K562r6r1",
  "MCF7r3r1", "MCF7r4r1"
)

# Data Pre-processing
for (v in varnames) {
  assign(v, process_df(get(v)))
}
'
# save data in csv format
for (v in varnames) {
  df <- get(v)
  write.csv(df, file = paste0("task2/", v, ".csv"), row.names = FALSE)
}
'

#  ======== Create Random Forest Pre-processing Features ======== 
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

# List all reads
files <- c("task2/A549r5r1.csv", "task2/A549r6r1.csv",
           "task2/Hct116r3r1.csv", "task2/Hct116r3r4.csv", "task2/Hct116r4r3.csv",
           "task2/HepG2r5r2.csv", "task2/HepG2r6r1.csv",
           "task2/K562r4r1.csv", "task2/K562r5r1.csv", "task2/K562r6r1.csv",
           "task2/MCF7r3r1.csv", "task2/MCF7r4r1.csv")

# Load Random Forest model
bundle <- readRDS("../models/final_rf_model_bundle.rds")
final_model   <- bundle$model
feature_cols  <- bundle$features
top_features  <- read.csv("top_features.csv")[[1]]

# Load & Predict each read individually
pred_list <- lapply(files, function(f) {
  message("Processing: ", f)
  
  # Read data efficiently
  df <- fread(f)
  names(df) <- make.names(names(df))
  df$cell_line <- sub("task2/|\\.csv", "", f)
  
  # Feature engineering (must match training step)
  df_feat <- create_engineered_features(df, top_features)
  
  # Predict using trained RF model
  preds <- predict(final_model, 
                   data = df_feat %>% select(all_of(feature_cols)))$predictions[, "1"]
  
  df_feat$preds <- preds
  df_feat$cell_line <- df$cell_line
  df_feat
})

# Combine and process
df_all <- rbindlist(pred_list, use.names = TRUE, fill = TRUE)

# Add identifiers
df_all[, site_id := paste(transcript, position, sep=":")]
df_all[, line := sub("r[0-9]+r[0-9]+", "", cell_line)]
df_all[, replicate := cell_line]

# Save combined predictions
# names(df_all) <- make.names(names(df_all)) 
# fwrite(df_all, "RF_predictions_task2_master.csv")

cat("\n Per-cell-line predictions complete. File saved as RF_predictions_task2_master.csv\n")

# Process combined model output file
new_data = fread("RF_predictions_task2_master.csv")
names(new_data) <- make.names(names(new_data))
new_data$cell_line <- sub("\\.csv$", "", new_data$cell_line)
new_data$site_id <- paste(new_data$transcript, new_data$position, sep=":")
new_data$line <- sub("r[0-9]+r[0-9]+", "", new_data$cell_line)  
new_data$replicate <- new_data$cell_line

# List all info files
info_files <- c(
  "task2/info/SGNex_A549_directRNA_replicate5_run1.info",
  "task2/info/SGNex_A549_directRNA_replicate6_run1.info",
  "task2/info/SGNex_Hct116_directRNA_replicate3_run1.info",
  "task2/info/SGNex_Hct116_directRNA_replicate3_run4.info",
  "task2/info/SGNex_Hct116_directRNA_replicate4_run3.info",
  "task2/info/SGNex_HepG2_directRNA_replicate5_run2.info",
  "task2/info/SGNex_HepG2_directRNA_replicate6_run1.info",
  "task2/info/SGNex_K562_directRNA_replicate4_run1.info",
  "task2/info/SGNex_K562_directRNA_replicate5_run1.info",
  "task2/info/SGNex_K562_directRNA_replicate6_run1.info",
  "task2/info/SGNex_MCF7_directRNA_replicate3_run1.info",
  "task2/info/SGNex_MCF7_directRNA_replicate4_run1.info"
)

# Read, combine & process all info files
info_all <- info_files %>%
  map_dfr(function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    df$source_file <- str_remove(basename(f), ".info$")
    
    parts <- str_split(df$source_file, "_", simplify = TRUE)
    cell_line <- parts[, 2]
    replicate <- str_extract(parts[, 4], "\\d+")
    run <- str_extract(parts[, 5], "\\d+")
    
    df$line <- cell_line
    df$replicate <- replicate
    df$run <- run
    df$line_id <- paste0(cell_line, "r", replicate, "r", run)
    
    df
  })

# Standardize column names for merging
info_all <- info_all %>%
  rename(transcript = transcript_id, position = transcript_position, cell_line = line_id)

# Merge with new_data
merged_data <- new_data %>%
  left_join(info_all, by = c("transcript", "position", "cell_line"))

head(merged_data)

# Save Merged Data
# names(merged_data) <- make.names(names(merged_data)) 
# fwrite(merged_data, "RF_predictions_task2_master_2.csv")
cat("\n Per-cell-line predictions complete. File saved as RF_predictions_task2_master_2.csv\n")

# Load Full dataset
new_data = fread("RF_predictions_task2_master_2.csv")

# ======== Visualisations ======== 
# distribution of predicted score by reads
p = ggplot(new_data, aes(x = preds, color = cell_line)) +
  geom_density() +
  facet_wrap(~ line.x, scales = "free_y") +
  labs(
    x = "Predicted m6A probability",
    y = "Density",
    title = "Distribution of predicted m6A probabilities per replicate/run"
  ) +
  theme_classic()
# Since the plots are very similar to one another, we can use mean.

ggsave("plots/plot1.png", plot = p, width = 6, height = 4, dpi = 300)

# Creating df_mean
df_mean <- new_data %>%
  group_by(line.x, transcript, position, site_id) %>%
  summarize(score = mean(preds), n_rep = n(), .groups = "drop")

# Summary stats

# Identifying top 5 transcripts per cell line
top_tx <- df_mean %>%
  group_by(line.x, transcript) %>%
  summarize(mean_tx_score = median(score, na.rm=TRUE), .groups="drop") %>%
  group_by(line.x) %>%
  slice_max(mean_tx_score, n = 5)

p = ggplot(top_tx, aes(x = reorder_within(transcript, mean_tx_score, line.x),
                       y = mean_tx_score, fill = line.x)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ line.x, scales="free_y") +
  scale_x_reordered() +
  labs(x="Transcript", y="Median m6A probability",
       title="Top 10 transcripts with highest predicted modification per cell line") +
  coord_flip() +
  theme_bw()
ggsave("plots/plot1_1.png", plot = p, width = 9, height = 6, dpi = 300)

# summary table
sum_tab <- df_mean %>%
  group_by(line.x) %>%
  summarize(
    n_sites = n(),
    n_high = sum(score >= 0.6),
    frac_high = n_high / n_sites,
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    p95 = quantile(score, 0.95, na.rm = TRUE),
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    .groups = "drop"
  )

print(sum_tab)

# density plot
p = ggplot(df_mean, aes(x = score, fill = line.x, color = line.x)) +
  geom_density(alpha = 0.3) +
  labs(
    x = "Predicted m6A probability (score)",
    y = "Density",
    title = "Distribution of m6A site probabilities across cell lines"
  ) +
  theme_classic()
ggsave("plots/plot2.png", plot = p, width = 6, height = 4, dpi = 300)

# Histogram
p = ggplot(df_mean, aes(x = score, fill = line.x)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
  facet_wrap(~ line.x, scales = "free_y") +
  labs(x = "Predicted score", y = "Count") +
  theme_bw()
ggsave("plots/plot3.png", plot = p, width = 6, height = 4, dpi = 300)

# Boxplot of per-transcript mean scores
tx_tab <- df_mean %>%
  group_by(line.x, transcript) %>%
  summarize(tx_score = mean(score, na.rm = TRUE), .groups = "drop")

p = ggplot(tx_tab, aes(x = line.x, y = tx_score, fill = line.x)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    x = "Cell line",
    y = "Mean transcript-level m6A probability",
    title = "Distribution of transcript-level m6A scores"
  ) +
  theme_classic()
ggsave("plots/plot4.png", plot = p, width = 6, height = 4, dpi = 300)

# Boxplot of per-transcript median scores
tx_tab <- df_mean %>%
  group_by(line.x, transcript) %>%
  summarize(tx_score = median(score, na.rm = TRUE), .groups = "drop")

p = ggplot(tx_tab, aes(x = line.x, y = tx_score, fill = line.x)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    x = "Cell line",
    y = "Median transcript-level m6A probability",
    title = "Distribution of transcript-level m6A scores"
  ) +
  theme_classic()
ggsave("plots/plot4_2.png", plot = p, width = 6, height = 4, dpi = 300)

# Scatter plot comparing mean & median
tx_stats <- df_mean %>%
  group_by(line.x, transcript) %>%
  summarise(mean_score = mean(score, na.rm = TRUE),
            median_score = median(score, na.rm = TRUE),
            .groups = "drop")
p = ggplot(tx_stats, aes(x = mean_score, y = median_score, color = line.x)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")
ggsave("plots/plot4_3.png", plot = p, width = 6, height = 4, dpi = 300)

# Checking if modifications overlap via heatmap
thresh <- 0.6 

calls <- df_mean %>%
  group_by(line.x) %>%
  filter(score >= thresh) %>%
  summarize(mod_sites = list(unique(site_id)), .groups = "drop")

cell_lines <- calls$line.x
jaccard <- matrix(NA, nrow = length(cell_lines), ncol = length(cell_lines),
                  dimnames = list(cell_lines, cell_lines))

for (i in seq_along(cell_lines)) {
  for (j in seq_along(cell_lines)) {
    A <- calls$mod_sites[[i]]
    B <- calls$mod_sites[[j]]
    jaccard[i, j] <- length(intersect(A, B)) / length(union(A, B))
  }
}

p = pheatmap(jaccard,
             main = "Jaccard similarity of predicted m6A sites",
             color = colorRampPalette(c("white","red"))(50))
ggsave("plots/plot5.png", plot = p, width = 6, height = 4, dpi = 300)

# Median probability scores of a typical transcript
df_rel <- new_data %>%
  group_by(transcript) %>%
  mutate(rel_pos = (position - min(position)) / (max(position) - min(position))) %>%
  ungroup() %>%
  filter(!is.na(rel_pos))

df_bins <- df_rel %>%
  mutate(bin = cut(rel_pos, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)) %>%
  group_by(line.x, bin) %>%
  summarize(mean_score = mean(preds, na.rm = TRUE),
            median_score = median(preds, na.rm = TRUE),
            .groups = "drop")

p = ggplot(df_bins, aes(x = as.numeric(bin), y = median_score, color = line.x)) +
  geom_line(linewidth = 1) +
  labs(x = "Relative transcript position",
       y = "Median predicted m6A probability",
       title = "Positional trend of predicted m6A sites in a typical transcripts") +
  theme_classic()
ggsave("plots/plot6.png", plot = p, width = 6, height = 4, dpi = 300)

# Percentage of high-confidence sites (score ≥ 0.6)
p = ggplot(sum_tab, aes(x = reorder(line.x, -frac_high), y = frac_high, fill = line.x)) +
  geom_col() +
  labs(x="Cell line", y="Fraction of high-confidence sites (≥0.6)",
       title="Relative abundance of predicted m6A sites") +
  theme_classic()
ggsave("plots/plot7.png", plot = p, width = 6, height = 4, dpi = 300)
