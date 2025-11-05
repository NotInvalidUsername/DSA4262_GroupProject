
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
new_data = fread("RF_predictions_task2_master_2.csv") %>%
  select(-line.y, -replicate.x) %>%
  rename(line = line.x, replicate = replicate.y)

# ======== Visualisations ======== 
## Getting an idea of the data
# distribution of predicted score by reads
p = ggplot(new_data, aes(x = preds, color = cell_line)) +
  geom_density() +
  facet_wrap(~ line, scales = "free_y") +
  labs(
    x = "Predicted m6A probability",
    y = "Density",
    title = "Distribution of predicted m6A probabilities per replicate/run"
  ) +
  theme_classic()

ggsave("plots/plot1.png", plot = p, width = 6, height = 4, dpi = 300)

# ======== Distribution of sites per transcript ========
p = new_data %>%
  group_by(line, transcript) %>%
  summarise(n_positions = n(), .groups = "drop") %>%
  ggplot(aes(x = n_positions, fill = n_positions < 10)) +
  geom_histogram(bins = 50, color = "black") +
  scale_fill_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "tomato"),
    labels = c("FALSE" = "≥ 10 sites", "TRUE" = "< 10 sites")
  ) +
  labs(
    x = "Number of sites per transcript",
    y = "Count",
    fill = "Category",
    title = "Distribution of sites per transcript across cell lines"
  ) +
  facet_wrap(~ line, scales = "free_y") +
  theme_classic()
ggsave("plots/plot2.png", plot = p, width = 8, height = 4, dpi = 300)

# Filter out transcripts with less than 10 data points 
new_data <- new_data %>% 
  group_by(line, transcript) %>%
  filter(n() >= 10) %>%
  ungroup()

# Distribution of data points according to position
p = ggplot(new_data, aes(x = position)) +
  geom_density(fill = "steelblue") +
  labs(
    x = "Position",
    y = "Count of data points",
    title = "Distribution of data point counts by position across cell lines"
  ) +
  facet_wrap(~ line, scales = "free_y") +
  theme_classic()
ggsave("plots/plot3.png", plot = p, width = 6, height = 4, dpi = 300)

# Scatter plot comparing mean & median
tx_stats <- new_data %>%
  group_by(line, position) %>%
  summarise(mean_score = mean(preds, na.rm = TRUE),
            median_score = median(preds, na.rm = TRUE),
            .groups = "drop")

p = ggplot(tx_stats, aes(x = mean_score, y = median_score, color = line)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")
ggsave("plots/plot4.png", plot = p, width = 6, height = 4, dpi = 300)


df_tx <- new_data %>%
  group_by(line, transcript) %>%
  summarise(mid_pred = median(preds, na.rm = TRUE), .groups = "drop")

p = ggplot(df_tx, aes(x = line, y = mid_pred, fill = line)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    x = "Cell line",
    y = "Median transcript-level m6A probability",
    title = "Distribution of transcript-level m6A scores"
  )
ggsave("plots/plot5.png", plot = p, width = 6, height = 4, dpi = 300)

threshold <- 0.6 

rel_abundance <- new_data %>%
  group_by(line) %>%
  summarise(
    total_sites = n(),
    high_conf_sites = sum(preds >= threshold, na.rm = TRUE)
  ) %>%
  mutate(
    rel_abundance = high_conf_sites / total_sites
  )

p = ggplot(rel_abundance, aes(x = line, y = rel_abundance, fill = line)) +
  geom_col(color = "black") +
  geom_text(aes(label = scales::percent(rel_abundance, accuracy = 0.1)),
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Cell line",
    y = "Relative abundance of high-confidence m6A sites",
    title = "Proportion of high-confidence predicted m6A sites per cell line"
  ) +
  theme_classic() +
  theme(legend.position = "none")
ggsave("plots/plot6.png", plot = p, width = 6, height = 4, dpi = 300)


## This is where results start
tx_median <- new_data %>%
  group_by(line, position) %>%
  summarise(mid_pred = median(preds, na.rm = TRUE), .groups = "drop")

tx_ranked <- tx_median %>%
  group_by(line) %>%
  mutate(
    top_cut = quantile(mid_pred, 0.95, na.rm = TRUE),
    bot_cut = quantile(mid_pred, 0.05, na.rm = TRUE),
    group = case_when(
      mid_pred >= top_cut ~ "Highly modified",
      mid_pred <= bot_cut ~ "Lowly modified",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  select(line, position, group)

df_binned <- new_data %>%
  inner_join(tx_ranked, by = c("line", "position")) %>%
  mutate(bin = cut(position, breaks = 100)) %>%
  group_by(line, group, bin) %>%
  summarise(
    mid_pred = median(preds, na.rm = TRUE),
    mid_pos = median(as.numeric(sub("\\((.+),.*", "\\1", bin))),
    .groups = "drop"
  )

df_summary_all <- df_binned %>%
  group_by(line, group, bin) %>%
  summarise(
    mid_pred = median(mid_pred, na.rm = TRUE),
    mid_pos = median(mid_pos, na.rm = TRUE),
    .groups = "drop"
  )

p = ggplot(df_summary_all, aes(x = mid_pos, y = mid_pred, color = group)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ line, ncol = 3, scales = "free_y") +
  labs(
    x = "Binned position",
    y = "Median predicted m6A probability",
    title = "Top vs. bottom 5% modification profiles by cell line"
  ) +
  theme_classic()
ggsave("plots/plot7.png", plot = p, width = 6, height = 4, dpi = 300)

# Heatmappie map
high_sites <- new_data %>%
  filter(preds >= threshold) %>%
  select(line, position) %>%
  distinct()

site_matrix <- high_sites %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = line, values_from = value, values_fill = 0)

lines <- colnames(site_matrix)[-1]
jaccard_mat <- matrix(NA, nrow = length(lines), ncol = length(lines),
                      dimnames = list(lines, lines))

for (i in lines) {
  for (j in lines) {
    a <- site_matrix[[i]]
    b <- site_matrix[[j]]
    intersection <- sum(a & b)
    union <- sum(a | b)
    jaccard_mat[i, j] <- ifelse(union == 0, NA, intersection / union)
  }
}

jaccard_df <- reshape2::melt(jaccard_mat, varnames = c("Line1", "Line2"), value.name = "Jaccard")

p = ggplot(jaccard_df, aes(x = Line1, y = Line2, fill = Jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 3.5, color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey90") +
  labs(
    x = "",
    y = "",
    fill = "Jaccard\nsimilarity",
    title = paste("Jaccard similarity of high-confidence (≥", threshold, ") m6A positions")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
ggsave("plots/plot8.png", plot = p, width = 6, height = 4, dpi = 300)
