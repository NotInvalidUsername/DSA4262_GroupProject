library(jsonlite)
library(dplyr)

##Info data
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

##Dataset1
dataset1 <- parse_m6anet_jsonl_full(
  "dataset1.json.gz"
)

### Aggregate dataset1 
dataset_1 <- dataset1 %>% #calculating average for each transcript and position
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




##Dataset2
dataset2 <- parse_m6anet_jsonl_full(
  "dataset2.json.gz"
)

### Aggregate dataset2
dataset_2 <- dataset2 %>% #calculating average for each transcript and position
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



## Export

# Dataset0
write.csv(train_df_full, 
          file = "train_df_full.csv", 
          row.names = FALSE)

cat("✓ Exported to: train_df_full.csv\n")

#Dataset1
write.csv(dataset_1, 
          file = "df1_full_no_labels.csv", 
          row.names = FALSE)

cat("✓ Exported to: df1_full.csv\n")

#Dataset2
write.csv(dataset_2, 
          file = "df2_full_no_labels.csv", 
          row.names = FALSE)

cat("✓ Exported to: df2_full.csv\n")

#Site features full
write.csv(site_features_full, 
          file = "site_features_full.csv", 
          row.names = FALSE)

cat("✓ Exported to: site_features_full.csv\n")

#Info
write.csv(info, 
          file = "info.csv", 
          row.names = FALSE)

cat("✓ Exported to: info.csv\n")

