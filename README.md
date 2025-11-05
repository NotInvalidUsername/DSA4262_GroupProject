# Project Overview

This repository contains all scripts, data, and models used for m6A site prediction and analysis across multiple cancer cell lines.
It includes two main tasks:
	1.	Task 1 – Model Training (Random Forest and XGBoost)
	2.	Task 2 – Applying Trained Models and Visualisation

# Schema:
```
project/
├── data/                       # **This is the main Working Directory**
│   ├── plots/                  # All visualisations for Task 2
│   ├── task2/                  # All datasets for Task 2
│   │   ├── index/              # Transcript index files
│   │   ├── info/               # Transcript information files
│   │   ├── readcount/          # Read count data per transcript
│   │   └── *.csv, *.json ...   # Processed datasets for Task 2
│   └── task1/                  # (If exists) Raw or processed data for Task 1
│
├── models/                     # Trained machine learning models
│   └── *.xgb, *.rds, *.pkl     # Saved model files
│
└── README.md                   # Project documentation
```

# Running the Code

There are two ways to run this project:
	1.	Full Run: Execute all scripts sequentially to reproduce model training and data processing (takes longer).
	2.	Quick Run: Use the preprocessed data and pretrained models available on SharePoint.

SharePoint link: https://nusu-my.sharepoint.com/:f:/g/personal/e0959087_u_nus_edu/Ekxtlk19TDBCs5FGJIc3s_0BNioPh-XZW0uqs6Stvlywfw
     

## Task 1: Training the model
**Goal**: Train and evaluate Random Forest and XGBoost models on parsed nanopore m6A data.

Prerequisites
	•	Language: R (≥ 4.3)
	•	Packages: dplyr, ggplot2, xgboost, randomForest, readr, tibble
	•	Working directory:
  
```
setwd("project/data")
```

### Input Files

Located in project/data:
```
dataset0.json.gz
dataset1.json.gz
dataset2.json.gz
dataset3.json.gz
```

### Run Instructions
#### Option 1 – Full Run (train models from scratch)
1.	Parse datasets:
```
source("Data_Preparation.R")
```
Generates:
	•	site_features_full.csv
	•	info.csv
	•	train_df_full.csv
	•	df1_full_no_labels.csv
	•	df2_full_no_labels.csv
	•	df3_full_no_labels.csv

2.	Train Model:
```
source("Final_Model_Train.R")
```
Outputs saved in project/models/:
	•	final_rf_model_bundle.rds

3. Apply Model:
```
source("Apply_Final_Model.R")
```
Outputs saved in project/data/:
	•	Predicted scores for df0–df3

#### Option 2 – Pretrained Models
Skip training and use the saved model on `Apply_Final_Model.R`
	•	../models/final_rf_model_bundle.rds


## Task 2: Applying the model on various cancer cell lines

**Goal**: Apply pretrained models to predict m6A modification probabilities and visualize transcriptomic distributions.

Steps
	1.	Run the shell script to download all necessary data:
```
bash data_download.sh
```


