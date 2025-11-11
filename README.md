# Project Overview

This repository contains all scripts, data, and models used for m6A site prediction and analysis across multiple cancer cell lines.
It includes two main tasks:
- Task 1 – Model Training (Random Forest and XGBoost)
- Task 2 – Applying Trained Models and Visualisation

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
- Full Run: Execute all scripts sequentially to reproduce model training and data processing (takes longer).
- Quick Run: Use the preprocessed data and pretrained models available on SharePoint.

- Please ensure to select Machine: **M5A.12XLARGE** (Ronin)
- **SharePoint link**: [Click here](https://nusu-my.sharepoint.com/:f:/g/personal/e0959087_u_nus_edu/Ekxtlk19TDBCs5FGJIc3s_0BNioPh-XZW0uqs6Stvlywfw)
- **Please refer to the demo in the video presentation if required, the demo has more step by step details if you are unsure of how to run the required methods**
     
## To install R
Upon entering your server, run these lines.

```
sudo apt update
```
```
sudo apt install r-base
```
It will take some time. If there is a pop-up, press Enter

## Task 1: Training the model
**Goal**: Train and evaluate Random Forest and XGBoost models on parsed nanopore m6A data.

Prerequisites
- Language: R (≥ 4.3)
- Packages: 
```
install.packages("jsonlite")
install.packages("dplyr")
install.packages("pROC")
install.packages("PRROC")
install.packages("ranger")
install.packages("randomForest")
install.packages("ggplot2")
install.packages("xgboost")
install.packages("readr")
install.packages("tibble")
```
- Working directory:
  
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
- site_features_full.csv
- info.csv
- train_df_full.csv
- df1_full_no_labels.csv
- df2_full_no_labels.csv
- df3_full_no_labels.csv

2.	Train Model:
```
source("Final_Model_Train.R")
```
Outputs saved in project/models/:
- final_rf_model_bundle.rds

3. Apply Model:
```
source("Apply_Final_Model.R")
```
Outputs saved in project/data/:
- Predicted scores for df0–df3

#### Option 2 – Pretrained Models
Skip training and use the saved model on `Apply_Final_Model.R`
- ../models/final_rf_model_bundle.rds


## Task 2: Applying the model on various cancer cell lines

**Goal**: Apply the trained Random Forest model to predict m6A modification probabilities across multiple cancer cell lines and generate visualisations.

Prerequisites
- Script: Task2.R
- Language: R
- Working directory:	
```
setwd("project/data")
```
	

### Run Instructions
#### Option 1 - Full Run
1.	Run the shell script to download all necessary data:
```
bash data_download.sh
```
(All required files are also on SharePoint.)

2.	Parse JSONs for all cell lines into CSV format using:
```
source("Task2.R")
```
This step:
- Reads each cell line JSON file (Uses a lot of ram. Might have to parse 1 by 1)
- Converts it into CSV format
- Uses the pretrained Random Forest model (final_rf_model_bundle.rds) for prediction on each dataset
3.	The script saves:
- Individual prediction outputs per cell line (*.csv)
- A concatenated master dataset (RF_predictions_task2_master_2.csv)
	
4.	The master dataset (RF_predictions_task2_master_2.csv) is then used for:
- Data analysis
- Generating visualisations in data/plots/ (distribution plots, median/mean prediction trends, summary statistics)

#### Option 2:
Skip json file parsing and skip to *line 260*

#### Outputs
	
- Predictions:
CSVs for each cancer cell line in data/task2/
- Master dataset:
data/RF_predictions_task2_master_2.csv
- Plots:
Generated in data/plots/

