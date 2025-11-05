This readme provides instructions on how to run our code. 
```
Schema:
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

There are two ways that you can run the code
  1) Running the code and running all functions organically which will take some time to train and run the model.
  2) Use the data available on sharepoint where the output of lengthy processes can be found.
  sharepoint link: https://nusu-my.sharepoint.com/:f:/g/personal/e0959087_u_nus_edu/Ekxtlk19TDBCs5FGJIc3s_0BNioPh-XZW0uqs6Stvlywfw
     

## Task 1: Training the model
.
.
.
.


## Task 2: Applying the model on various cancer cell lines
1) run the `data_download.sh` to download all the necessary data files to your current working directory. (All data are available in CSV format in the sharepoint folder)
