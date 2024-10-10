# brGMM_BAE model

BrGMM_BAE  is a clustering model to predict lake pH and salinity based on brGDGT compounds.

## Project Overview

This project utilizes brGDGT compounds found in lakes to predict bacterial assemblages using a GMM (Gaussian Mixture Model).

## Features

- Predicts lake bacterial assemblages using brGDGT compounds.
- Estimates lake pH and salinity both in modern and paleo.

## Installation Guide

1. **Clone the repository**

   ```bash
   git clone https://github.com/yourusername/brGMM_BAE.git
   cd brGMM_BAE

2. **Install dependencies**

   Ensure you have `joblib` and `pandas` installed. If not, install them using the following command:

   ```bash
   pip install joblib pandas
   ```

3. **Load the model and make predictions (Python)**

   Run the following code in your Python environment:

   ```python
   import joblib
   import pandas as pd
   
   # Load the GMM model
   loaded_gmm = joblib.load('brGMM-BAE.pkl')
   
   # Read the data
   data = pd.read_excel('Example.xlsx', sheet_name='Sheet1')
   data_X = data[["MBT'5me", 'IR']]
   
   # Predict
   predicted_labels = loaded_gmm.predict(data_X)
   data['Predicted_Labels'] = predicted_labels + 1
   
   data['Bacterial Cluster'] = ''
   data.loc[TP_data['Predicted_Labels'] == 1, 'Bacterial Cluster'] = 'Halo-alkalophylic Species'
   data.loc[data['Predicted_Labels'] == 2, 'Bacterial Cluster'] = 'Freshwater Species'
   ```

## GUI Version

For those unfamiliar with Python, we create a graphical user interface (GUI) that allows you to make predictions without needing to install or run Python code.

![alt 属性文本](https://github.com/jieliangbio/brGMM_BAE/blob/main/Images/brGMMapp.png)

## Usage Instructions

1. If using Python, follow the installation guide to set up all dependencies and use the provided code snippet to load the model and data, then perform predictions.
2. If using the Windows GUI, simply download and run the GUI to make predictions without needing any additional setup.

## Contribution Guidelines

Contributions are welcome! If you want to contribute to this project, please fork the repository, create a new branch for your changes, and submit a Pull Request.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
```

Please replace the repository URL `https://github.com/yourusername/brLCP.git` and the Excel file name `XXX.xlsx` with the actual values. Let me know if there are any other details you would like to add or modify!
