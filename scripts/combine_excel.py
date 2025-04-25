import pandas as pd
import os
from glob import glob

# Define the path pattern to match the uploaded CSV files
path_pattern = "../out/*gap_ot.csv.sets.csv"

# List all matching CSV files
csv_files = glob(path_pattern)

# Dictionary to hold DataFrames with modified column order
dataframes = {}

for file_path in csv_files:
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Move 'PAM_coordinate' to the third column if it exists
    if 'PAM_coordinate' in df.columns:
        cols = list(df.columns)
        cols.remove('PAM_coordinate')
        df = df[cols[:2] + ['PAM_coordinate'] + cols[2:]]
    
    # Clean up the file name to use as sheet name
    base_name = os.path.basename(file_path)
    sheet_name = base_name.replace("gap_ot.csvgap_ot.csv", "").replace(".sets.csv", "")
    
    # Ensure the sheet name is Excel-safe (max 31 characters, no special characters)
    sheet_name = sheet_name[:31]
    
    # Store the modified DataFrame
    dataframes[sheet_name] = df

# Write all DataFrames to a single Excel file
output_excel_path = "SWOfftarget_output.xlsx"
with pd.ExcelWriter(output_excel_path) as writer:
    for sheet_name, df in dataframes.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)

output_excel_path
