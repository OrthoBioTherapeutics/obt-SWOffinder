import pandas as pd
import os

# Get a list of all CSV files in the current directory
csv_files = [f for f in os.listdir('.') if f.endswith('.csv')]

for file in csv_files:
    try:
        df = pd.read_csv(file)
        if 'score' not in df.columns:
            continue
        score_sum = df['score'].sum()
        result = 100/(score_sum + 100)
        print(f"'{file}': {result}")
    
    except Exception as e:
        print(f"An error occurred with file '{file}': {e}")
