import pandas as pd
import sys

# File path to input csv file
ot_df = pd.read_csv(sys.argv[1])
human_diseases_df = pd.read_csv('../misc/human_diseases_final.csv')
tumours_df = pd.read_csv('../misc/tumours_COSMICv100_21052024.csv')

# First merge on the 'gene' column
merged_df = pd.merge(ot_df, human_diseases_df, on='gene', how='left')
merged_df = pd.merge(merged_df, tumours_df, on='gene', how='left')

# Now merge again based on 'nearest_gene' to account for that condition
merged_df = pd.merge(merged_df, human_diseases_df[['gene', 'human_disease_association']],
                     left_on='nearest_gene', right_on='gene', how='left', suffixes=('', '_nearest'))
merged_df = pd.merge(merged_df, tumours_df[['gene', 'tumour_cosmic100_2024']],
                     left_on='nearest_gene', right_on='gene', how='left', suffixes=('', '_nearest'))

# Combine the human_disease_association columns into one column (wherever data exists)
merged_df['human_disease_association'] = merged_df['human_disease_association'].combine_first(
    merged_df['human_disease_association_nearest']
)
merged_df['tumour_cosmic100_2024'] = merged_df['tumour_cosmic100_2024'].combine_first(
    merged_df['tumour_cosmic100_2024_nearest']
)
# Drop the extra 'gene' and 'human_disease_association_nearest' columns
merged_df = merged_df.drop(columns=['human_disease_association_nearest', 'gene_nearest'])
merged_df = merged_df.drop(columns=['tumour_cosmic100_2024_nearest'])

# Save the final dataframe to a new CSV file
merged_df.to_csv(sys.argv[2], index=False)
