import pandas as pd
import sys

original_file = sys.argv[1]
ot_df_file = sys.argv[2]
original_df = pd.read_csv(original_file)
ot_df = pd.read_csv(ot_df_file)

human_diseases_df = pd.read_csv('../misc/human_diseases_final.csv')
tumours_df = pd.read_csv('../misc/tumours_COSMICv100_21052024.csv')
canine_osteosarcoma_df = pd.read_csv('../misc/canine_osteosarcoma_genes_Gardner2019ComBiol.csv')

# Convert 'gene' column to uppercase
for df in [ot_df, human_diseases_df, tumours_df, canine_osteosarcoma_df]:
    df['gene'] = df['gene'].str.upper()
for df in [ot_df]:
    df['nearest_gene'] = df['nearest_gene'].str.upper()

# Merge based on the 'gene' column
merged_df = pd.merge(ot_df, human_diseases_df, on='gene', how='left')
merged_df = pd.merge(merged_df, tumours_df, on='gene', how='left')
merged_df = pd.merge(merged_df, canine_osteosarcoma_df, on='gene', how='left')

# Merge again based on 'nearest_gene' to account for that condition
merged_df = pd.merge(merged_df, human_diseases_df[['gene', 'human_disease_association']],
                     left_on='nearest_gene', right_on='gene', how='left', suffixes=('', '_nearest'))
merged_df = pd.merge(merged_df, tumours_df[['gene', 'tumour_cosmic100_2024']],
                     left_on='nearest_gene', right_on='gene', how='left', suffixes=('', '_nearest'))
merged_df = pd.merge(merged_df, canine_osteosarcoma_df[['gene', 'canine_osteosarcoma_gardner2019']],
                     left_on='nearest_gene', right_on='gene', how='left', suffixes=('', '_nearest'))

# Combine disease association columns where data exists
merged_df['human_disease_association'] = merged_df['human_disease_association'].combine_first(
    merged_df['human_disease_association_nearest']
)
merged_df['tumour_cosmic100_2024'] = merged_df['tumour_cosmic100_2024'].combine_first(
    merged_df['tumour_cosmic100_2024_nearest']
)
merged_df['canine_osteosarcoma_gardner2019'] = merged_df['canine_osteosarcoma_gardner2019'].combine_first(
    merged_df['canine_osteosarcoma_gardner2019_nearest']
)

# Drop unnecessary columns
merged_df = merged_df.drop(columns=['human_disease_association_nearest', 'tumour_cosmic100_2024_nearest',
                                     'canine_osteosarcoma_gardner2019_nearest'])
merged_df.columns = merged_df.columns.str.replace(r"^mcols\.", "", regex=True)
merged_df.columns = merged_df.columns.str.replace(r"^X\.", "", regex=True)

# Account for difference in coords
original_df['EndPosition'] = original_df.apply(
    lambda row: row['EndPosition'] - 5 if row['Strand'] == '+' 
    else row['EndPosition'] - 17 if row['Strand'] == '-' 
    else row['EndPosition'], axis=1
)

merged_df = merged_df.drop(columns=['gene_nearest', 'start', 'width', 'strand', 'SiteSeqPlusMaxEditsBefore'])
original_df = original_df.drop(columns=['X.Bulges', 'X.Mismatches', 'X.Edit', 'SiteSeqPlusMaxEditsBefore'])

merged_df = merged_df.rename(columns={
    'seqnames': 'Chromosome',
    'end': 'Cut_coordinate',
    'Strand': 'Protospacer_strand'
})

original_df = original_df.rename(columns={
    'Strand': 'Protospacer_strand',
    'EndPosition': 'Cut_coordinate'
})

final_df = pd.merge(original_df, merged_df, on=['Chromosome', 'Cut_coordinate', 'Protospacer_strand', 'AlignedTarget', 'AlignedText', 'score'])

final_df['PAM_coordinate'] = final_df.apply(
    lambda row: row['Cut_coordinate'] + 3 if row['Protospacer_strand'] == '+' else row['Cut_coordinate'] - 3,
    axis=1
)


# Rename 'score' to 'Gap_aware_CFD'
final_df = final_df.rename(columns={'score': 'Gap_aware_CFD'})

# Reorder columns: move 'PAM_coordinate' to the third column
cols = list(final_df.columns)
pam_idx = cols.index('PAM_coordinate')
# Remove PAM_coordinate and insert at index 2
cols.insert(2, cols.pop(pam_idx))
final_df = final_df[cols]

# Save the final dataframe
final_df.to_csv(ot_df_file + ".sets.csv", index=False)

