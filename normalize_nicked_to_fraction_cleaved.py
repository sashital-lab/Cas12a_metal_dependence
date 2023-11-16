# This script will take nicked values that were normalized to control and normalize against fraction nicked
# Be sure to specify the directories and files below!

import pandas as pd
import numpy as np
import os

# Define the directory and file names
input_dir = '' # specify input directory
input_file = '' # specify input file name, this should be output containing all replicates for all Mg2+ concnetrations from normalize_nicked_to_control.py
norm_file = '' # specify normalization factors file name i.e. average fraction nicked at each Mg concentration for the given library/Cas12a ortholog

# Load the input file into the pandas df
input_path = os.path.join(input_dir, input_file)
df = pd.read_csv(input_path, sep='\t', header=None)

# Split the input filename based on underscores
filename_parts = input_file.split('_')

# Remove the file extension
input_name = os.path.splitext(input_file)[0]

# Load the normalization factors into a separate dataframe
norm_path = os.path.join(input_dir, norm_file)
norm_df = pd.read_csv(norm_path, sep='\t', header=None)

# Find the row corresponding to the input file name in norm_df
norm_factors = norm_df.loc[norm_df[0] == input_name].iloc[0, 1:]

# Multiply the values in df by the corresponding normalization factors
for col in range(1, 13):
    df.loc[:, col] = df[col] * norm_factors[col]

# Add column names to df
df.columns = ['Sequence', 1, 1, 1, 2, 2, 2, 5, 5, 5, 10, 10, 10]

# Get sequences where there are at least 1 value
valid_seqs = df.notnull().any(axis=1)

# Get the valid subset of the dataframe
df_valid = df.loc[valid_seqs].copy()

# Calculate the average values for each Mg concentration, ignoring NaN cells
df_avg = pd.DataFrame(columns=['Sequence', 1, 2, 5, 10])
for index, row in df_valid.iterrows():
    avg_vals = []
    for i in [1, 2, 5, 10]:
        vals = [x for x in row.values[1:] if not pd.isna(x) and int(df_valid.columns.values[row.values.tolist().index(x)]) == i]
        if vals:
            avg = sum(vals)/len(vals)
        else:
            avg = np.nan
        avg_vals.append(avg)
    df_new_row = pd.DataFrame({'Sequence': row['Sequence'], 1: avg_vals[0], 2: avg_vals[1], 5: avg_vals[2], 10: avg_vals[3]}, index=[0])
    df_avg = pd.concat([df_avg, df_new_row], ignore_index=True)

# Define the output file paths
output_norm_file = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + '_norm.txt')
output_all_file = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + '_all.txt')
output_means_file = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + '_means.txt')
output_1mM_file = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + '_1mM.txt')
output_2mM_file = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + '_2mM.txt')
output_5mM_file = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + '_5mM.txt')
output_10mM_file = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + '_10mM.txt')

# Output the final dataframes to tab-delimited files
df_valid[['Sequence', 1, 2, 5, 10]].to_csv(output_norm_file, sep='\t', index=False, header=False)
df_avg.to_csv(output_all_file, sep='\t', index=False, header=False)
df_avg.to_csv(output_means_file, sep='\t', index=False, header=False)
df_avg[['Sequence', 1]].to_csv(output_1mM_file, sep='\t', index=False, header=False)
df_avg[['Sequence', 2]].to_csv(output_2mM_file, sep='\t', index=False, header=False)
df_avg[['Sequence', 5]].to_csv(output_5mM_file, sep='\t', index=False, header=False)
df_avg[['Sequence', 10]].to_csv(output_10mM_file, sep='\t', index=False, header=False)
