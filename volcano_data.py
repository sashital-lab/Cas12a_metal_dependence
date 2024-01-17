import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import os

# Define the directory and input file name
directory_in = 'C:/Users/sashital/Box/Manuscripts/Cas12a_metal-dependent_specificity/figure_workshop/231220_libraries/plasmid_library/nicked_input'
directory_out = 'C:/Users/sashital/Box/Manuscripts/Cas12a_metal-dependent_specificity/figure_workshop/231220_libraries/plasmid_library/nicked_input'
input_file = 'Lb_W_30_nicked_norm.txt'
data_file_path = os.path.join(directory_in, input_file)
data = pd.read_csv(data_file_path, delimiter='\t', header=None)

# Read the file containing mutant names
mutant_names_file_path = 'W_mutant_names.txt'  # Replace with your mutant names file path
mutant_names = pd.read_csv(mutant_names_file_path, header=None)

# Check if mutant_names has only one column
if mutant_names.shape[1] == 1:
    # Replace the first column of the 'data' DataFrame with the values from 'mutant_names'
    data.iloc[:, 0] = mutant_names.iloc[:, 0]

# Convert non-numeric values to NaN for specific columns
columns_to_convert = list(range(1, 4)) + list(range(10, 13))
data.iloc[:, columns_to_convert] = data.iloc[:, columns_to_convert].apply(pd.to_numeric, errors='coerce')

# Calculate fold change between columns 1-3 and columns 10-12 for each row
data['fold_change'] = data.iloc[:, 1:4].mean(axis=1)/data.iloc[:, 10:13].mean(axis=1)

# Perform a t-test for each row between the two sets of columns
p_values = []
problematic_rows = []

for index, row in data.iterrows():
    first = row.iloc[1:4]
    second = row.iloc[10:13]
    
    try:
        first_numeric = pd.to_numeric(first.dropna())
        second_numeric = pd.to_numeric(second.dropna())

        if len(first_numeric) > 1 and len(second_numeric) > 1:
            t_stat, p_value = ttest_ind(first_numeric, second_numeric)
            p_values.append(p_value)
        else:
            p_values.append(float('nan'))  # Append NaN if there are insufficient numeric values for t-test
    except ValueError as e:
        problematic_rows.append(index)
        print(f"For row {index}: ValueError -", e)
        p_values.append(float('nan'))  # Append NaN if t-test is not possible

           
data['p_values'] = p_values

log2_fold_change = []
negative_log_p = []

# Iterate through rows to calculate log2 fold change and -log10 p-values row-wise
for index, row in data.iterrows():
      
    # Calculate log2 fold change for the current row
    log2_fold_change.append(np.log2(row['fold_change']))
   
    # Calculate -log10(p-value) for the current row
    negative_log_p.append(-np.log10(row['p_values']))

# Create a DataFrame with results
result_df = pd.DataFrame({
    'mutant': mutant_names[0],
    'log2 fold change': log2_fold_change,
    '-log P': negative_log_p
})


# Create the output file name based on the input file name
output_file_name = os.path.splitext(os.path.basename(data_file_path))[0].split('_')[:4]
output_file_name = "_".join(output_file_name) + '_volcano.csv'

# Save the selected columns to a new file
output_file_path = os.path.join(directory_out, output_file_name)  # Output file path with the new name
result_df.to_csv(output_file_path, index=False)
