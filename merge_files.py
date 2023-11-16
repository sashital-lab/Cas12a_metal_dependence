# This script is to merge the output files containing all replicates for each individual Mg2+ concentration

import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) < 5:
    print("Usage: python script.py input_file1 input_file2 input_file3 input_file4 output_file")
    sys.exit(1)

# Get input file paths and output file path from command-line arguments
input_file_paths = sys.argv[1:5]
output_file_path = sys.argv[5]

# Initialize an empty list to store merged data
merged_data = []

# Read data from the first file and keep only the first column
with open(input_file_paths[0], 'r') as file:
    for line in file:
        columns = line.strip().split('\t')
        merged_data.append([columns[0]])  # Keep only the first column

# Read data from the remaining files and merge columns based on the length of columns
for input_file_path in input_file_paths[0:]:
    with open(input_file_path, 'r') as file:
        for idx, line in enumerate(file):
            columns = line.strip().split('\t')
            # Determine how to merge columns based on the length of columns
            if len(columns) == 4:
                merged_columns = columns[-3:]
            elif len(columns) == 3:
                merged_columns = columns[-2:] + ['']  # Add an empty string for the missing column
            elif len(columns) == 2:
                merged_columns = columns[-1:] + [''] * 2  # Add empty strings for the missing columns
            elif len(columns) == 1:
                merged_columns = [''] * 3  # Add three empty strings for missing columns
            else:
                merged_columns = [''] * 3  # Default to three empty strings if len(columns) > 4
            # Add the merged columns to the corresponding row in merged data
            merged_data[idx].extend(merged_columns)

# Write merged data to the specified output file
with open(output_file_path, 'w') as output_file:
    for row in merged_data:
        # Join columns with tabs and write to the output file
        output_file.write('\t'.join(map(str, row)) + '\n')

print(f"Merged data has been written to {output_file_path}")
