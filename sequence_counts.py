# This script is to merge the output of normalize_to_counts.py, which is used to normalize the nicked data to the controls

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

[sashital@biocrunch6 230119_libraries_DNAfacility]$ more sequence_counts.py
from collections import defaultdict
import sys
import os

# Define the reference sequences
ref_sequences = {
    'L': "TTTGTACTGTCCACGCCGACGGAA",  # L reference sequence
    'W': "TTTACGGCCACTTCCGTGTCTGAC"  # W reference sequence
}

# Read input file names from the command line arguments
input_files = sys.argv[1:]

# Process each input file
for input_file in input_files:
    file_parts = input_file.split('_')
    file_prefix = file_parts[1]
    if file_prefix in ref_sequences:
        reference_sequence = ref_sequences[file_prefix]

# Generate potential sequences with 1 or 2 mutations
potential_sequences = set()

# Single nucleotide mutations
for i in range(len(reference_sequence)):
    for base in "ACGT":
        mutated_sequence = list(reference_sequence)
        mutated_sequence[i] = base
        potential_sequences.add("".join(mutated_sequence))

# Single nucleotide deletions
for i in range(len(reference_sequence)):
    deleted_sequence = list(reference_sequence)
    del deleted_sequence[i]
    potential_sequences.add("".join(deleted_sequence))

# Double nucleotide mutations
for i in range(len(reference_sequence)):
    for j in range(i + 1, len(reference_sequence)):
        for base1 in "ACGT":
            for base2 in "ACGT":
                mutated_sequence = list(reference_sequence)
                mutated_sequence[i] = base1
                mutated_sequence[j] = base2
                potential_sequences.add("".join(mutated_sequence))

# Add spike-in sequences
spikein_sequences = [
    "GCATTCTTATTAATACATTTGAAA",
    "CGCGCCCAACTGACGCTAGGCAAG",
    "TAAGGGTAAACATACAAGTCGATA",
    "TCAGTGCAGGCTCCCGTGTTAGGA",
    "CATCCAGCACTCTACGGTTCCTCC"
]
potential_sequences.update(spikein_sequences)

# Create a dictionary to store counts for each potential sequence
sequence_counts = {seq: [0, 0, 0] for seq in potential_sequences}  # 3 columns for 3 files


# Iterate through each input fastq file
for file_idx, file in enumerate(input_files):
    with open(file, "r") as handle:
        for line in handle:
            sequence = line.strip()  # Assuming each line contains a sequence
            for potential_seq in potential_sequences:
                if potential_seq in sequence:
                    sequence_counts[potential_seq][file_idx] += 1

# Write the counts to a tab-delimited output file

parts = input_files[0].split('_')
output_parts = parts[:4] + [parts[-1].split('.')[0] + "_counts.txt"]
output_file = '_'.join(output_parts)

with open(output_file, "w") as out_handle:
    # Write counts for each sequence
    for sequence, counts in sequence_counts.items():
        counts_str = "\t".join(map(str, counts))
        out_handle.write(f"{sequence}\t{counts_str}\n")

print(f"Counts written to {output_file}")
