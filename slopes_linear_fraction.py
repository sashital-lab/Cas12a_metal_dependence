# This script takes the mean value abundance scores for nicked and SC fractions and determines the fraction of DNA that was not nicked or supercoiled (i.e. linearized)
# It then determines a slope based on the fraction linear at the four Mg2+ concentrations
# It outputs both the fraction linear values at each Mg2+ concentration and the slopes. The latter can be used to generate the slope heatmaps
# Be sure to entier the directory and filenames for the two means files.

import os
from scipy.stats import linregress

# Set directory and filenames
directory = '' # enter the name of the input/output directory
file1_name = '' # enter the name for the mean abundance file for the nicked fraction
file2_name = '' # enter the name for the mean abundance file for the SC faction

# Extract prefix from input filenames
prefix = '_'.join(file1_name.split('_')[:3])

# Initialize dictionaries
file1_dict = {}
file2_dict = {}

# Read data from file 1
with open(os.path.join(directory, file1_name), 'r') as file1:
    for line in file1:
        line = line.strip().split('\t')
        values = [float(x) if x else 0 for x in line[1:]]
        values += [0] * (4 - len(values))  # Pad with 0s if necessary
        file1_dict[line[0]] = values

# Read data from file 2
with open(os.path.join(directory, file2_name), 'r') as file2:
    header = file2.readline().strip()  # Read the header line
    headers = header.split('\t')[1:]  # Extract column headers
    for line in file2:
        line = line.strip().split('\t')
        values = [float(x) if x else 0 for x in line[1:]]
        values += [0] * (4 - len(values))  # Pad with 0s if necessary
        file2_dict[line[0]] = values


# Create a set of all unique sequences
sequences = set(file1_dict.keys()) | set(file2_dict.keys())

# Combine the two dictionaries
combined_values = []
for seq in sequences:
    value1 = file1_dict.get(seq, [0]*4)  # Get values for sequence from file 1
    value2 = file2_dict.get(seq, [0]*4)  # Get values for sequence from file 2
    combined_value = [max(0, 1 - (value1[i] + value2[i])) for i in range(4)] # Calculate combined value for each column
    combined_values.append([seq] + combined_value)


        
# assign x values
x = [1, 2, 5, 10]

# Calculate slopes and write to output file 1
slopes = []
for row in combined_values:
    y = [row[col] for col in range(1, 5)]
    slope, _, r_value, _, _ = linregress(x, y)
    slopes.append(slope)
    
with open(os.path.join(directory, f"{prefix}_lin_slopes.txt"), 'w') as outfile:
    for sequence, slope in zip(sequences, slopes):
        outfile.write(f"{sequence}\t{slope}\n")

# Write the combined data to a new file
with open(os.path.join(directory, f"{prefix}_lin_means.txt"), 'w') as outfile:
    # Write data
    for row in combined_values:
        outfile.write('\t'.join([str(x) for x in row]) + '\n')
