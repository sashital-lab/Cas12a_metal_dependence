# This script is for normalizing fraction total reads for an experimental sample to the control samples.
# This is for analyzing nicked data.
# It ignores any sequences that don't show up in the controls.
# The triplicate data should be compiled with other Mg concentrations using merge_nicked.py for an input file for slopes_individual_replicates_fracnorm_nicked.py

import sys
import os

# Define input/output file names
input_file = sys.argv[1]
base_filename = os.path.splitext(os.path.basename(input_file))[0]
output_file = base_filename.replace("_counts", "") + "_3reps.txt"
output_file_ave = base_filename.replace("_counts", "_ave.txt")

# Read data from the experimental sample input file
experimental_data = {}
with open(sys.argv[1], "r") as exp_file:
    for line in exp_file:
        fields = line.strip().split("\t")
        sequence = fields[0]
        values = [float(field) if field != '' else 0.0 for field in fields[1:]]
        experimental_data[sequence] = values

# Read data from the control input files (two replicates)
control_data = [{}, {}]
for control_idx in range(2):
    with open(sys.argv[2 + control_idx], "r") as control_file:
        for line in control_file:
            fields = line.strip().split("\t")
            sequence = fields[0]
            values = [float(field) if field != '' else 0.0 for field in fields[1:]]
            control_data[control_idx][sequence] = values

# Calculate the sum of total reads for each replicate
total_reads_exp = [sum(values) for values in zip(*experimental_data.values())]
total_reads_control = [[sum(values) for values in zip(*control.values())] for control in control_data]

# Calculate the fraction of total reads (FTR) for each sequence in each replicate
ftr_exp = {sequence: [value / total_reads_exp[i] for i, value in enumerate(values)] for sequence, values in experimental_data.items()}
ftr_control = [{sequence: [value / total_reads_control[control_idx][i] for i, value in enumerate(values)] for sequence, values in control.items()} for control_idx, control in enumerate(control_data)]

# Calculate average FTR values for each sequence in the control data
average_ftr_control = {}
for sequence in ftr_control[0]:
    ftr_values = [ftr_control[control_idx][sequence] for control_idx in range(2)]
    average_ftr = sum(sum(values) for values in zip(*ftr_values)) / (2 * len(ftr_values[0]))
    average_ftr_control[sequence] = average_ftr

# Normalize experimental data using average control FTR values
normalize_to_control = {}
for sequence, values in ftr_exp.items():
    average_control_ftr = average_ftr_control.get(sequence, 0)  # Default to 1.0 if sequence not found in control data

    # Ignore sequences not found in control data or with average_control_ftr value of 0
    if average_control_ftr != 0:
        normalized_values = [value / average_control_ftr for value in values]
    else:
        normalized_values = [""] * len(values)  # Set normalized values to NaN for ignored sequences

    normalize_to_control[sequence] = normalized_values

# Write normalized data to an output file
with open(output_file, "w") as out_file:
    sorted_sequences = sorted(normalize_to_control.keys())  # Sort sequences alphabetically
    for sequence in sorted_sequences:
        values = normalize_to_control[sequence]
        formatted_values = '\t'.join(['' if value == 0.0 else str(value) for value in values])
        out_file.write("{}\t{}\n".format(sequence, formatted_values))

# Calculate average values for the three replicates
#average_values = {}
#for sequence in sorted(normalize_to_control.keys()):
#    values = normalize_to_control[sequence]
#    avg_value = sum(values) / len(values)
#    average_values[sequence] = avg_value

# Write average values to a separate output file
#with open(output_file_ave, "w") as out_file:
#    sorted_sequences = sorted(average_values.keys())  # Sort sequences alphabetically
#    for sequence in sorted_sequences:
#        avg_value = average_values[sequence]
#
#        # Check if the absolute difference is very small, treat it as 0
#        if abs(avg_value) < 1e-10:
#            formatted_avg_value = ''
#        else:
#            formatted_avg_value = str(avg_value)
#
#        out_file.write("{}\t{}\n".format(sequence, formatted_avg_value))

print(f"Normalized data written to {output_file}")
