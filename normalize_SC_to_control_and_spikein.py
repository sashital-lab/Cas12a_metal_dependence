# This script is for normalizing both to the control and to the spike-in samples for the SC fraction.
# It also puts out an average file that can be used to generate heatmaps.
# It is agnostic to how many replicates there were in each experimental or control file, in case there were reps with low reads.


import sys
import os

# Define input/output file names
input_file = sys.argv[1]
base_filename = os.path.splitext(os.path.basename(input_file))[0]
output_file = base_filename.replace("_counts", "") + ".txt"
output_file_ave = base_filename.replace("_counts", "_ave") + ".txt"
output_label = base_filename.replace("_counts", "")

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
    average_control_ftr = average_ftr_control.get(sequence, 1.0)  # Default to 1.0 if sequence not found in control data
    normalized_values = [value / average_control_ftr for value in values]
    normalize_to_control[sequence] = normalized_values

# Determine number of usable replicates
num_replicates_exp = len(ftr_exp[next(iter(ftr_exp))])  # Number of replicates in experimental data
num_replicates_con1 = len(ftr_control[0][next(iter(ftr_control[0]))])  # Number of replicates in control data
num_replicates_con2 = len(ftr_control[1][next(iter(ftr_control[1]))])

# Calculate average FTR for the spike-in sequences in each experimental replicate and in control
spikein_sequences = [
    "GCATTCTTATTAATACATTTGAAA",
    "CGCGCCCAACTGACGCTAGGCAAG",
    "TAAGGGTAAACATACAAGTCGATA",
    "TCAGTGCAGGCTCCCGTGTTAGGA",
    "CATCCAGCACTCTACGGTTCCTCC"
]

num_spikein_sequences = len(spikein_sequences)

average_ftr_spikein = [
    sum([ftr_exp[sequence][replicate] for sequence in spikein_sequences]) / num_spikein_sequences
    for replicate in range(num_replicates_exp)
]

average_ftr_control_spikein = (sum(
    ftr_control[0][sequence][replicate]
    for replicate in range(num_replicates_con1)
    for sequence in spikein_sequences
    ) + sum(
    ftr_control[1][sequence][replicate]
    for replicate in range(num_replicates_con2)
    for sequence in spikein_sequences
    )) / ((num_replicates_con1+num_replicates_con2) * num_spikein_sequences )

# Calculate normalization factors for each experimental replicate
normalization_factors = [average_ftr_spikein[replicate] / average_ftr_control_spikein for replicate in range(num_replicates_exp)]

# Divide values in normalize_to_control by normalization factors
for sequence, values in normalize_to_control.items():
    for replicate in range(num_replicates_exp):
        values[replicate] /= normalization_factors[replicate]

# Write normalized data to an output file
with open(output_file, "w") as out_file:
    sorted_sequences = sorted(normalize_to_control.keys())  # Sort sequences alphabetically
    for sequence in sorted_sequences:
        values = normalize_to_control[sequence]
        formatted_values = '\t'.join(['' if value == 0 else str(value) for value in values])
        out_file.write("{}\t{}\n".format(sequence, formatted_values))

# Calculate average values for the three replicates
average_values = {}
for sequence in sorted(normalize_to_control.keys()):
    values = normalize_to_control[sequence]
    avg_value = sum(values) / len(values)
    average_values[sequence] = avg_value

# Write average values to a separate output file
with open(output_file_ave, "w") as out_file:
    sorted_sequences = sorted(average_values.keys())  # Sort sequences alphabetically
    for sequence in sorted_sequences:
        avg_value = average_values[sequence]

        # Check if the absolute difference is very small, treat it as 0
        if abs(avg_value) < 1e-10:
            formatted_avg_value = ''
        else:
            formatted_avg_value = str(avg_value)

        out_file.write("{}\t{}\n".format(sequence, formatted_avg_value))

# print(f"Normalized data written to {output_file}")
# print(f"Average values written to {output_file_ave}")
print(output_label + "\t" + "\t".join(map(str, normalization_factors)))
