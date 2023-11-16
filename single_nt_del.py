# This script extracts sequences containing a single nt deletion. 
# This is somewhat quick and dirty and could be improved and incorporated directly to the heatmap scripts
# The order_file is a list of all the single nt deletion sequences in order. 
# The input_file can be an output of one of the scripts that puts out normalized values or slopes

import sys

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 4:
    print("Usage: python single_nt_del.py <order_file> <input_file> <output_file>")
    sys.exit(1)

# Get input and output file names from command-line arguments
order_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

# Read sequences in the desired order from the order file into a list
with open(order_file, 'r') as order_file:
    order_list = [line.strip() for line in order_file]

# Read input sequences and store them in a dictionary with sequence as key and full line as value
sequences_dict = {}
with open(input_file, 'r') as infile:
    for line in infile:
        sequence = line.strip().split('\t')[0]
        sequences_dict[sequence] = line

# Write the reordered sequences to the output file
with open(output_file, 'w') as outfile:
    for sequence in order_list:
        if sequence in sequences_dict:
            outfile.write(sequences_dict[sequence])
