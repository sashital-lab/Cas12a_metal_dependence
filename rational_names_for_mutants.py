def generate_mutant_names(reference_sequence, mutant_sequences):
    # Check if the reference sequence length is 24
    if len(reference_sequence) != 24:
        print("Reference sequence length must be 24 nucleotides.")
        return []

    mutant_names = []

    # Compare each mutant sequence to the reference sequence
    for mutant_seq in mutant_sequences:
        if len(mutant_seq) == 24:
            mutations = []
            if mutant_seq == reference_seq:
                mutations.append("perfect")
            for pos in range(24):
                if mutant_seq[pos] != reference_sequence[pos]:
                    if pos < 4:
                        mutations.append(f"{reference_sequence[pos]}{pos - 4}{mutant_seq[pos]}")
                    else:
                        mutations.append(f"{reference_sequence[pos]}{pos - 3}{mutant_seq[pos]}")


            mutant_names.append(" ".join(mutations))
        elif len(mutant_seq) == 23:
            # Find the missing nucleotide in the shorter sequence
            for pos in range(24):
                if reference_sequence[:pos] + reference_sequence[pos + 1 :] == mutant_seq:
                    if pos < 4:
                        mutant_names.append(f"del {reference_sequence[pos]}{pos - 3}")
                        break
                    else:
                        mutant_names.append(f"del {reference_sequence[pos]}{pos - 2}")
                        break
        else:
            print(f"Skipping mutant sequence '{mutant_seq}' - length must be 23 or 24 nucleotides.")

    return mutant_names


def read_sequences_from_file(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            sequences.append(line.strip().split()[0])  # Assuming sequences are in the first column
    return sequences

# Example usage:
reference_seq = "TTTACGGCCACTTCCGTGTCTGAC"  # Replace with your reference sequence
input_file_path = "Fn_W_1_SC_1mM.txt"  # Replace with your input file path
output_file_path = "W_mutant_names.txt"  # Path for the output file

mutant_sequences = read_sequences_from_file(input_file_path)
mutant_names = generate_mutant_names(reference_seq, mutant_sequences)

# Write the mutant names to a file
with open(output_file_path, 'w') as output_file:
    for name in mutant_names:
        output_file.write(name + '\n')
