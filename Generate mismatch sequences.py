
import itertools
import xlsxwriter

perfect_target = 'TTTGTACTGTCCACGCCGACGGAA'
mismatch_count = [1,2] #include sequences with all of these number of mismatches
deletion_count = [1] #include sequences with all of these number of deletions
include_perfect = True
def gen_mismatch_sequences(base_sequence,num_mismatches):
    '''
    returns a list of sequences with some
    number of mismatches relative to
    the base sequence
    '''
    possible_bases = ['A', 'T', 'C', 'G']
    position_combos = list(itertools.combinations(range(len(base_sequence)),num_mismatches))
    base_combos = list(itertools.product(possible_bases,repeat = num_mismatches))
    sequence_list = []
    for positions in position_combos:
        new_mm_seq = base_sequence
        for bases in base_combos:
            for pair in range(num_mismatches):
                new_mm_seq = new_mm_seq[:positions[pair]] + bases[pair] + new_mm_seq[positions[pair]+1:]
            mismatches = 0
            for pos in range(len(base_sequence)):
                if new_mm_seq[pos] != base_sequence[pos]:
                    mismatches += 1
            if mismatches == num_mismatches:
                sequence_list.append(new_mm_seq)
    return  sequence_list

def gen_deletion_sequences(base_sequence,num_deletions):
    '''
    returns a list of sequences with some
    number of deletions relative to the base sequence.
    note: same deletions caused by repeated bases in the base
    sequence are removed
    '''

    position_combos = list(itertools.combinations(range(len(base_sequence)),num_deletions))
    sequence_list = []
    for positions in position_combos:
        new_seq = ''
        for pos in range(len(base_sequence)):
            if pos not in positions:
                new_seq += base_sequence[pos]
        sequence_list.append(new_seq)
    sequence_list = list(set(sequence_list)) # removing duplicates
    return  sequence_list

#generating and adding all sequences to a list
all_sequences = []
for pos in mismatch_count:
    all_sequences += gen_mismatch_sequences(perfect_target,pos)
for pos in deletion_count:
    all_sequences += gen_deletion_sequences(perfect_target,pos)
if include_perfect == True:
    all_sequences.append(perfect_target)

#writing sequences to an excel sheet all in the first column
workbook = xlsxwriter.Workbook('mm_del_sequences.xlsx')
worksheet = workbook.add_worksheet('sheet1')
column = 0
row = 0
for sequence in all_sequences:
        worksheet.write(row, column, sequence)
        row += 1
workbook.close()

