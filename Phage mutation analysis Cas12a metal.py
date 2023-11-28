
import xlsxwriter
import glob


def deletion_fill(in_seq, perfect):
    length_dif = len(perfect) - len(in_seq)
    out_seq = ''
    mm_found = False
    if length_dif > 0:
        for position in range(len(in_seq)):
            if perfect[position] != in_seq[position]:
                out_seq = perfect[:position] + (length_dif * 'X') + perfect[position + length_dif:]
                mm_found = True
                break
        if mm_found == False:  # happens when deletion is taken from the end
            out_seq = in_seq + (length_dif * 'X')
        return out_seq
    else:
        return in_seq

def reverse_comp(in_sequence):  # faster way lifted from stackoverflow
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq = in_sequence
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def combine_file_dicts(list_of_dicts):
    output_dict = {}
    for dict in list_of_dicts:
        for item in dict:
            if item not in output_dict:
                output_dict[item] = dict[item]
            elif item in ['valid_sequences', 'wt_sequences', 'mutated_sequences', 'sequence_lines']:
                output_dict[item] += dict[item]
            else:
                output_dict[item]['count'] += dict[item]['count']
    return output_dict

def calculate_diversity(file_dict): # need to prevent matching of sequences twice
    diversity = 0
    counted_list = []
    for item1 in file_dict:
        counted_list.append(item1)
        for item2 in file_dict:
            if item1 != perfect_target and item2 != perfect_target and item1 != item2 and item2 not in counted_list:
                if item1 not in ['valid_sequences', 'wt_sequences', 'mutated_sequences', 'sequence_lines']:
                    if item2 not in ['valid_sequences', 'wt_sequences', 'mutated_sequences', 'sequence_lines']:
                        div_mismatch_count = 0  # number of mismatches comparing item1 and item2
                        if len(deletion_fill(item1, perfect_target)) == 24 and len(deletion_fill(item2, perfect_target)) == 24:
                            for h in range(24):
                                if deletion_fill(item1, perfect_target)[h] != deletion_fill(item2, perfect_target)[h]:
                                    div_mismatch_count += 1
                        else:
                            print("ERROR", perfect_target, len(deletion_fill(item1, perfect_target)),
                                  len(deletion_fill(item2, perfect_target)))

                        diversity += (2 * (file_dict[item1]['count'] / file_dict['valid_sequences']) * (file_dict[item2]['count'] / file_dict['valid_sequences']) * (div_mismatch_count / 24))
    return diversity

def target_information_lib(csv_file,target_name):
    file = open(csv_file)
    lines = file.readlines()
    headers = ''
    info = ''
    for line in lines:
        if 'target' in line.split(','):
            headers = line.split(',')
        if target_name in line.split(','):
            info = line.split(',')
    headers[-1] = headers[-1][:-1]
    info[-1] = info[-1][:-1]
    info_dict = dict(zip(headers, info))
    return info_dict


specific_target_name = 'Gene_A'
target_info = target_information_lib("Gene target flanking sequences.csv",specific_target_name)
perfect_target = target_info['perfect_target']  # FORWARD orientation
to_find1_R1 = target_info['R1_left'] # extract the target sequence between these two sequences from the R1 file
to_find2_R1 = target_info['R1_right']
to_find1_R2 = target_info['R2_left'] # extract the target sequence between these two sequences from the R2 file
to_find2_R2 = target_info['R2_right']

workbook_new = xlsxwriter.Workbook(target_info['workbook_name'])  # naming the output excel sheet

# selecting files to analyze, here it's all the Cas12a ones and the control
total_file_list = []
total_file_list += glob.glob('*As_A_*')
total_file_list += glob.glob('*Fn_A_*')
total_file_list += glob.glob('*Lb_A_*')
total_file_list += glob.glob('*NT_A*')



print(total_file_list)
print(len(total_file_list), 'files gathered')

file_list_length_half = int(len(total_file_list) / 2)
print(file_list_length_half)

total_file_list_pairs = []  # this list will contain pairs of files ordered by their order in total file list
# total file list must be in order or the wrong files will be matched!!

for x in range(file_list_length_half):
    total_file_list_pairs.append([total_file_list[2 * x]] + [total_file_list[2 * x + 1]])

print(total_file_list_pairs)


# for file_pair in total_file_list_pairs:
def make_file_dict(file_pair):
    output_dict = {'valid_sequences': 0, 'wt_sequences': 0, 'mutated_sequences': 0, 'sequence_lines': 0}

    R1_file = open(file_pair[0])
    R2_file = open(file_pair[1])
    lines_R1 = R1_file.readlines()
    lines_R2 = R2_file.readlines()
    line_count_r1 = 0
    line_count_r2 = 0
    output_dict['sequence_lines'] += min(len(lines_R1), len(lines_R2))
    for line_num in range(len(lines_R1)):  # lines_R1 and R2 SHOULD have the same length
        temp_sequence_R1 = ''
        temp_sequence_R2 = ''
        if lines_R1[line_num].find(to_find1_R1) != -1 and lines_R1[line_num].find(to_find2_R1) != -1:
            line_count_r1 += 1
            temp_sequence_R1 = lines_R1[line_num][
                               lines_R1[line_num].find(to_find1_R1) + len(to_find1_R1):   lines_R1[line_num].find(
                                   to_find2_R1)]

        if lines_R2[line_num].find(to_find1_R2) != -1 and lines_R2[line_num].find(to_find2_R2) != -1:
            line_count_r2 += 1
            temp_sequence_R2 = lines_R2[line_num][
                               lines_R2[line_num].find(to_find1_R2) + len(to_find1_R2): lines_R2[line_num].find(
                                   to_find2_R2)]

        if len(temp_sequence_R1) == len(temp_sequence_R2) and len(temp_sequence_R1) > 20 and len(temp_sequence_R2) > 20:
            # if they are same length, rev complement R2, match with R1
            # same_length_sequences += 1

            if reverse_comp(temp_sequence_R2) == temp_sequence_R1:
                output_dict['valid_sequences'] += 1
                temp_sequence = temp_sequence_R1
                if temp_sequence not in output_dict:
                    output_dict[temp_sequence] = {'mismatch_positions': [], 'deletion': False}
                    output_dict[temp_sequence]['count'] = 0
                else:
                    output_dict[temp_sequence]['count'] += 1

    for sequence in output_dict:
        if sequence not in ['valid_sequences', 'wt_sequences', 'mutated_sequences', 'sequence_lines']:
            if len(sequence) < len(perfect_target):  # accounting for deletion sequences first
                difference = len(perfect_target) - len(sequence)
                deletion_pos_found = False
                for x in range(len(sequence)):
                    if sequence[x] != perfect_target[x]:
                        deletion_pos_found = True
                        for mm in range(difference):
                            output_dict[sequence]['mismatch_positions'].append(x + mm)
                        break
                if deletion_pos_found == False:  # accounting for deletions at the end of the sequence
                    for mm in range(difference):
                        output_dict[sequence]['mismatch_positions'].append(len(perfect_target) - mm - 1)
                output_dict[sequence]['deletion'] = True
            if len(sequence) == len(perfect_target):  # non-deletion sequences
                for position in range(len(sequence)):
                    if sequence[position] != perfect_target[position]:
                        output_dict[sequence]['mismatch_positions'].append(position)

        if sequence != perfect_target and sequence not in ['valid_sequences', 'wt_sequences', 'mutated_sequences',
                                                           'sequence_lines']:
            output_dict['mutated_sequences'] += output_dict[sequence]['count']

    if perfect_target in output_dict:
        output_dict['wt_sequences'] += output_dict[perfect_target]['count']

    return output_dict


def make_mismatch_lib(in_file_dict):
    mismatch_library = {}  # like former mismatch list for writing to excel sheet for bar graph and stuff
    categories = ['A', 'T', 'C', 'G', 'deletion', 'multi', 'position_sum']
    for cat in categories:
        mismatch_library[cat] = [0] * 24

    for seq in in_file_dict:  # output of make_file_dict
        if seq != perfect_target and seq not in ['valid_sequences', 'wt_sequences', 'mutated_sequences',
                                                 'sequence_lines']:
            if len(in_file_dict[seq]['mismatch_positions']) == 1 and in_file_dict[seq][
                'deletion'] == False:  # single mm
                for pos in in_file_dict[seq]['mismatch_positions']:
                    base = seq[pos]
                    mismatch_library[base][pos] += in_file_dict[seq]['count']
            if len(in_file_dict[seq]['mismatch_positions']) > 1 and in_file_dict[seq]['deletion'] == False:  # multi mm
                for pos in in_file_dict[seq]['mismatch_positions']:
                    mismatch_library['multi'][pos] += in_file_dict[seq]['count']
            if in_file_dict[seq]['deletion'] == True:  # deletion
                for pos in in_file_dict[seq]['mismatch_positions']:
                    mismatch_library['deletion'][pos] += in_file_dict[seq]['count']
            if len(seq) < 25:  # basically all not crazy sequences
                for pos in in_file_dict[seq]['mismatch_positions']:
                    mismatch_library['position_sum'][pos] += in_file_dict[seq]['count']
    return mismatch_library


def create_file_worksheet(workbook,temp_file_dict, mismatch_library,file_name):  # directly adds worksheet to workbook, does not return
    # for file_pair in total_file_list_pairs:


    temp_mismatch_library = mismatch_library

    if 'NT' in file_name:
        file_name_short = file_name[0:file_name.find('hr_') + 2]
    else:
        file_name_short = file_name[0:file_name.find('hr_') + 4]
    worksheet = workbook.add_worksheet(file_name_short)

    perfect_target_list = []
    for x in range(len(perfect_target)):
        perfect_target_list += perfect_target[x]
    worksheet.write_column('B3', perfect_target_list)
    categories_ordered = ['A', 'T', 'C', 'G', 'deletion', 'multi', 'position_sum']  # order to write data in excel sheet
    row = 1
    column = 2
    for mut_category in categories_ordered:
        worksheet.write(row, column, mut_category)
        row += 1
        for position in range(len(temp_mismatch_library[mut_category])):
            worksheet.write(row + position, column, temp_mismatch_library[mut_category][position])
        column += 1
        row = 1

    # Creating chart object
    chart = workbook.add_chart({'type': 'column', 'subtype': 'stacked'})
    for type in range(len(categories_ordered)):
        if categories_ordered[type] != 'position_sum':
            chart.add_series({'values': [file_name_short, 2, 2 + type, 2 + len(perfect_target) - 1, 2 + type],
                              # [sheetname, first_row, first_col, last_row, last_col]
                              'categories': [file_name_short, 2, 1, 25, 1],
                              'name': categories_ordered[type]
                              })  # values strings specifies what values in the excel sheet to make the chart from

    chart.set_x_axis({'name': 'mutation position'})
    chart.set_y_axis({'name': 'sequence counts'})
    worksheet.insert_chart('N3', chart)

    row = 2
    column = 10
    for item in ['sequence_lines', 'valid_sequences', 'wt_sequences', 'mutated_sequences']:
        worksheet.write(row, column, temp_file_dict[item])
        worksheet.write(row, column + 1, item)
        row += 1
    worksheet.write(row, column, temp_file_dict['mutated_sequences'] / temp_file_dict['valid_sequences'])
    worksheet.write(row, column + 1, 'fraction_mutated')


def create_all_info_worksheet(workbook,dict_collection):
    all_file_info = workbook.add_worksheet('info from all files')
    row = 0
    column = 0
    categories = ['file_name', 'fraction_mutated', 'diversity']
    for item in categories:
        all_file_info.write(row, column, item)
        column += 1
    column = 0
    triple_count = 1
    triple_dict_list = []
    for item in dict_collection:
        all_file_info.write(row+1, column, item)
        all_file_info.write(row+1, column+1, dict_collection[item]['file_dict']['mutated_sequences'] /
                            dict_collection[item]['file_dict']['valid_sequences'] )
        #writing position_sum values for heat maps
        for position in range(len(dict_collection[item]['mismatch_lib']['position_sum'])):
            all_file_info.write(row+1, column+3+position,dict_collection[item]['mismatch_lib']['position_sum'][position]/
                                dict_collection[item]['file_dict']['valid_sequences'])
        triple_dict_list.append(dict_collection[item]['file_dict'])
        if triple_count == 3:
            diversity = calculate_diversity(combine_file_dicts(triple_dict_list))
            all_file_info.write(row+1, column+2, diversity)
            triple_count = 0
            triple_dict_list = []
        triple_count += 1
        row += 1


# execution step for processing all files that were pooled
file_dict_collection = {}  # library with file dictionaries and mismatch libs for each file pair
#make once for each file to use for further processing
for file_pair in total_file_list_pairs:
    file_dict_collection[file_pair[0]] = {'file_dict': make_file_dict(file_pair), 'mismatch_lib': 'NA'}
    file_dict_collection[file_pair[0]]['mismatch_lib'] = make_mismatch_lib( file_dict_collection[file_pair[0]]['file_dict'])

create_all_info_worksheet(workbook_new,file_dict_collection)
for item in file_dict_collection:
    create_file_worksheet(workbook_new,file_dict_collection[item]['file_dict'], file_dict_collection[item]['mismatch_lib'] ,item)

workbook_new.close()

