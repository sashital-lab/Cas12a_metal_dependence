
import xlsxwriter
import glob


def deletion_fill(in_seq,perfect):
    length_dif = len(perfect) - len(in_seq)
    out_seq = ''
    mm_found = False
    if length_dif > 0:
        for position in range(len(in_seq)):
            if perfect[position] != in_seq[position]:
                out_seq = perfect[:position] + (length_dif * 'X') + perfect[position+length_dif:]
                mm_found = True
                break
        if mm_found == False: #happens when deletion is taken from the end
            out_seq = in_seq + (length_dif * 'X')
        return out_seq
    else:
        return in_seq


def reverse_comp(in_sequence): #faster way lifted from stackoverflow
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq = in_sequence
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return  reverse_complement


complement_list = [['A', 'T'], ['T', 'A'], ['C', 'G'], ['G', 'C']]




perfect_target = 'TTTGATGATGATATTGAACAGGAA' #FORWARD orientation
perfect_target_list = []
for x in range(len(perfect_target)):
    perfect_target_list += perfect_target[x]
print(perfect_target_list)

to_find1_R1 = 'TTGCTGCT'  #extract the target sequence between these two sequences from the R1 file
to_find2_R1 = 'GGCTCTCC'

to_find1_R2 = 'GGAGAGCC'  #extract the target sequence between these two sequences from the R2 file
to_find2_R2 = 'AGCAGCAA'

workbook = xlsxwriter.Workbook('(top)cas12a_escapers A R1R2.xlsx')  # naming the output excel sheet

# selecting files to analyze, here it's all the Cas12a ones and the control
total_file_list = []
total_file_list += glob.glob('*As_A_*')
total_file_list += glob.glob('*Fn_A_*')
total_file_list += glob.glob('*Lb_A_*')
total_file_list += glob.glob('*NT_A*')


print(total_file_list)
print(len(total_file_list) , 'files gathered')

file_list_length_half = int(len(total_file_list) / 2)
print(file_list_length_half)

total_file_list_pairs = [] #this list will contain pairs of files ordered by their order in total file list
#total file list must be in order or the wrong files will be matched!!

for x in range(file_list_length_half):
    total_file_list_pairs.append([total_file_list[2*x]] + [total_file_list[2*x + 1]])

print(total_file_list_pairs)

all_sample_info = [] # this will be written on the last sheet containing info for all of the samples processed
worksheet_all = workbook.add_worksheet('All info')

triple_count = 0 # iterated to keep track of triplicate samples for nucleotide diversity calculations
triple_sequence_lines = 0
seq_count_list = []  # used for collecting sequences from three replicates for nucleotide diversity calculations, need to reset after triplicate is complete

nucleotide_diversity = 0

for file_pair in total_file_list_pairs:

    triple_count += 1 #dealing with combined triplicates for nuc diversity

    info_list_one_sample = [] # contains the info for only one sample: sample name, % mutated, number of samples over 1% of all reads

    mismatch_list = []  # list for making the mismatch distribution for each file
    for a in range(len(perfect_target)):
        mismatch_list.append([0, 0, 0, 0, 0, 0])  # mutation/deletion types: A,T,C,G,deletion,multiple mismatches where A means previous base mutated to A

    sequence_lines = 0    #this represents the total number of valid sequences found, valid = R1R2 match
    wild_type_sequences = 0
    single_mm_sequences = 0
    indel_sequences = 0
    multi_mm_sequences = 0




    # make a list with values for each file
    # total counts, WT, 1MM, 1deletion, locations of mm/del, 0,1,2,3,4,5... ,23



    R1_file = open(file_pair[0])
    R2_file = open(file_pair[1])

    lines_R1 = R1_file.readlines()
    lines_R2 = R2_file.readlines()

    print(len(lines_R1))
    print(len(lines_R2))

    line_count_r1 = 0
    line_count_r2 = 0

    same_length_sequences = 0
    matching_sequences = 0

    for line_num in range(len(lines_R1)):  # lines_R1 and R2 SHOULD have the same length
        temp_sequence_R1 = ''
        temp_sequence_R2 = ''
        if lines_R1[line_num].find(to_find1_R1) != -1 and lines_R1[line_num].find(to_find2_R1) != -1:
            temp_sequence_R1 = lines_R1[line_num][lines_R1[line_num].find(to_find1_R1) + len(to_find1_R1): lines_R1[line_num].find(to_find2_R1)]   #72 here because there are many occcurences of "to_find_R1" in the line
            line_count_r1 += 1
            # print('true')
            # print(temp_sequence_R1)
            # print(lines_R1[line_num].find(to_find1_R1))
            # print(lines_R1[line_num].find(to_find2_R1, 72) )

        if lines_R2[line_num].find(to_find1_R2) != -1 and lines_R2[line_num].find(to_find2_R2) != -1:
            temp_sequence_R2 = lines_R2[line_num][lines_R2[line_num].find(to_find1_R2) + len(to_find1_R2): lines_R2[line_num].find(to_find2_R2)]
            # print(temp_sequence_R2)
            line_count_r2 += 1



        if len(temp_sequence_R1) == len(temp_sequence_R2) and len(temp_sequence_R1) > 20 and len(temp_sequence_R2) > 20 :  # if they are same length, rev complement R2, match with R1
            same_length_sequences += 1


            if reverse_comp(temp_sequence_R2) == temp_sequence_R1:

                matching_sequences += 1
                sequence_lines += 1
                triple_sequence_lines += 1
                # print(temp_sequence_R1)
                temp_sequence = temp_sequence_R1  # the final sequence used for analysis and comparison
                # print(temp_sequence)

                if len(temp_sequence) < len(perfect_target):  # counting single indels, multiple deletions are ignored
                    difference = len(perfect_target) - len(temp_sequence)
                    for x in range(len(temp_sequence)):
                        if temp_sequence[x] != perfect_target[x]:
                            for index in range(difference):
                                mismatch_list[x + index][4] += 1  # 4th position is number of deletions in "mismatch list"

                            break
                    indel_sequences += 1


                all_matching = True
                complement_base_pos_dict = {'T': 1, 'A': 0, 'G': 3, 'C': 2}

                if len(temp_sequence) == len(perfect_target):
                    mm_count_temp = 0
                    for x in range(len(temp_sequence)):  #counting total dmismatches
                        if temp_sequence[x] != perfect_target[x]:
                            mm_count_temp += 1

                    for x in range(len(temp_sequence)):
                        if temp_sequence[x] != perfect_target[x]:
                            if mm_count_temp == 1: #single mismatches
                                mismatch_list[x][complement_base_pos_dict[temp_sequence[x]]] += 1
                                single_mm_sequences += 1
                                break
                            if mm_count_temp > 1: #multiple mismatches
                                mismatch_list[x][5] += 1


                    if mm_count_temp == 0 :
                        wild_type_sequences += 1
                    if mm_count_temp > 1:
                        multi_mm_sequences += 1

                # making seq count list - to compare to control file to see if mutated sequence was there already
                in_list = False
                for x in range(len(seq_count_list)):
                    if temp_sequence in seq_count_list[x][0]:
                        seq_count_list[x][1] += 1
                        in_list = True
                        break
                if in_list == False:
                    seq_count_list.append([temp_sequence,0,0,0,0,0,0])  #sequence, count, control count, status, base before, base after, position

    print(sequence_lines, 'sequence lines')
    seq_count_list = [item for item in seq_count_list if item[1] > 0 ]  # WT Cas12a protospacer
    print('seq_count_list has length' ,len(seq_count_list))

    if triple_count == 3:

        for sequence4 in seq_count_list:
            short_length_match = min(len(perfect_target), len(sequence4[0])) #avoiding out of index errors for deletions or extensions
            for x in range(short_length_match):
                if sequence4[0][x] != perfect_target[x]:
                    sequence4[6] = x + 1  #sixth position is mismatch/indel location
                    sequence4[4] = perfect_target[x]
                    sequence4[5] = sequence4[0][x]
                    break

            mm_count_temp_2 = 0

            for x in range(short_length_match):
                if sequence4[0][x] != perfect_target[x]:
                    mm_count_temp_2 += 1
            if len(sequence4[0]) < len(perfect_target):
                sequence4[3] = 'deletion'
            if len(sequence4[0]) > len(perfect_target):
                sequence4[3] = 'extension'
            if len(sequence4[0]) == len(perfect_target) and mm_count_temp_2 == 1:
                sequence4[3] = 'SNP'
            if len(sequence4[0]) == len(perfect_target) and mm_count_temp_2 > 1:
                sequence4[3] = 'multi-mismatch'
            if sequence4[0] == perfect_target:
                sequence4[3],sequence4[4],sequence4[5],sequence4[6] = 'perfect','perfect','perfect','perfect'


        #NEW SECTION FOR nucleotide diversity

        for index1 in range(len(seq_count_list) -1):
            for index2 in range(index1+1, len(seq_count_list)):
                temp1 = seq_count_list[index1][0]
                temp2 = seq_count_list[index2][0]
                if temp1 != perfect_target and temp2 != perfect_target: #ignoring WT sequence contribution to diversity in this version
                    div_mismatch_count = 0 # number of mismatches comparing temp1 and temp2
                    if len(deletion_fill(temp1,perfect_target)) == 24 and len(deletion_fill(temp2,perfect_target)) == 24:
                        for h in range(24):
                            if deletion_fill(temp1,perfect_target)[h] != deletion_fill(temp2,perfect_target)[h]:
                                div_mismatch_count += 1
                    else:
                        print("ERROR", perfect_target, len(deletion_fill(temp1,perfect_target)),len(deletion_fill(temp2,perfect_target)))

                    nucleotide_diversity += (  2 * (seq_count_list[index1][1]/triple_sequence_lines) * (seq_count_list[index2][1]/triple_sequence_lines) * (div_mismatch_count/24)  )
                        # print(div_mismatch_count, 'mismatch count', 2 * (seq_count_list[index1][1]/sequence_lines) * (seq_count_list[index2][1]/sequence_lines) * (div_mismatch_count/20) )

        print(seq_count_list)
        print(nucleotide_diversity, 'nucleotide diversity')

    for x in range(len(mismatch_list)): #making a fraction of mutated sequences column for heat maps
        sum = 0
        for y in range(len(mismatch_list[x])):
            sum += mismatch_list[x][y]

        mismatch_list[x].append(sum/sequence_lines)
    print(mismatch_list)

    if sequence_lines != 0 :  # shouldn't actually need this line if you glob files correctly above

        if 'NT' in file_pair[0]:
            file_name = file_pair[0][0:file_pair[0].find('hr_') + 2]
        else:
            file_name = file_pair[0][0:file_pair[0].find('hr_') + 4]
        worksheet = workbook.add_worksheet(file_name)

        info_list_one_sample.append(file_name)
        #every three times we will actually add nucleotide diversity for the triplicate
        if triple_count == 3:
            info_list_one_sample.append(nucleotide_diversity)
        else:
            info_list_one_sample.append('triplicate not complete')
        info_list_one_sample.append(round(((indel_sequences + single_mm_sequences + multi_mm_sequences) / sequence_lines), 2))
        numerous_sequences = 0

        for sequence in range(len(seq_count_list)):
            if seq_count_list[sequence][1] > (sequence_lines/100):
                numerous_sequences += 1
        info_list_one_sample.append(numerous_sequences)
        add_list = []
        for item in mismatch_list:
            add_list.append(item[6])
        info_list_one_sample.append(add_list)
        print(info_list_one_sample)
        all_sample_info.append(info_list_one_sample)

        # Creating chart object
        chart = workbook.add_chart({'type': 'column', 'subtype': 'stacked'})
        # worksheet.write_column('A1', mismatch_list)

        row = 2
        column = 2
        for A0, T1, C2, G3, Deletion4, multi_mm_5, fractions in mismatch_list:
            worksheet.write(row, column, A0)
            worksheet.write(row, column + 1, T1)
            worksheet.write(row, column + 2, C2)
            worksheet.write(row, column + 3, G3)
            worksheet.write(row, column + 4, Deletion4)
            worksheet.write(row, column + 5, multi_mm_5)
            worksheet.write(row, column + 6, fractions)
            row += 1

        worksheet.write_column('B1', perfect_target_list)
        name_list = ['A', 'T', 'C', 'G', 'Deletion', 'multi_mm']
        for type in range(len(name_list)):
            chart.add_series({'values': [file_name, 2, 2 + type, 2 + len(perfect_target) - 1, 2 + type],
                              # [sheetname, first_row, first_col, last_row, last_col]
                              'categories': [file_name, 0, 1, 23, 1],
                              'name': name_list[type]
                              })  # values strings specifies what values in the excel sheet to make the chart from

        chart.set_x_axis({'name': 'mutation position'})
        chart.set_y_axis({'name': 'sequence counts'})

        worksheet.insert_chart('N3', chart)

        worksheet.write(39, 2, 'mutated sequence')
        worksheet.write(39, 3, 'count')
        worksheet.write(39, 4, 'count in control ')

        if seq_count_list != []:  # writing counts of sequences found in control file
            row = 40
            column = 2
            for sequence1, count, control_count, status, base_before, base_after, position in seq_count_list:
                worksheet.write(row, column, sequence1)
                worksheet.write(row, column + 1, count)
                worksheet.write(row, column + 2, control_count)
                worksheet.write(row, column + 3, status)
                worksheet.write(row, column + 4, position)
                worksheet.write(row, column + 5, base_before)
                worksheet.write(row, column + 6, 'to')
                worksheet.write(row, column + 7, base_after)
                row += 1

        # writing checks to make sure all sequences are accounted for
        worksheet.write(2, 9, str(sequence_lines) )
        worksheet.write(2, 10, str("total sequences") )
        worksheet.write(3, 9, str(single_mm_sequences + indel_sequences) + " single mutation or indel sequences")
        worksheet.write(4, 9, str(wild_type_sequences) + " wild-type sequences")
        worksheet.write(5, 9, str(sequence_lines - (indel_sequences + single_mm_sequences + wild_type_sequences + multi_mm_sequences)) + " other sequences")
        worksheet.write(6, 9, "fraction mutated is " )
        worksheet.write(6, 10, str(round(((indel_sequences + single_mm_sequences + multi_mm_sequences) / sequence_lines), 2)))

        # values_string = "=" + file_name + "!$A$1:$A$24"
        # print(values_string)
        options = {'fill': {'none': True}, 'line': {'none': True}, 'font': {'color': 'black',
                                                                            'size': 14},
                   'font': {'underline': True}}
        worksheet.insert_textbox('R5', 'Mutated Base', options)

        #resetting triplicate things for nucleotide diversity
        #reset after going through the third replicate

        if triple_count == 3:
            triple_count = 0
            seq_count_list = []
            triple_sequence_lines = 0
            nucleotide_diversity = 0

row = 0
column = 0

for name, diversity, percent, numerous,fractions  in all_sample_info:
    worksheet_all.write(row, column, name)
    worksheet_all.write(row, column + 1, diversity)
    worksheet_all.write(row, column + 2, percent)
    worksheet_all.write(row, column + 3, numerous)

    column_add = 0
    for position in fractions:
        worksheet_all.write(row, column + column_add + 4, position)
        column_add += 1

    row +=1










workbook.close()

