# This script generates volcano plots using input csv files generated with volcano_data.py
# For double-mutant sequences, the dots in the plot can be colored based on the position of the first (lines 52-62) or second mutation (64-74)

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
import os

# Read the data from the CSV file# Define the directory and input file name
directory_in = 'C:/Users/sashital/Box/Manuscripts/Cas12a_metal-dependent_specificity/figure_workshop/231220_libraries/plasmid_library/nicked_input'
directory_out = 'C:/Users/sashital/Box/Manuscripts/Cas12a_metal-dependent_specificity/figure_workshop/231220_libraries/plasmid_library/nicked_volcano'
input_file = 'Fn_L_30_nicked_volcano.csv'
file_path = os.path.join(directory_in, input_file)
data = pd.read_csv(file_path)

# Extracting columns from the data
mutants = data['mutant']
log2_fold_change = data['log2 fold change']
neg_log10_p_value = data['-log P']

def extract_numbers(mutant):
    if mutant == 'perfect':  # Ignore perfect sequence
        return []

    # Extracting all numbers and hyphen-separated values using regular expression
    numbers = re.findall(r'[-+]?\d*\.*\d+|\w+-\d+', mutant)
    # Converting extracted strings to integers
    numbers = [int(re.search(r'[-+]?\d+', num).group()) if re.search(r'[-+]?\d+', num) else 0 for num in numbers]
    return numbers

# Specify an RGB color using a tuple
PAM = (1, 0.859, 0.090)
seed = (0.710, 0.871, 0.969)
mid = (0.816, 1, 0.831)
distal = (0.945, 0.816, 1)

# Assigning colors based on the number of values in numbers
def assign_color(numbers):
    if len(numbers) == 1:
        number = numbers[0]
        if number in [-4, -3, -2, -1, 0]:
            return PAM
        elif 1 <= number <= 8:
            return seed
        elif 9 <= number <= 13:
            return mid
        elif 14 <= number <= 20:
            return distal
    # Use this for uncleaved fraction
    # elif len(numbers) == 2:
    #     if numbers:  # Check if numbers array is not empty
    #         first_number = numbers[0]
    #         if first_number in [-4, -3, -2, -1]:
    #             return PAM
    #         elif 1 <= first_number <= 8:
    #             return seed
    #         elif 9 <= first_number <= 13:
    #             return mid
    #         elif 14 <= first_number <= 20:
    #             return distal
    # Use this for nicked fraction
    elif len(numbers) == 2:
        if numbers:  # Check if numbers array is not empty
            second_number = numbers[1]
            if second_number in [-4, -3, -2, -1, 0]:
                return PAM
            elif 1 <= second_number <= 8:
                return seed
            elif 9 <= second_number <= 13:
                return mid
            elif 14 <= second_number <= 20:
                return distal    
    return 'grey'  # Default color for empty array or other values

# Extracting numbers and assigning colors
numbers = [extract_numbers(mutant) for mutant in mutants]
colors = [assign_color(number) for number in numbers]

# Assign colors for significant and non-significant points based on p-value threshold
threshold_colors = ['lightgrey' if (-1 <= fold_change <= 1) or p_val < 1.3 else color
                    for fold_change, p_val, color in zip(log2_fold_change, neg_log10_p_value, colors)]

# Set the font properties using rcParams for consistency
plt.rcParams.update({
    'font.family': 'Arial',  # Change the font family
    'font.size': 6,           # Change the font size
})

# Create the volcano plot with standardized settings
plt.figure(figsize=(1.25, 1.75), dpi=300)  # Set a standardized figure size and DPI
plt.scatter(log2_fold_change, neg_log10_p_value, color=threshold_colors, edgecolors='grey', linewidths=0.25, s=10)
# plt.scatter([], [], color=PAM, label='PAM', edgecolor='grey')
# plt.scatter([], [], color=seed, label='seed', edgecolor='grey')
# plt.scatter([], [], color=mid, label='mid', edgecolor='grey')
# plt.scatter([], [], color=distal, label='distal', edgecolor='grey')
plt.axhline(y=1.3, color='black', linestyle='--', linewidth=0.5)
plt.axvline(x=-1, color='black', linestyle='--', linewidth=0.5)
plt.axvline(x=1, color='black', linestyle='--', linewidth=0.5)
plt.xlabel('log2 fold change',fontname='Arial',fontsize=8)
plt.ylabel('-log10 p-value',fontname='Arial',fontsize=8)
#plt.title('Uncleaved sequences in 1 mM Mg2+ versus 10 mM Mg2+',fontname='Arial')
#plt.legend()
plt.grid(False)
plt.tight_layout()

# Get the current axes
ax = plt.gca()

# Set the linewidth for the axis lines
ax.spines['top'].set_linewidth(0.5)    # Top axis line
ax.spines['bottom'].set_linewidth(0.5) # Bottom axis line
ax.spines['left'].set_linewidth(0.5)   # Left axis line
ax.spines['right'].set_linewidth(0.5)  # Right axis line

# Set the linewidth for the tick marks
ax.tick_params(axis='both', width=0.5, length=2, pad = 2) # Width of tick marks on x and y axes

# Setting custom axis ranges
plt.xticks(range(-6, 7,2))
plt.xlim(-6, 6)  # Set the x-axis range from -5 to 5
plt.ylim(0, 7)  # Set the y-axis range from 0 to 30 (adjust as needed)



# Save the current figure instance as a PDF in the output directory
filename = os.path.splitext(os.path.basename(input_file))[0]



# Adjusting the plot and saving as a PDF with font embedding
mpl.rcParams['pdf.fonttype'] = 42  # Embed fonts as TrueType (editable)
mpl.rcParams['font.family'] = 'Arial'  # Set the default font for the entire plot

# Save the plot with standardized settings and a tight bounding box
plt.tight_layout(pad=0.1)
plt.savefig(
    os.path.join(directory_out, filename + '.pdf'),
    dpi=300,
    bbox_inches='tight',
    format='pdf',
    transparent=True,
)

# Show plot
plt.show()
