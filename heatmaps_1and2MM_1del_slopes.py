# This script generates the heatmaps using slope files as input files. The slope files can be generated using "slopes_linear_fraction.py"
# The 1 nt deletion sequences are extracted separately using "single_nt_del.py" and should be the second input file
# Be sure to specify input/output directory and the input files below
# Note that the first input filename must be something like Fn_W_1_lin_slopes.txt, i.e. must start with ortholog_gene_time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.colors as colors
import re
import os


# Define the directory and input file name
directory_in = '' # specify input directory
directory_out = '' # specify output directory
input_file1 = '' # specify input file containing slopes for all sequences
input_file2 = '' # specify input file containing slopes for 1 nt deletion sequences
input_path1 = os.path.join(directory_in, input_file1)
input_path2 = os.path.join(directory_in, input_file2)
filename_parts = input_file1.split("_")
gene_name = filename_parts[0]
target = filename_parts[1]
time = filename_parts[2]
concentration = filename_parts[4].replace("mM.txt", " mM")

# Create the directory if it doesn't exist
if not os.path.exists(directory_out):
    os.makedirs(directory_out)


# Define the reference sequence
W_ref_seq = "TTTACGGCCACTTCCGTGTCTGAC" #W reference sequence
L_ref_seq = "TTTGTACTGTCCACGCCGACGGAA" #L reference sequence
if input_file1.split('_')[1] == 'L':
    ref_seq = L_ref_seq
elif input_file1.split('_')[1] == 'W':
    ref_seq = W_ref_seq
else:
    raise ValueError("Invalid file name: second identifier should be 'L' or 'W'")

# Read the input file into a dictionary
seq_dict = {}
with open(input_path1, 'r') as f:
    for line in f:
        line = line.strip()
        if line:
            split_line = line.split('\t')
            if len(split_line) == 2:
                seq, val = split_line
                if val == '':
                    seq_dict[seq] = np.nan
                else:
                    seq_dict[seq] = float(val)
            elif len(split_line) == 1:
                seq = split_line[0]
                seq_dict[seq] = np.nan
            else:
                raise ValueError("Invalid line format: expected 1 or 2 values separated by a tab")

data = np.zeros((len(ref_seq), len(ref_seq), 3, 3))

# Read the 1 deletion values from the second column of input file2
values = []
with open(input_path2, 'r') as infile:
    for line in infile:
        columns = line.strip().split('\t')
        if len(columns) >= 2:
            if columns[1] == '':
                value = np.nan
            else:
                value = float(columns[1])
            values.append(value)

# Fill the 3D array with the slope values from the dictionary
for i in range(len(ref_seq)):
    # 1 mutation sequences
    for k, nuc1 in enumerate([x for x in ["A", "C", "G", "T"] if x != ref_seq[i]]):
        for l in range(3):
            data[i, 0, k, l] = np.nan
            seq_1mm = ref_seq[:i] + nuc1 + ref_seq[i+1:]
            slope_list = []
            for seq, val in seq_dict.items():
                if re.search(seq_1mm, seq):
                    slope_list.append(val)
            if slope_list:
                slope = slope_list[-1]
                if slope == 0:
                    slope = np.nan
            else:
                slope = np.nan
            data[i, 0, k, k] = slope

    # 2 mutation sequences
    for j in range(i+1, len(ref_seq)):
        for k, nuc1 in enumerate([x for x in ["A", "C", "G", "T"] if x != ref_seq[i]]):
            for l, nuc2 in enumerate([x for x in ["A", "C", "G", "T"] if x != ref_seq[j]]):
                seq_2mm = ref_seq[:i] + nuc1 + ref_seq[i+1:j] + nuc2 + ref_seq[j+1:]
                slope_list = []
                for seq, val in seq_dict.items():
                    if re.search(seq_2mm, seq):
                        slope_list.append(val)
                if slope_list:
                    slope = slope_list[-1]
                    if slope == 0:
                        slope = np.nan
                else:
                    slope = np.nan
                data[i, j, k, l] = slope
                   
# Create a figure and axis object
plt.rcParams['font.family'] = 'Arial'
fig, ax = plt.subplots(figsize=(15,15), sharex=True)
ax.set_aspect('equal')

# Create a color map
cmap = plt.get_cmap('coolwarm')
cmap.set_bad('white')

# Define the color scale
norm = colors.Normalize(vmin=-0.1, vmax=0.1)

# Create a color map for the nucleotide colors
nucleotide_colors = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}

# Define variables for creating space between boxes in heatmap
offset_x = 0
increment = 0.2

# Loop through the 3x3 sub-arrays and plot each one as a heatmap cell
# 2 mutation sequences
for i in range(23):
    offset_y = offset_x
    for j in range(24-i):
        if j > 0:
            for k in range(3):
                for l in range(3):
                    x_inner = i + (k-1) / 3 + 1/3 + 0.025 + offset_x # Add a small offset to create space
                    y_inner = (i + j) + (l-1) / 3 + 1/3 + 0.025 + offset_y # Add a small offset to create space
                    value = data[i, (i+j), k, l]
                    color = cmap(norm(value))
                    rect_inner = plt.Rectangle((x_inner, y_inner), 1/3 - 0.05, 1/3 - 0.05, facecolor=color)
                    ax.add_patch(rect_inner)                                 
        offset_y += increment
    offset_x += increment

# 1 mutation sequences
offset = 0
for i in range(24):
    for k in range(3):
        for l in range(3):
            if k == l:
                x_inner = i + (k-1) / 3 + 1/3 + 0.025 + offset # Add a small offset to create space
                y_inner = i + (l-1) / 3 + 1/3 + 0.025 + offset # Add a small offset to create space
                value = data[i, 0, k, l]
                color = cmap(norm(value))
                rect_inner = plt.Rectangle((x_inner, y_inner), 1/3 - 0.05, 1/3 - 0.05, facecolor=color) 
                ax.add_patch(rect_inner)
    offset += increment

# Create 1D heatmap for the 1 nt deletions
offset=0
for i in range(24):
    value = values[i]
    color = cmap(norm(value))
    rect = patches.Rectangle((i + offset, 29.3), 1, 1, facecolor=color)
    ax.add_patch(rect)
    offset += 0.2
    
# # Add label for 1D heatmaps for the 1 nt deletions
# ax.text(-4,30.3,"1 nt Δ",color='black', fontname='Arial', fontsize = 32)
    
# Add squares on the x-axis indicating the mutation type
offset = 0
for i in range(24):
    for k, nuc in enumerate([x for x in ["A", "C", "G", "T"] if x != ref_seq[i]]):
        x_square = i + (k-1) / 3 + 1/3 + 0.025 + offset
        rect_x = patches.Rectangle((x_square, 28.825), 1/3 - 0.05, 1/3 - 0.05, facecolor=nucleotide_colors[nuc])
        ax.add_patch(rect_x)
    offset += increment

# Add squares on the y-axis indicating the mutation type
offset = 0
for i in range(24):
    for l, nuc in enumerate([x for x in ["A", "C", "G", "T"] if x != ref_seq[i]]):
        y_square = i + (l-1) / 3 + 1/3 + 0.025 + offset
        rect_y = patches.Rectangle((-1/3-increment, y_square), 1/3 - 0.05, 1/3 - 0.05, facecolor=nucleotide_colors[nuc])
        ax.add_patch(rect_y)
    offset += increment
        
plt.xlim(-2,32)
plt.ylim(-2,32)

# Invert the y-axis to plot the arrays from top to bottom
ax.invert_yaxis()

# Remove spines
for spine in ax.spines.values():
    spine.set_visible(False)

# Remove ticks and labels
ax.set_xticks([])
ax.set_yticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])

# Add all labels to x-axis
# Set the x-axis tick labels to be centered below the 3x3 arrays
offset_x = 0
xticks = []
xticklabels = []
for i in range(24):
    xticks.append(i + offset_x + 0.425)
    xticklabels.append(ref_seq[i])
    offset_x += increment

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, fontname='Courier New', weight='bold', fontsize=36, bbox = None)

# Get the current tick positions and labels
ticks = ax.get_xticks()
ticklabels = ax.get_xticklabels()

# Set the desired vertical position for the tick marks (e.g., 0.1 for 10% from the bottom)
vertical_position = 0.043

# Adjust the vertical position of the tick marks
for tick, label in zip(ticks, ticklabels):
    label.set_y(vertical_position)
    
# Add rectangles under x-axis tick labels
pam_rect = patches.Rectangle((0, 30.5),4.75,1.5, facecolor = (1, 0.859, 0.090))
ax.add_patch(pam_rect)
seed_rect = patches.Rectangle((4.65, 30.5),9.6,1.5, facecolor = (0.710, 0.871, 0.969))
ax.add_patch(seed_rect)
mid_rect = patches.Rectangle((14.2, 30.5),6.2,1.5, facecolor = (0.816, 1, 0.831))
ax.add_patch(mid_rect)
distal_rect = patches.Rectangle((20.28, 30.5),8.3,1.5, facecolor = (0.945, 0.816, 1))
ax.add_patch(distal_rect)
    
# Add numbers under x-axis tick labels
offset_x = 0.4
for i in range(-4, 21):
    if i != 0:
        ax.text(offset_x,32.5,i,ha='center',va='center',color='black', fontname='Arial', fontsize = 18)
        offset_x += 1.2
        
# Add labels under x-axis    
ax.text(2.25,33.7,"PAM",ha='center',va='center',color='black', fontname='Arial', fontsize = 36)
ax.text(9.4,33.7,"seed",ha='center',va='center',color='black', fontname='Arial', fontsize = 36)
ax.text(17.2,33.7,"mid",ha='center',va='center',color='black', fontname='Arial', fontsize = 36)
ax.text(24.5,33.7,"distal",ha='center',va='center',color='black', fontname='Arial', fontsize = 36)
ax.text(14.3,35.25,"Δ or first mismatch position",ha='center',va='center',color='black', fontname='Arial', fontsize = 36)

# Add all labels to y-axis
# Set the y-axis tick labels to be centered below the 3x3 arrays
offset_y = 0
yticks = []
yticklabels = []
for i in range(24):
    yticks.append(i + offset_y + 0.425)
    yticklabels.append(ref_seq[i])
    offset_y += increment

ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels, fontname='Courier New', weight='bold', fontsize=36)

# Get the current tick positions and labels
ticks = ax.get_yticks()
ticklabels = ax.get_yticklabels()

# Set the desired vertical position for the tick marks (e.g., 0.1 for 10% from the bottom)
horizontal_position = 0.04

# Adjust the vertical position of the tick marks
for tick, label in zip(ticks, ticklabels):
    label.set_x(horizontal_position)

# Remove tick marks
ax.tick_params(bottom=False, top=False, left=False, right=False)

# Add rectangles under y-axis tick labels
pam_rect = patches.Rectangle((-2.25, 0), 1.5, 4.95, facecolor = (1, 0.859, 0.090))
ax.add_patch(pam_rect)
seed_rect = patches.Rectangle((-2.25, 4.7),1.5, 9.8, facecolor = (0.710, 0.871, 0.969))
ax.add_patch(seed_rect)
mid_rect = patches.Rectangle((-2.25, 14.3),1.5, 6,facecolor = (0.816, 1, 0.831))
ax.add_patch(mid_rect)
distal_rect = patches.Rectangle((-2.25, 20.25),1.5, 8.4,facecolor = (0.945, 0.816, 1))
ax.add_patch(distal_rect)

# Add numbers next to y-axis tick labels
offset_y = 0.5
for i in range(-4, 21):
    if i != 0:
        ax.text(-2.5,offset_y,i,ha='center',va='center',color='black', fontname='Arial', fontsize = 18)
        offset_y += 1.2

# Add labels next to y-axis    
ax.text(-3.75,2.25,"PAM",ha='center',va='center',color='black', fontname='Arial', fontsize = 36, rotation = 90)
ax.text(-3.75,9.5,"seed",ha='center',va='center',color='black', fontname='Arial', fontsize = 36, rotation = 90)
ax.text(-3.75,17.25,"mid",ha='center',va='center',color='black', fontname='Arial', fontsize = 36, rotation = 90)
ax.text(-3.75,24.25,"distal",ha='center',va='center',color='black', fontname='Arial', fontsize = 36, rotation = 90)
ax.text(-5.5,14.5,"second mismatch position",ha='center',va='center',color='black', fontname='Arial', fontsize = 36, rotation = 90)



# Create the upper and lower labels
upper_label = f"{gene_name}Cas12a gene {target} target {time} min"

# Add the labels to the heatmap
ax.text(15, 0, upper_label, ha='center', va='bottom', fontsize=36)

# Adjust the plot layout to accommodate the labels
fig.subplots_adjust(top=0.85, bottom=0.15)  # Adjust the values as needed

# Set the tick parameters
ax.tick_params(axis='both', which='both', length=0)
ax.tick_params(axis='x', which='major', pad=2)


# Set the color bar and its label
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-0.1, vmax=0.1))
sm.set_array([])
cbar_ax = fig.add_axes([0.5, 0.65, 0.25, 0.02])  # Adjust the coordinates as needed
cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
cbar.outline.set_linewidth(2)

# Set the ticks and labels for the color bar
ticks = [-0.1, 0, 0.1]
labels = [str(t) for t in ticks]
cbar.set_ticks(ticks)
cbar.set_ticklabels(labels, fontsize = 36)
cbar.ax.tick_params(width=2, length = 10)

# Set the label for the color bar
label_text = 'slope'
cbar.ax.text(0, 1.5, label_text, ha='center', va='bottom', fontsize=36)


# Create the legend patches and labels
legend_patches = []
legend_labels = []

box_size = 1.25
x_position = 16
y_position = 11
nuc_x_position = x_position
for nucleotide, color in nucleotide_colors.items():
    # Create a square patch for each nucleotide color
    patch = patches.Rectangle((nuc_x_position, y_position), box_size, box_size, facecolor=color)
    ax.add_patch(patch)

    # Add the label at the center of the box
    label_x = nuc_x_position + box_size / 2
    label_y = y_position + box_size / 1.75
    ax.text(label_x, label_y, nucleotide, ha='center', va='center', color='white', fontname='Courier New', weight='bold', fontsize=36)  
    nuc_x_position += box_size + 0.25
    
# Create a square patch for missing sequences
patch = patches.Rectangle((x_position, y_position+1.5), box_size, box_size, facecolor='white')
ax.add_patch(patch)

# Add the label at the center of the box
ax.text(x_position + 1.5, y_position + 2.65, 'missing', color='black', fontname='Arial', fontsize=36)  

# Show the plot
plt.show()

# Save the plot as png and svg in the input_dir
filename = os.path.splitext(os.path.basename(input_file1))[0]
fig.savefig(os.path.join(directory_out, filename + '_HM.pdf'))
# output_file = os.path.join(directory_out, filename + '_HM.png')
# plt.savefig(output_file, dpi=300, bbox_inches='tight')
