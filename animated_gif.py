# This script can be used to compile animated gifs from the individual pngs created by heatmaps_1and2MM_1del_individual.py

from PIL import Image
import glob
import os

# Specify the input and output directories
input_dir = '' # specify input directory
output_dir = '' # specify output directory

# Specify the prefix of the frame filenames
frame_prefix = 'As_L_1_nicked' # specify prefix of all pngs to be compiled into gif 

# Get the file paths of the frames
frame_paths = glob.glob(os.path.join(input_dir, f'{frame_prefix}*.png'))

# Sort the frame paths based on the Mg concentration designator
sorted_frame_paths = sorted(frame_paths, key=lambda path: int(path.split('_')[-2].split('mM')[0]))

# Create a list to store the image frames
frames = []

# Load each frame, crop it, convert to RGB mode, and append it to the frames list
for frame_path in sorted_frame_paths:
    frame = Image.open(frame_path)
    rgb_frame = frame.convert('RGB')
    frames.append(rgb_frame)

# Create a reversed copy of the frames list (excluding the first and last frames)
#reversed_frames = frames[-2:0:-1] # for forward and backward
reversed_frames = frames[::-1] # for backward only

# Concatenate the original frames and the reversed frames
#all_frames = frames + reversed_frames # for forward and backward
#all_frames = frames # for forward only
all_frames = reversed_frames # for backward only

# Specify the output file path using the frame_prefix
output_file = os.path.join(output_dir, f'{frame_prefix}.gif')

# Specify the frame rate (in milliseconds)
frame_rate = 1000

# Save the frames as an animated GIF
all_frames[0].save(output_file, format='GIF', append_images=all_frames[1:], save_all=True, 
                   duration=frame_rate, loop=0, dither=Image.NONE, optimize=False, quality=95)
