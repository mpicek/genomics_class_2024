import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Function to load depth data from a file
def load_depth_data(file_path):
    data_frame = pd.read_csv(file_path, delimiter='\t', header=None, names=['chromosome', 'position', 'depth'])
    return data_frame

# Function to compute average depth across specified base pair intervals
def compute_average_depth(data_frame):
    interval = 1000
    data_frame['interval'] = data_frame['position'] // interval
    average_depth = data_frame.groupby(['chromosome', 'interval'])['depth'].mean().reset_index()
    average_depth['center_position'] = average_depth['interval'] * interval + interval // 2
    return average_depth

# Function to compute the log2 ratio of tumor to wild-type depths
def log2_ratio(tumor_data, normal_data):
    combined_data = pd.merge(tumor_data, normal_data, on=['chromosome', 'center_position'], suffixes=('_tumor', '_normal'))
    combined_data['log2_ratio'] = np.log2(combined_data['depth_tumor'] / combined_data['depth_normal'])
    return combined_data

# Function to plot the log2 ratio and save to a file
def plot_log2_ratio(data_frame):
    plt.figure(figsize=(10, 6))
    plt.plot(data_frame['center_position'], data_frame['log2_ratio'], color='red', linestyle='-', marker='o', markersize=2)
    plt.xlabel('Genomic Position')
    plt.ylabel('Log2 Ratio (Tumor vs. Normal)')
    plt.title('Depth of Read Plot')
    plt.grid(True)
    plt.savefig('output/read_depth_plot.png', bbox_inches='tight')  # Save the plot as a PNG file
    plt.close()

# Main variables for file paths
tumor_file = 'output/tu_read_depth.csv'
normal_file = 'output/wt_read_depth.csv'

# Load, process, and plot data
tumor_df = load_depth_data(tumor_file)
normal_df = load_depth_data(normal_file)

average_tumor_depth = compute_average_depth(tumor_df)
average_normal_depth = compute_average_depth(normal_df)

ratio_data = log2_ratio(average_tumor_depth, average_normal_depth)

plot_log2_ratio(ratio_data)