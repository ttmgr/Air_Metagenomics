import matplotlib.pyplot as plt
import numpy as np

# Function to extract read lengths from a FASTQ file
def extract_read_lengths(fastq_file_path):
    read_lengths = []
    with open(fastq_file_path, 'r') as file:
        for i, line in enumerate(file):
            if i % 4 == 1:  # Sequence lines are every 4th line starting from the second line (0-based indexing)
                read_lengths.append(len(line.strip()))
    return read_lengths

# Path to your FASTQ file
fastq_file_path = '/path/to/readfile/reads.fastq'  # Update this path to your FASTQ file

# Extract read lengths
read_lengths = extract_read_lengths(fastq_file_path)

# Convert read lengths to numpy array for efficient computation
read_lengths_array = np.array(read_lengths)

# Plotting
plt.figure(figsize=(10, 6))

# Histogram with light green color
plt.hist(read_lengths_array, bins=np.logspace(np.log10(min(read_lengths_array)), np.log10(max(read_lengths_array)), 80), color='lightgreen', edgecolor='black')

# Setting the x-axis to log scale
plt.xscale('log')

# Calculate and plot the median read length
median_read_length = np.median(read_lengths_array)
plt.axvline(median_read_length, color='black', linestyle='--', label=f'Median: {median_read_length:.2f} bp')

# Adding labels, title, and legend
plt.xlabel('Read Length (bp, log scale)')
plt.ylabel('Number of Reads')
plt.title('Distribution of Read Lengths')
plt.legend()

# Display the plot
plt.show()
