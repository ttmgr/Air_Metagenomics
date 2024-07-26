#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 09:53:20 2024

Generic script for plotting read counts from a CSV file.

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(file_path, delimiter=';'):
    """
    Load dataset from a CSV file.

    Parameters:
    file_path (str): Path to the CSV file.
    delimiter (str): Delimiter used in the CSV file.

    Returns:
    DataFrame: Loaded data.
    """
    data = pd.read_csv(file_path, delimiter=delimiter)
    data['Location'] = data['Sample ID'].str.extract(r'(^.*?)(?: \d)')[0]
    return data

def create_plots(data, ylim_max):
    """
    Create bar plots for each location.

    Parameters:
    data (DataFrame): Data to plot.
    ylim_max (int): Maximum limit for y-axis.
    """
    sns.set(style="whitegrid")
    for location in data['Location'].unique():
        location_data = data[data['Location'] == location]
        location_data['Short Sample ID'] = location_data['Sample ID'].apply(lambda x: x.replace(location + ' ', ''))
        
        plt.figure(figsize=(14, 10))
        ax = location_data.set_index('Short Sample ID')[['Reads after filtering', 'Reads mapped to phylum', 'Reads mapped to genus', 'Reads mapped to top 20 genera']].plot(kind='bar')
        
        ax.set_title(f'{location} Read Counts')
        ax.set_xlabel(location)
        ax.set_ylabel('Number of Reads')
        ax.set_ylim(0, ylim_max)
        ax.axhline(y=25000, color=sns.color_palette()[1], linestyle='--', label='Downsampling Threshold (25k reads)')
        ax.legend(title='Metrics')
        ax.grid(False)
        
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.show()

def main(data_file_path):
    """
    Main function to load data and create plots.

    Parameters:
    data_file_path (str): Path to the CSV file containing the data.
    """
    # Load the dataset
    data = load_data(data_file_path)
    
    # Determine the maximum y-axis limit
    ylim_max = data[['Reads after filtering', 'Reads mapped to phylum', 'Reads mapped to genus', 'Reads mapped to top 20 genera']].max().max()
    
    # Create plots
    create_plots(data, ylim_max)

if __name__ == "__main__":
    # Example usage: update this path with your actual file path
    data_file_path = 'path_to_csv_file.csv'
    main(data_file_path)
