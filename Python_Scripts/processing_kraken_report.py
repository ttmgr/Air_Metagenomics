import pandas as pd
import os

def process_kraken2_file(file_path):
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None, names=["Percentage", "Reads", "Reads_at_Level", "Rank", "TaxID", "Name"])
    df['Name'] = df['Name'].str.strip()

    # Filter out unclassified reads and reads matching the genus "Homo"
    df_filtered = df[(df['Rank'] != 'U') & (df['Name'] != 'Homo')]

    # Filtering for genera and phyla (Rank codes 'G' and 'P')
    df_taxa = df_filtered[df_filtered['Rank'].isin(['G', 'P'])]

    # Summing up the reads at the genus and phylum levels
    taxa_reads = df_taxa.groupby(['Rank', 'Name'])['Reads'].sum()

    # Calculate total reads for calculating relative abundance
    total_taxa_reads = taxa_reads.sum(level='Rank')

    # Calculate relative abundances
    relative_abundances = taxa_reads / total_taxa_reads * 100

    # Creating a DataFrame for the relative abundances
    df_relative_abundances = relative_abundances.reset_index()
    df_relative_abundances.columns = ['Rank.code', 'Name', 'Relative_Abundance']

    # Filtering out taxa with relative abundance below 1%
    df_relative_abundances_filtered = df_relative_abundances[df_relative_abundances['Relative_Abundance'] >= 1]

    return df_relative_abundances_filtered

def process_all_files_in_directory(directory_path):
    combined_df = pd.DataFrame()
    for filename in os.listdir(directory_path):
        if filename.endswith(".txt"):  # Assuming all files are .txt files
            file_path = os.path.join(directory_path, filename)
            df = process_kraken2_file(file_path)
            df['Barcode'] = filename.split('.')[0]  # Extracting barcode name from filename
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    return combined_df

# Usage
directory_path = '/path/to/kraken_report_files/'  # Replace with your directory path
combined_results = process_all_files_in_directory(directory_path)

# Saving the combined results to a new file
combined_results.to_csv('/path/to/output_dir/combined_kraken2_results.csv', index=False)
