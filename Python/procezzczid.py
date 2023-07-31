```
import pandas as pd
import os

# Function to process a single file
def process_file(filename):
    # Load the data
    data = pd.read_csv(filename)

    # Retain only the specified columns
    data_filtered = data[['tax_id', 'tax_level', 'nt_base_count', 'phylum', 'class', 'order', 'family', 'genus']]

    # Calculate the total nt_base_count for tax_level 1 and 2
    total_nt_base_count_level1 = data_filtered[data_filtered['tax_level'] == 1]['nt_base_count'].sum()
    total_nt_base_count_level2 = data_filtered[data_filtered['tax_level'] == 2]['nt_base_count'].sum()

    # Add a new column "relative_abundance"
    data_filtered.loc[data_filtered['tax_level'] == 1, 'relative_abundance'] = data_filtered.loc[data_filtered['tax_level'] == 1, 'nt_base_count'] / total_nt_base_count_level1
    data_filtered.loc[data_filtered['tax_level'] == 2, 'relative_abundance'] = data_filtered.loc[data_filtered['tax_level'] == 2, 'nt_base_count'] / total_nt_base_count_level2

    # Save the processed data to a new file
    output_filename = filename.replace('.csv', '_processed.csv')
    data_filtered.to_csv(output_filename, index=False)

    print(f"Processed {filename} and saved to {output_filename}")


# Define the folder where your files are located
folder = '/path/to/your/folder'

# Get a list of all files in the folder
files = os.listdir(folder)

# Process each file
for file in files:
    # Ensure we're only processing .csv files
    if file.endswith('.csv'):
        process_file(os.path.join(folder, file))
