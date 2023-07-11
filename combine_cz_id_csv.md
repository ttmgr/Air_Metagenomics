```
import pandas as pd
import os

# Define the folder where your processed files are located
folder = '.'

# Get a list of all files in the folder
files = os.listdir(folder)

# Create an empty list to hold dataframes
dfs = []

# Process each file
for file in files:
    # Ensure we're only processing the processed .csv files
    if file.endswith('_processed.csv'):
        # Generate the sample name from the file name
        sample_name = file.replace('_filtered.concat_42198', '').replace('_taxon_report_processed.csv', '')
        
        # Load the data
        df = pd.read_csv(os.path.join(folder, file))
        
        # Add a new column with the sample name
        df['sample'] = sample_name
        
        # Append the dataframe to the list
        dfs.append(df)

# Combine all dataframes
combined_df = pd.concat(dfs)

# Save the combined data to a new file
combined_df.to_csv('./combined_data.csv', index=False)
