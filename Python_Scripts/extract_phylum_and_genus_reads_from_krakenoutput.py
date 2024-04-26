import os
import pandas as pd

# Specify the directory containing the report files
directory = '/path/to/kraken/report_files/'

# Initialize an empty DataFrame to collect data
data_list = []

# Loop through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        # Extract barcode correctly from the filename
        barcode = filename.split('_')[2]  # Assuming format 'report_filtered_barcode01_passed.txt'
        
        # Initialize counters for phylum and genus
        phylum_sum = 0
        genus_sum = 0
        
        # Read the file
        with open(os.path.join(directory, filename), 'r') as file:
            for line in file:
                parts = line.split('\t')
                if len(parts) >= 6:
                    rank = parts[3].strip()
                    count = int(parts[1].strip())
                    
                    if rank == 'P':  # Phylum level
                        phylum_sum += count
                    elif rank == 'G':  # Genus level
                        genus_sum += count
        
        # Append the results for phylum and genus sums to the data list
        data_list.append({'Barcode': barcode, 'Phylum_Reads': phylum_sum, 'Genus_Reads': genus_sum})

# Create a DataFrame from the list of dictionaries
data_df = pd.DataFrame(data_list)

# Save the collected data to a CSV file
data_df.to_csv('/path/to/results.csv', index=False)
