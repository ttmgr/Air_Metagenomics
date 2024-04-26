import os
import pandas as pd

# Specify the directory containing the nanostat report files
directory = '/path/to/nanostat_results/'

# Initialize a list to collect data
data_list = []

# Loop through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        # Extract barcode from filename
        barcode = filename.split('_')[1]  # Assuming format 'filtered_barcode01_passed_nanostats.txt'
        
        # Initialize stats for this file
        file_stats = {
            'Barcode': barcode,
            'Mean_Read_Length': 0,
            'Median_Read_Length': 0,
            'Read_Length_N50': 0,
            'Number_of_Reads': 0
        }
        
        # Read the file
        with open(os.path.join(directory, filename), 'r') as file:
            for line in file:
                if 'Mean read length:' in line:
                    file_stats['Mean_Read_Length'] = float(line.split(':')[1].strip().replace(',', ''))
                elif 'Median read length:' in line:
                    file_stats['Median_Read_Length'] = float(line.split(':')[1].strip().replace(',', ''))
                elif 'Read length N50:' in line:
                    file_stats['Read_Length_N50'] = float(line.split(':')[1].strip().replace(',', ''))
                elif 'Number of reads:' in line:
                    file_stats['Number_of_Reads'] = float(line.split(':')[1].strip().replace(',', ''))

        # Append the stats for this file to the list
        data_list.append(file_stats)

# Convert the list of dictionaries to a DataFrame
data_df = pd.DataFrame(data_list)

# Save the DataFrame to a CSV file
data_df.to_csv('/path/to/results/results_summary.csv', index=False)
