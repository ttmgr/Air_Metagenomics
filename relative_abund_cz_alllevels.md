```
import pandas as pd

# Load the combined data
data = pd.read_csv('/path/to/your/folder/combined_data.csv')

# Define the number of top entries to select for each level
top_entries = {'genus': 20, 'order': 10, 'class': 10, 'family': 10, 'phylum': 5}

# For each level, group the data by both sample and level, then write the top most abundant taxa to a CSV file
for level, top in top_entries.items():
    result = []
    # Process each sample
    for sample, sample_data in data.groupby('sample'):
        # Group by level, then calculate the sum of relative_abundance
        grouped_data = sample_data.groupby(level)['relative_abundance'].sum().reset_index()

        # Sort by relative_abundance in descending order
        grouped_data = grouped_data.sort_values(by='relative_abundance', ascending=False)

        # Add a new column for adjusted level, where entries not in the top are labeled as "Rare"
        grouped_data[f'{level}_adjusted'] = grouped_data[level].where(grouped_data[level].isin(grouped_data[level].head(top)), 'Rare')

        # Group by the adjusted level, then calculate the sum of relative_abundance
        grouped_data_adjusted = grouped_data.groupby(f'{level}_adjusted')['relative_abundance'].sum().reset_index()

        # Normalize relative_abundance to sum up to 1
        grouped_data_adjusted['relative_abundance'] = grouped_data_adjusted['relative_abundance'] / grouped_data_adjusted['relative_abundance'].sum()

        # Add the sample name to the dataframe
        grouped_data_adjusted['sample'] = sample

        # Append the dataframe to the result list
        result.append(grouped_data_adjusted)

    # Concatenate all dataframes in the result list
    result_df = pd.concat(result)

    # Write the adjusted data to a CSV file
    result_df.to_csv(f'/path/to/your/folder/top_{level}_by_sample.csv', index=False)
