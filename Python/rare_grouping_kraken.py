import pandas as pd

def compute_rare_stats(input_file, output_file):
    # Load the data
    df = pd.read_csv(input_file)

    # Identify the top 20 taxa by relative abundance across all samples
    top_taxa = df.nlargest(20, 'Relative Abundance')['5'].tolist()

    # Define a function to categorize taxa
    def categorize_taxa(row):
        if row['5'] in top_taxa or row['Relative Abundance'] >= 0.01:
            return 'Top'
        else:
            return 'Rare'

    # Apply the categorization function
    df['Category'] = df.apply(categorize_taxa, axis=1)

    # Filter for the "Rare" group and group by sample and rank to compute the sum of relative abundances
    rare_group = df[df['Category'] == 'Rare'].groupby(['Sample Name', 'Rank'])['Relative Abundance'].sum().reset_index()

    # Calculate median and standard deviation across all samples for each rank
    stats = rare_group.groupby('Rank')['Relative Abundance'].agg(['median', 'std']).reset_index()

    # Save the results
    stats.to_csv(output_file, index=False)
  
  # Compute stats for combined read abundances
compute_rare_stats("./processing/kraken2_stats/combined_read_abundances.csv", 
                   "./processing/kraken2_stats/rare_read_stats.csv")

# Compute stats for combined assembly abundances
compute_rare_stats("./processing/kraken2_stats/combined_assembly_abundances.csv", 
                   "./processing/kraken2_stats/rare_assembly_stats.csv")
