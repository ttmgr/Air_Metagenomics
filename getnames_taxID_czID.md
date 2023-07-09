from ete3 import NCBITaxa
import pandas as pd
import os

ncbi = NCBITaxa()

# Define the directory that contains the CSV files
csv_dir = '.'

# Define the taxonomic ranks you're interested in
ranks = ['phylum', 'class', 'order', 'family', 'genus']

# Loop through all the CSV files in the directory
for filename in os.listdir(csv_dir):
    if filename.endswith(".csv"):
        filepath = os.path.join(csv_dir, filename)
        
        # Load the data
        df = pd.read_csv(filepath)

        # Initialize columns for each rank
        for rank in ranks:
            df[rank] = ""

        # Fill in the taxonomic ranks for each taxon
        for i, taxon in df['tax_id'].items():
            try:
                lineage = ncbi.get_lineage(taxon)
                names = ncbi.get_taxid_translator(lineage)
                lineage2ranks = ncbi.get_rank(names.keys())
                
                for taxid in lineage:
                    if lineage2ranks[taxid] in ranks:
                        df.at[i, lineage2ranks[taxid]] = names[taxid]
            except:
                print(f"Unable to fetch data for taxID: {taxon}")

        # Save the updated DataFrame to a new CSV file
        df.to_csv(filepath, index=False)
