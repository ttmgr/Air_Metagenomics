```
from ete3 import NCBITaxa
import pandas as pd

ncbi = NCBITaxa()

# Load the data
df = pd.read_csv('classification_wimp_v2-v1(2).csv')

# Define the taxonomic ranks you're interested in
ranks = ['phylum', 'class', 'order', 'family', 'genus']

# Initialize columns for each rank
for rank in ranks:
    df[rank] = ""

# Fill in the taxonomic ranks for each taxon
for i, taxon in df['taxID'].items():
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
df.to_csv('updated_classification.csv', index=False)
