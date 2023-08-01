import pandas as pd
from Bio import Entrez
import time

# Always tell NCBI who you are
Entrez.email = "your_email@example.com"

def get_taxonomy(protein_id):
    """Fetch taxonomy for a given protein accession number"""
    search = Entrez.efetch(db="protein", id=protein_id, retmode="xml")
    record = Entrez.read(search)
    return record[0]["GBSeq_organism"], record[0]["GBSeq_taxonomy"]

# Load the TSV file
data = pd.read_csv("diamond_output.tsv", sep='\t', header=None)

# Define the column names
column_names = ["query_seq_id", "subject_seq_id", "pct_identical_matches", "alignment_length", "num_mismatches", 
                "num_gap_opens", "start_alignment_query", "end_alignment_query", "start_alignment_subject", 
                "end_alignment_subject", "e_value", "bit_score"]

# Add column names to the dataframe
data.columns = column_names

# Create new columns for organism and taxonomy
data["organism"] = ""
data["taxonomy"] = ""

# Get the total number of rows
total_rows = len(data)

# Go through each row in the dataframe
for i, row in data.iterrows():
    # Get the subject sequence ID
    subject_seq_id = row["subject_seq_id"]

    print(f"Processing row {i+1} of {total_rows}: {subject_seq_id}")

    # Get the taxonomy
    try:
        organism, taxonomy = get_taxonomy(subject_seq_id)
        data.at[i, "organism"] = organism
        data.at[i, "taxonomy"] = taxonomy
        print(f"  Fetched taxonomy for {subject_seq_id}")
    except Exception as e:
        print(f"  Couldn't fetch taxonomy for {subject_seq_id}: {e}")

    # To prevent overwhelming the NCBI servers, sleep for a little bit between each request
    time.sleep(0.2)

# Save the dataframe to a new TSV file
data.to_csv("diamond_output_with_taxonomy.tsv", sep='\t', index=False)

print("Finished processing all rows.")
