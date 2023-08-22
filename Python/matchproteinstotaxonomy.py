import os
import pandas as pd
from Bio import Entrez
import time
import socket
import numpy as np

Entrez.email = "your@email"

def is_connected():
    """Check if there is an internet connection"""
    try:
        socket.create_connection(("www.google.com", 80))
        return True
    except OSError:
        pass
    return False

def get_taxonomy(protein_id):
    """Fetch taxonomy for a given protein accession number"""
    search = Entrez.efetch(db="protein", id=protein_id, retmode="xml")
    record = Entrez.read(search)
    return record[0]["GBSeq_organism"], record[0]["GBSeq_taxonomy"]

# Define sample name
sample_name = "diamondoutput_samplename"

# Load data from the original file
data = pd.read_csv(f"{sample_name}.tsv", sep='\t', header=None)
data.columns = ["query_seq_id", "subject_seq_id", "pct_identical_matches", "alignment_length", "num_mismatches", 
                "num_gap_opens", "start_alignment_query", "end_alignment_query", "start_alignment_subject", 
                "end_alignment_subject", "e_value", "bit_score"]
data["organism"] = np.nan
data["taxonomy"] = np.nan

# Load the checkpoint file if it exists
if os.path.exists(f"{sample_name}_checkpoint.tsv"):
    checkpoint_data = pd.read_csv(f"{sample_name}_checkpoint.tsv", sep='\t')
else:
    checkpoint_data = data.copy()

total_rows = len(data)

for i, row in checkpoint_data.iterrows():
    subject_seq_id = row["subject_seq_id"]

    print(f"Processing row {i+1} of {total_rows}: {subject_seq_id}")

    # Skip rows that have already been processed
    if pd.notna(row["organism"]) and pd.notna(row["taxonomy"]):
        continue

    while not is_connected():
        print("  Waiting for internet connection...")
        time.sleep(5)

    try:
        organism, taxonomy = get_taxonomy(subject_seq_id)
        checkpoint_data.at[i, "organism"] = organism
        checkpoint_data.at[i, "taxonomy"] = taxonomy
        print(f"  Fetched taxonomy for {subject_seq_id}")

        # Save progress after each successful fetch
        checkpoint_data.to_csv(f"{sample_name}_checkpoint.tsv", sep='\t', index=False)
    except Exception as e:
        print(f"  Couldn't fetch taxonomy for {subject_seq_id}: {e}")

# Save the final result
checkpoint_data.to_csv(f"{sample_name}_with_taxonomy.tsv", sep='\t', index=False)

print("Finished processing all rows.")

