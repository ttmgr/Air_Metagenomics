import os
import csv
import time
from Bio import Entrez

# Define the function to fetch genus and phylum for a given TaxID
def fetch_tax_info(tax_id, email, api_key):
    """ Fetch genus and phylum for a given TaxID from NCBI. """
    Entrez.email = email
    Entrez.api_key = api_key
    try:
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        genus = phylum = None
        for lineage in records[0]['LineageEx']:
            if lineage['Rank'] == 'genus':
                genus = lineage['ScientificName']
            elif lineage['Rank'] == 'phylum':
                phylum = lineage['ScientificName']
        return genus, phylum
    except Exception as e:
        print(f"Error fetching data for TaxID {tax_id}: {e}")
        return None, None

# Define the function to save a checkpoint
def save_checkpoint(file_path, line_num):
    """ Save the current line number to a checkpoint file. """
    with open(file_path, 'w') as f:
        f.write(str(line_num))

# Define the function to load a checkpoint
def load_checkpoint(file_path):
    """ Load the last line number from a checkpoint file. """
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            return int(f.read().strip())
    return 0

# Define the function to process each TaxID and fetch taxonomic information
def process_tax_ids(input_file, output_folder, checkpoints_folder, email, api_key):
    """ Process each TaxID from the input file to fetch genus and phylum. """
    checkpoint_file = os.path.join(checkpoints_folder, 'tax_id_processing.checkpoint')
    start_line = load_checkpoint(checkpoint_file)

    with open(input_file, 'r') as infile:
        tax_ids = infile.read().strip().split('\n')
        outfile_path = os.path.join(output_folder, 'tax_id_taxonomy.csv')
        
        if start_line == 0 and os.path.exists(outfile_path):
            os.remove(outfile_path)  # Remove existing output file if starting from scratch

        with open(outfile_path, 'a', newline='') as outfile:
            writer = csv.writer(outfile)
            if start_line == 0:
                writer.writerow(['TaxID', 'Genus', 'Phylum'])

            for i, tax_id in enumerate(tax_ids, start=start_line):
                print(f"Processing TaxID {tax_id}")
                genus, phylum = fetch_tax_info(tax_id, email, api_key)
                writer.writerow([tax_id, genus, phylum])

                if (i + 1) % 100 == 0:
                    save_checkpoint(checkpoint_file, i + 1)
                    print(f"Checkpoint saved: {i + 1} TaxIDs processed")
                    time.sleep(5)  # Pause to avoid rate limiting

# Define the main function to process TaxIDs from the input file
def main(input_file, output_folder, checkpoints_folder, email, api_key):
    """ Main function to process TaxIDs from the input file. """
    print("Processing TaxIDs...")
    process_tax_ids(input_file, output_folder, checkpoints_folder, email, api_key)
    print("All TaxIDs processed.")

# Set your parameters here
input_file = "/input/path/unique_tax_ids.txt"
output_folder = "/output/path"
checkpoints_folder = "/output/path"  # Create this directory for storing checkpoint files
email = "insert_email_registered_with_ncbi"
api_key = "insert_api_key"

main(input_file, output_folder, checkpoints_folder, email, api_key)
