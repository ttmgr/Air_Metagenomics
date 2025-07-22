import os
import pandas as pd
import argparse
import time
from Bio import Entrez

def extract_nanostat_metrics(directory):
    """Extracts metrics from NanoStat report files."""
    data_list = []
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            barcode = filename.split('_')[1]
            file_stats = {'Barcode': barcode, 'Mean_Read_Length': 0, 'Median_Read_Length': 0, 'Read_Length_N50': 0, 'Number_of_Reads': 0}
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
            data_list.append(file_stats)
    return pd.DataFrame(data_list)

def extract_kraken_reads(directory):
    """Extracts phylum and genus read counts from Kraken2 report files."""
    data_list = []
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            barcode = filename.split('_')[2]
            phylum_sum = 0
            genus_sum = 0
            with open(os.path.join(directory, filename), 'r') as file:
                for line in file:
                    parts = line.split('\t')
                    if len(parts) >= 6:
                        rank = parts[3].strip()
                        count = int(parts[1].strip())
                        if rank == 'P':
                            phylum_sum += count
                        elif rank == 'G':
                            genus_sum += count
            data_list.append({'Barcode': barcode, 'Phylum_Reads': phylum_sum, 'Genus_Reads': genus_sum})
    return pd.DataFrame(data_list)

def process_kraken_report(file_path):
    """Processes a single Kraken2 report file."""
    df = pd.read_csv(file_path, sep='\t', header=None, names=["Percentage", "Reads", "Reads_at_Level", "Rank", "TaxID", "Name"])
    df['Name'] = df['Name'].str.strip()
    df_filtered = df[(df['Rank'] != 'U') & (df['Name'] != 'Homo')]
    df_taxa = df_filtered[df_filtered['Rank'].isin(['G', 'P'])]
    taxa_reads = df_taxa.groupby(['Rank', 'Name'])['Reads'].sum()
    total_taxa_reads = taxa_reads.sum(level='Rank')
    relative_abundances = taxa_reads / total_taxa_reads * 100
    df_relative_abundances = relative_abundances.reset_index()
    df_relative_abundances.columns = ['Rank.code', 'Name', 'Relative_Abundance']
    return df_relative_abundances[df_relative_abundances['Relative_Abundance'] >= 1]

def process_all_kraken_reports(directory_path):
    """Processes all Kraken2 report files in a directory."""
    combined_df = pd.DataFrame()
    for filename in os.listdir(directory_path):
        if filename.endswith(".txt"):
            file_path = os.path.join(directory_path, filename)
            df = process_kraken_report(file_path)
            df['Barcode'] = filename.split('.')[0]
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    return combined_df

def fetch_tax_info(tax_id, email, api_key):
    """Fetches genus and phylum for a given TaxID from NCBI."""
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

def process_tax_ids(input_file, output_folder, checkpoints_folder, email, api_key):
    """Processes each TaxID from the input file to fetch genus and phylum."""
    # ... (rest of the function from fetch_taxonomy_for_diamond.py)

def main():
    parser = argparse.ArgumentParser(description="Process metagenomic data.")
    parser.add_argument('--task', required=True, choices=['nanostat', 'kraken_reads', 'kraken_reports', 'fetch_taxonomy'], help="The task to perform.")
    parser.add_argument('--input_dir', help="Input directory for nanostat, kraken_reads, and kraken_reports tasks.")
    parser.add_argument('--output_dir', required=True, help="Output directory for the results.")
    # Add other arguments as needed

    args = parser.parse_args()

    if args.task == 'nanostat':
        if not args.input_dir:
            parser.error("--input_dir is required for nanostat task.")
        df = extract_nanostat_metrics(args.input_dir)
        df.to_csv(os.path.join(args.output_dir, 'nanostat_summary.csv'), index=False)
        print("NanoStat metrics extracted successfully.")
    # Add other tasks here

if __name__ == '__main__':
    main()
