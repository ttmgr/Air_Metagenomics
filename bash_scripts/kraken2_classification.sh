#!/bin/bash

# Create directories for output
mkdir -p processing/kraken2_read_classification
mkdir -p processing/kraken2_contig_classification

# Specify Kraken2 DB path
kraken2_db_path="/path/to/kraken2/database"

# Loop over all .fastq files in ./processing/nanofilt
for fastq_file in ./processing/nanofilt/*.fastq; do
    # Extract the base name of the file without the directory and extension
    base_name=$(basename -- "$fastq_file")
    base_name_no_ext="${base_name%.*}"

    # Run Kraken2 on each .fastq file
    if [ -f "$fastq_file" ]; then
        # Kraken2 is a tool for assigning taxonomic labels to short DNA sequences
        # --db: specifies the path to the Kraken2 database
        # --use-names: outputs scientific names instead of taxonomic IDs
        # --report: specifies the path to save the report file
        # --output: specifies the path to save the classification output
        # --memory-mapping: enables memory mapping to handle large databases
        # --threads: specifies the number of threads to use
        kraken2 --db "${kraken2_db_path}" --use-names --report "./processing/kraken2_read_classification/report_${base_name_no_ext}.txt" --output "./processing/kraken2_read_classification/output_${base_name_no_ext}.txt" "$fastq_file" --memory-mapping --threads 28
    else
        # Print a message if the file does not exist
        echo "File $fastq_file does not exist. Skipping."
    fi
done

# Loop over all contig files generated by Flye in ./processing/flye
for contig_file in ./processing/flye/*/assembly.fasta; do
    # Extract the directory name and base name of the file without the directory and extension
    dir_name=$(dirname -- "$contig_file")
    base_name=$(basename -- "$dir_name")
    
    # Run Kraken2 on each contig file
    if [ -f "$contig_file" ]; then
        # Kraken2 classification for contig files
        kraken2 --db "${kraken2_db_path}" --use-names --report "./processing/kraken2_contig_classification/report_${base_name}_contigs.txt" --output "./processing/kraken2_contig_classification/output_${base_name}_contigs.txt" "$contig_file" --memory-mapping --threads 28
    else
        # Print a message if the file does not exist
        echo "File $contig_file does not exist. Skipping."
    fi
done