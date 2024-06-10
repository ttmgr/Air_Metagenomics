#!/bin/bash

# Create a directory for Prokka output
mkdir -p processing/prokka_annotations

# Base directory for input files
input_dir="./processing/racon"

# Output base directory
output_dir="./processing/prokka_annotations"

# Loop through all polished FASTA files in the racon directory
for fasta_file in ${input_dir}/polished_barcode*.fasta; do
    # Check if the input FASTA file exists
    if [ -f "$fasta_file" ]; then
        # Extract the barcode number (or identifier) from the filename
        barcode=$(basename "$fasta_file" .fasta | sed 's/polished_barcode//')
        
        # Define Prokka output directory uniquely for each input file
        prokka_outdir="${output_dir}/prokka_barcode${barcode}"
        
        # Create a directory for Prokka output for each barcode
        mkdir -p "$prokka_outdir"
        
        # Log the processing of the current input file
        echo "Running Prokka on $fasta_file..."
        
        # Run Prokka for annotation
        # --outdir: specify the output directory
        # --prefix: prefix for output files
        prokka --outdir "$prokka_outdir" --prefix "barcode${barcode}" "$fasta_file"
    else
        # Log if the input file does not exist
        echo "File $fasta_file does not exist. Skipping."
    fi
done

# Log the completion of the annotation process
echo "Prokka annotation completed."
