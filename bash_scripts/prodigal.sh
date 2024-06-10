#!/bin/bash

# Make directories for output
mkdir -p processing/prodigal_metagenomic

# Run Prodigal on the metagenomic assemblies (FASTA files)
for i in $(seq -f "%02g" 1 24); do
    # Input FASTA file for the current barcode
    input_file="./processing/racon/polished_barcode${i}.fasta"
    
    # Output GFF file for Prodigal predictions
    output_gff="./processing/prodigal_metagenomic/prodigal_output_barcode${i}.gff"
    
    # Output FAA file for Prodigal protein translations
    output_faa="./processing/prodigal_metagenomic/prodigal_output_barcode${i}.faa"
    
    # Check if the input FASTA file exists
    if [ -f "$input_file" ]; then
        # Log the processing of the current input file
        echo "Processing $input_file..."
        
        # Run Prodigal for gene prediction
        # -i: input FASTA file
        # -o: output GFF file for predicted genes
        # -a: output FAA file for translated proteins
        # -p meta: use metagenomic mode
        prodigal -i "$input_file" -o "$output_gff" -a "$output_faa" -p meta
    else
        # Log if the input file does not exist
        echo "File $input_file does not exist. Skipping."
    fi
done
