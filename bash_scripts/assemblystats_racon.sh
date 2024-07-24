#!/bin/bash

# Define the directories
FASTA_DIR="./processing/racon"
OUTPUT_DIR="./processing/assemblystats"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over each .fasta file in the directory
for fasta_file in "$FASTA_DIR"/*.fasta; do
    # Check if the file exists
    if [[ -f "$fasta_file" ]]; then
        # Extract the base name of the fasta file (without the directory and extension)
        base_name=$(basename "$fasta_file" .fasta)
        # Run assembly-stats and save the result to a .txt file in the output directory
        assembly-stats "$fasta_file" > "$OUTPUT_DIR/${base_name}_assemblystats.txt"
    fi
done
