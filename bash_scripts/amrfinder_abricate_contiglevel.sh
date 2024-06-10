#!/bin/bash

# Define paths
input_folder="./processing/racon"    # Directory containing input contig FASTA files
amrfinder_output_dir="./processing/amrfinder_results"  # Directory to store amrfinderplus results
abricate_output_dir="./processing/abricate_contig"  # Directory to store abricate results
amrfinder_db_path="/path/to/amrfinderplus/data"  # Path to the amrfinderplus database

# Create output directories if they don't exist
mkdir -p "$amrfinder_output_dir"
mkdir -p "$abricate_output_dir"

# Run amrfinderplus on contig FASTA files
for fasta_file in "$input_folder"/*.fasta; do
    # Get the filename without the path and extension
    filename=$(basename -- "$fasta_file" .fasta)
    
    # Define the output file for amrfinderplus results
    amrfinder_output="${amrfinder_output_dir}/${filename}_amrfinder.txt"
    
    # Log the amrfinderplus processing
    echo "Running amrfinderplus on $fasta_file..."
    
    # Run amrfinderplus
    amrfinder --threads 10 -n "$fasta_file" -d "$amrfinder_db_path" > "$amrfinder_output"
done

# Run abricate on contig FASTA files
for fasta_file in "$input_folder"/*.fasta; do
    # Get the filename without the path and extension
    filename=$(basename -- "$fasta_file" .fasta)
    
    # Define the output file for abricate results
    abricate_output="${abricate_output_dir}/${filename}_abricate.txt"
    
    # Log the abricate processing
    echo "Running abricate on $fasta_file..."
    
    # Run abricate
    abricate "$fasta_file" > "$abricate_output"
done

# Log the completion of all processing
echo "All processing complete."
