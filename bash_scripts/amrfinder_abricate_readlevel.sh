#!/bin/bash

# Define paths
input_folder="./processing_downsampled/nanofilt"    # Directory containing input FASTQ files
output_folder="./processing_downsampled/fasta_output"  # Directory to store output FASTA files
database_path="/path/to/amrfinderplus/data"  # Path to the amrfinderplus database
abricate_output_dir="./processing/abricate"  # Directory to store abricate results

# Create output directories if they don't exist
mkdir -p "$output_folder"
mkdir -p "$abricate_output_dir"

# Convert FASTQ to FASTA using seqkit
for fastq_file in "$input_folder"/*.fastq; do
    # Define the output FASTA file path
    fasta_file="${fastq_file%.fastq}.fasta"
    fasta_file="${output_folder}/$(basename "$fasta_file")"
    
    # Log the conversion process
    echo "Converting $fastq_file to $fasta_file"
    
    # Convert FASTQ to FASTA
    seqkit fq2fa "$fastq_file" -o "$fasta_file"
done

# Run amrfinderplus on converted FASTA files
for fasta_file in "$output_folder"/*.fasta; do
    # Define the output file for amrfinderplus results
    output_file="${fasta_file%.fasta}_amrfinder.txt"
    
    # Log the amrfinderplus processing
    echo "Running amrfinderplus on $fasta_file..."
    
    # Run amrfinderplus
    amrfinder --threads 10 -n "$fasta_file" -d "$database_path" > "$output_file"
done

# Run abricate on converted FASTA files
for fasta_file in "$output_folder"/*.fasta; do
    # Get the filename without the path and extension
    filename=$(basename -- "$fasta_file" .fasta)
    
    # Define the output file path for abricate results
    abricate_output="${abricate_output_dir}/${filename}_abricate.txt"
    
    # Log the abricate processing
    echo "Running abricate on $fasta_file..."
    
    # Run abricate
    abricate "$fasta_file" > "$abricate_output"
done

# Log the completion of all processing
echo "All processing complete."
