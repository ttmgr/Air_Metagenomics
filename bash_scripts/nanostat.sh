#!/bin/bash

# Define directories
fastq_dir="./processing/nanofilt/"  # Directory containing input FASTQ files
output_dir="./processing/nanostat/"  # Directory to store NanoStat results

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Log the start of the NanoStat processing
echo "Running NanoStat on .fastq files..."

# Loop through each FASTQ file in the input directory
for fastq_file in "${fastq_dir}"*.fastq; do
    # Check if the FASTQ file exists
    if [ -f "$fastq_file" ]; then
        # Define the output file path
        output_file="${output_dir}$(basename "${fastq_file%.fastq}")_nanostats.txt"
        
        # Log the processing of the current file
        echo "Processing $fastq_file..."
        
        # Run NanoStat on the FASTQ file and save the output to the output directory
        NanoStat --fastq "$fastq_file" > "$output_file"
    else
        # Log if the FASTQ file does not exist
        echo "File $fastq_file does not exist. Skipping."
    fi
done

# Log the completion of the NanoStat processing
echo "NanoStat processing completed."
