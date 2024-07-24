#!/bin/bash

# Define the path to the Bakta database
BAKTA_DB="/path/to/bakta_db"

# Create the output directory if it doesn't exist
RACON_OUTPUT_DIR="./processing/bakta_racon"
mkdir -p $RACON_OUTPUT_DIR

# Function to process files
process_files() {
    local input_dir=$1
    local output_dir=$2
    local file_pattern=$3

    for input_file in ${input_dir}/${file_pattern}; do
        if [ -f "$input_file" ]; then
            local base_name=$(basename "$input_file" .fasta)
            local barcode_output_dir="${output_dir}/${base_name}"
            mkdir -p $barcode_output_dir
            bakta --db $BAKTA_DB --output $barcode_output_dir $input_file --threads 15 --force --skip-plot
        else
            echo "File $input_file does not exist."
        fi
    done
}

# Process files from racon
process_files "./processing/racon" $RACON_OUTPUT_DIR "*.fasta"
