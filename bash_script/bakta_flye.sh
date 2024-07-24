#!/bin/bash

# Define the path to the Bakta database
BAKTA_DB="/path/to/bakta_db/"
# Define the output directory
OUTPUT_DIR="./processing/bakta_flye"
# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop over all the barcode directories
for i in $(seq -w 1 24); do
    # Define the input FASTA file path
    INPUT_FILE="processing/flye/assembly_filtered_barcode${i}/assembly.fasta"
    
    # Check if the input file exists
    if [ -f "$INPUT_FILE" ]; then
        # Define the output directory for this barcode
        BARCODE_OUTPUT_DIR="${OUTPUT_DIR}/barcode${i}"
        # Create the output directory if it doesn't exist
        mkdir -p $BARCODE_OUTPUT_DIR
        
        # Run Bakta
        bakta --db $BAKTA_DB --output $BARCODE_OUTPUT_DIR $INPUT_FILE --threads 15 --force --skip-plot
    else
        echo "File $INPUT_FILE does not exist."
    fi
done
