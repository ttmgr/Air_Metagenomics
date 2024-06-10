#!/bin/bash

# Define the absolute path to the dorado binary
DORADO_BIN="/path/to/dorado/exe"

# Define the absolute path to the configuration file
CONFIG_FILE="/path/to/dorado/mmodel/file/dna_r10.4.1_e8.2_400bps_hac@v4.3.0"

# Step 1: Run dorado basecaller
echo "Running dorado basecaller..."
# dorado basecaller: runs the basecalling process to convert raw nanopore signal data to nucleotide sequences
# --emit-fastq: specifies that the output should be in FASTQ format
# -r: specifies the directory containing the raw data
# --kit-name: specifies the sequencing kit used
# --no-trim: disables the trimming of adapter sequences
${DORADO_BIN} basecaller --emit-fastq ${CONFIG_FILE} -r 20231220_1339_MC-113930_FAU98948_259115e5/ > basecalled.fastq --kit-name SQK-RBK114-24 --no-trim

# Check if basecalling was successful
if [ $? -ne 0 ]; then
    # $? checks the exit status of the last command executed
    # A non-zero exit status indicates an error
    echo "Error in basecalling step."
    exit 1
fi

# Step 2: Run dorado demux
echo "Running dorado demux..."
# dorado demux: performs demultiplexing to separate reads based on their barcodes
# --output-dir: specifies the directory where the demultiplexed files should be saved
# --emit-fastq: specifies that the output should be in FASTQ format
# --kit-name: specifies the sequencing kit used
${DORADO_BIN} demux --output-dir basecalled/ --emit-fastq --kit-name SQK-RBK114-24 basecalled.fastq

# Check if demultiplexing was successful
if [ $? -ne 0 ]; then
    # A non-zero exit status indicates an error
    echo "Error in demultiplexing step."
    exit 1
fi

echo "Process completed successfully."

