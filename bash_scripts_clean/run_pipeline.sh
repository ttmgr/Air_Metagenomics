#!/bin/bash
# Main orchestrator for the metagenomics pipeline.
# This script executes each stage of the analysis in a sequential order.

set -euo pipefail

# --- User Configuration ---
# Set the base directory for all analysis outputs
OUTPUT_BASE_DIR="./processing"

# Set the path to your demultiplexed FASTQ files (e.g., from Guppy or Dorado)
# The pipeline assumes files are named like 'barcode01.fastq', 'barcode02.fastq', etc.
INPUT_FASTQ_DIR="./raw_fastq"

# Set the number of threads to use for multithreaded tasks
THREADS=24

# Set paths to your downloaded databases
# These should match the paths in download_databases.sh
KRAKEN2_DB_PATH="/path/to/kraken2/database"
AMRFINDER_DB_PATH="/path/to/amrfinderplus/data"
BAKTA_DB_PATH="/path/to/bakta/db"
EGGNOG_DATA_DIR="/path/to/eggnog/data"

# Export variables to be used by subscripts
export OUTPUT_BASE_DIR
export INPUT_FASTQ_DIR
export THREADS
export KRAKEN2_DB_PATH
export AMRFINDER_DB_PATH
export BAKTA_DB_PATH
export EGGNOG_DATA_DIR

# --- Pipeline Execution ---

echo "--- Starting Metagenomics Pipeline ---"

# Create the base output directory
mkdir -p "$OUTPUT_BASE_DIR"
echo "Output will be written to: $OUTPUT_BASE_DIR"

# Get the path of the currently executing script
SCRIPT_DIR=$(dirname "$0")

# --- Stage 1: Read Processing ---
echo "--- Stage 1: Running Read Processing ---"
bash "${SCRIPT_DIR}/01_read_processing.sh"
echo "--- Stage 1: Completed ---"

# --- Stage 2: Assembly and Polishing ---
echo "--- Stage 2: Running Assembly and Polishing ---"
bash "${SCRIPT_DIR}/02_assembly_and_polishing.sh"
echo "--- Stage 2: Completed ---"

# --- Stage 3: Metagenome Binning ---
echo "--- Stage 3: Running Metagenome Binning ---"
bash "${SCRIPT_DIR}/03_binning.sh"
echo "--- Stage 3: Completed ---"

# --- Stage 4: Functional Annotation ---
echo "--- Stage 4: Running Functional Annotation ---"
bash "${SCRIPT_DIR}/04_annotation.sh"
echo "--- Stage 4: Completed ---"

echo "--- Metagenomics Pipeline Finished Successfully ---"
