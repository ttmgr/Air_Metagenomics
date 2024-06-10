#!/bin/bash

# Number of CPUs to use
cpu_count=INT

# Define the eggNOG-mapper data directory where your databases are located
eggnog_data_dir="/path/to/eggnog/data"

# Define the output base directory
output_base_dir="./processing/eggnog_mappedprodigalcontigs"

# Ensure the output directory exists
mkdir -p "$output_base_dir"

# Run eggNOG-mapper on the Prodigal output proteins
for i in $(seq -w 1 24); do
  # Define the input file for the current barcode
  input_file="./processing/prodigal_metagenomic/prodigal_output_barcode${i}.faa"
  
  # Define the base name for the output files
  output_base="${output_base_dir}/eggnog_output_barcode${i}"

  # Check if the input file exists
  if [ -f "$input_file" ]; then
    # Log the processing of the current input file
    echo "Processing $input_file..."
    
    # Run eggNOG-mapper
    # -m diamond: use DIAMOND for fast sequence alignment
    # --data_dir: specify the directory containing eggNOG data
    # -d bactNOG: specify the database to use (bactNOG for bacterial sequences)
    # -i: input file
    # --output: base name for output files
    # --override: overwrite existing output files
    # --target_orthologs all: map to all orthologous groups
    # --query-cover 20: minimum query coverage
    # --subject-cover 20: minimum subject coverage
    # --tax_scope auto: automatically detect taxonomic scope
    # --cpu: number of CPUs to use
    emapper.py -m diamond \
      --data_dir "$eggnog_data_dir" \
      -d bactNOG \
      -i "$input_file" \
      --output "$output_base" \
      --override \
      --target_orthologs all \
      --query-cover 20 \
      --subject-cover 20 \
      --tax_scope auto \
      --cpu "$cpu_count"
  else
    # Log if the input file does not exist
    echo "File $input_file does not exist. Skipping."
  fi
done

# Log the completion of the eggNOG-mapper processing
echo "eggNOG-mapper processing completed."
