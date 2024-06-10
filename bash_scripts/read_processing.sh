#!/bin/bash

# Iterate over barcode files numbered from 01 to 24
for i in {01..24}
do
   # Construct the filename for the current barcode
   file="SQK-RBK114-24_barcode${i}.fastq"
   
   # Check if the file exists
   if [ -f "$file" ]; then
      # Run Porechop on the current barcode file
      # Porechop is a tool for removing adapter sequences from nanopore reads
      # The trimmed output is stored in the "porechop" directory
      porechop -i "$file" -o "./processing/porechop/trimmed_barcode${i}_passed.fastq"
   else
      # Print a message if the file does not exist
      echo "File $file does not exist. Skipping."
   fi
done

# Iterate over the trimmed files numbered from 01 to 24
for i in {01..24}
do
   # Construct the filename for the trimmed file
   trimmed_file="./processing/porechop/trimmed_barcode${i}_passed.fastq"
   
   # Check if the trimmed file exists
   if [ -f "$trimmed_file" ]; then
      # Run NanoFilt on the trimmed file to filter reads longer than 100 bp and store the output in the "nanofilt" directory
      cat "$trimmed_file" | NanoFilt -l 100 > "./processing/nanofilt/filtered_barcode${i}_passed.fastq"
   else
      # Print a message if the trimmed file does not exist
      echo "File $trimmed_file does not exist. Skipping."
   fi
done
