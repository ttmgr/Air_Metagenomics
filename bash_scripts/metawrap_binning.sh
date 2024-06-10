#!/bin/bash

# Directories for input and output
nanofilt_dir="./processing/nanofilt"
racon_dir="./processing/racon"
binning_output_base="./processing/binning_output"

# Ensure minimap2, samtools, and metawrap are installed and accessible in your environment

# Loop through each FASTQ file in the nanofilt directory
for fastq_file in ${nanofilt_dir}/filtered_barcode*.fastq; do
    # Extract barcode number from the file name
    barcode=$(echo "$fastq_file" | grep -oP 'barcode\d+' | head -1)
    
    # Corresponding FASTA file in the racon directory
    fasta_file="${racon_dir}/polished_${barcode}.fasta"
    
    # Check if the corresponding FASTA file exists
    if [ ! -f "$fasta_file" ]; then
        echo "Corresponding FASTA file for $fastq_file not found: $fasta_file"
        continue
    fi

    # Create binning output directory for this barcode
    binning_output="${binning_output_base}/${barcode}"
    mkdir -p "$binning_output/work_files"

    # Align reads to the assembly with Minimap2 and process with Samtools
    # Minimap2: align reads to the assembly
    # -ax map-ont: specify Oxford Nanopore read mapping
    # -L: log level
    # -t: number of threads
    minimap2 -ax map-ont -L -t 20 "$fasta_file" "$fastq_file" > "$binning_output/work_files/${barcode}.sam"
    
    # Samtools: sort and index the SAM file
    # sort: sort the BAM file
    # -T: temporary file prefix
    # -@: number of threads
    # -O: output format (BAM)
    samtools sort -T tmp-samtools -@ 20 -O BAM -o "$binning_output/work_files/${barcode}.bam" "$binning_output/work_files/${barcode}.sam"
    
    # Index the BAM file
    samtools index -@ 20 -b "$binning_output/work_files/${barcode}.bam"

    # Run the binning module with MetaWRAP
    # metawrap binning: binning of metagenomic contigs
    # --metabat2, --maxbin2, --concoct: specify binning tools
    # -t: number of threads
    # -m: memory in GB
    # --single-end: single-end reads
    # --universal: universal mode
    # --run-checkm: run CheckM for bin quality estimation
    # -l: minimum contig length
    # -a: assembly file
    # -o: output directory
    metawrap binning --metabat2 --maxbin2 --concoct -t 20 -m 64 --single-end --universal --run-checkm -l 10000 -a "$fasta_file" -o "$binning_output" "$fastq_file"

    # Refine bins with MetaWRAP
    # metawrap bin_refinement: refine bins
    # -o: output directory
    # -t: number of threads
    # -A, -B, -C: directories of bins from different tools
    # -c: minimum bin completion percentage
    # -x: maximum bin contamination percentage
    metawrap bin_refinement -o "$binning_output/BIN_REFINEMENT" -t 20 -A "$binning_output/metabat2_bins/" -B "$binning_output/maxbin2_bins/" -C "$binning_output/concoct_bins/" -c 50 -x 10

done

echo "Processing complete."
