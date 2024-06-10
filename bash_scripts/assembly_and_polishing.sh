#!/bin/bash

# Run Flye on the filtered files in the nanofilt directory
for fq in processing/nanofilt/*.fastq; do
    # Extract the base name of the file without the "_passed.fastq" suffix
    base_name=$(basename "$fq" "_passed.fastq")
    
    # Check if the filtered file exists
    if [ -f "$fq" ]; then
        # Flye is a de novo assembler for long reads
        # --meta: metagenome assembly mode
        # --nano-hq: specify that input is high-quality nanopore reads
        # -o: specify the output directory
        # --threads: number of threads to use
        flye --meta --nano-hq "$fq" -o "processing/flye/assembly_${base_name}" --threads 15
    else
        # Print a message if the file does not exist
        echo "File $fq does not exist. Skipping."
    fi
done

# Run minimap2 on the filtered reads against the respective assemblies
for fq in processing/nanofilt/*.fastq; do
    # Extract the base name of the file without the "_passed.fastq" suffix
    base_name=$(basename "$fq" "_passed.fastq")
    
    # Check if both the filtered file and the corresponding assembly exist
    if [ -f "$fq" ] && [ -f "processing/flye/assembly_${base_name}/assembly.fasta" ]; then
        # minimap2 is a fast aligner for long nucleotide sequences
        # -ax: specify the preset for mapping Oxford Nanopore reads
        # map-ont: preset for Oxford Nanopore reads
        # -t: number of threads to use
        minimap2 -ax map-ont -t 15 "processing/flye/assembly_${base_name}/assembly.fasta" "$fq" > "processing/minimap2/aligned_${base_name}.sam"
    else
        # Print a message if one or more required files do not exist
        echo "One or more required files for minimap2 do not exist. Skipping ${base_name}."
    fi
done

# Convert SAM to sorted BAM using samtools
for sam_file in processing/minimap2/*.sam; do
    # Extract the base name of the file without the ".sam" suffix
    base_name=$(basename "$sam_file" ".sam")
    
    # Check if the SAM file exists
    if [ -f "$sam_file" ]; then
        # samtools is a suite of programs for interacting with high-throughput sequencing data
        # view: convert SAM to BAM
        # -@: number of threads to use
        # -b: output in BAM format
        # sort: sort BAM file
        samtools view -@ 15 -b "$sam_file" | samtools sort -@ 15 -o "processing/samtools/sorted_${base_name}.bam"
    else
        # Print a message if the SAM file does not exist
        echo "File $sam_file does not exist. Skipping."
    fi
done

# Run Racon on the sorted BAM files, the assemblies, and the filtered reads
# Note: Racon by default uses all available threads, but you can specify it with the -t option if needed.
for fq in processing/nanofilt/*.fastq; do
    # Extract the base name of the file without the "_passed.fastq" suffix
    base_name=$(basename "$fq" "_passed.fastq")
    
    # Check if all required files for Racon exist
    if [ -f "$fq" ] && [ -f "processing/minimap2/aligned_${base_name}.sam" ] && [ -f "processing/flye/assembly_${base_name}/assembly.fasta" ]; then
        # Racon is a tool for consensus sequence generation
        # -t: number of threads to use
        racon -t 15 "$fq" "processing/minimap2/aligned_${base_name}.sam" "processing/flye/assembly_${base_name}/assembly.fasta" > "processing/racon/polished_${base_name}.fasta"
    else
        # Print a message if one or more required files do not exist
        echo "One or more required files for Racon do not exist. Skipping ${base_name}."
    fi
done
