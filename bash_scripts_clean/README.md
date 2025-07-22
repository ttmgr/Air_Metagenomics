# Metagenomics Analysis Bash Pipeline

This directory contains a streamlined, modular bash pipeline for processing nanopore metagenomics data. The workflow is orchestrated by the main `run_pipeline.sh` script, which executes a series of numbered scripts in a logical order.

## Pipeline Workflow

The pipeline is divided into four main stages:

1.  **Read Processing (`01_read_processing.sh`)**: This script takes the demultiplexed FASTQ files as input and performs adapter trimming, quality/length filtering, generates read-level statistics, and performs initial taxonomic classification on the reads.
2.  **Assembly & Polishing (`02_assembly_and_polishing.sh`)**: This script takes the filtered reads from Stage 1, assembles them into contigs using Flye, polishes the assembly with Racon, and generates final assembly statistics.
3.  **Metagenome Binning (`03_binning.sh`)**: This script uses the polished contigs from Stage 2 and the filtered reads from Stage 1 to group contigs into high-quality Metagenome-Assembled Genomes (MAGs) using MetaWRAP.
4.  **Functional Annotation (`04_annotation.sh`)**: The final stage annotates both the polished assemblies from Stage 2 and the filtered reads from Stage 1. It predicts genes, assigns functions, and screens for antimicrobial resistance (AMR) and virulence factors.

## How to Run

1.  **Setup**: First, ensure all required databases are downloaded by running the `download_databases.sh` script. Make sure to edit the script to point to your desired database locations.
2.  **Configure**: Edit the `run_pipeline.sh` script to set the correct paths for your input data, output directory, and database locations.
3.  **Execute**: Run the main pipeline script from the root of the repository:
    ```bash
    bash bash_scripts/run_pipeline.sh
    ```

The main script will handle the creation of directories and execute each stage in the correct order.
