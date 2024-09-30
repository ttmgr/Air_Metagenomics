# Air Monitoring by Nanopore Sequencing

## Project Overview

This project utilizes metagenomic analysis of bioaerosols to monitor air quality, employing the latest advancements in nanopore sequencing technology. Our workflow is designed to extract, process, and analyze genetic material from air samples to understand microbial compositions.

## Key Stages

1. Base Calling with Guppy and Dorado
2. Processing of Reads
3. Assembly of Metagenomes
4. Binning
5. Quality Control
6. Taxonomic Classification
7. Protein Classification and Functional Annotation

## Pipeline Components

Our analysis pipeline integrates various bioinformatics tools for comprehensive metagenomic analysis. The main components are:

- Sequencing Data Processing
- Metagenomic Assembly
- Genome Binning
- Quality Assessment
- Taxonomic Classification
- Functional Annotation
- Antimicrobial Resistance Gene Detection

## Repository Structure

- `/Functionality`: Contains detailed information about functional annotation tools.
- `/Taxonomy`: Provides explanations and usage guides for taxonomic classification tools.
- `Installation_tutorial.md`: Step-by-step guide for installing all required tools and databases.

## Tools Used

1. Guppy Basecaller
2. Porechop
3. NanoFilt
4. Flye
5. Minimap2
6. Racon
7. MetaWRAP
8. CheckM
9. Kraken 2
10. Prodigal
11. DIAMOND
12. Prokka
13. Seqkit
14. AMRFinderPlus
15. ABRicate

For detailed explanations of functionality and taxonomy tools, please refer to the respective directories.

## Installation

For comprehensive installation instructions for all tools and required databases, please refer to the `Installation_tutorial.md` file in the main directory.

## Usage

This pipeline is designed for the analysis of metagenomic data from air samples sequenced using Oxford Nanopore technology. Follow these steps to use the pipeline:

1. **Base Calling**:
   Use Guppy or Dorado to convert raw electrical signals to DNA sequences.
   ```
   guppy_basecaller -i /path/to/raw_data -s /path/to/output --config dna_r9.4.1_450bps_hac.cfg
   ```

   **Dorado Base Calling** (alternative to Guppy):
    If using Dorado instead of Guppy for base calling.
    ```
    dorado basecaller --model /path/to/model.json /path/to/pod5_files > basecalled_reads.bam
    ```

2. **Read Processing**:
   Trim adapters and filter reads based on quality and length.
   ```
   porechop -i input.fastq -o trimmed.fastq
   NanoFilt -q 10 -l 1000 < trimmed.fastq > filtered.fastq
   ```

3. **Metagenome Assembly**:
   Assemble the filtered reads into contigs.
   ```
   flye --nano-raw filtered.fastq --out-dir assembly_output --meta
   ```

4. **Read Alignment and Assembly Polishing**:
   Align reads to the assembly and polish it.
   ```
   minimap2 -ax map-ont assembly.fasta filtered.fastq | samtools sort > aligned.bam
   racon filtered.fastq aligned.bam assembly.fasta > polished_assembly.fasta
   ```

5. **Binning**:
   Use MetaWRAP to bin the assembled contigs.
   ```
   metawrap binning -o BINNING_OUT -t 16 -a polished_assembly.fasta --metabat2 --maxbin2 --concoct filtered.fastq
   ```

6. **Quality Control**:
   Assess the quality of the bins using CheckM.
   ```
   checkm lineage_wf -t 16 -x fasta BINNING_OUT/bin_refinement/metawrap_70_10_bins ./checkm_output
   ```

7. **Taxonomic Classification**:
   Classify the sequences using Kraken2.
   ```
   kraken2 --db /path/to/kraken_db --threads 16 --output kraken_output.txt --report kraken_report.txt filtered.fastq
   ```

8. **Gene Prediction**:
   Predict genes using Prodigal.
   ```
   prodigal -i polished_assembly.fasta -a proteins.faa -d genes.fna -o genes.gff
   ```

9. **Protein Alignment**:
   Align predicted proteins against a reference database using DIAMOND.
   ```
   diamond blastp --db /path/to/nr_database --query proteins.faa --out diamond_output.txt --outfmt 6 --evalue 1e-5 --threads 16
   ```

10. **Genome Annotation**:
    Annotate the assembled genomes using Prokka and Bakta.
    ```
    prokka --outdir prokka_output --prefix sample_name polished_assembly.fasta
    bakta --db /path/to/bakta_db --output bakta_output polished_assembly.fasta
    ```

11. **Sequence Manipulation**:
    Use Seqkit for various sequence operations.
    ```
    seqkit stats *.fastq > sequencing_stats.txt
    seqkit seq -m 1000 filtered.fastq > filtered_1k.fastq
    ```

12. **AMR Gene Detection**:
    Screen for antimicrobial resistance genes.
    ```
    abricate --db card polished_assembly.fasta > abricate_output.txt
    amrfinder -n polished_assembly.fasta -o amrfinder_output.txt
    ```

13. **Functional Annotation**:
    Annotate genes with eggNOG-mapper.
    ```
    emapper.py -i proteins.faa -o eggnog_output --cpu 16
    ```

14. **Protein Family Analysis**:
    Identify protein families using PfamScan.
    ```
    pfam_scan.pl -fasta proteins.faa -dir /path/to/pfam/database -outfile pfam_output.txt
    ```

15. **Generate Sequencing Statistics**:
    Produce summary statistics for your sequencing data.
    ```
    NanoStat --fastq filtered.fastq --outdir nanostat_output
    ```

16. **Generate Assembly Statistics**:
    Produce summary statistics for your assembly.
    ```
    assembly-stats polished_assembly.fasta > assembly_stats.txt
    ```

Note: 
- Replace `/path/to/` with the actual paths on your system.
- Adjust thread counts and other parameters based on your available computational resources.
- Refer to the individual tool documentation in the Functionality and Taxonomy directories for more detailed usage instructions and parameter optimization.

For a more detailed walkthrough of the entire pipeline, including input preparation and output interpretation, please refer to our comprehensive user guide [link_to_user_guide].

[Remaining sections stay the same]
