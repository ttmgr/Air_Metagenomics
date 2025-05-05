# Air Monitoring by Nanopore Sequencing (PRJEB76446)

## Project Overview

This project utilizes metagenomic analysis of bioaerosols collected via air sampling to monitor microbial communities, employing Oxford Nanopore Technologies (ONT) sequencing. The repository contains a workflow designed to process raw sequencing data, assemble metagenomes, perform binning, and conduct taxonomic and functional analysis.

## ENA Files: https://www.ebi.ac.uk/ena/browser/view/PRJEB76446 
---

**❗ Important First Step ❗**

Before using the general pipeline described below, you **must** consult the specific data guides for the dataset you are analyzing. These guides contain critical information about:

* **Data Access:** Where to find the raw (FAST5/POD5) and/or processed (FASTQ) data.
* **Sample Mapping:** Which barcode corresponds to which sample.
* **Basecalling/Demultiplexing:** Specific commands and configurations used for the initial conversion of raw data to FASTQ, including required demultiplexing based on barcodes.
* **Preprocessing:** Any necessary steps like concatenating paired FASTQ files for certain samples.

**Find the essential guides here:**

* **Greenhouse vs. Natural Environment Study (Pilot Study):** [`pilot_study_guide.md`](https://github.com/ttmgr/Air_Metagenomics/blob/main/pilot_study_guide.md)
    * Describes samples from Greenhouse and Natural environments with different collection times.
    * Uses **FAST5** data format and **Guppy** for basecalling/demultiplexing.
    * Requires **concatenation** of FASTQ pairs for Natural Environment samples.
* **Urban Air Study:** [`urban_study_sample_guide.md`](https://github.com/ttmgr/Air_Metagenomics/blob/main/urban_study_sample_guide.md)
    * Describes samples collected from various urban locations (City Center, Greenbelt, etc.).
    * Uses **POD5** data format and **Dorado** for basecalling/demultiplexing.

**The pipeline described below assumes you have already completed the necessary steps from the relevant guide and have demultiplexed, preprocessed FASTQ files ready for analysis (one file or concatenated file per sample).**

---

## Analysis Pipeline Overview

This pipeline processes individual FASTQ sample files through several key bioinformatics stages:

1.  **Read Processing:** Adapter trimming and quality/length filtering.
2.  **Metagenome Assembly:** Assembling reads into contiguous sequences (contigs).
3.  **Assembly Polishing:** Improving assembly accuracy using raw reads.
4.  **Metagenome Binning:** Grouping contigs into putative genome bins (MAGs).
5.  **Quality Control (QC):** Assessing the completeness and contamination of bins.
6.  **Taxonomic Classification:** Identifying the microbial composition of reads and/or bins.
7.  **Gene Prediction & Functional Annotation:** Finding genes and predicting their functions.
8.  **Antimicrobial Resistance (AMR) Gene Detection:** Screening for AMR genes.
9.  **Statistics Generation:** Calculating metrics for reads and assemblies.

## Tools Used

This pipeline integrates the following key bioinformatics tools:

1.  **Guppy / Dorado:** Basecallers (ONT) - *Initial processing covered in linked guides.*
2.  **Porechop:** Adapter trimming for ONT reads.
3.  **NanoFilt:** Quality and length filtering for ONT reads.
4.  **Flye:** Long-read assembler (meta-assembly mode).
5.  **Minimap2:** Long-read alignment.
6.  **Racon:** Consensus correction/polishing for assemblies.
7.  **MetaWRAP:** Meta-pipeline for binning (using MetaBAT2, MaxBin2, CONCOCT) and refinement.
8.  **CheckM:** Assesses the quality (completeness, contamination) of MAGs.
9.  **Kraken 2:** K-mer based taxonomic classification.
10. **Prodigal:** Gene prediction for prokaryotic genomes.
11. **DIAMOND:** High-throughput protein sequence alignment (BLASTp equivalent).
12. **Prokka / Bakta:** Prokaryotic genome annotation pipelines.
13. **Seqkit:** Toolkit for FASTA/Q sequence manipulation and statistics.
14. **AMRFinderPlus / ABRicate:** Tools for detecting acquired AMR genes.
15. **eggNOG-mapper:** Functional annotation based on orthologous groups.
16. **NanoStat:** Generates statistics for Nanopore sequencing data.
17. **Assembly-stats:** Calculates basic assembly metrics.

## Repository Structure

-   `/Functionality`: Contains detailed information about functional annotation tools (e.g., DIAMOND, eggNOG-mapper). *(Consider adding details here)*
-   `/Taxonomy`: Provides explanations and usage guides for taxonomic classification tools (e.g., Kraken 2). *(Consider adding details here)*
-   `Installation_tutorial.md`: Step-by-step guide for installing all required tools and databases.
-   `urban_study_sample_guide.md`: Essential guide for accessing and preprocessing the Urban Air study data.
-   `pilot_study_guide.md`: Essential guide for accessing and preprocessing the Greenhouse/Natural Environment study data.

## Installation

For comprehensive installation instructions for all pipeline tools and required databases, please refer to the [`Installation_tutorial.md`](Installation_tutorial.md) file.

## Usage Workflow (General Pipeline)

This section outlines the commands for processing **one sample** (i.e., one demultiplexed FASTQ file obtained via the steps in the linked guides). You will typically need to run these steps for each sample, potentially using loops or a workflow manager.

**Input:** A single FASTQ file per sample (e.g., `sample_X.fastq`). Remember Natural Environment samples require concatenation first (see `pilot_study_guide.md`).

1.  **Read Processing:**
    * Trim adapters (if not done during basecalling/demux).
    * Filter by quality (e.g., Q8) and length (e.g., min 100 bp).
    ```bash
    # Input: sample_X.fastq
    porechop -i sample_X.fastq -o sample_X.trimmed.fastq
    NanoFilt -q 8 -l 100 < sample_X.trimmed.fastq > sample_X.filtered.fastq
    # Output: sample_X.filtered.fastq
    ```

2.  **Metagenome Assembly:**
    * Assemble the filtered reads using Flye in meta mode.
    ```bash
    # Input: sample_X.filtered.fastq
    flye --nano-hq sample_X.filtered.fastq --out-dir sample_X_assembly --meta --threads 16
    # Output: sample_X_assembly/assembly.fasta
    ```

3.  **Assembly Polishing (Optional but Recommended):**
    * Align filtered reads back to the assembly.
    * Use Racon for polishing (can be run multiple rounds).
    ```bash
    # Input: sample_X_assembly/assembly.fasta, sample_X.filtered.fastq
    minimap2 -ax map-ont sample_X_assembly/assembly.fasta sample_X.filtered.fastq | samtools sort -@ 8 -o sample_X.aligned.bam
    samtools index sample_X.aligned.bam
    racon -t 16 sample_X.filtered.fastq sample_X.aligned.bam sample_X_assembly/assembly.fasta > sample_X.polished.fasta
    # Output: sample_X.polished.fasta (Repeat alignment and racon for more rounds if desired)
    ```
    *(Note: Subsequent steps use the polished assembly if available, otherwise use `sample_X_assembly/assembly.fasta`)*

4.  **Metagenome Binning:**
    * Use MetaWRAP's binning module with multiple binners.
    ```bash
    # Input: sample_X.polished.fasta (or assembly.fasta), sample_X.filtered.fastq
    ASSEMBLY_FA=sample_X.polished.fasta # Or sample_X_assembly/assembly.fasta
    FILTERED_FQ=sample_X.filtered.fastq
    metawrap binning -o sample_X_BINNING -t 16 -a $ASSEMBLY_FA --metabat2 --maxbin2 --concoct $FILTERED_FQ
    # Output: Bins in sample_X_BINNING/metabat2_bins/, maxbin2_bins/, concoct_bins/
    ```

5.  **Bin Consolidation & Refinement (MetaWRAP):**
    * Combine results from different binners and refine bins based on CheckM scores (e.g., Completeness > 70%, Contamination < 10%).
    ```bash
    # Input: Bins from previous step
    metawrap bin_refinement -o sample_X_BIN_REFINEMENT -t 16 -A sample_X_BINNING/metabat2_bins/ -B sample_X_BINNING/maxbin2_bins/ -C sample_X_BINNING/concoct_bins/ -c 70 -x 10
    # Output: Refined bins in sample_X_BIN_REFINEMENT/metawrap_70_10_bins/
    ```

6.  **Bin Quality Control:**
    * Assess the quality of the final, refined bins using CheckM.
    ```bash
    # Input: Refined bins
    REFINED_BINS_DIR=sample_X_BIN_REFINEMENT/metawrap_70_10_bins
    checkm lineage_wf -t 16 -x fasta $REFINED_BINS_DIR ./sample_X_checkm_output
    # Output: Quality report in sample_X_checkm_output/
    ```

7.  **Taxonomic Classification (Reads):**
    * Classify the processed reads using Kraken 2.
    ```bash
    # Input: sample_X.filtered.fastq
    kraken2 --db /path/to/kraken_db --threads 16 --output sample_X.kraken_output.txt --report sample_X.kraken_report.txt sample_X.filtered.fastq
    # Output: sample_X.kraken_output.txt, sample_X.kraken_report.txt
    ```
    *(Note: You can also run Kraken 2 on the assembled contigs or individual bins if needed)*

8.  **Gene Prediction (Assembly/Bins):**
    * Predict protein-coding genes on the assembly or high-quality bins.
    ```bash
    # Input: sample_X.polished.fasta (or individual bin .fa files)
    ASSEMBLY_FA=sample_X.polished.fasta
    prodigal -i $ASSEMBLY_FA -a sample_X.proteins.faa -d sample_X.genes.fna -o sample_X.genes.gff -p meta
    # Output: sample_X.proteins.faa, sample_X.genes.fna, sample_X.genes.gff
    ```

9.  **Protein Alignment (Functional Annotation):**
    * Align predicted proteins against a reference database (e.g., NCBI nr) using DIAMOND for initial functional hints.
    ```bash
    # Input: sample_X.proteins.faa
    diamond blastp --db /path/to/nr_database.dmnd --query sample_X.proteins.faa --out sample_X.diamond_output.txt --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --evalue 1e-5 --threads 16 --sensitive
    # Output: sample_X.diamond_output.txt
    ```

10. **Genome Annotation (Assembly/Bins):**
    * Annotate the assembly or bins using Prokka or Bakta (Bakta often preferred for MAGs).
    ```bash
    # Input: sample_X.polished.fasta (or individual bin .fa files)
    ASSEMBLY_FA=sample_X.polished.fasta
    # Using Prokka (Example)
    prokka --outdir sample_X_prokka_output --prefix sample_X --cpus 16 $ASSEMBLY_FA
    # Using Bakta (Example - replace --db path)
    # bakta --db /path/to/bakta_db --output sample_X_bakta_output --threads 16 $ASSEMBLY_FA
    # Output: Annotation files in respective output directories
    ```

11. **Sequence Manipulation (Example):**
    * Use Seqkit for various operations, e.g., getting stats or filtering.
    ```bash
    # Input: sample_X.filtered.fastq
    seqkit stats sample_X.filtered.fastq >> all_samples_filtered_stats.txt
    # Filter sequences longer than 5kbp
    # seqkit seq -m 5000 sample_X.filtered.fastq > sample_X.filtered_5k.fastq
    ```

12. **AMR Gene Detection:**
    * Screen the assembly or bins for antimicrobial resistance genes using ABRicate and/or AMRFinderPlus.
    ```bash
    # Input: sample_X.polished.fasta (or individual bin .fa files)
    ASSEMBLY_FA=sample_X.polished.fasta
    # Using ABRicate (Example with CARD db)
    abricate --db card --threads 16 $ASSEMBLY_FA > sample_X.abricate_card.txt
    # Using AMRFinderPlus (Example for contigs)
    amrfinder -n $ASSEMBLY_FA -o sample_X.amrfinder.txt --threads 16
    # Output: AMR gene reports
    ```

13. **Functional Annotation (eggNOG):**
    * Annotate predicted proteins using eggNOG-mapper. Requires downloaded database.
    ```bash
    # Input: sample_X.proteins.faa
    # Ensure EGGNOG_DATA_DIR environment variable is set or use --data_dir
    emapper.py -i sample_X.proteins.faa -o sample_X_eggnog --output_basename sample_X --cpu 16
    # Output: Annotation files in sample_X_eggnog/
    ```

14. **Generate Sequencing Statistics:**
    * Produce summary statistics for the processed reads.
    ```bash
    # Input: sample_X.filtered.fastq
    NanoStat --fastq sample_X.filtered.fastq --outdir sample_X_nanostat_output -n sample_X_nanostat_report.txt
    # Output: Reports in sample_X_nanostat_output/
    ```

15. **Generate Assembly Statistics:**
    * Produce summary statistics for the final assembly.
    ```bash
    # Input: sample_X.polished.fasta
    assembly-stats sample_X.polished.fasta > sample_X_assembly_stats.txt
    # Output: sample_X_assembly_stats.txt
    ```

**Notes:**

* Replace placeholders like `/path/to/`, `sample_X`, `<N>`, `<M>`, database paths, etc., with your actual paths and desired parameters.
* Adjust thread counts (`-t`, `--threads`, `--cpu`, `-@`) based on your system's resources.
* This workflow provides a template. You may need to adapt steps, parameters, or tools based on specific research questions and data characteristics.
* Refer to the individual tool documentation (and potentially guides in `/Functionality` and `/Taxonomy` directories) for detailed usage and optimization.
