# Wetland Metagenomics & Viromics by Nanopore Sequencing (PRJEBXXXXX)

## Project Overview

This project utilizes metagenomic (DNA) and viromic (RNA) analysis of environmental water and air samples to characterize microbial communities, antimicrobial resistance (AMR) genes, and Avian Influenza Viruses (AIV), employing Oxford Nanopore Technologies (ONT) sequencing[cite: 36]. This repository contains workflows designed to process raw sequencing data for these analyses.

## ENA Files: [https://www.ebi.ac.uk/ena/browser/view/PRJEBXXXXX](https://www.ebi.ac.uk/ena/browser/view/PRJEBXXXXX) ---

**❗ Important First Step ❗**

Before using the general pipelines described below, you **must** consult the specific data guides. These guides contain critical information about:

* **Data Access:** Where to find the raw (POD5) and/or processed (FASTQ) data.
* **Sample Mapping & Barcoding:** Which barcode corresponds to which sample for different sequencing runs (DNA shotgun, AIV). Details on library prep (e.g., RBK114-24 for DNA[cite: 78], SQK-RBK114.24 for AIV [cite: 108]) and pooling strategies[cite: 79, 80, 81, 82].
* **Basecalling/Demultiplexing:** Specific Dorado commands and configurations used for the initial conversion of raw POD5 data to FASTQ files, including required demultiplexing based on barcodes[cite: 83, 109].

**Find the essential guides here:**

* **DNA Shotgun Data Guide:** [`dna_shotgun_data_guide.md`](./dna_shotgun_data_guide.md)
    * Describes samples for DNA metagenomic analysis (active water, passive water, air)[cite: 117].
    * Uses **POD5** data format and **Dorado** (v5.0.0, model dna_r10.4.1_e8.2_400bps_sup@v5.0.0) for basecalling/demultiplexing[cite: 83].
* **AIV (RNA) Data Guide:** [`aiv_rna_data_guide.md`](./aiv_rna_data_guide.md)
    * Describes samples positive for AIV and processed for RNA sequencing[cite: 104].
    * Uses **POD5** data format and **Dorado** for basecalling/demultiplexing (primers and adapters removed by Dorado)[cite: 109].

**The pipelines described below assume you have already completed the necessary steps from the relevant guide and have demultiplexed FASTQ files ready for analysis.**

---

## Analysis Pipeline Overview

This repository outlines two main analysis pipelines:

1.  **DNA Shotgun Metagenomics:** Processes FASTQ files from environmental DNA samples.
    * Read Processing: Adapter trimming and quality/length filtering[cite: 84, 85].
    * Taxonomic Classification: Assigning taxonomy to reads[cite: 88].
    * Metagenome Assembly: Assembling reads into contigs using two different assemblers[cite: 91, 92].
    * Assembly Polishing: Improving assembly accuracy[cite: 91, 92].
    * AMR Gene Detection: Identifying antimicrobial resistance genes from reads and contigs[cite: 95, 96, 98].
    * Taxonomic Origin of AMR: Determining the host of AMR genes on contigs[cite: 99, 100, 101, 102].
2.  **AIV (RNA) Analysis:** Processes FASTQ files from samples amplified for AIV[cite: 104].
    * Read Processing: Quality and length filtering[cite: 110].
    * Alignment: Aligning reads to AIV reference genomes[cite: 111, 112].
    * Consensus Sequence Generation: Creating a consensus AIV genome for each sample[cite: 116].

## Tools Used

This project integrates the following key bioinformatics tools:

* **Dorado:** Basecaller for ONT data (v5.0.0 for DNA[cite: 83], specified version for AIV [cite: 109]).
* **Porechop:** Adapter and barcode trimming (v0.2.4)[cite: 84].
* **NanoFilt:** Quality and length filtering for ONT reads (v2.8.0)[cite: 85].
* **Filtlong:** Quality and length filtering for AIV reads[cite: 110].
* **Kraken2:** K-mer based taxonomic classification (v2.1.2)[cite: 88, 94, 101].
* **metaFlye:** Long-read assembler for metagenomes (v2.9.6)[cite: 91].
* **nanoMDBG:** Long-read assembler for metagenomes (v1.1)[cite: 92].
* **Minimap2:** Long-read alignment (v2.28 for polishing DNA assemblies[cite: 91], v2.26 for AIV alignment [cite: 111]).
* **Racon:** Consensus correction/polishing for assemblies (v1.5)[cite: 91].
* **Medaka:** Consensus correction/polishing for assemblies (v2.0.1)[cite: 92].
* **AMRFinderPlus:** Detection of AMR genes (v3.12.8)[cite: 95].
* **DIAMOND:** Protein sequence alignment for taxonomic assignment of AMR-carrying contigs[cite: 101].
* **Seqkit:** Toolkit for FASTA/Q sequence manipulation (v2.10.0)[cite: 86, 96, 97].
* **SAMtools:** Utilities for SAM/BAM alignment files (v1.17)[cite: 113].
* **BCFtools:** Utilities for variant calling and consensus generation (v1.17)[cite: 116].
* **(Python Libraries for PCoA):** scikit-bio v0.6.3, Matplotlib v3.10.0, Pandas v2.2.3, NumPy v1.26.4[cite: 90].

## Repository Structure

* `dna_shotgun_data_guide.md`: Essential guide for accessing and preprocessing DNA shotgun data.
* `aiv_rna_data_guide.md`: Essential guide for accessing and preprocessing AIV (RNA) data.
* `Installation_tutorial.md`: Step-by-step guide for installing all required tools and databases.
* `/scripts` (example directory): May contain helper scripts or detailed command sequences.

## Installation

For comprehensive installation instructions for all pipeline tools and required databases, please refer to the [`Installation_tutorial.md`](./Installation_tutorial.md) file.

## Usage Workflow

This section outlines the commands for processing samples. You will typically need to run these steps for each sample, potentially using loops or a workflow manager.

### DNA Shotgun Metagenomics Workflow

**Input:** Demultiplexed FASTQ files per sample (e.g., `sample_X.fastq`), obtained via steps in `dna_shotgun_data_guide.md`.

1.  **Read Processing (Adapters, Barcodes, Quality/Length):**
    * Dorado (during basecalling) demultiplexes based on barcodes[cite: 83].
    * Porechop removes sequencing adapters and barcodes[cite: 84].
    * NanoFilt filters reads by length (minimum 100 bp)[cite: 85].
    ```bash
    # Input: sample_X_raw.fastq (example name post-demux)
    porechop -i sample_X_raw.fastq -o sample_X.trimmed.fastq # [cite: 84]
    NanoFilt --length 100 sample_X.trimmed.fastq > sample_X.filtered.fastq # [cite: 85]
    # Output: sample_X.filtered.fastq
    ```

2.  **Taxonomic Classification (Reads):**
    * Classify filtered reads using Kraken2[cite: 88]. (Downsample for PCoA as per methods [cite: 86]).
    ```bash
    # Input: sample_X.filtered.fastq
    # For PCoA: First downsample to 14,000 reads (if sample has more) [cite: 86]
    # seqkit sample -n 14000 sample_X.filtered.fastq -o sample_X.14k.fastq # [cite: 86]
    # KRAKEN_DB_PATH=/path/to/nt_core_database_May2025
    kraken2 --db $KRAKEN_DB_PATH --threads <N> --output sample_X.kraken_output.txt --report sample_X.kraken_report.txt sample_X.filtered.fastq # Or sample_X.14k.fastq for PCoA set [cite: 88]
    # Output: sample_X.kraken_output.txt, sample_X.kraken_report.txt
    ```

3.  **Metagenome Assembly & Polishing:**
    * Option A: metaFlye + Racon polishing [cite: 91]
        ```bash
        # Input: sample_X.filtered.fastq
        flye --nano-raw sample_X.filtered.fastq --out-dir sample_X_metaflye_assembly --meta --threads <N> # [cite: 91]
        # Polishing with Racon (example: 3 rounds) [cite: 91]
        ASSEMBLY_FA=sample_X_metaflye_assembly/assembly.fasta
        POLISHED_FA_RAC=$ASSEMBLY_FA
        for i in {1..3}; do
          # PDF uses Minimap2 then Racon [cite: 91]
          minimap2 -ax map-ont $POLISHED_FA_RAC sample_X.filtered.fastq | samtools sort -@ <N> -o sample_X.aligned_for_racon.bam #
          samtools index sample_X.aligned_for_racon.bam #
          racon -t <N> sample_X.filtered.fastq sample_X.aligned_for_racon.bam $POLISHED_FA_RAC > sample_X.racon_round${i}.fasta # [cite: 91]
          POLISHED_FA_RAC=sample_X.racon_round${i}.fasta
        done
        # Final metaFlye polished assembly: $POLISHED_FA_RAC (e.g., sample_X.racon_round3.fasta)
        # Further polish with Medaka v2.0.1 (Applied to both assemblers in paper) [cite: 92]
        # medaka_consensus -i sample_X.filtered.fastq -d $POLISHED_FA_RAC -o sample_X_metaflye_medaka_polished -t <N> -m <model_appropriate_for_r10.4.1> # [cite: 92]
        # Output: sample_X_metaflye_medaka_polished/consensus.fasta
        ```
    * Option B: nanoMDBG + Medaka polishing [cite: 92]
        ```bash
        # Input: sample_X.filtered.fastq
        # nanoMDBG v1.1 command (refer to tool's documentation for specific parameters) [cite: 92]
        # nanomdbg --reads sample_X.filtered.fastq --out sample_X_nanomdbg_assembly # (placeholder command)
        # ASSEMBLY_FA_NMDBG=sample_X_nanomdbg_assembly/contigs.fasta (example output)
        # Polish with Medaka v2.0.1 [cite: 92]
        # medaka_consensus -i sample_X.filtered.fastq -d $ASSEMBLY_FA_NMDBG -o sample_X_nanomdbg_medaka_polished -t <N> -m <model_appropriate_for_r10.4.1> # [cite: 92]
        # Output: sample_X_nanomdbg_medaka_polished/consensus.fasta (Chosen for AMR analysis in paper [cite: 139])
        ```

4.  **AMR Gene Detection:**
    * On downsampled reads and nanoMDBG contigs[cite: 96, 98].
    ```bash
    # Input: sample_X.filtered.fastq, sample_X_nanomdbg_medaka_polished/consensus.fasta
    # On reads (downsample first based on sample type: 87k for Passive Water, 93k for Active Water, 14k for Air) [cite: 96]
    # seqkit sample -n <read_count_threshold> sample_X.filtered.fastq > sample_X.downsampled.fastq # [cite: 96]
    # seqkit fq2fa sample_X.downsampled.fastq -o sample_X.downsampled.fasta # [cite: 97]
    AMRFINDER_DB_PATH=/path/to/amrfinder_db
    amrfinder -n sample_X.downsampled.fasta -o sample_X.amr_reads.txt --threads <N> --plus --database $AMRFINDER_DB_PATH # [cite: 95, 98]

    # On nanoMDBG contigs
    CONTIGS_FA=sample_X_nanomdbg_medaka_polished/consensus.fasta
    amrfinder -n $CONTIGS_FA -o sample_X.amr_contigs.txt --threads <N> --plus --database $AMRFINDER_DB_PATH # [cite: 95, 98]
    # Output: sample_X.amr_reads.txt, sample_X.amr_contigs.txt
    ```

5.  **Taxonomic Origin of AMR on Contigs:**
    * Use DIAMOND and Kraken2 on contigs with AMR genes[cite: 99, 100, 101]. An AMR gene was assigned to a specific species only if both DIAMOND and Kraken2 showed the same species-level classification[cite: 102].
    ```bash
    # Input: Contigs identified by AMRFinderPlus as having AMR genes (e.g., from sample_X.amr_contigs.txt) [cite: 100]
    # For each contig (CONTIG_WITH_AMR.fa):
    # DIAMOND_DB_PATH=/path/to/nr_database_May2025
    # KRAKEN_DB_PATH=/path/to/nt_core_database_May2025
    # diamond blastx -d $DIAMOND_DB_PATH -q CONTIG_WITH_AMR.fa -o contig.diamond.txt --threads <N> --outfmt 6 # [cite: 101]
    # kraken2 --db $KRAKEN_DB_PATH --threads <N> --output contig.kraken.txt CONTIG_WITH_AMR.fa # [cite: 101]
    # Assign species if both DIAMOND and Kraken2 agree [cite: 102]
    ```

### AIV (RNA) Analysis Workflow

**Input:** Demultiplexed FASTQ files per sample (e.g., `aiv_sample_Y.fastq`), obtained via steps in `aiv_rna_data_guide.md`. Dorado removed primers/adapters[cite: 109].

1.  **Read Processing (Quality/Length Filtering):**
    * Filter reads using Filtlong (minimum Phred score >8, min length >150 bp)[cite: 110].
    ```bash
    # Input: aiv_sample_Y.fastq
    filtlong --min_mean_q 9 --min_length 150 aiv_sample_Y.fastq > aiv_sample_Y.filtered.fastq # PDF states ">8", commonly Q9 [cite: 110]
    # Output: aiv_sample_Y.filtered.fastq
    ```

2.  **Alignment to Reference Genomes:**
    * Align filtered reads to AIV reference database (all AIV nucleotide sequences from Europe as of 04/03/2023) using Minimap2 (v2.26) with -ax map-ont setting[cite: 111, 112].
    ```bash
    # Input: aiv_sample_Y.filtered.fastq
    AIV_REF_DB_PATH=/path/to/aiv_segment_references_europe_04032023.fasta
    minimap2 -ax map-ont $AIV_REF_DB_PATH aiv_sample_Y.filtered.fastq > aiv_sample_Y.aligned.sam # [cite: 111]
    # Output: aiv_sample_Y.aligned.sam
    ```

3.  **Consensus Sequence Generation:**
    * Convert SAM to BAM, sort, and index using SAMtools (v1.17)[cite: 113].
    * Select best reference for each segment using `samtools idxstats` to find the reference to which most reads mapped across every segment[cite: 113, 114].
    * Map all reads to the best reference for each of the eight segments[cite: 115].
    * Generate consensus using BCFtools (v1.17)[cite: 116].
    ```bash
    # Input: aiv_sample_Y.aligned.sam (or aiv_sample_Y.filtered.fastq if re-mapping to best refs)
    samtools view -bS aiv_sample_Y.aligned.sam | samtools sort -o aiv_sample_Y.sorted.bam # [cite: 113]
    samtools index aiv_sample_Y.sorted.bam # [cite: 113]

    # For each of the 8 AIV segments (example for one chosen best_ref_segment.fa):
    # minimap2 -ax map-ont best_ref_chosen_segment.fa aiv_sample_Y.filtered.fastq | samtools sort -@ <N> -o aiv_sample_Y.segment_aligned.bam # [cite: 115]
    # samtools index aiv_sample_Y.segment_aligned.bam #
    # bcftools mpileup -Ou -f best_ref_chosen_segment.fa aiv_sample_Y.segment_aligned.bam | bcftools call -mv -Oz -o aiv_sample_Y.segment.vcf.gz # [cite: 116]
    # bcftools index aiv_sample_Y.segment.vcf.gz #
    # bcftools consensus -f best_ref_chosen_segment.fa aiv_sample_Y.segment.vcf.gz -o aiv_sample_Y.segment_consensus.fasta # [cite: 116]
    # Repeat for all 8 segments
    # Output: Consensus fasta files for each AIV segment
    ```

**Notes:**

* Replace placeholders like `/path/to/`, `sample_X`, `aiv_sample_Y`, `<N>`, `<read_count_threshold>`, database paths, and Medaka models with your actual paths, names, desired parameters, and appropriate models for R10.4.1 data.
* Adjust thread counts (`-t`, `--threads`, `-@`) based on your system's resources.
* This workflow provides a template based on the provided PDF. Refer to individual tool documentation for detailed usage and optimization.
* The PDF specifies "Minimap2 v2.28 [ref] and three rounds of Racon v1.5 [ref]" for metaFlye polishing[cite: 91]. Racon typically takes alignments in SAM/BAM format for correction.
* The PDF specifies that assemblies from metaFlye and nanoMDBG were polished with Medaka v2.0.1[cite: 92].
