# Air Monitoring by Nanopore Sequencing

## Detailed Overview

This project utilizes metagenomic analysis of bioaerosols to monitor air quality, employing the latest advancements in nanopore sequencing technology. Our workflow is designed to extract, process, and analyze genetic material from air samples to understand microbial compositions. Key stages include:

### Processing of Reads
Initial processing involves quality control and preparation of raw sequencing reads, including error correction and adapter trimming, to ensure data accuracy.

### Base Calling with Guppy and Dorado
Base calling translates the electrical signals from nanopore sequencing into nucleotide sequences using Guppy and Dorado. This critical step converts raw data at a sequencing frequency of 4 or 5 kilohertz, enabling precise identification of nucleotides from the electrical signal generated by the DNA passing through nanopores.

### Assembly of Metagenomes
We assemble metagenomes from high-quality reads, reconstructing the genomic sequences of organisms present in bioaerosol samples. This process creates comprehensive genomic profiles for in-depth analysis.

### Binning
Sequences are sorted into bins based on genomic content, separating different organisms' genomes for targeted analysis.

### Quality Control
Quality control is integral at all stages, ensuring the integrity and reliability of data from sequencing to assembly and bining.

### Taxonomic Classification
Our analysis includes taxonomic classification at the read, contig, and metagenome-assembled genome (MAG) levels, providing insights into the diversity and structure of microbial communities.

### Protein Classification and Functional Annotation
We classify proteins and annotate genes and genomes to understand the biological functions and potential environmental impacts of identified organisms.

By integrating high-frequency nanopore sequencing with detailed bioinformatic analysis, this project aims to enhance our understanding of air microbial compositions.

### KEGG-Pathway Analysis

Prokka > KOBAS


# Nanopore Sequencing Data Analysis Pipeline

This repository contains a pipeline for the analysis of Oxford Nanopore sequencing data. The pipeline uses the following tools:

1. **Guppy Basecaller**: This software tool, provided by Oxford Nanopore Technologies (ONT), is used to convert the raw electrical signal data from nanopore sequencing into DNA sequences, a process known as basecalling.

2. **Porechop**: This tool is developed for Oxford Nanopore sequencing data. It finds and removes adapters from Oxford Nanopore reads.

3. **NanoFilt**: This is a simple tool to filter Oxford Nanopore sequencing data. It reads in a FASTQ file (or stdin), filters reads based on a minimum quality and/or a minimum length, and writes out the filtered reads to stdout.

4. **Flye**: This is a de novo assembler for single-molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies. It can be run in 'meta' mode for metagenomic projects.

5. **Minimap2**: This tool serves as an aligner for bio-sequences. In the context of nanopore sequencing, it aligns the raw nanopore reads to the assembled sequences produced by Flye. The `-ax map-ont` option specifies that the input data is Oxford Nanopore reads. The output is a BAM file sorted by samtools.

6. **Racon**: It is a tool that performs consensus sequencing. It takes as input the raw nanopore reads, and the alignment file produced by Minimap2, and the assembly produced by Flye, and then it generates a consensus sequence that represents an 'improved' version of the initial assembly. This step is also known as 'polishing' the assembly.

7. **MetaWRAP**: MetaWRAP is a comprehensive pipeline designed for the efficient analysis of metagenomic data. It integrates several bioinformatics tools and workflows to streamline the process of quality control, assembly, binning, and annotation of metagenomic samples. MetaWRAP simplifies the complex task of analyzing metagenomic data, making it more accessible while maintaining high standards of robustness and accuracy. The input for MetaWRAP typically consists of raw sequencing reads from microbial communities, which are then processed through its various modules to yield meaningful insights into the microbial composition and function.

9. **CheckM**: CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes. It uses lineage-specific marker sets to assess genome completeness and contamination. It provides robust estimates of genome completeness and contamination by using collocated sets of genes that are ubiquitous and single-copy within a lineage.

9. **Kraken 2**: This is a system for assigning taxonomic labels to short DNA sequences. It's generally used in metagenomics projects to identify the species present in a sample. The input to Kraken 2 are the filtered reads, the contigs, and the bins.

10. **Prodigal**: Prodigal (PROkaryotic DYnamic programming Gene-finding ALgorithm) is a microbial gene prediction program. It predicts protein-coding genes in the contigs. The input to Prodigal is the binned contigs from Vamb.

11. **DIAMOND**: Diamond (Double-Index Alignment of Next-Generation Sequencing Data) is a high-performance bioinformatics software for rapid alignment of DNA sequences to protein databases. It is particularly useful for large-scale genomic analysis, such as metagenomic studies. Diamond operates by translating nucleotide queries into protein sequences and aligning them against a protein reference database, offering speeds significantly faster than traditional BLAST. The typical input for Diamond includes short DNA sequences, often from metagenomic or other genomic datasets.

12. **MetaWRAP**: MetaWRAP is a comprehensive pipeline designed for the efficient analysis of metagenomic data. It integrates several bioinformatics tools and workflows to streamline the process of quality control, assembly, binning, and annotation of metagenomic samples. MetaWRAP simplifies the complex task of analyzing metagenomic data, making it more accessible while maintaining high standards of robustness and accuracy. The input for MetaWRAP typically consists of raw sequencing reads from microbial communities, which are then processed through its various modules to yield meaningful insights into the microbial composition and function.

13. **Prokka**: Prokka (Prokaryotic Genome Annotation System) is a comprehensive software tool designed for the rapid annotation of prokaryotic genomes. It is widely used in the microbiology research community for annotating bacterial, archaeal, and viral genomes. Prokka compiles a suite of software tools to predict coding sequences, rRNA genes, tRNA genes, and other features. It also assigns functional annotations to the predicted gene models using a combination of reference databases and heuristic methods. The output of Prokka is a richly annotated genome that can be used for further biological analysis and comparative genomics. Prokka's efficiency and ease of use make it an indispensable tool for microbial genomics studies, providing a comprehensive overview of functional genomics in a matter of hours.

14. **Seqkit**: Seqkit is a powerful and versatile toolkit for processing biological sequence data in FASTA/Q formats. It provides a collection of commands for various operations on sequencing reads, such as filtering, trimming, converting, and manipulating sequences. Seqkit is designed to be fast, efficient, and user-friendly, making it suitable for processing large-scale sequencing datasets. One of its key features is the ability to convert between FASTA and FASTQ formats, which is essential for compatibility with downstream analysis tools. Seqkit's modular design and command-line interface make it easily integrable into bioinformatics pipelines, enabling seamless data preprocessing and format conversion. Its performance and flexibility have made it a popular choice among bioinformaticians for handling sequencing data.

15. **AMRFinderPlus**: AMRFinderPlus is a comprehensive tool for the identification and annotation of antimicrobial resistance (AMR) genes in bacterial genome sequences. It utilizes a curated database of AMR gene sequences and a combination of sequence alignment and machine learning algorithms to predict the presence of AMR genes in assembled contigs or reads. AMRFinderPlus provides detailed information about the identified AMR genes, including the gene name, antibiotic class, and mechanism of resistance. By accurately identifying AMR genes, AMRFinderPlus enables researchers to assess the resistance profile of bacterial isolates and understand the potential for antibiotic resistance. Its ability to handle both contigs and reads as input makes it flexible for various stages of genome analysis. AMRFinderPlus has become an essential tool in the fight against antimicrobial resistance, aiding in the surveillance, monitoring, and understanding of AMR in bacterial populations.

16. **ABRicate**: ABRicate is a tool for mass screening of contigs or reads for antimicrobial resistance genes. It allows rapid identification of AMR genes in large datasets, such as those generated from whole-genome sequencing of bacterial isolates. ABRicate uses a database of known AMR gene sequences and performs sequence alignment to identify potential matches in the input contigs or reads. It provides a summary report of the identified AMR genes, along with their respective coverage and identity percentages. ABRicate's speed and efficiency make it suitable for high-throughput analysis of bacterial genomes for the presence of AMR genes. Its compatibility with both contigs and reads as input enables its use at different stages of the genome assembly and analysis pipeline. ABRicate has become a go-to tool for researchers studying antimicrobial resistance, facilitating the rapid screening and characterization of AMR in bacterial populations.

Note: To use AMRFinderPlus and ABRicate with sequencing reads as input, the reads typically need to be in FASTA format. Seqkit can be used to convert FASTQ reads to FASTA format, enabling their compatibility with these AMR analysis tools. This conversion step is crucial for ensuring the proper input format and allowing the accurate identification of AMR genes from the sequencing data.

