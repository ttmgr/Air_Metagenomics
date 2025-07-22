# Guide to Taxonomic, Functional, and AMR Annotation Tools

This document provides a detailed explanation of the primary tools used for annotation in this pipeline, combining information on taxonomic classification, functional annotation, and AMR/virulence gene detection.

## 1. Taxonomic Classification Tools

### Kraken2
- **What it is**: An ultrafast, k-mer-based taxonomic classification tool.
- **How it Works**: It maps k-mers (short, fixed-length DNA sequences) from each read to the lowest common ancestor (LCA) of all genomes that contain that k-mer in a pre-built database. This makes it extremely fast.
- **Input**: Raw or assembled nucleotide sequences (FASTA/FASTQ).
- **Output**: A text file with a classification for each read, and a summary report file with the abundance of each detected taxon.
- **Best for**: Rapidly assigning taxonomy to millions of reads, providing a quick overview of community composition.

### Kaiju
- **What it is**: A fast, protein-level taxonomic classifier.
- **How it Works**: It translates nucleotide reads into all six reading frames and searches for maximum exact matches (MEMs) against a protein reference database (like NCBI nr). This approach is more tolerant of sequencing errors than k-mer methods.
- **Input**: Raw nucleotide sequences (FASTA/FASTQ).
- **Output**: A classification file similar to Kraken2's output.
- **Best for**: Analyzing data with higher error rates (like uncorrected nanopore reads) or for identifying novel or divergent organisms that might be missed by nucleotide-based methods.

## 2. Functional Annotation Tools

### Prodigal
- **What it is**: A highly accurate gene prediction tool for prokaryotic genomes.
- **How it Works**: It uses a dynamic programming algorithm to identify protein-coding genes. Its "meta" mode is specifically optimized for short and incomplete fragments of DNA from metagenomic samples.
- **Input**: Assembled contigs (FASTA).
- **Output**: Gene coordinates (GFF), nucleotide sequences of genes (FASTA), and translated protein sequences (FASTA).
- **Best for**: The essential first step of finding genes in your assembly before you can annotate their function.

### Prokka & Bakta
- **What it is**: Rapid, all-in-one prokaryotic genome annotation pipelines.
- **How it Works**: These tools automate the process of gene finding (using Prodigal), rRNA/tRNA detection, and assigning functions by searching against multiple databases of known proteins (e.g., UniProt) and protein domains.
- **Input**: Assembled contigs (FASTA).
- **Output**: Comprehensive annotation files in standard formats like GFF, GenBank, and summary tables.
- **Best for**: Quickly generating high-quality, standardized annotations for assembled genomes or high-completeness MAGs.

### eggNOG-mapper
- **What it is**: A tool for fast functional annotation based on orthology.
- **How it Works**: It maps your predicted protein sequences to the eggNOG database of orthologous groups. By identifying the group a protein belongs to, it can transfer a wealth of pre-computed functional information, including Gene Ontology (GO) terms, KEGG pathways, and COG categories.
- **Input**: Predicted protein sequences (FASTA).
- **Output**: A detailed tab-delimited file with functional annotations for each protein.
- **Best for**: Getting broad functional profiling and metabolic pathway information for your entire metagenome.

## 3. AMR & Virulence Gene Detection Tools

### ABRicate
- **What it is**: A tool for mass screening of contigs for AMR and virulence genes.
- **How it Works**: It takes your assembled contigs and performs BLASTn searches against multiple curated databases simultaneously, including CARD (AMR), ResFinder (AMR), and VFDB (virulence factors). It then generates a convenient summary report.
- **Input**: Assembled contigs (FASTA).
- **Output**: A single, easy-to-read summary file tabulating all hits across the different databases.
- **Best for**: A rapid and broad initial screen for a wide variety of known AMR and virulence genes.

### AMRFinderPlus
- **What it is**: A tool developed by NCBI for identifying acquired AMR genes.
- **How it Works**: It uses a highly curated database of AMR proteins and hidden Markov models (HMMs), making it very precise. It is the standard tool used by NCBI for annotating resistance genes in submitted genomes.
- **Input**: Can take either nucleotide contigs or predicted protein sequences (FASTA).
- **Output**: A detailed tab-delimited report on found AMR genes, including information on resistance mechanisms.
- **Best for**: A highly accurate, targeted search for AMR genes that aligns with NCBI standards.
