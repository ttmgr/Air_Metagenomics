# Methods for Functional Annotation of Metagenomic Data Derived from Nanopore Sequencing

Functional annotation of metagenomic data involves assigning biological functions to sequences obtained from environmental samples. Nanopore sequencing provides long-read data, which is advantageous for resolving complex genomic regions. This document outlines tools commonly used for the functional annotation of metagenomic datasets derived from nanopore sequencing.

## Table of Contents

- [DIAMOND](#diamond)
- [MEGAN](#megan)
- [eggNOG-mapper](#eggnog-mapper)
- [InterProScan](#interproscan)
- [Prokka](#prokka)
- [MetaErg](#metaerg)
- [DRAM](#dram)
- [MetaGeneMark](#metagenemark)
- [HUMAnN](#humann)
- [GhostKOALA and BlastKOALA](#ghostkoala-and-blastkoala)
- [SqueezeMeta](#squeezemetameta)
- [Pfam Scan](#pfam-scan)
- [COG Annotation](#cog-annotation)
- [PANNZER2](#pannzer2)
- [DeepARG](#deeparg)
- [antiSMASH](#antismash)
- [Kaiju](#kaiju)
- [MG-RAST](#mg-rast)
- [Databases Used](#databases-used)
- [General Workflow](#general-workflow)
- [Considerations for Nanopore Data](#considerations-for-nanopore-data)

## DIAMOND

DIAMOND is a fast sequence aligner for protein and translated DNA searches.

| Aspect | Description |
|--------|-------------|
| How it Works | Uses double-indexed search algorithms for fast alignment against large protein databases |
| Input | - Nucleotide sequences (translated in all six frames) or protein sequences in FASTA format<br>- Protein databases (e.g., NCBI NR) |
| Output | Alignment files in BLAST tabular format, including match statistics like e-values and bit scores |

## MEGAN

MEGAN (MEtaGenome ANalyzer) is a tool for analyzing metagenomic data, allowing taxonomic and functional interpretation.

| Aspect | Description |
|--------|-------------|
| How it Works | Parses alignment files to assign sequences to taxonomic and functional categories |
| Input | - BLAST or DIAMOND alignment files<br>- Reference databases (e.g., NCBI taxonomy, SEED, KEGG) |
| Output | - Interactive charts and graphs showing taxonomic and functional distributions<br>- Exportable summary tables |

## eggNOG-mapper

eggNOG-mapper provides fast functional annotation of sequences using the eggNOG database of orthologous groups.

| Aspect | Description |
|--------|-------------|
| How it Works | Maps query sequences to precomputed orthologous groups and transfers functional information based on homology |
| Input | Protein or nucleotide sequences in FASTA format |
| Output | - Functional annotations including Gene Ontology (GO) terms, KEGG pathways, and COG categories<br>- Tab-delimited annotation files |

## InterProScan

InterProScan integrates multiple protein signature recognition methods into one resource.

| Aspect | Description |
|--------|-------------|
| How it Works | Scans protein sequences against InterPro's member databases to identify domains, motifs, and sites |
| Input | Protein sequences in FASTA format |
| Output | - Annotated sequences with identified protein families and domains<br>- Various output formats including XML and TSV |

## Prokka

Prokka is a rapid prokaryotic genome annotation tool.

| Aspect | Description |
|--------|-------------|
| How it Works | Predicts coding sequences and RNA features, then assigns functions based on similarity to known proteins |
| Input | Assembled genomes or contigs in FASTA format |
| Output | - Annotated genome files in GFF, GenBank, and other formats<br>- Functional annotation tables |

## MetaErg

MetaErg is an automated pipeline for functional annotation of metagenomic assemblies.

| Aspect | Description |
|--------|-------------|
| How it Works | Combines gene prediction with functional annotation using multiple databases |
| Input | Assembled metagenomic contigs in FASTA format |
| Output | - Annotated genes with functional assignments<br>- Summary reports and visualization files |

## DRAM

DRAM (Distilled and Refined Annotation of Metabolism) annotates microbial genomes and metagenomes with metabolic functions.

| Aspect | Description |
|--------|-------------|
| How it Works | Uses a combination of databases and heuristics to assign metabolic functions, focusing on carbohydrate metabolism and energy conservation |
| Input | Assembled genomes or metagenomes in FASTA format |
| Output | - Annotated metabolic functions<br>- Excel-compatible summary tables |

## MetaGeneMark

MetaGeneMark predicts protein-coding genes in metagenomic sequences.

| Aspect | Description |
|--------|-------------|
| How it Works | Applies species-unspecific gene-finding algorithms suitable for metagenomic data |
| Input | Nucleotide sequences in FASTA format |
| Output | - Predicted gene coordinates<br>- Translated protein sequences |

## HUMAnN

HUMAnN (The HMP Unified Metabolic Analysis Network) profiles the presence, absence, and abundance of microbial pathways.

| Aspect | Description |
|--------|-------------|
| How it Works | Maps metagenomic reads to a reference database to quantify gene and pathway abundances |
| Input | Quality-controlled nucleotide sequences |
| Output | - Tables of gene family and pathway abundances<br>- Stratified and unstratified functional profiles |

## GhostKOALA and BlastKOALA

Online tools for KEGG Orthology (KO) assignment.

| Aspect | Description |
|--------|-------------|
| How it Works | Performs homology searches against the KEGG database to assign KO numbers and reconstruct metabolic pathways |
| Input | Protein sequences in FASTA format |
| Output | - KO assignments<br>- Pathway maps highlighting detected functions |

## SqueezeMeta

An automatic pipeline for metagenomic analysis, integrating assembly, binning, and annotation.

| Aspect | Description |
|--------|-------------|
| How it Works | Combines various software tools to process metagenomic data end-to-end |
| Input | Raw sequencing reads or assembled contigs |
| Output | - Annotated genes and genomes<br>- Functional and taxonomic profiles |

## Pfam Scan

Identifies Pfam protein domains within sequences.

| Aspect | Description |
|--------|-------------|
| How it Works | Uses HMMER to search protein sequences against Pfam HMM profiles |
| Input | Protein sequences in FASTA format |
| Output | Detected Pfam domains with statistical significance scores |

## COG Annotation

Assigns genes to Clusters of Orthologous Groups (COGs).

| Aspect | Description |
|--------|-------------|
| How it Works | Performs sequence similarity searches against the COG database to categorize proteins |
| Input | Protein sequences |
| Output | COG assignments and associated functional categories |

## PANNZER2

Automated functional annotation for prokaryotic and eukaryotic proteins.

| Aspect | Description |
|--------|-------------|
| How it Works | Uses sequence similarity and machine learning to assign functions |
| Input | Protein sequences in FASTA format |
| Output | - Functional descriptions<br>- GO terms and Enzyme Commission (EC) numbers |

## DeepARG

Predicts antibiotic resistance genes (ARGs) using deep learning.

| Aspect | Description |
|--------|-------------|
| How it Works | Utilizes neural networks trained on known ARGs to classify new sequences |
| Input | DNA or protein sequences in FASTA format |
| Output | ARG predictions with confidence scores |

## antiSMASH

Identifies secondary metabolite biosynthesis gene clusters.

| Aspect | Description |
|--------|-------------|
| How it Works | Detects and analyzes biosynthetic pathways in genomic data |
| Input | Assembled genomes or contigs |
| Output | - Annotated gene clusters<br>- Visualization of biosynthetic pathways |

## Kaiju

A tool for taxonomic classification and functional annotation.

| Aspect | Description |
|--------|-------------|
| How it Works | Translates DNA reads into protein sequences and searches them against a protein database |
| Input | Nucleotide sequences |
| Output | - Taxonomic assignments<br>- Functional annotations based on protein matches |

## MG-RAST

An online metagenomics service providing automated analysis.

| Aspect | Description |
|--------|-------------|
| How it Works | Processes data through quality control, gene prediction, and functional annotation pipelines |
| Input | Raw sequencing reads |
| Output | - Taxonomic and functional profiles<br>- Comparative analysis tools |

## Databases Used

- [x] NCBI NR: Non-redundant protein database
- [x] KEGG: Pathway and enzyme information
- [x] Pfam: Protein families and domains
- [x] eggNOG: Orthologous groups and functional annotation
- [x] COG: Clusters of Orthologous Groups for function prediction
- [x] UniProt: Protein sequence and functional information

---

For more information on each tool, please refer to their respective documentation and publications.
