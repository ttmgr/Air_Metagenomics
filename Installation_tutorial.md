# Installation Guide for Air Monitoring Nanopore Sequencing Pipeline Tools

This guide provides instructions for installing all the tools used in the Air Monitoring Nanopore Sequencing Pipeline. We'll be using Mamba to create separate environments for each tool, ensuring clean and isolated installations.

## Prerequisites

- Mamba (a faster alternative to Conda) should be installed on your system. If not, you can install it using:

```bash
conda install mamba -n base -c conda-forge
```

## General Installation Process

For each tool, we'll follow these steps:
1. Create a new Mamba environment
2. Activate the environment
3. Install the tool using Mamba
4. (Optional) Download required databases

## Tool Installation Instructions

### 1. Guppy Basecaller

Guppy is typically provided by Oxford Nanopore Technologies and may require a specific installation process. Please refer to the ONT community for the latest installation instructions.

### 2. Porechop

```bash
mamba create -n porechop -c bioconda porechop
mamba activate porechop
```

### 3. NanoFilt

```bash
mamba create -n nanofilt -c bioconda nanofilt
mamba activate nanofilt
```

### 4. Flye

```bash
mamba create -n flye -c bioconda flye
mamba activate flye
```

### 5. Minimap2

```bash
mamba create -n minimap2 -c bioconda minimap2
mamba activate minimap2
```

### 6. Racon

```bash
mamba create -n racon -c bioconda racon
mamba activate racon
```

### 7. MetaWRAP

```bash
mamba create -n metawrap -c ursky metawrap-mg
mamba activate metawrap
```

### 8. Prodigal

```bash
mamba create -n prodigal -c bioconda prodigal
mamba activate prodigal
```

### 9. EggNOG-mapper

```bash
mamba create -n eggnog-mapper -c bioconda eggnog-mapper
mamba activate eggnog-mapper
```

Download the EggNOG database:
```bash
download_eggnog_data.py
```

### 10. CheckM

```bash
mamba create -n checkm -c bioconda checkm-genome
mamba activate checkm
```

Download the CheckM database:
```bash
checkm data setRoot /path/to/checkm_data
checkm download_data -t /path/to/checkm_data
```

### 11. Kraken2

```bash
mamba create -n kraken2 -c bioconda kraken2
mamba activate kraken2
```

Download the Kraken2 database:
```bash
kraken2-build --standard --threads 4 --db /path/to/kraken2_db
```

### 12. DIAMOND

```bash
mamba create -n diamond -c bioconda diamond
mamba activate diamond
```

Download and format the NCBI nr database:
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
diamond makedb --in nr -d nr
```

### 13. ABRicate

```bash
mamba create -n abricate -c bioconda abricate
mamba activate abricate
```

Update the ABRicate databases:
```bash
abricate-get_db --setupdb
```

### 14. NCBI-AMRFinderPlus

```bash
mamba create -n amrfinder -c bioconda ncbi-amrfinderplus
mamba activate amrfinder
```

Update the AMRFinder database:
```bash
amrfinder_update --force_update
```

### 15. PfamScan

```bash
mamba create -n pfamscan -c bioconda pfam_scan
mamba activate pfamscan
```

Download the Pfam database:
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

### 16. Prokka

```bash
mamba create -n prokka -c bioconda prokka
mamba activate prokka
```

### 17. Bakta

```bash
mamba create -n bakta -c conda-forge -c bioconda bakta
mamba activate bakta
```

Download the Bakta database:
```bash
bakta_db download --output /path/to/bakta_db
```

### 18. NanoStat

```bash
mamba create -n nanostat -c bioconda nanostat
mamba activate nanostat
```

### 19. Assembly-Stats

```bash
mamba create -n assembly-stats -c bioconda assembly-stats
mamba activate assembly-stats
```

## Usage

To use a specific tool, activate its environment before running:

```bash
mamba activate tool_name
# Run the tool
mamba deactivate
```

Replace `tool_name` with the name of the environment you created for that tool.

## Note

- Always ensure you're using the correct environment for each tool to avoid conflicts.
- Regularly update your tools and databases to benefit from the latest improvements and data.
- Some databases are large and may take significant time and storage space to download and process.
- Adjust the paths in the database download commands according to your system's directory structure.
