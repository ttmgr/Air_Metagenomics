# Nanopore Metagenomics Pipeline

A comprehensive bioinformatics pipeline for analyzing environmental microbiome data from air, water, and soil samples using Oxford Nanopore sequencing technology.

## 🔬 Overview

This pipeline provides an end-to-end solution for processing nanopore sequencing data from environmental samples, including:
- Basecalling and demultiplexing
- Quality control and filtering
- Taxonomic classification
- Functional annotation
- Antimicrobial resistance (AMR) and virulence factor detection
- Metagenomic assembly and binning

## 📋 Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Structure](#pipeline-structure)
- [Detailed Documentation](#detailed-documentation)
- [Citation](#citation)
- [License](#license)

## ✨ Features

- **Long-read optimized**: Designed specifically for Oxford Nanopore long-read data
- **Comprehensive analysis**: From raw signals to functional insights
- **AMR detection**: Multiple tools for antimicrobial resistance gene identification
- **Taxonomic profiling**: High-resolution microbial community composition
- **Functional annotation**: Detailed metabolic pathway analysis
- **Quality control**: Rigorous filtering and validation at each step
- **Modular design**: Run individual components or the full pipeline

## 🛠️ Requirements

### Core Tools
- **Dorado** (≥0.5.0) - Basecalling
- **Porechop** (≥0.2.4) - Adapter removal
- **NanoFilt** (≥2.8.0) - Read filtering
- **Flye** (≥2.9) - Metagenomic assembly
- **Kraken2** (≥2.1.2) - Taxonomic classification
- **DIAMOND** (≥2.0.15) - Protein alignment
- **ABRicate** (≥1.0.1) - AMR/virulence screening
- **AMRFinderPlus** (≥3.11) - AMR gene detection
- **Prokka** (≥1.14.6) - Genome annotation
- **eggNOG-mapper** (≥2.1.9) - Functional annotation
- **PlasmodFinder** - Plasmid identification

### Additional Tools
- Minimap2 (≥2.24)
- Samtools (≥1.16)
- Racon (≥1.5.0)
- MetaWRAP (≥1.3.2)
- Prodigal (≥2.6.3)
- CheckM (≥1.2.2)
- Bakta (≥1.8.1)
- NanoStat (≥1.6.0)

### System Requirements
- Linux operating system (Ubuntu 20.04+ or CentOS 7+)
- Minimum 64GB RAM (128GB recommended for large datasets)
- 500GB+ free disk space
- 16+ CPU cores recommended

## 📦 Installation

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/nanopore-metagenomics-pipeline.git
cd nanopore-metagenomics-pipeline
```

### 2. Set up conda environment
```bash
conda env create -f environment.yml
conda activate nanopore-metagenomics
```

### 3. Download databases
```bash
# Download and setup required databases
bash scripts/setup/download_databases.sh

# This will download:
# - Kraken2 database
# - DIAMOND nr database
# - AMR databases (CARD, ResFinder)
# - eggNOG database
# - Bakta database
```

### 4. Configure paths
```bash
cp config/config_template.yaml config/config.yaml
# Edit config.yaml with your specific paths
```

## 🚀 Quick Start

### Basic usage
```bash
# Run the complete pipeline
bash run_pipeline.sh -i input_folder/ -o output_folder/ -t 24

# Run with custom configuration
bash run_pipeline.sh -i input_folder/ -o output_folder/ -c config/custom_config.yaml
```

### Example with test data
```bash
# Download test dataset
bash scripts/download_test_data.sh

# Run pipeline on test data
bash run_pipeline.sh -i test_data/ -o test_output/ -t 8
```

## 📁 Pipeline Structure

```
nanopore-metagenomics-pipeline/
├── README.md
├── LICENSE
├── environment.yml
├── run_pipeline.sh
├── config/
│   ├── config_template.yaml
│   └── database_paths.yaml
├── docs/
│   ├── installation.md
│   ├── usage.md
│   ├── troubleshooting.md
│   └── tools/
│       ├── AMR_detection.md
│       ├── functional_annotation.md
│       └── taxonomic_classification.md
├── scripts/
│   ├── setup/
│   │   └── download_databases.sh
│   ├── preprocessing/
│   │   ├── dorado_basecalling.sh
│   │   ├── read_processing.sh
│   │   └── quality_control.sh
│   ├── assembly/
│   │   ├── assembly_and_polishing.sh
│   │   └── metagenomic_binning.sh
│   ├── annotation/
│   │   ├── prokka_annotation.sh
│   │   ├── bakta_annotation.sh
│   │   ├── prodigal_prediction.sh
│   │   └── eggnog_mapping.sh
│   ├── classification/
│   │   └── kraken2_classification.sh
│   ├── amr_detection/
│   │   ├── abricate_screening.sh
│   │   └── amrfinder_analysis.sh
│   └── analysis/
│       ├── extract_metrics.py
│       ├── generate_plots.py
│       └── create_report.py
├── workflows/
│   ├── preprocessing.nf
│   ├── assembly.nf
│   └── annotation.nf
└── test_data/
    └── README.md
```

## 📚 Detailed Documentation

### 1. [Preprocessing](docs/preprocessing.md)
- Basecalling with Dorado
- Demultiplexing
- Adapter removal with Porechop
- Quality filtering with NanoFilt

### 2. [Assembly](docs/assembly.md)
- Metagenomic assembly with Flye
- Read mapping with Minimap2
- Polishing with Racon
- Assembly statistics

### 3. [Taxonomic Classification](docs/tools/taxonomic_classification.md)
- Read-level classification with Kraken2
- Contig-level classification
- Abundance estimation
- Visualization

### 4. [Functional Annotation](docs/tools/functional_annotation.md)
- Gene prediction with Prodigal
- Annotation with Prokka/Bakta
- Functional mapping with eggNOG-mapper
- Pathway analysis

### 5. [AMR & Virulence Detection](docs/tools/AMR_detection.md)
- ABRicate screening
- AMRFinderPlus analysis
- Resistance gene quantification
- Virulence factor identification

### 6. [Metagenomic Binning](docs/binning.md)
- MetaWRAP binning
- Bin refinement
- Quality assessment with CheckM

## 🔧 Configuration

The pipeline behavior can be customized through the `config.yaml` file:

```yaml
# General settings
threads: 24
memory: 128G

# Quality thresholds
min_read_length: 100
min_read_quality: 9
min_contig_length: 1000

# Database paths
kraken2_db: /path/to/kraken2/db
diamond_db: /path/to/diamond/db
eggnog_db: /path/to/eggnog/db
bakta_db: /path/to/bakta/db

# Tool-specific parameters
flye:
  mode: meta
  min_overlap: 1000

kraken2:
  confidence: 0.1
  minimum_hit_groups: 2
```

## 📊 Output Structure

```
output_folder/
├── 01_preprocessing/
│   ├── basecalled/
│   ├── demultiplexed/
│   ├── trimmed/
│   └── filtered/
├── 02_assembly/
│   ├── flye/
│   ├── minimap2/
│   └── racon/
├── 03_classification/
│   ├── kraken2_reads/
│   └── kraken2_contigs/
├── 04_annotation/
│   ├── prokka/
│   ├── bakta/
│   └── eggnog/
├── 05_amr_detection/
│   ├── abricate/
│   └── amrfinder/
├── 06_binning/
│   └── metawrap/
├── 07_analysis/
│   ├── statistics/
│   ├── plots/
│   └── reports/
└── logs/
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

## 📄 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{nanopore_metagenomics_pipeline,
  title = {Nanopore Metagenomics Pipeline},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/nanopore-metagenomics-pipeline}
}
```

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Oxford Nanopore Technologies for sequencing technology
- All tool developers whose software is integrated in this pipeline
- The metagenomics community for continuous feedback and improvements

## 📞 Support

- 📧 Email: your.email@example.com
- 🐛 Issues: [GitHub Issues](https://github.com/yourusername/nanopore-metagenomics-pipeline/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/yourusername/nanopore-metagenomics-pipeline/discussions)
