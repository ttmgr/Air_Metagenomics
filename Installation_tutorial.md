# Installation Guide for the Nanopore Metagenomics Pipeline

This guide provides streamlined instructions for setting up the necessary software environment and databases to run the entire Air Monitoring Nanopore Sequencing Pipeline.

## 🔬 Overview of the Installation Process

The setup process consists of three main steps:

1.  **Create the Conda Environment**: Use the provided `environment.yaml` file to create a single, unified Conda environment with all the required pipeline tools. This is the recommended and most straightforward method.
2.  **Install ONT Basecaller (Manual)**: Manually install the appropriate Oxford Nanopore basecaller (Guppy or Dorado), as these are not available through Conda.
3.  **Download Databases**: Run the provided helper script to download the large databases required for taxonomic classification, functional annotation, and AMR gene detection.

---

## ✅ Prerequisites

Before you begin, ensure you have **Mamba** installed, which is a much faster alternative to Conda for creating environments. If you don't have it, you can install it into your base Conda environment:

```bash
conda install mamba -n base -c conda-forge
```

---

## Step 1: Create the Conda Environment

This is the primary installation step. The `environment.yaml` file will automatically install all the necessary pipeline tools into a single, isolated environment named `nanopore-metagenomics`.

```bash
# Navigate to the repository's root directory
# Create the environment using the provided file
mamba env create -f env/environment.yaml
```

This process will take some time as it downloads and installs numerous bioinformatics packages. Once complete, you can activate the environment at any time using:

```bash
mamba activate nanopore-metagenomics
```

---

## Step 2: Install ONT Basecallers (Manual)

The Oxford Nanopore basecallers, **Guppy** (for FAST5 data) and **Dorado** (for POD5 data), must be installed manually. They are not available on Conda.

Please refer to the **ONT Community** for instructions on downloading and installing the correct version for your system's architecture (e.g., Linux/Windows, CPU/GPU).

---

## Step 3: Download Databases

The pipeline requires several large databases. A helper script, `download_databases.sh`, is provided to automate this process.

1.  **❗️ IMPORTANT ❗️:** First, open the script **`bash_scripts/download_databases.sh`** and edit the `DB_BASE_DIR` variable to the full path where you want to store the databases. This location requires **~400-500 GB** of free space.

2.  Run the script. Make sure you have activated the `nanopore-metagenomics` environment first, as it contains some of the required helper tools.
    ```bash
    mamba activate nanopore-metagenomics
    bash bash_scripts/download_databases.sh
    ```

3.  The script will download and set up databases for:
    * Kraken2
    * AMRFinderPlus
    * Bakta
    * eggNOG-mapper
    * ABRicate

4.  After the script finishes, it will print the final paths for each database. Copy these paths and update the configuration section at the top of **`bash_scripts/run_pipeline.sh`** to ensure the pipeline can find them.

---

## 🧰 Key Tools Included in the Environment

Creating the environment with the `environment.yaml` file provides all the necessary tools for the automated pipeline, including:

* **Read Processing**: `Porechop`, `NanoFilt`
* **Quality Control**: `NanoStat`, `Assembly-stats`
* **Assembly & Polishing**: `Flye`, `Minimap2`, `Racon`
* **Binning**: `MetaWRAP`, `MetaBAT2`, `MaxBin2`, `CONCOCT`
* **Annotation**: `Kraken2`, `Prokka`, `Bakta`, `Prodigal`, `eggNOG-mapper`
* **AMR Detection**: `ABRicate`, `NCBI-AMRFinderPlus`
* **Utilities**: `Seqkit`, `Samtools`

With the environment activated and databases configured, you are now ready to run the main pipeline.
