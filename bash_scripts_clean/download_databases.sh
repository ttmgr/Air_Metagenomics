#!/bin/bash
# Script to download and set up all required databases for the pipeline.
#
# IMPORTANT:
# 1. Edit the DB_BASE_DIR to the location where you want to store ~400GB of databases.
# 2. Ensure you have the required tools (wget, tar, kraken2, diamond, etc.) installed.
#    It's best to run this after setting up the conda environment from environment.yaml.

set -euo pipefail

# --- User Configuration ---
# The main directory to store all databases.
# IMPORTANT: This requires a large amount of disk space (~400-500 GB).
DB_BASE_DIR="/path/to/your/databases"

# Number of threads for database building
THREADS=16

# --- Script ---
echo "--- Database Download and Setup Script ---"
echo "Databases will be installed in: $DB_BASE_DIR"
echo "This will download several hundred gigabytes of data and may take many hours."
read -p "Are you sure you want to continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Setup cancelled."
    exit 1
fi

mkdir -p "$DB_BASE_DIR"

# 1. Kraken2 Standard Database (~100 GB)
KRAKEN2_DB_PATH="${DB_BASE_DIR}/kraken2_standard"
if [ ! -d "$KRAKEN2_DB_PATH" ]; then
    echo "--- Setting up Kraken2 Database ---"
    mkdir -p "$KRAKEN2_DB_PATH"
    kraken2-build --standard --threads "$THREADS" --db "$KRAKEN2_DB_PATH"
else
    echo "--- Kraken2 Database already exists. Skipping. ---"
fi

# 2. AMRFinderPlus Database (~2 GB)
AMRFINDER_DB_PATH="${DB_BASE_DIR}/amrfinderplus"
if [ ! -d "$AMRFINDER_DB_PATH" ]; then
    echo "--- Setting up AMRFinderPlus Database ---"
    mkdir -p "$AMRFINDER_DB_PATH"
    amrfinder_update --database "$AMRFINDER_DB_PATH"
else
    echo "--- AMRFinderPlus Database already exists. Updating... ---"
    amrfinder_update --database "$AMRFINDER_DB_PATH"
fi

# 3. Bakta Database (~30 GB)
BAKTA_DB_PATH="${DB_BASE_DIR}/bakta_db"
if [ ! -d "$BAKTA_DB_PATH" ]; then
    echo "--- Setting up Bakta Database ---"
    bakta_db download --output "$BAKTA_DB_PATH"
else
    echo "--- Bakta Database already exists. Skipping. ---"
fi

# 4. eggNOG-mapper Database (~50 GB)
EGGNOG_DATA_DIR="${DB_BASE_DIR}/eggnog_data"
if [ ! -d "$EGGNOG_DATA_DIR" ]; then
    echo "--- Setting up eggNOG-mapper Database ---"
    mkdir -p "$EGGNOG_DATA_DIR"
    download_eggnog_data.py --data_dir "$EGGNOG_DATA_DIR" -y
else
    echo "--- eggNOG-mapper Database already exists. Skipping. ---"
fi

# 5. ABRicate Databases (~1 GB)
echo "--- Setting up ABRicate Databases ---"
abricate-get_db

# --- Final Instructions ---
echo "--- Database Setup Complete ---"
echo "Please update the database paths in 'bash_scripts/run_pipeline.sh' to the following:"
echo "KRAKEN2_DB_PATH=\"$KRAKEN2_DB_PATH\""
echo "AMRFINDER_DB_PATH=\"$AMRFINDER_DB_PATH\""
echo "BAKTA_DB_PATH=\"${BAKTA_DB_PATH}/db\""
echo "EGGNOG_DATA_DIR=\"$EGGNOG_DATA_DIR\""
