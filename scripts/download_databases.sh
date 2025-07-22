#!/bin/bash

#############################################
# Database Download and Setup Script
# Downloads and configures all required databases
#############################################

set -euo pipefail

# Default parameters
DB_DIR="${HOME}/metagenomics_databases"
THREADS=8
DOWNLOAD_ALL=true

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Options:
    -d DB_DIR       Database directory (default: ~/metagenomics_databases)
    -t THREADS      Number of threads (default: 8)
    --skip-kraken2  Skip Kraken2 database
    --skip-diamond  Skip DIAMOND database
    --skip-eggnog   Skip eggNOG database
    --skip-bakta    Skip Bakta database
    --skip-amr      Skip AMR databases
    -h, --help      Display this help message

Example:
    $0 -d /data/databases -t 16
    $0 --skip-diamond --skip-eggnog
EOF
    exit 1
}

# Logging function
log() {
    local level=$1
    shift
    local message="$@"
    
    case $level in
        "INFO")
            echo -e "${BLUE}[INFO]${NC} $message"
            ;;
        "SUCCESS")
            echo -e "${GREEN}[SUCCESS]${NC} $message"
            ;;
        "WARNING")
            echo -e "${YELLOW}[WARNING]${NC} $message"
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} $message"
            ;;
    esac
}

# Function to check disk space
check_disk_space() {
    local required_gb=$1
    local target_dir=$2
    
    local available_kb=$(df -k "$target_dir" | awk 'NR==2 {print $4}')
    local available_gb=$((available_kb / 1024 / 1024))
    
    if [ "$available_gb" -lt "$required_gb" ]; then
        log "ERROR" "Insufficient disk space. Required: ${required_gb}GB, Available: ${available_gb}GB"
        exit 1
    fi
    
    log "INFO" "Disk space check passed. Available: ${available_gb}GB"
}

# Function to download and setup Kraken2 database
setup_kraken2() {
    log "INFO" "Setting up Kraken2 database..."
    
    local kraken_dir="${DB_DIR}/kraken2"
    mkdir -p "$kraken_dir"
    
    # Check disk space (standard database ~100GB)
    check_disk_space 150 "$kraken_dir"
    
    cd "$kraken_dir"
    
    # Download prebuilt standard database
    log "INFO" "Downloading Kraken2 standard database (~100GB)..."
    wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz
    
    log "INFO" "Extracting database..."
    tar -xzf k2_standard_20230605.tar.gz
    
    # Clean up
    rm -f k2_standard_20230605.tar.gz
    
    log "SUCCESS" "Kraken2 database setup complete"
}

# Function to download and setup DIAMOND database
setup_diamond() {
    log "INFO" "Setting up DIAMOND database..."
    
    local diamond_dir="${DB_DIR}/diamond"
    mkdir -p "$diamond_dir"
    
    # Check disk space (NR database ~200GB)
    check_disk_space 250 "$diamond_dir"
    
    cd "$diamond_dir"
    
    # Download NCBI NR database
    log "INFO" "Downloading NCBI NR database (~200GB)..."
    wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
    
    log "INFO" "Extracting database..."
    gunzip nr.gz
    
    log "INFO" "Building DIAMOND database..."
    diamond makedb --in nr --db nr --threads "$THREADS"
    
    # Clean up
    rm -f nr
    
    log "SUCCESS" "DIAMOND database setup complete"
}

# Function to download and setup eggNOG database
setup_eggnog() {
    log "INFO" "Setting up eggNOG database..."
    
    local eggnog_dir="${DB_DIR}/eggnog"
    mkdir -p "$eggnog_dir"
    
    # Check disk space (~50GB)
    check_disk_space 60 "$eggnog_dir"
    
    cd "$eggnog_dir"
    
    # Download eggNOG data
    log "INFO" "Downloading eggNOG database..."
    download_eggnog_data.py -y --data_dir .
    
    log "SUCCESS" "eggNOG database setup complete"
}

# Function to download and setup Bakta database
setup_bakta() {
    log "INFO" "Setting up Bakta database..."
    
    local bakta_dir="${DB_DIR}/bakta"
    mkdir -p "$bakta_dir"
    
    # Check disk space (~30GB)
    check_disk_space 40 "$bakta_dir"
    
    cd "$bakta_dir"
    
    # Download Bakta database
    log "INFO" "Downloading Bakta database..."
    bakta_db download --output .
    
    log "SUCCESS" "Bakta database setup complete"
}

# Function to download and setup AMR databases
setup_amr_databases() {
    log "INFO" "Setting up AMR databases..."
    
    local amr_dir="${DB_DIR}/amr"
    mkdir -p "$amr_dir"
    
    # AMRFinderPlus database
    log "INFO" "Setting up AMRFinderPlus database..."
    local amrfinder_dir="${amr_dir}/amrfinderplus"
    mkdir -p "$amrfinder_dir"
    cd "$amrfinder_dir"
    amrfinder_update --database .
    
    # ABRicate databases
    log "INFO" "Updating ABRicate databases..."
    abricate-get_db --db ncbi
    abricate-get_db --db card
    abricate-get_db --db resfinder
    abricate-get_db --db argannot
    abricate-get_db --db vfdb
    
    log "SUCCESS" "AMR databases setup complete"
}

# Function to download and setup CheckM database
setup_checkm() {
    log "INFO" "Setting up CheckM database..."
    
    local checkm_dir="${DB_DIR}/checkm"
    mkdir -p "$checkm_dir"
    
    # Check disk space (~1.4GB)
    check_disk_space 2 "$checkm_dir"
    
    cd "$checkm_dir"
    
    # Download CheckM data
    log "INFO" "Downloading CheckM database..."
    wget -c https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    
    log "INFO" "Extracting database..."
    tar -xzf checkm_data_2015_01_16.tar.gz
    
    # Set CheckM data path
    checkm data setRoot .
    
    # Clean up
    rm -f checkm_data_2015_01_16.tar.gz
    
    log "SUCCESS" "CheckM database setup complete"
}

# Function to create database configuration file
create_db_config() {
    log "INFO" "Creating database configuration file..."
    
    cat > "${DB_DIR}/database_paths.yaml" << EOF
# Database paths for Nanopore Metagenomics Pipeline
# Generated on $(date)

databases:
  kraken2: "${DB_DIR}/kraken2"
  diamond: "${DB_DIR}/diamond/nr.dmnd"
  eggnog: "${DB_DIR}/eggnog"
  bakta: "${DB_DIR}/bakta/db"
  amrfinder: "${DB_DIR}/amr/amrfinderplus"
  checkm: "${DB_DIR}/checkm"
  
# ABRicate databases are managed by ABRicate itself
# Use 'abricate --list' to see available databases
EOF
    
    log "SUCCESS" "Database configuration saved to: ${DB_DIR}/database_paths.yaml"
}

# Parse command line arguments
SKIP_KRAKEN2=false
SKIP_DIAMOND=false
SKIP_EGGNOG=false
SKIP_BAKTA=false
SKIP_AMR=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -d)
            DB_DIR="$2"
            shift 2
            ;;
        -t)
            THREADS="$2"
            shift 2
            ;;
        --skip-kraken2)
            SKIP_KRAKEN2=true
            shift
            ;;
        --skip-diamond)
            SKIP_DIAMOND=true
            shift
            ;;
        --skip-eggnog)
            SKIP_EGGNOG=true
            shift
            ;;
        --skip-bakta)
            SKIP_BAKTA=true
            shift
            ;;
        --skip-amr)
            SKIP_AMR=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Main execution
cat << EOF
╔══════════════════════════════════════════════════════╗
║        Database Setup for Metagenomics Pipeline      ║
╚══════════════════════════════════════════════════════╝

Database directory: $DB_DIR
Threads: $THREADS

This script will download and configure the following databases:
- Kraken2 (~100GB)
- DIAMOND NR (~200GB)
- eggNOG (~50GB)
- Bakta (~30GB)
- AMR databases (~10GB)
- CheckM (~1.4GB)

Total space required: ~400GB

EOF

read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    log "INFO" "Setup cancelled"
    exit 1
fi

# Create main database directory
mkdir -p "$DB_DIR"

# Check total disk space
check_disk_space 400 "$DB_DIR"

# Setup databases
if [ "$SKIP_KRAKEN2" = false ]; then
    setup_kraken2
fi

if [ "$SKIP_DIAMOND" = false ]; then
    setup_diamond
fi

if [ "$SKIP_EGGNOG" = false ]; then
    setup_eggnog
fi

if [ "$SKIP_BAKTA" = false ]; then
    setup_bakta
fi

if [ "$SKIP_AMR" = false ]; then
    setup_amr_databases
fi

# Always setup CheckM (small database)
setup_checkm

# Create configuration file
create_db_config

# Final summary
cat << EOF

╔══════════════════════════════════════════════════════╗
║              Database Setup Complete!                 ║
╚══════════════════════════════════════════════════════╝

All databases have been successfully downloaded and configured.

Database location: $DB_DIR
Configuration file: ${DB_DIR}/database_paths.yaml

To use these databases with the pipeline:
1. Copy the database paths to your config file:
   cp ${DB_DIR}/database_paths.yaml config/database_paths.yaml

2. Or specify the database directory when running the pipeline:
   ./run_pipeline.sh -i input/ -o output/ --db-dir $DB_DIR

EOF

log "SUCCESS" "Setup complete!"
