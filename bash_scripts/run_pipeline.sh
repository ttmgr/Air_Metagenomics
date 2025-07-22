#!/bin/bash

#############################################
# Nanopore Metagenomics Pipeline Runner
# Main orchestration script
#############################################

set -euo pipefail

# Default parameters
THREADS=24
CONFIG_FILE="config/config.yaml"
SKIP_BASECALLING=false
SKIP_ASSEMBLY=false
SKIP_BINNING=false
DUPLEX_MODE=false
VERBOSE=false

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to display usage
usage() {
    cat << EOF
Usage: $0 -i INPUT_DIR -o OUTPUT_DIR [OPTIONS]

Required:
    -i INPUT_DIR     Input directory containing raw nanopore data
    -o OUTPUT_DIR    Output directory for results

Optional:
    -t THREADS       Number of threads (default: 24)
    -c CONFIG        Configuration file (default: config/config.yaml)
    --skip-basecalling   Skip basecalling step (if already done)
    --skip-assembly      Skip assembly step
    --skip-binning       Skip binning step
    --duplex             Enable duplex basecalling mode
    -v, --verbose        Verbose output
    -h, --help           Display this help message

Example:
    $0 -i raw_data/ -o results/ -t 32
    $0 -i raw_data/ -o results/ --skip-basecalling --duplex
EOF
    exit 1
}

# Function to log messages
log() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo -e "${BLUE}[${timestamp}]${NC} INFO: $message"
            ;;
        "SUCCESS")
            echo -e "${GREEN}[${timestamp}]${NC} SUCCESS: $message"
            ;;
        "WARNING")
            echo -e "${YELLOW}[${timestamp}]${NC} WARNING: $message"
            ;;
        "ERROR")
            echo -e "${RED}[${timestamp}]${NC} ERROR: $message"
            ;;
    esac
    
    # Also log to file
    echo "[${timestamp}] ${level}: $message" >> "${OUTPUT_DIR}/logs/pipeline.log"
}

# Function to check if required tools are installed
check_dependencies() {
    log "INFO" "Checking dependencies..."
    
    local required_tools=(
        "dorado"
        "porechop"
        "NanoFilt"
        "flye"
        "minimap2"
        "samtools"
        "racon"
        "kraken2"
        "prokka"
        "bakta"
        "emapper.py"
        "abricate"
        "amrfinder"
        "prodigal"
        "NanoStat"
        "metawrap"
    )
    
    local missing_tools=()
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        log "ERROR" "Missing required tools: ${missing_tools[*]}"
        log "ERROR" "Please install missing tools or activate the conda environment"
        exit 1
    fi
    
    log "SUCCESS" "All required tools are installed"
}

# Function to create directory structure
setup_directories() {
    log "INFO" "Setting up directory structure..."
    
    local dirs=(
        "01_preprocessing/basecalled"
        "01_preprocessing/demultiplexed"
        "01_preprocessing/trimmed"
        "01_preprocessing/filtered"
        "01_preprocessing/qc"
        "02_assembly/flye"
        "02_assembly/minimap2"
        "02_assembly/racon"
        "02_assembly/stats"
        "03_classification/kraken2_reads"
        "03_classification/kraken2_contigs"
        "03_classification/diamond"
        "04_annotation/prokka"
        "04_annotation/bakta"
        "04_annotation/prodigal"
        "04_annotation/eggnog"
        "05_amr_detection/abricate"
        "05_amr_detection/amrfinder"
        "06_binning/metawrap"
        "06_binning/checkm"
        "07_analysis/statistics"
        "07_analysis/plots"
        "07_analysis/reports"
        "logs"
        "temp"
    )
    
    for dir in "${dirs[@]}"; do
        mkdir -p "${OUTPUT_DIR}/${dir}"
    done
    
    log "SUCCESS" "Directory structure created"
}

# Function to run preprocessing
run_preprocessing() {
    log "INFO" "Starting preprocessing pipeline..."
    
    # Basecalling and demultiplexing
    if [ "$SKIP_BASECALLING" = false ]; then
        log "INFO" "Running basecalling with Dorado..."
        bash scripts/preprocessing/dorado_basecalling.sh \
            -i "$INPUT_DIR" \
            -o "${OUTPUT_DIR}/01_preprocessing" \
            -t "$THREADS" \
            --duplex "$DUPLEX_MODE" 2>&1 | tee -a "${OUTPUT_DIR}/logs/basecalling.log"
    fi
    
    # Read processing (adapter removal and filtering)
    log "INFO" "Running read processing..."
    bash scripts/preprocessing/read_processing.sh \
        -i "${OUTPUT_DIR}/01_preprocessing/demultiplexed" \
        -o "${OUTPUT_DIR}/01_preprocessing" \
        -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/read_processing.log"
    
    # Quality control
    log "INFO" "Running quality control..."
    bash scripts/preprocessing/quality_control.sh \
        -i "${OUTPUT_DIR}/01_preprocessing/filtered" \
        -o "${OUTPUT_DIR}/01_preprocessing/qc" 2>&1 | tee -a "${OUTPUT_DIR}/logs/qc.log"
    
    log "SUCCESS" "Preprocessing completed"
}

# Function to run assembly
run_assembly() {
    if [ "$SKIP_ASSEMBLY" = true ]; then
        log "INFO" "Skipping assembly step as requested"
        return
    fi
    
    log "INFO" "Starting assembly pipeline..."
    
    bash scripts/assembly/assembly_and_polishing.sh \
        -i "${OUTPUT_DIR}/01_preprocessing/filtered" \
        -o "${OUTPUT_DIR}/02_assembly" \
        -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/assembly.log"
    
    log "SUCCESS" "Assembly completed"
}

# Function to run classification
run_classification() {
    log "INFO" "Starting taxonomic classification..."
    
    bash scripts/classification/kraken2_classification.sh \
        -i "${OUTPUT_DIR}/01_preprocessing/filtered" \
        -a "${OUTPUT_DIR}/02_assembly/racon" \
        -o "${OUTPUT_DIR}/03_classification" \
        -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/classification.log"
    
    log "SUCCESS" "Classification completed"
}

# Function to run annotation
run_annotation() {
    log "INFO" "Starting functional annotation..."
    
    # Run multiple annotation tools in parallel
    (
        bash scripts/annotation/prokka_annotation.sh \
            -i "${OUTPUT_DIR}/02_assembly/racon" \
            -o "${OUTPUT_DIR}/04_annotation/prokka" \
            -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/prokka.log"
    ) &
    
    (
        bash scripts/annotation/bakta_annotation.sh \
            -i "${OUTPUT_DIR}/02_assembly/racon" \
            -o "${OUTPUT_DIR}/04_annotation/bakta" \
            -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/bakta.log"
    ) &
    
    (
        bash scripts/annotation/prodigal_prediction.sh \
            -i "${OUTPUT_DIR}/02_assembly/racon" \
            -o "${OUTPUT_DIR}/04_annotation/prodigal" 2>&1 | tee -a "${OUTPUT_DIR}/logs/prodigal.log"
    ) &
    
    wait  # Wait for all annotation jobs to complete
    
    # Run eggNOG-mapper on Prodigal results
    bash scripts/annotation/eggnog_mapping.sh \
        -i "${OUTPUT_DIR}/04_annotation/prodigal" \
        -o "${OUTPUT_DIR}/04_annotation/eggnog" \
        -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/eggnog.log"
    
    log "SUCCESS" "Annotation completed"
}

# Function to run AMR detection
run_amr_detection() {
    log "INFO" "Starting AMR and virulence detection..."
    
    # Run on reads
    bash scripts/amr_detection/amrfinder_abricate_readlevel.sh \
        -i "${OUTPUT_DIR}/01_preprocessing/filtered" \
        -o "${OUTPUT_DIR}/05_amr_detection" \
        -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/amr_reads.log"
    
    # Run on contigs
    bash scripts/amr_detection/amrfinder_abricate_contiglevel.sh \
        -i "${OUTPUT_DIR}/02_assembly/racon" \
        -o "${OUTPUT_DIR}/05_amr_detection" \
        -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/amr_contigs.log"
    
    log "SUCCESS" "AMR detection completed"
}

# Function to run binning
run_binning() {
    if [ "$SKIP_BINNING" = true ]; then
        log "INFO" "Skipping binning step as requested"
        return
    fi
    
    log "INFO" "Starting metagenomic binning..."
    
    bash scripts/assembly/metagenomic_binning.sh \
        -r "${OUTPUT_DIR}/01_preprocessing/filtered" \
        -a "${OUTPUT_DIR}/02_assembly/racon" \
        -o "${OUTPUT_DIR}/06_binning" \
        -t "$THREADS" 2>&1 | tee -a "${OUTPUT_DIR}/logs/binning.log"
    
    log "SUCCESS" "Binning completed"
}

# Function to run analysis and generate reports
run_analysis() {
    log "INFO" "Starting analysis and report generation..."
    
    # Extract metrics
    python3 scripts/analysis/extract_metrics.py \
        --input-dir "$OUTPUT_DIR" \
        --output-dir "${OUTPUT_DIR}/07_analysis/statistics" \
        2>&1 | tee -a "${OUTPUT_DIR}/logs/analysis.log"
    
    # Generate plots
    python3 scripts/analysis/generate_plots.py \
        --input-dir "${OUTPUT_DIR}/07_analysis/statistics" \
        --output-dir "${OUTPUT_DIR}/07_analysis/plots" \
        2>&1 | tee -a "${OUTPUT_DIR}/logs/plots.log"
    
    # Create final report
    python3 scripts/analysis/create_report.py \
        --input-dir "$OUTPUT_DIR" \
        --output-file "${OUTPUT_DIR}/07_analysis/reports/final_report.html" \
        2>&1 | tee -a "${OUTPUT_DIR}/logs/report.log"
    
    log "SUCCESS" "Analysis and reporting completed"
}

# Function to clean up temporary files
cleanup() {
    log "INFO" "Cleaning up temporary files..."
    
    if [ -d "${OUTPUT_DIR}/temp" ]; then
        rm -rf "${OUTPUT_DIR}/temp"/*
    fi
    
    # Compress log files
    if command -v gzip &> /dev/null; then
        gzip "${OUTPUT_DIR}/logs"/*.log
    fi
    
    log "SUCCESS" "Cleanup completed"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t)
            THREADS="$2"
            shift 2
            ;;
        -c)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --skip-basecalling)
            SKIP_BASECALLING=true
            shift
            ;;
        --skip-assembly)
            SKIP_ASSEMBLY=true
            shift
            ;;
        --skip-binning)
            SKIP_BINNING=true
            shift
            ;;
        --duplex)
            DUPLEX_MODE=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
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

# Check required arguments
if [ -z "${INPUT_DIR:-}" ] || [ -z "${OUTPUT_DIR:-}" ]; then
    echo "Error: Input and output directories are required"
    usage
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Set up logging
mkdir -p "${OUTPUT_DIR}/logs"
exec 1> >(tee -a "${OUTPUT_DIR}/logs/pipeline.log")
exec 2>&1

# Print pipeline header
cat << EOF
╔══════════════════════════════════════════════════════╗
║       Nanopore Metagenomics Pipeline v1.0.0          ║
║                                                      ║
║  Processing environmental microbiome samples         ║
╚══════════════════════════════════════════════════════╝

Input directory:  $INPUT_DIR
Output directory: $OUTPUT_DIR
Threads:          $THREADS
Config file:      $CONFIG_FILE
Duplex mode:      $DUPLEX_MODE

EOF

# Start pipeline
START_TIME=$(date +%s)
log "INFO" "Pipeline started at $(date)"

# Check dependencies
check_dependencies

# Setup directories
setup_directories

# Run pipeline stages
run_preprocessing
run_assembly
run_classification
run_annotation
run_amr_detection
run_binning
run_analysis
cleanup

# Calculate runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
HOURS=$((RUNTIME / 3600))
MINUTES=$(((RUNTIME % 3600) / 60))
SECONDS=$((RUNTIME % 60))

log "SUCCESS" "Pipeline completed successfully!"
log "INFO" "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
log "INFO" "Results available in: $OUTPUT_DIR"
log "INFO" "Final report: ${OUTPUT_DIR}/07_analysis/reports/final_report.html"

# Print summary
cat << EOF

╔══════════════════════════════════════════════════════╗
║                  Pipeline Summary                     ║
╚══════════════════════════════════════════════════════╝

✓ Preprocessing completed
✓ Assembly completed
✓ Taxonomic classification completed
✓ Functional annotation completed
✓ AMR detection completed
✓ Binning completed
✓ Analysis and reporting completed

View the complete results in:
${OUTPUT_DIR}/07_analysis/reports/final_report.html

EOF
