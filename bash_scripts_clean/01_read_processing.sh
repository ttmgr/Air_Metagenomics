#!/bin/bash
# Stage 1: Read Processing
# Performs adapter trimming, quality filtering, QC stats, and read-based taxonomic classification.

set -euo pipefail

# --- Directory Setup ---
PORECHOP_DIR="${OUTPUT_BASE_DIR}/01_porechop"
NANOFILT_DIR="${OUTPUT_BASE_DIR}/02_nanofilt"
NANOSTAT_DIR="${OUTPUT_BASE_DIR}/03_nanostat"
KRAKEN2_READS_DIR="${OUTPUT_BASE_DIR}/04_kraken2_reads"

mkdir -p "$PORECHOP_DIR" "$NANOFILT_DIR" "$NANOSTAT_DIR" "$KRAKEN2_READS_DIR"

echo "Stage 1: Processing reads from $INPUT_FASTQ_DIR"

# Loop through all input FASTQ files
for fq_file in ${INPUT_FASTQ_DIR}/barcode*.fastq; do
    if [ -f "$fq_file" ]; then
        base_name=$(basename "$fq_file" .fastq)
        echo "Processing sample: $base_name"

        # 1. Adapter Removal with Porechop
        trimmed_fq="${PORECHOP_DIR}/trimmed_${base_name}.fastq"
        echo "  -> Running Porechop..."
        porechop -i "$fq_file" -o "$trimmed_fq" --threads "$THREADS"

        # 2. Quality and Length Filtering with NanoFilt
        filtered_fq="${NANOFILT_DIR}/filtered_${base_name}.fastq"
        echo "  -> Running NanoFilt..."
        cat "$trimmed_fq" | NanoFilt -q 9 -l 500 > "$filtered_fq"

        # 3. Quality Control with NanoStat
        nanostat_out="${NANOSTAT_DIR}/${base_name}_nanostats.txt"
        echo "  -> Running NanoStat..."
        NanoStat --fastq "$filtered_fq" > "$nanostat_out"

        # 4. Taxonomic Classification of Reads with Kraken2
        kraken_report="${KRAKEN2_READS_DIR}/report_${base_name}.txt"
        kraken_output="${KRAKEN2_READS_DIR}/output_${base_name}.txt"
        echo "  -> Running Kraken2 on reads..."
        kraken2 --db "${KRAKEN2_DB_PATH}" \
                --threads "$THREADS" \
                --use-names \
                --report "$kraken_report" \
                --output "$kraken_output" \
                "$filtered_fq"
    else
        echo "Warning: No FASTQ files found in $INPUT_FASTQ_DIR matching 'barcode*.fastq'"
        break
    fi
done

echo "Stage 1: Read Processing complete."
