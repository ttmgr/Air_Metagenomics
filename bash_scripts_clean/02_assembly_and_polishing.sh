#!/bin/bash
# Stage 2: Assembly and Polishing
# Assembles filtered reads, polishes the assembly, and generates stats.

set -euo pipefail

# --- Directory Setup ---
NANOFILT_DIR="${OUTPUT_BASE_DIR}/02_nanofilt"
FLYE_DIR="${OUTPUT_BASE_DIR}/05_flye"
MINIMAP2_DIR="${OUTPUT_BASE_DIR}/06_minimap2"
RACON_DIR="${OUTPUT_BASE_DIR}/07_racon"
ASSEMBLYSTATS_DIR="${OUTPUT_BASE_DIR}/08_assembly_stats"

mkdir -p "$FLYE_DIR" "$MINIMAP2_DIR" "$RACON_DIR" "$ASSEMBLYSTATS_DIR"

echo "Stage 2: Assembling filtered reads from $NANOFILT_DIR"

# Loop through all filtered FASTQ files
for filtered_fq in ${NANOFILT_DIR}/filtered_barcode*.fastq; do
    if [ -f "$filtered_fq" ]; then
        base_name=$(basename "$filtered_fq" .fastq | sed 's/filtered_//')
        echo "Assembling sample: $base_name"

        # 1. De novo assembly with Flye
        flye_out_dir="${FLYE_DIR}/flye_${base_name}"
        assembly_fasta="${flye_out_dir}/assembly.fasta"
        echo "  -> Running Flye..."
        flye --nano-hq "$filtered_fq" -o "$flye_out_dir" --threads "$THREADS" --meta

        # 2. Polishing with Racon (1 round)
        # Align reads to the assembly
        sam_file="${MINIMAP2_DIR}/${base_name}.sam"
        echo "  -> Running Minimap2 for polishing..."
        minimap2 -ax map-ont -t "$THREADS" "$assembly_fasta" "$filtered_fq" > "$sam_file"

        # Run Racon
        polished_fasta="${RACON_DIR}/polished_${base_name}.fasta"
        echo "  -> Running Racon for polishing..."
        racon -t "$THREADS" "$filtered_fq" "$sam_file" "$assembly_fasta" > "$polished_fasta"

        # 3. Generate Assembly Statistics on the polished assembly
        assemblystats_out="${ASSEMBLYSTATS_DIR}/${base_name}_assemblystats.txt"
        echo "  -> Running assembly-stats..."
        assembly-stats "$polished_fasta" > "$assemblystats_out"
    else
        echo "Warning: No filtered FASTQ files found in $NANOFILT_DIR"
        break
    fi
done

echo "Stage 2: Assembly and Polishing complete."
