#!/bin/bash
# Stage 3: Metagenome Binning
# Uses MetaWRAP to bin contigs into Metagenome-Assembled Genomes (MAGs).

set -euo pipefail

# --- Directory Setup ---
NANOFILT_DIR="${OUTPUT_BASE_DIR}/02_nanofilt"
RACON_DIR="${OUTPUT_BASE_DIR}/07_racon"
BINNING_DIR="${OUTPUT_BASE_DIR}/09_metawrap_binning"

mkdir -p "$BINNING_DIR"

echo "Stage 3: Binning polished contigs from $RACON_DIR"

# Loop through all polished assemblies
for polished_fasta in ${RACON_DIR}/polished_barcode*.fasta; do
    if [ -f "$polished_fasta" ]; then
        base_name=$(basename "$polished_fasta" .fasta | sed 's/polished_//')
        echo "Binning sample: $base_name"

        # Find the corresponding filtered reads
        filtered_fq="${NANOFILT_DIR}/filtered_${base_name}.fastq"

        if [ ! -f "$filtered_fq" ]; then
            echo "Error: Cannot find filtered reads for $base_name. Skipping."
            continue
        fi

        # Run MetaWRAP binning
        binning_out_dir="${BINNING_DIR}/${base_name}"
        echo "  -> Running MetaWRAP binning..."
        metawrap binning -o "$binning_out_dir" \
                         -t "$THREADS" \
                         -a "$polished_fasta" \
                         --metabat2 --maxbin2 --concoct \
                         "$filtered_fq"

        # Run MetaWRAP bin refinement
        echo "  -> Running MetaWRAP bin refinement..."
        metawrap bin_refinement -o "${binning_out_dir}/BIN_REFINEMENT" \
                                -t "$THREADS" \
                                -A "${binning_out_dir}/metabat2_bins/" \
                                -B "${binning_out_dir}/maxbin2_bins/" \
                                -C "${binning_out_dir}/concoct_bins/" \
                                -c 50 -x 10
    else
        echo "Warning: No polished FASTA files found in $RACON_DIR"
        break
    fi
done

echo "Stage 3: Metagenome Binning complete."
