#!/bin/bash
# Stage 4: Functional Annotation
# Annotates both reads and assemblies for function, AMR, and virulence factors.

set -euo pipefail

# --- Directory Setup ---
NANOFILT_DIR="${OUTPUT_BASE_DIR}/02_nanofilt"
RACON_DIR="${OUTPUT_BASE_DIR}/07_racon"

ANNOTATION_DIR="${OUTPUT_BASE_DIR}/10_annotation"
READ_FASTA_DIR="${ANNOTATION_DIR}/read_fasta"
READ_AMR_DIR="${ANNOTATION_DIR}/read_amr"
CONTIG_KRAKEN_DIR="${ANNOTATION_DIR}/contig_kraken"
CONTIG_PROKKA_DIR="${ANNOTATION_DIR}/contig_prokka"
CONTIG_PRODIGAL_DIR="${ANNOTATION_DIR}/contig_prodigal"
CONTIG_EGGNOG_DIR="${ANNOTATION_DIR}/contig_eggnog"
CONTIG_BAKTA_DIR="${ANNOTATION_DIR}/contig_bakta"
CONTIG_AMR_DIR="${ANNOTATION_DIR}/contig_amr"

mkdir -p "$READ_FASTA_DIR" "$READ_AMR_DIR" "$CONTIG_KRAKEN_DIR" "$CONTIG_PROKKA_DIR" "$CONTIG_PRODIGAL_DIR" "$CONTIG_EGGNOG_DIR" "$CONTIG_BAKTA_DIR" "$CONTIG_AMR_DIR"

echo "Stage 4: Running annotation on reads and assemblies."

# --- Annotation of Filtered Reads ---
echo "--- Annotating Filtered Reads ---"
for filtered_fq in ${NANOFILT_DIR}/filtered_barcode*.fastq; do
    if [ -f "$filtered_fq" ]; then
        base_name=$(basename "$filtered_fq" .fastq | sed 's/filtered_//')
        echo "Annotating reads for sample: $base_name"

        # 1. Convert FASTQ to FASTA for AMR tools
        read_fasta="${READ_FASTA_DIR}/${base_name}.fasta"
        echo "  -> Converting FASTQ to FASTA..."
        seqkit fq2fa "$filtered_fq" -o "$read_fasta"

        # 2. Run Abricate and AMRFinder+ on reads
        abricate_out="${READ_AMR_DIR}/${base_name}_abricate.txt"
        amrfinder_out="${READ_AMR_DIR}/${base_name}_amrfinder.txt"
        echo "  -> Running Abricate on reads..."
        abricate --threads "$THREADS" "$read_fasta" > "$abricate_out"
        echo "  -> Running AMRFinder+ on reads..."
        amrfinder --threads "$THREADS" -n "$read_fasta" -d "$AMRFINDER_DB_PATH" > "$amrfinder_out"
    fi
done

# --- Annotation of Polished Assemblies ---
echo "--- Annotating Polished Assemblies ---"
for polished_fasta in ${RACON_DIR}/polished_barcode*.fasta; do
    if [ -f "$polished_fasta" ]; then
        base_name=$(basename "$polished_fasta" .fasta | sed 's/polished_//')
        echo "Annotating assembly for sample: $base_name"

        # 1. Kraken2 on Contigs
        kraken_report_contig="${CONTIG_KRAKEN_DIR}/report_${base_name}.txt"
        echo "  -> Running Kraken2 on contigs..."
        kraken2 --db "${KRAKEN2_DB_PATH}" --threads "$THREADS" --use-names --report "$kraken_report_contig" "$polished_fasta" > /dev/null

        # 2. Prokka Annotation
        prokka_out_dir="${CONTIG_PROKKA_DIR}/prokka_${base_name}"
        echo "  -> Running Prokka..."
        prokka --outdir "$prokka_out_dir" --prefix "$base_name" --cpus "$THREADS" "$polished_fasta"

        # 3. Bakta Annotation
        bakta_out_dir="${CONTIG_BAKTA_DIR}/bakta_${base_name}"
        echo "  -> Running Bakta..."
        bakta --db "$BAKTA_DB_PATH" --output "$bakta_out_dir" "$polished_fasta" --threads "$THREADS"

        # 4. Prodigal + eggNOG-mapper
        proteins_faa="${CONTIG_PRODIGAL_DIR}/${base_name}.faa"
        echo "  -> Running Prodigal for gene prediction..."
        prodigal -i "$polished_fasta" -a "$proteins_faa" -p meta -q
        
        eggnog_out="${CONTIG_EGGNOG_DIR}/${base_name}"
        echo "  -> Running eggNOG-mapper..."
        emapper.py -i "$proteins_faa" --output "$eggnog_out" -m diamond --cpu "$THREADS" --data_dir "$EGGNOG_DATA_DIR" --override

        # 5. Abricate and AMRFinder+ on Contigs
        abricate_contig_out="${CONTIG_AMR_DIR}/${base_name}_abricate.txt"
        amrfinder_contig_out="${CONTIG_AMR_DIR}/${base_name}_amrfinder.txt"
        echo "  -> Running Abricate on contigs..."
        abricate --threads "$THREADS" "$polished_fasta" > "$abricate_contig_out"
        echo "  -> Running AMRFinder+ on contigs..."
        amrfinder --threads "$THREADS" -n "$polished_fasta" -d "$AMRFINDER_DB_PATH" > "$amrfinder_contig_out"
    fi
done

echo "Stage 4: Annotation complete."
