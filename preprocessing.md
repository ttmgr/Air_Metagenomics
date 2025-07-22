# Preprocessing Pipeline Documentation

This document provides detailed information about the preprocessing steps in the nanopore metagenomics pipeline.

## Overview

The preprocessing pipeline transforms raw nanopore sequencing data into high-quality reads ready for downstream analysis. The main steps include:

1. **Basecalling**: Converting raw electrical signals to DNA sequences
2. **Demultiplexing**: Separating reads by barcode
3. **Adapter Removal**: Removing sequencing adapters
4. **Quality Filtering**: Removing low-quality reads
5. **Quality Control**: Generating QC reports

## Detailed Steps

### 1. Basecalling with Dorado

Dorado is Oxford Nanopore's high-accuracy basecaller that uses deep learning models to convert raw electrical signals (in FAST5 or POD5 format) into DNA sequences.

#### Standard Basecalling
```bash
dorado basecaller \
    --emit-fastq \
    dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
    -r input_folder/ \
    > basecalled.fastq \
    --kit-name SQK-RBK114-24 \
    --no-trim
```

#### Duplex Basecalling
For higher accuracy, duplex basecalling can be enabled:
```bash
dorado duplex \
    dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
    -r input_folder/ \
    > basecalled_duplex.bam \
    -t 2
```

**Key Parameters:**
- `--emit-fastq`: Output in FASTQ format (default is BAM)
- `--kit-name`: Specifies the sequencing kit used
- `--no-trim`: Disables automatic adapter trimming (we use Porechop later)
- Model selection depends on your flowcell and kit

**Duplex Read Identification:**
- `dx:i:1`: Duplex reads
- `dx:i:0`: Simplex reads without duplex offspring
- `dx:i:-1`: Simplex reads with duplex offspring

### 2. Demultiplexing

Separates multiplexed reads into individual barcode bins:

```bash
dorado demux \
    --output-dir basecalled/ \
    --emit-fastq \
    --kit-name SQK-RBK114-24 \
    basecalled.fastq
```

**Output:**
- Individual FASTQ files for each barcode
- `unclassified.fastq` for reads without clear barcode assignment

### 3. Adapter Removal with Porechop

Porechop identifies and removes adapter sequences that may remain after basecalling:

```bash
porechop \
    -i barcode01.fastq \
    -o trimmed_barcode01.fastq \
    --threads 8 \
    --verbosity 2
```

**Key Features:**
- Searches for adapters at both ends of reads
- Can split chimeric reads
- Removes middle adapters
- Quality-aware trimming

**Important Options:**
- `--no_split`: Disable splitting of chimeric reads
- `--discard_middle`: Remove reads with middle adapters
- `--min_split_read_size`: Minimum size of split reads to keep

### 4. Quality Filtering with NanoFilt

NanoFilt filters reads based on quality and length criteria:

```bash
cat trimmed_barcode01.fastq | \
NanoFilt \
    -q 9 \
    -l 100 \
    --maxlength 50000 \
    > filtered_barcode01.fastq
```

**Filtering Parameters:**
- `-q 9`: Minimum average read quality score
- `-l 100`: Minimum read length
- `--maxlength 50000`: Maximum read length
- `--headcrop 50`: Remove first 50 bases
- `--tailcrop 50`: Remove last 50 bases

### 5. Quality Control with NanoStat

Generate comprehensive statistics about your reads:

```bash
NanoStat \
    --fastq filtered_barcode01.fastq \
    --outdir qc_reports/ \
    --name barcode01 \
    --threads 8
```

**Output Metrics:**
- Number of reads
- Total bases
- Read length distribution (mean, median, N50)
- Quality score distribution
- Read length vs quality plots

## Quality Control Checkpoints

### 1. Post-Basecalling QC
- Check basecalling quality scores
- Verify expected number of reads
- Check for failed basecalling runs

### 2. Post-Demultiplexing QC
- Verify barcode distribution
- Check percentage of unclassified reads
- Ensure balanced representation across samples

### 3. Post-Filtering QC
- Compare read counts before/after filtering
- Check read length distribution
- Verify quality score improvement

## Best Practices

### 1. Basecalling Model Selection
Choose the appropriate model based on:
- Flowcell type (R9.4.1, R10.3, R10.4.1)
- Sequencing kit
- Desired accuracy vs speed trade-off

### 2. Quality Thresholds
Recommended settings for environmental samples:
- Minimum quality: Q9-Q10
- Minimum length: 100-500 bp
- Maximum length: 30,000-50,000 bp

### 3. Adapter Removal
- Always verify adapter removal with spot checks
- Consider manual inspection of a subset of reads
- Be cautious with aggressive trimming

### 4. Data Retention
- Keep raw data (FAST5/POD5) for potential re-basecalling
- Archive trimmed but unfiltered reads
- Document all preprocessing parameters

## Troubleshooting

### Low Read Yield
- Check sequencing run metrics
- Verify correct basecalling model
- Adjust quality filtering thresholds

### High Adapter Content
- Verify correct kit specification
- Check Porechop parameters
- Consider manual adapter sequences

### Unbalanced Barcodes
- Check library preparation
- Verify demultiplexing parameters
- Consider adjusting barcode detection threshold

## Output Structure

```
01_preprocessing/
├── basecalled/
│   ├── basecalled.fastq
│   └── sequencing_summary.txt
├── demultiplexed/
│   ├── barcode01.fastq
│   ├── barcode02.fastq
│   └── unclassified.fastq
├── trimmed/
│   ├── trimmed_barcode01.fastq
│   └── trimmed_barcode02.fastq
├── filtered/
│   ├── filtered_barcode01.fastq
│   └── filtered_barcode02.fastq
└── qc/
    ├── barcode01_nanostats.txt
    ├── barcode02_nanostats.txt
    └── preprocessing_summary.html
```

## Advanced Options

### Custom Basecalling Parameters
```bash
dorado basecaller \
    --modified-bases 5mCG_5hmCG \
    --min-qscore 8 \
    --chunk-size 8000 \
    --overlap 500 \
    ...
```

### Parallel Processing
For large datasets, consider parallel processing:
```bash
# Split FAST5 files into batches
# Run multiple basecalling jobs
# Merge results
```

### Real-time Analysis
Enable real-time basecalling and analysis:
```bash
dorado basecaller \
    --watch \
    --recursive \
    ...
```

## References

1. [Dorado Documentation](https://github.com/nanoporetech/dorado)
2. [Porechop GitHub](https://github.com/rrwick/Porechop)
3. [NanoFilt Documentation](https://github.com/wdecoster/nanofilt)
4. [NanoStat Documentation](https://github.com/wdecoster/nanostat)
