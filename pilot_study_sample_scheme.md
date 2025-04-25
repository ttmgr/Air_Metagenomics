# Greenhouse vs. Natural Environment Air Sample Analysis (PRJEBXXXXX)

This project analyzed air samples collected from a controlled Greenhouse environment and a nearby Natural Environment using Oxford Nanopore Technologies (ONT) sequencing. Samples were collected at different time points for each environment and are grouped into separate archives based on environment and collection duration.

* **Greenhouse:** 1-hour (1h) and 3-hour (3h) collections.
* **Natural Environment:** 3-hour (3h) and 6-hour (6h) collections.

Samples were prepared and barcoded using the ONT Rapid Barcoding Kit (SQK-RBK114-24), allowing for multiplexed sequencing.

The raw sequencing data is provided as individual FASTQ (`.fastq`) files grouped into **four separate TAR archives**:

1.  `greenhouse_1h_samples.tar.gz` (Assumed name)
2.  `greenhouse_3h_samples.tar.gz` (Assumed name)
3.  `natural_environment_3h_samples.tar.gz` (Assumed name)
4.  `natural_environment_6h_samples.tar.gz` (Assumed name)

Users will need to unpack the relevant archive(s) (e.g., `tar -xzvf <archive_name>.tar.gz`) to access the files listed below.

**Important Note for Natural Environment Samples:** Each sample from the natural environment (both 3h and 6h collections) was sequenced using **two different barcodes**. To obtain the complete data for a single natural environment replicate, users **must concatenate** the corresponding pair of FASTQ files from the *same* archive.

---

## Greenhouse Samples

### Greenhouse 1h Samples

Samples collected from the greenhouse environment after 1 hour. These files are contained within the `greenhouse_1h_samples.tar.gz` archive (assumed filename). Each file represents a complete dataset for the specified replicate.

| Sample Name        | Collection Time | Replicate | Barcode   | Filename (within archive)                       |
| :----------------- | :-------------- | :-------- | :-------- | :---------------------------------------------- |
| Greenhouse_1h_Rep1 | 1h              | 1         | barcode03 | `SQK-RBK114-24_greenhouse_1h_1_barcode03.fastq` |
| Greenhouse_1h_Rep2 | 1h              | 2         | barcode04 | `SQK-RBK114-24_greenhouse_1h_2_barcode04.fastq` |
| Greenhouse_1h_Rep3 | 1h              | 3         | barcode05 | `SQK-RBK114-24_greenhouse_1h_3_barcode05.fastq` |
| Greenhouse_1h_Rep4 | 1h              | 4         | barcode08 | `SQK-RBK114-24_greenhouse_1h_4_barcode08.fastq` |
| Greenhouse_1h_Rep5 | 1h              | 5         | barcode09 | `SQK-RBK114-24_greenhouse_1h_5_barcode09.fastq` |
| Greenhouse_1h_Rep6 | 1h              | 6         | barcode10 | `SQK-RBK114-24_greenhouse_1h_6_barcode10.fastq` |
| Greenhouse_1h_Rep7 | 1h              | 7         | barcode13 | `SQK-RBK114-24_greenhouse_1h_7_barcode13.fastq` |
| Greenhouse_1h_Rep8 | 1h              | 8         | barcode14 | `SQK-RBK114-24_greenhouse_1h_8_barcode14.fastq` |
| Greenhouse_1h_Rep9 | 1h              | 9         | barcode15 | `SQK-RBK114-24_greenhouse_1h_9_barcode15.fastq` |

### Greenhouse 3h Samples

Samples collected from the greenhouse environment after 3 hours. These files are contained within the `greenhouse_3h_samples.tar.gz` archive (assumed filename). Each file represents a complete dataset for the specified replicate.

| Sample Name        | Collection Time | Replicate | Barcode   | Filename (within archive)                       |
| :----------------- | :-------------- | :-------- | :-------- | :---------------------------------------------- |
| Greenhouse_3h_Rep1 | 3h              | 1         | barcode24 | `SQK-RBK114-24_greenhouse_3h_1_barcode24.fastq` |
| Greenhouse_3h_Rep2 | 3h              | 2         | barcode22 | `SQK-RBK114-24_greenhouse_3h_2_barcode22.fastq` |
| Greenhouse_3h_Rep3 | 3h              | 3         | barcode01 | `SQK-RBK114-24_greenhouse_3h_3_barcode01.fastq` |

---

## Natural Environment Samples

**ACTION REQUIRED FOR ALL NATURAL ENVIRONMENT SAMPLES:** Remember that each sample replicate listed below consists of **two** FASTQ files found within the *same* archive. You must combine the pair of files for each replicate before analysis.

Example concatenation command using `cat`:
```bash
cat <file_barcode_A>.fastq <file_barcode_B>.fastq > <combined_sample_name>.fastq
