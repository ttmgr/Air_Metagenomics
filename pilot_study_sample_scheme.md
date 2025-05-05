# Greenhouse vs. Natural Environment Air Sample Analysis (PRJEB76446)

This project analyzed air samples collected from a controlled Greenhouse environment and a nearby Natural Environment using Oxford Nanopore Technologies (ONT) sequencing. Samples were collected at different time points for each environment.

* **Greenhouse:** 1-hour (1h) and 3-hour (3h) collections.
* **Natural Environment:** 3-hour (3h) and 6-hour (6h) collections.

Samples were prepared and barcoded using the ONT Rapid Barcoding Kit (SQK-RBK114-24), allowing for multiplexed sequencing.

The raw sequencing data is provided as individual FASTQ (`.fastq`) files grouped into **four separate TAR archives**:

1.  **Greenhouse 1h:** `greenhouse_1h_samples.tar.gz` (Assumed name)
2.  **Greenhouse 3h:** `greenhouse_3h_samples.tar.gz` (Assumed name)
3.  **Natural Environment 3h:** `natural_environment_3h_samples.tar.gz` (Assumed name)
4.  **Natural Environment 6h:** `natural_environment_6h_samples.tar.gz` (Assumed name)

Users will need to unpack the relevant archive(s) (e.g., `tar -xzvf <archive_name>.tar.gz`) to access the files listed below.

**Important Note for Natural Environment Samples:** Each sample from the natural environment (both 3h and 6h collections) was sequenced using **two different barcodes**. To obtain the complete data for a single natural environment replicate, users **must concatenate** the corresponding pair of FASTQ files retrieved from the *same* archive (either the 3h or 6h archive).

---

## Greenhouse Samples

Samples collected from the greenhouse environment. The 1h samples are in `greenhouse_1h_samples.tar.gz` and the 3h samples are in `greenhouse_3h_samples.tar.gz` (assumed filenames). Each file represents a complete dataset for the specified replicate and time point.

| Sample Name        | Collection Time | Replicate | Barcode   | Filename (within respective archive)            |
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
| Greenhouse_3h_Rep1 | 3h              | 1         | barcode24 | `SQK-RBK114-24_greenhouse_3h_1_barcode24.fastq` |
| Greenhouse_3h_Rep2 | 3h              | 2         | barcode22 | `SQK-RBK114-24_greenhouse_3h_2_barcode22.fastq` |
| Greenhouse_3h_Rep3 | 3h              | 3         | barcode01 | `SQK-RBK114-24_greenhouse_3h_3_barcode01.fastq` |

---

## Natural Environment Samples

Samples collected from the natural environment. The 3h samples are in `natural_environment_3h_samples.tar.gz` and the 6h samples are in `natural_environment_6h_samples.tar.gz` (assumed filenames).

| Sample Name        | Collection Time | Replicate | Barcodes (Pair)     | Filenames to Concatenate (within respective archive)                                                                                  |
| :----------------- | :-------------- | :-------- | :------------------ | :-------------------------------------------------------------------------------------------------------------------------------------- |
| NaturalEnv_3h_Rep1 | 3h              | 1         | barcode01, barcode02 | `SQK-RBK114-24_naturalenvironment_3h_1_barcode01.fastq`<br>`SQK-RBK114-24_naturalenvironment_3h_1_barcode02.fastq`                    |
| NaturalEnv_3h_Rep2 | 3h              | 2         | barcode03, barcode04 | `SQK-RBK114-24_naturalenvironment_3h_2_barcode03.fastq`<br>`SQK-RBK114-24_naturalenvironment_3h_2_barcode04.fastq`                    |
| NaturalEnv_3h_Rep3 | 3h              | 3         | barcode05, barcode06 | `SQK-RBK114-24_naturalenvironment_3h_3_barcode05.fastq`<br>`SQK-RBK114-24_naturalenvironment_3h_3_barcode06.fastq`                    |
| NaturalEnv_6h_Rep1 | 6h              | 1         | barcode11, barcode12 | `SQK-RBK114-24_naturalenvironment_6h_1_barcode11.fastq`<br>`SQK-RBK114-24_naturalenvironment_6h_1_barcode12.fastq`                    |
| NaturalEnv_6h_Rep2 | 6h              | 2         | barcode13, barcode14 | `SQK-RBK114-24_naturalenvironment_6h_2_barcode13.fastq`<br>`SQK-RBK114-24_naturalenvironment_6h_2_barcode14.fastq`                    |
| NaturalEnv_6h_Rep3 | 6h              | 3         | barcode15, barcode16 | `SQK-RBK114-24_naturalenvironment_6h_3_barcode15.fastq`<br>`SQK-RBK114-24_naturalenvironment_6h_3_barcode16.fastq`                    |

**ACTION REQUIRED FOR ALL NATURAL ENVIRONMENT SAMPLES:** Remember that each sample replicate listed below consists of **two** FASTQ files found within the *same* archive (either the 3h or 6h archive). You must combine the pair of files for each replicate before analysis.

Example concatenation command using `cat`:
```bash
# For a 3h sample:
cat SQK-RBK114-24_naturalenvironment_3h_1_barcode01.fastq SQK-RBK114-24_naturalenvironment_3h_1_barcode02.fastq > NaturalEnv_3h_Rep1_combined.fastq

# For a 6h sample:
cat SQK-RBK114-24_naturalenvironment_6h_1_barcode11.fastq SQK-RBK114-24_naturalenvironment_6h_1_barcode12.fastq > NaturalEnv_6h_Rep1_combined.fastq
