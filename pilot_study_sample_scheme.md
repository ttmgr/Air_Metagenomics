# Raw Data Access, Basecalling, and Sample Mapping (Greenhouse vs. Natural Environment - PRJEB76446)

This project analyzed air samples collected from a controlled Greenhouse environment and a nearby Natural Environment using Oxford Nanopore Technologies (ONT) sequencing. Samples were collected at different time points: 1-hour (1h) and 3-hour (3h) for Greenhouse; 3-hour (3h) and 6-hour (6h) for Natural Environment.

Two forms of data are available:

1.  **Raw Signal Data:** Unprocessed `.fast5` files from the sequencing runs, grouped and archived by run/condition.
2.  **Processed Sequence Data:** Basecalled and demultiplexed FASTQ files (`.fastq`), as described in the tables below. These may have been generated using ONT Guppy basecaller v6.3.2.

This document first describes how to access and process the raw `.fast5` files using Guppy to generate FASTQ data. Subsequently, it provides the mapping between the resulting FASTQ files and the corresponding sample information, including specific instructions for handling the Natural Environment samples.

## Raw FAST5 Data and Basecalling Instructions

The raw sequencing signal data is available as multi-part `.fast5` files. These files are grouped by experimental condition (Greenhouse 1h, Greenhouse 3h, Natural Environment 3h, Natural Environment 6h) and archived into potentially large, multi-part TAR archives.

**File Naming Convention:**
The archives likely follow a pattern like:
* `greenhouse_1h_partXY.tar.gz`
* `greenhouse_3h_partXY.tar.gz`
* `natural_environment_3h_partXY.tar.gz`
* `natural_environment_6h_partXY.tar.gz`

**Important:**
* You need to download **all parts** associated with a specific experimental condition (e.g., all `greenhouse_1h_part*.tar.gz` files) to reconstruct the complete raw dataset for that condition.
* Each sequencing run likely contains multiplexed samples prepared using the ONT Rapid Barcoding Kit (SQK-RBK114-24).

**Processing Steps using Guppy:**

1.  **Download & Prepare:** Download all TAR parts for the experimental condition of interest. Extract them so that all `.fast5` files for that condition reside within a single input directory.
2.  **Basecalling & Demultiplexing:** Use ONT's Guppy basecaller (specifically version 6.3.2 recommended) to convert the raw signal data into basecalled, demultiplexed FASTQ files. Guppy can perform basecalling and demultiplexing simultaneously.
    * Specify the input directory containing the `.fast5` files (`-i`).
    * Specify a save directory for the output (`-s`). Guppy will create subdirectories here for each barcode.
    * Specify the appropriate configuration file for the chemistry and flowcell used (e.g., `-c dna_r10.4.1_e8.2_400bps_hac.cfg` for R10.4.1 flowcells, E8.2 chemistry, 400bps speed, high accuracy model).
    * Specify the barcoding kit used (`--barcode_kits "SQK-RBK114-24"`).
    * Enable barcode trimming (`--trim_barcodes`).
    * Ensure FASTQ output (often default, or use `--compress_fastq` for gzipped FASTQ).
    * Allocate sufficient compute resources (CPU: `--num_callers`, `--cpu_threads_per_caller`; or GPU: `--device "cuda:all"` or similar).

    ```bash
    # Example Guppy command (adjust config, paths, resources, and kit as needed)
    # Basecalls FAST5, demultiplexes using SQK-RBK114-24, trims barcodes, and saves FASTQ files
    guppy_basecaller \
      -i <path_to_fast5_directory> \
      -s <output_save_directory> \
      -c dna_r10.4.1_e8.2_400bps_hac.cfg \
      --barcode_kits "SQK-RBK114-24" \
      --trim_barcodes \
      --num_callers <N> \
      --cpu_threads_per_caller <M>
      # Or use GPU options e.g.: --device "cuda:0"
    ```
3.  **Output:** Guppy will create a directory structure within your `<output_save_directory>`, such as `barcode01`, `barcode02`, ..., `unclassified`. Inside each `barcodeXX` folder, you will find the corresponding basecalled FASTQ file(s).

Once you have generated these demultiplexed FASTQ files, you can use the mapping tables below to identify which sample corresponds to each barcode file.

## Processed Data Mapping and Sample Information

This section details the mapping between barcodes and the original Greenhouse or Natural Environment samples. This mapping applies to the FASTQ files generated from the raw FAST5 data using Guppy (as described above).

The processed FASTQ files (`.fastq`) may originally have been distributed within **four separate TAR archives** corresponding to the experimental conditions (assumed names provided for reference):

1.  **Greenhouse 1h:** `greenhouse_1h_samples.tar.gz`
2.  **Greenhouse 3h:** `greenhouse_3h_samples.tar.gz`
3.  **Natural Environment 3h:** `natural_environment_3h_samples.tar.gz`
4.  **Natural Environment 6h:** `natural_environment_6h_samples.tar.gz`

If using pre-generated FASTQ files, unpack the relevant archive(s) (e.g., `tar -xzvf <archive_name>.tar.gz`) to access the files listed below.

**Crucial Note for Natural Environment Samples:** Each sample replicate from the natural environment (both 3h and 6h collections) was sequenced using **two different barcodes**. To obtain the complete data for a single natural environment replicate, users **must concatenate** the corresponding pair of FASTQ files (whether generated by Guppy or downloaded). These pairs correspond to the barcodes listed in the "Natural Environment Samples" table below.

---

## Greenhouse Samples

Samples collected from the greenhouse environment. The 1h samples correspond to files generated under `barcodeXX` folders from the `greenhouse_1h` FAST5 data (or found in `greenhouse_1h_samples.tar.gz`). The 3h samples correspond to files generated from the `greenhouse_3h` FAST5 data (or found in `greenhouse_3h_samples.tar.gz`). Each barcode listed represents the complete dataset for that specific replicate.

| Sample Name       | Collection Time | Replicate | Barcode   | Expected FASTQ Filename Pattern (from Guppy or Archive) |
| :---------------- | :-------------- | :-------- | :-------- | :------------------------------------------------------ |
| Greenhouse_1h_Rep1 | 1h              | 1         | barcode03 | `.../barcode03/SQK-RBK114-24_greenhouse_1h_1_barcode03*.fastq.gz` (or similar) |
| Greenhouse_1h_Rep2 | 1h              | 2         | barcode04 | `.../barcode04/SQK-RBK114-24_greenhouse_1h_2_barcode04*.fastq.gz` |
| Greenhouse_1h_Rep3 | 1h              | 3         | barcode05 | `.../barcode05/SQK-RBK114-24_greenhouse_1h_3_barcode05*.fastq.gz` |
| Greenhouse_1h_Rep4 | 1h              | 4         | barcode08 | `.../barcode08/SQK-RBK114-24_greenhouse_1h_4_barcode08*.fastq.gz` |
| Greenhouse_1h_Rep5 | 1h              | 5         | barcode09 | `.../barcode09/SQK-RBK114-24_greenhouse_1h_5_barcode09*.fastq.gz` |
| Greenhouse_1h_Rep6 | 1h              | 6         | barcode10 | `.../barcode10/SQK-RBK114-24_greenhouse_1h_6_barcode10*.fastq.gz` |
| Greenhouse_1h_Rep7 | 1h              | 7         | barcode13 | `.../barcode13/SQK-RBK114-24_greenhouse_1h_7_barcode13*.fastq.gz` |
| Greenhouse_1h_Rep8 | 1h              | 8         | barcode14 | `.../barcode14/SQK-RBK114-24_greenhouse_1h_8_barcode14*.fastq.gz` |
| Greenhouse_1h_Rep9 | 1h              | 9         | barcode15 | `.../barcode15/SQK-RBK114-24_greenhouse_1h_9_barcode15*.fastq.gz` |
| Greenhouse_3h_Rep1 | 3h              | 1         | barcode24 | `.../barcode24/SQK-RBK114-24_greenhouse_3h_1_barcode24*.fastq.gz` |
| Greenhouse_3h_Rep2 | 3h              | 2         | barcode22 | `.../barcode22/SQK-RBK114-24_greenhouse_3h_2_barcode22*.fastq.gz` |
| Greenhouse_3h_Rep3 | 3h              | 3         | barcode01 | `.../barcode01/SQK-RBK114-24_greenhouse_3h_3_barcode01*.fastq.gz` |

---

## Natural Environment Samples

Samples collected from the natural environment. The 3h samples correspond to files generated under `barcodeXX` folders from the `natural_environment_3h` FAST5 data (or found in `natural_environment_3h_samples.tar.gz`). The 6h samples correspond to files generated from the `natural_environment_6h` FAST5 data (or found in `natural_environment_6h_samples.tar.gz`).

**ACTION REQUIRED FOR ALL NATURAL ENVIRONMENT SAMPLES:** Remember that each sample replicate listed below consists of **two** FASTQ files, identified by the barcode pair. You must combine the pair of files for each replicate before analysis.

| Sample Name       | Collection Time | Replicate | Barcodes (Pair)     | Expected FASTQ Files to Concatenate (from Guppy output folders or Archive)                                                                                                                               |
| :---------------- | :-------------- | :-------- | :------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| NaturalEnv_3h_Rep1 | 3h              | 1         | barcode01, barcode02 | `.../barcode01/SQK-RBK114-24_naturalenvironment_3h_1_barcode01*.fastq.gz`<br>`.../barcode02/SQK-RBK114-24_naturalenvironment_3h_1_barcode02*.fastq.gz` |
| NaturalEnv_3h_Rep2 | 3h              | 2         | barcode03, barcode04 | `.../barcode03/SQK-RBK114-24_naturalenvironment_3h_2_barcode03*.fastq.gz`<br>`.../barcode04/SQK-RBK114-24_naturalenvironment_3h_2_barcode04*.fastq.gz` |
| NaturalEnv_3h_Rep3 | 3h              | 3         | barcode05, barcode06 | `.../barcode05/SQK-RBK114-24_naturalenvironment_3h_3_barcode05*.fastq.gz`<br>`.../barcode06/SQK-RBK114-24_naturalenvironment_3h_3_barcode06*.fastq.gz` |
| NaturalEnv_6h_Rep1 | 6h              | 1         | barcode11, barcode12 | `.../barcode11/SQK-RBK114-24_naturalenvironment_6h_1_barcode11*.fastq.gz`<br>`.../barcode12/SQK-RBK114-24_naturalenvironment_6h_1_barcode12*.fastq.gz` |
| NaturalEnv_6h_Rep2 | 6h              | 2         | barcode13, barcode14 | `.../barcode13/SQK-RBK114-24_naturalenvironment_6h_2_barcode13*.fastq.gz`<br>`.../barcode14/SQK-RBK114-24_naturalenvironment_6h_2_barcode14*.fastq.gz` |
| NaturalEnv_6h_Rep3 | 6h              | 3         | barcode15, barcode16 | `.../barcode15/SQK-RBK114-24_naturalenvironment_6h_3_barcode15*.fastq.gz`<br>`.../barcode16/SQK-RBK114-24_naturalenvironment_6h_3_barcode16*.fastq.gz` |


**Example Concatenation Commands:**

Retrieve the FASTQ files from the Guppy output directories (e.g., `<output_save_directory>/barcodeXX/`).

```bash
# Example for NaturalEnv_3h_Rep1 (assuming gzipped FASTQ output from Guppy)
zcat <output_save_directory>/barcode01/SQK-RBK114-24_naturalenvironment_3h_1_barcode01*.fastq.gz \
     <output_save_directory>/barcode02/SQK-RBK114-24_naturalenvironment_3h_1_barcode02*.fastq.gz \
     > NaturalEnv_3h_Rep1_combined.fastq

# Example for NaturalEnv_6h_Rep1 (assuming gzipped FASTQ output from Guppy)
zcat <output_save_directory>/barcode11/SQK-RBK114-24_naturalenvironment_6h_1_barcode11*.fastq.gz \
     <output_save_directory>/barcode12/SQK-RBK114-24_naturalenvironment_6h_1_barcode12*.fastq.gz \
     > NaturalEnv_6h_Rep1_combined.fastq

# If Guppy produced uncompressed FASTQ, use 'cat' instead of 'zcat'
# cat <output_save_directory>/barcode01/....fastq <output_save_directory>/barcode02/....fastq > NaturalEnv_3h_Rep1_combined.fastq
