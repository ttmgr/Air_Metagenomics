# Raw Data Access, Basecalling, and Sample Mapping (Greenhouse vs. Natural Environment - PRJEB76446)

This project analyzed air samples collected from a controlled Greenhouse environment and a nearby Natural Environment using Oxford Nanopore Technologies (ONT) sequencing. Samples were collected at different time points: 1-hour (1h) and 3-hour (3h) for Greenhouse; 3-hour (3h) and 6-hour (6h) for Natural Environment.

Two forms of data related to ENA project PRJEB76446 are available:

1.  **Raw Signal Data:** Unprocessed `.fast5` files from the sequencing runs, grouped and archived by run/condition (as described below). Processing these requires basecalling.
2.  **Processed Sequence Data:** Basecalled and demultiplexed FASTQ files (`.fastq`), available directly from the ENA submission or generated via basecalling. These are grouped into specific archives.

This document first describes how to access and process the raw `.fast5` files using Guppy to generate FASTQ data. Subsequently, it provides the mapping between the FASTQ files (whether downloaded or generated) and the corresponding sample information, referencing the specific archives and filenames confirmed by the project's sample sheet. It also includes crucial instructions for handling the Natural Environment samples.

## Raw FAST5 Data and Basecalling Instructions

The raw sequencing signal data is available as multi-part `.fast5` files. These files are grouped by experimental condition (Greenhouse 1h, Greenhouse 3h, Natural Environment 3h, Natural Environment 6h) and archived into potentially large, multi-part TAR archives.

**File Naming Convention for Raw Data Archives:**
The raw data archives likely follow a pattern like (XY = numerical):
* `greenhouse_1h_partXY.tar.gz`
* `greenhouse_3h_partXY.tar.gz`
* `natural_environment_3h_partXY.tar.gz`
* `natural_environment_6h_partXY.tar.gz`

**Important:**
* You need to download **all parts** associated with a specific experimental condition (e.g., all `greenhouse_1h_part*.tar.gz` files) to reconstruct the complete raw dataset for that condition.
* Each sequencing run contains multiplexed samples prepared using the ONT Rapid Barcoding Kit (SQK-RBK114-24).

**Processing Steps using Guppy:**

1.  **Download & Prepare:** Download all raw data TAR parts for the experimental condition of interest. Extract them so that all `.fast5` files for that condition reside within a single input directory.
2.  **Basecalling & Demultiplexing:** Use ONT's Guppy basecaller (version 6.3.2 was used for generating the processed data) to convert the raw signal data into basecalled, demultiplexed FASTQ files. Guppy performs basecalling and demultiplexing simultaneously.
    * Specify the input directory containing the `.fast5` files (`-i`).
    * Specify a save directory for the output (`-s`). Guppy will create subdirectories here for each barcode (e.g., `barcode01`, `barcode02`, etc.).
    * Specify the appropriate configuration file: `-c dna_r10.4.1_e8.2_400bps_hac.cfg`.
    * Specify the barcoding kit used: `--barcode_kits "SQK-RBK114-24"`.
    * Allocate sufficient compute resources (CPU: `--num_callers`, `--cpu_threads_per_caller`; or GPU: `--device "cuda:all"` or similar).

    ```bash
    # Example Guppy command (adjust paths and resources as needed)
    # Basecalls FAST5, demultiplexes using SQK-RBK114-24, trims barcodes, and saves gzipped FASTQ files
    guppy_basecaller \
      -i <path_to_fast5_directory> \
      -s <output_save_directory> \
      -c dna_r10.4.1_e8.2_400bps_hac.cfg \
      --barcode_kits "SQK-RBK114-24" \
      --num_callers <N> \
      --cpu_threads_per_caller <M>
      # Or use GPU options e.g.: --device "cuda:0"
    ```
3.  **Output:** Guppy will create a directory structure within your `<output_save_directory>`, such as `barcode01`, `barcode02`, ..., `unclassified`. Inside each `barcodeXX` folder, you will find the corresponding basecalled FASTQ file(s). The naming should align with the `Filename` column in the tables below.

Once you have generated these demultiplexed FASTQ files, you can use the mapping tables below to identify which sample corresponds to each barcode file. **Remember the special handling required for Natural Environment samples.**

---

## Processed Data Mapping and Sample Information

This section details the mapping between barcodes and the original Greenhouse or Natural Environment samples, based on the official sample sheet (`ENA_PRJEB76446_pilot_study_sample_sheet.xlsx - Barcode_Mapping.csv`). This mapping applies both to the FASTQ files generated from raw data using Guppy (as described above) and to the pre-generated FASTQ files available from the ENA submission.

The processed FASTQ files (`.fastq`) submitted to ENA are grouped into **four specific TAR archives**:

1.  **Greenhouse 1h:** `greenhouse_1h_samples.tar.gz`
2.  **Greenhouse 3h:** `greenhouse_3h_samples.tar.gz`
3.  **Natural Environment 3h:** `natural_environment_3h_samples.tar.gz`
4.  **Natural Environment 6h:** `natural_environment_6h_samples.tar.gz`

To use the pre-generated data, download the relevant archive(s) and unpack (e.g., `tar -xzvf <archive_name>.tar.gz`) to access the individual FASTQ files listed in the tables below.

**Crucial Note for Natural Environment Samples:** Each sample replicate from the natural environment (both 3h and 6h collections) was sequenced using **two different barcodes**. To obtain the complete data for a single natural environment replicate, users **must concatenate** the corresponding pair of FASTQ files (whether generated by Guppy or extracted from the archives). These pairs are explicitly listed in the "Natural Environment Samples" table below.

---

## Greenhouse Samples

Samples collected from the greenhouse environment. Each FASTQ file listed represents the complete dataset for that specific replicate and time point. Files are located within the specified archive or generated by Guppy into the corresponding barcode folder.

| Sample Name       | Collection Time | Replicate | Barcode   | Filename (within Archive / Guppy Output)        | Archive Name                 |
| :---------------- | :-------------- | :-------- | :-------- | :---------------------------------------------- | :--------------------------- |
| Greenhouse_1h_Rep1 | 1h              | 1         | barcode03 | `SQK-RBK114-24_greenhouse_1h_1_barcode03.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_1h_Rep2 | 1h              | 2         | barcode04 | `SQK-RBK114-24_greenhouse_1h_2_barcode04.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_1h_Rep3 | 1h              | 3         | barcode05 | `SQK-RBK114-24_greenhouse_1h_3_barcode05.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_1h_Rep4 | 1h              | 4         | barcode08 | `SQK-RBK114-24_greenhouse_1h_4_barcode08.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_1h_Rep5 | 1h              | 5         | barcode09 | `SQK-RBK114-24_greenhouse_1h_5_barcode09.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_1h_Rep6 | 1h              | 6         | barcode10 | `SQK-RBK114-24_greenhouse_1h_6_barcode10.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_1h_Rep7 | 1h              | 7         | barcode13 | `SQK-RBK114-24_greenhouse_1h_7_barcode13.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_1h_Rep8 | 1h              | 8         | barcode14 | `SQK-RBK114-24_greenhouse_1h_8_barcode14.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_1h_Rep9 | 1h              | 9         | barcode15 | `SQK-RBK114-24_greenhouse_1h_9_barcode15.fastq` | greenhouse_1h_samples.tar.gz |
| Greenhouse_3h_Rep1 | 3h              | 1         | barcode24 | `SQK-RBK114-24_greenhouse_3h_1_barcode24.fastq` | greenhouse_3h_samples.tar.gz |
| Greenhouse_3h_Rep2 | 3h              | 2         | barcode22 | `SQK-RBK114-24_greenhouse_3h_2_barcode22.fastq` | greenhouse_3h_samples.tar.gz |
| Greenhouse_3h_Rep3 | 3h              | 3         | barcode01 | `SQK-RBK114-24_greenhouse_3h_3_barcode01.fastq` | greenhouse_3h_samples.tar.gz |

---

## Natural Environment Samples

Samples collected from the natural environment. Files are located within the specified archive or generated by Guppy into the corresponding barcode folders.

**ACTION REQUIRED FOR ALL NATURAL ENVIRONMENT SAMPLES:** Each sample replicate listed below consists of **two** FASTQ files found within the *same* archive (or corresponding Guppy barcode output folders). You must combine the pair of files for each replicate before analysis using the `cat` (or `zcat` if gzipped) command.

| Sample Name       | Collection Time | Replicate | Barcodes (Pair)     | Filenames to Concatenate (within Archive / Guppy Output)                                                                                             | Archive Name                        |
| :---------------- | :-------------- | :-------- | :------------------ | :--------------------------------------------------------------------------------------------------------------------------------------------------- | :---------------------------------- |
| NaturalEnv_3h_Rep1 | 3h              | 1         | barcode01, barcode02 | `SQK-RBK114-24_naturalenvironment_3h_1_barcode01.fastq`<br>`SQK-RBK114-24_naturalenvironment_3h_1_barcode02.fastq` | natural_environment_3h_samples.tar.gz |
| NaturalEnv_3h_Rep2 | 3h              | 2         | barcode03, barcode04 | `SQK-RBK114-24_naturalenvironment_3h_2_barcode03.fastq`<br>`SQK-RBK114-24_naturalenvironment_3h_2_barcode04.fastq` | natural_environment_3h_samples.tar.gz |
| NaturalEnv_3h_Rep3 | 3h              | 3         | barcode05, barcode06 | `SQK-RBK114-24_naturalenvironment_3h_3_barcode05.fastq`<br>`SQK-RBK114-24_naturalenvironment_3h_3_barcode06.fastq` | natural_environment_3h_samples.tar.gz |
| NaturalEnv_6h_Rep1 | 6h              | 1         | barcode11, barcode12 | `SQK-RBK114-24_naturalenvironment_6h_1_barcode11.fastq`<br>`SQK-RBK114-24_naturalenvironment_6h_1_barcode12.fastq` | natural_environment_6h_samples.tar.gz |
| NaturalEnv_6h_Rep2 | 6h              | 2         | barcode13, barcode14 | `SQK-RBK114-24_naturalenvironment_6h_2_barcode13.fastq`<br>`SQK-RBK114-24_naturalenvironment_6h_2_barcode14.fastq` | natural_environment_6h_samples.tar.gz |
| NaturalEnv_6h_Rep3 | 6h              | 3         | barcode15, barcode16 | `SQK-RBK114-24_naturalenvironment_6h_3_barcode15.fastq`<br>`SQK-RBK114-24_naturalenvironment_6h_3_barcode16.fastq` | natural_environment_6h_samples.tar.gz |


**Example Concatenation Commands:**

These examples assume you have extracted the archives or run Guppy (producing uncompressed FASTQ). If your FASTQ files are gzipped (e.g., `.fastq.gz`), use `zcat` instead of `cat` and adjust filenames accordingly.

```bash
# Example for NaturalEnv_3h_Rep1 (after extracting natural_environment_3h_samples.tar.gz)
cat SQK-RBK114-24_naturalenvironment_3h_1_barcode01.fastq \
    SQK-RBK114-24_naturalenvironment_3h_1_barcode02.fastq \
    > NaturalEnv_3h_Rep1_combined.fastq

# Example for NaturalEnv_6h_Rep1 (after extracting natural_environment_6h_samples.tar.gz)
cat SQK-RBK114-24_naturalenvironment_6h_1_barcode11.fastq \
    SQK-RBK114-24_naturalenvironment_6h_1_barcode12.fastq \
    > NaturalEnv_6h_Rep1_combined.fastq

# Example using files from Guppy output (assuming output in 'guppy_output' directory)
cat guppy_output/barcode01/SQK-RBK114-24_naturalenvironment_3h_1_barcode01.fastq \
    guppy_output/barcode02/SQK-RBK114-24_naturalenvironment_3h_1_barcode02.fastq \
    > NaturalEnv_3h_Rep1_combined.fastq
