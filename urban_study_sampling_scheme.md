# Raw Data Access, Basecalling, and Sample Mapping (ENA: PRJEB76446)

This project analyzed air samples collected from various locations using Oxford Nanopore Technologies (ONT) sequencing. Two forms of data are available:

1.  **Raw Signal Data:** Unprocessed `.pod5` (pod5.tar.gz) files from the sequencing runs, grouped by run/location.
2.  **Processed Sequence Data:** Basecalled and demultiplexed FASTQ files (`.fastq.tar.gz`), submitted to the European Nucleotide Archive (ENA) under accession PRJEB76446.

This document first describes how to access and process the raw `.pod5` files to generate FASTQ data. Subsequently, it provides the mapping between the resulting FASTQ files (or the pre-generated ENA files) and the corresponding sample information.

## Raw POD5 Data and Basecalling Instructions

The raw sequencing signal data is available as `.pod5` files. These files are typically large and may be split into multiple parts, grouped by the sequencing run. Based on the ENA submission structure (see image `grafik.png`), runs seem to be associated with locations like `outer_belt`, `gracia_eixample`, and `greenbelt`.

**Important:**
* A single sequencing run often contains multiplexed samples from different replicates or even locations. For example, the `greenbelt` run files also contain the data for `urbanbeach` samples.
* You need to download **all parts** associated with a specific sequencing run (e.g., `outer_belt_part*.tar.gz`, `greenbelt_1_part*.tar.gz`, `greenbelt_2_part*.tar.gz`, etc., assuming these archives contain the POD5 files or the POD5 files follow a similar naming pattern) to reconstruct the complete dataset for that run.

**Processing Steps using Dorado:**

1.  **Download & Prepare:** Download all parts for the sequencing run of interest. If they are archived (e.g., in `.tar.gz` files), extract them. Ensure all `.pod5` files for the run are in a single directory.
2.  **Basecalling:** Use ONT's Dorado basecaller to convert the raw signal data in the `.pod5` files into sequence reads (FASTQ format). Use the `--emit-fastq` flag to output FASTQ directly and `--no-trim` to skip adapter trimming at this stage. You will need to select an appropriate basecalling model matching the sequencing kit and flowcell version used (e.g., a model compatible with Kit 114 chemistry like `dna_r10.4.1_e8.2_400bps_sup@v4.2.0`).
    ```bash
    # Example Dorado basecalling command (adjust model and paths)
    # Outputs a single, un-trimmed FASTQ file containing all reads
    dorado basecaller --emit-fastq --no-trim <path_to_model> <path_to_pod5_directory> > combined_reads.fastq
    ```
3.  **Demultiplexing:** Use Dorado to demultiplex the basecalled reads from the combined FASTQ file based on the barcodes used. Samples were prepared using the ONT Rapid Barcoding Kit (SQK-RBK114-24). Use the `--emit-fastq` flag again. Adapter trimming is typically performed by default during demultiplexing when `--no-trim` is *not* specified.
    ```bash
    # Example Dorado demultiplexing command (adjust kit name and paths)
    # Takes the combined FASTQ as input, performs trimming, and outputs per-barcode FASTQ files
    dorado demux --kit-name SQK-RBK114-24 --emit-fastq --output-dir demultiplexed_fastq combined_reads.fastq
    ```
4.  **Output:** This process will generate trimmed, demultiplexed FASTQ files (likely named `barcodeXX.fastq` or similar) in the specified output directory (`demultiplexed_fastq` in the example).

Once you have generated these demultiplexed FASTQ files, you can use the mapping tables below to identify which sample corresponds to each barcode file.

---

## Sample Information and Barcode Mapping (ENA: PRJEB76446)

This section details the mapping between barcodes and the original environmental or blank samples. This mapping applies to both the FASTQ files generated from the raw POD5 data (as described above) and the pre-generated FASTQ files available directly from ENA.

Samples were prepared and barcoded using the ONT Rapid Barcoding Kit (SQK-RBK114-24), allowing for multiplexed sequencing. Blank samples were also collected and processed.

The pre-generated sequencing data submitted to ENA is provided as TAR archives containing gzipped FASTQ files (`.fastq.tar.gz`). The data is presented below in two sections: the main environmental samples and the associated blank samples.

### Urban Study Samples

These are the main environmental samples collected for the study. The original filenames (prior to archival for ENA) encode information about the sample origin:

* **Location**: e.g., `citycenter`, `greenbelt`, `outerbelt`, `residentialarea`, `urbanbeach`
* **Replicate Number**: The first number after the location (e.g., `_1_`, `_2_`, `_3_`) likely indicates the biological or site replicate visit.
* **Sub-sample/Technical Replicate**: The second number (e.g., `_1`, `_2`) might indicate a technical replicate or fraction, denoted here as Sub1 or Sub2.
* **Barcode**: The assigned ONT barcode (e.g., `barcode01`, `barcode19`).

For the ENA submission, each sample's FASTQ data was archived individually into its own `.fastq.tar.gz` file due to file size considerations.

| Sample name                 | Location         | Sampling Date | ENA File Name                                               | Barcode   |
| :-------------------------- | :--------------- | :------------ | :---------------------------------------------------------- | :-------- |
| CityCenter_Rep1_Sub1        | City center      | 19/10/2023    | `SQK-RBK114-24_citycenter_1_1_barcode19.fastq.tar.gz`       | barcode19 |
| CityCenter_Rep1_Sub2        | City center      | 19/10/2023    | `SQK-RBK114-24_citycenter_1_2_barcode20.fastq.tar.gz`       | barcode20 |
| CityCenter_Rep2_Sub1        | City center      | 31/10/2023    | `SQK-RBK114-24_citycenter_2_1_barcode21.fastq.tar.gz`       | barcode21 |
| CityCenter_Rep2_Sub2        | City center      | 31/10/2023    | `SQK-RBK114-24_citycenter_2_2_barcode22.fastq.tar.gz`       | barcode22 |
| CityCenter_Rep3_Sub1        | City center      | 03/11/2023    | `SQK-RBK114-24_citycenter_3_1_barcode23.fastq.tar.gz`       | barcode23 |
| CityCenter_Rep3_Sub2        | City center      | 03/11/2023    | `SQK-RBK114-24_citycenter_3_2_barcode24.fastq.tar.gz`       | barcode24 |
| Greenbelt_Rep1_Sub1         | Green belt       | 18/10/2023    | `SQK-RBK114-24_greenbelt_1_1_barcode01.fastq.tar.gz`        | barcode01 |
| Greenbelt_Rep1_Sub2         | Green belt       | 18/10/2023    | `SQK-RBK114-24_greenbelt_1_2_barcode02.fastq.tar.gz`        | barcode02 |
| Greenbelt_Rep2_Sub1         | Green belt       | 30/10/2023    | `SQK-RBK114-24_greenbelt_2_1_barcode03.fastq.tar.gz`        | barcode03 |
| Greenbelt_Rep2_Sub2         | Green belt       | 30/10/2023    | `SQK-RBK114-24_greenbelt_2_2_barcode04.fastq.tar.gz`        | barcode04 |
| Greenbelt_Rep3_Sub1         | Green belt       | 03/11/2023    | `SQK-RBK114-24_greenbelt_3_1_barcode05.fastq.tar.gz`        | barcode05 |
| Greenbelt_Rep3_Sub2         | Green belt       | 03/11/2023    | `SQK-RBK114-24_greenbelt_3_2_barcode06.fastq.tar.gz`        | barcode06 |
| Outerbelt_Rep1_Sub1         | Outer belt       | 17/10/2023    | `SQK-RBK114-24_outerbelt_1_1_barcode01.fastq.tar.gz`        | barcode01 |
| Outerbelt_Rep1_Sub2         | Outer belt       | 17/10/2023    | `SQK-RBK114-24_outerbelt_1_2_barcode02.fastq.tar.gz`        | barcode02 |
| Outerbelt_Rep2_Sub1         | Outer belt       | 20/10/2023    | `SQK-RBK114-24_outerbelt_2_1_barcode03.fastq.tar.gz`        | barcode03 |
| Outerbelt_Rep2_Sub2         | Outer belt       | 20/10/2023    | `SQK-RBK114-24_outerbelt_2_2_barcode04.fastq.tar.gz`        | barcode04 |
| Outerbelt_Rep3_Sub1         | Outer belt       | 02/11/2023    | `SQK-RBK114-24_outerbelt_3_1_barcode05.fastq.tar.gz`        | barcode05 |
| Outerbelt_Rep3_Sub2         | Outer belt       | 02/11/2023    | `SQK-RBK114-24_outerbelt_3_2_barcode06.fastq.tar.gz`        | barcode06 |
| ResidentialArea_Rep1_Sub1   | Residential area | 17/10/2023    | `SQK-RBK114-24_residentialarea_1_1_barcode13.fastq.tar.gz`  | barcode13 |
| ResidentialArea_Rep1_Sub2   | Residential area | 17/10/2023    | `SQK-RBK114-24_residentialarea_1_2_barcode14.fastq.tar.gz`  | barcode14 |
| ResidentialArea_Rep2_Sub1   | Residential area | 19/10/2023    | `SQK-RBK114-24_residentialarea_2_1_barcode15.fastq.tar.gz`  | barcode15 |
| ResidentialArea_Rep2_Sub2   | Residential area | 19/10/2023    | `SQK-RBK114-24_residentialarea_2_2_barcode16.fastq.tar.gz`  | barcode16 |
| ResidentialArea_Rep3_Sub1   | Residential area | 30/10/2023    | `SQK-RBK114-24_residentialarea_3_1_barcode17.fastq.tar.gz`  | barcode17 |
| ResidentialArea_Rep3_Sub2   | Residential area | 30/10/2023    | `SQK-RBK114-24_residentialarea_3_2_barcode18.fastq.tar.gz`  | barcode18 |
| UrbanBeach_Rep1_Sub1        | Urban beach      | 18/10/2023    | `SQK-RBK114-24_urbanbeach_1_1_barcode07.fastq.tar.gz`       | barcode07 |
| UrbanBeach_Rep1_Sub2        | Urban beach      | 18/10/2023    | `SQK-RBK114-24_urbanbeach_1_2_barcode08.fastq.tar.gz`       | barcode08 |
| UrbanBeach_Rep2_Sub1        | Urban beach      | 20/10/2023    | `SQK-RBK114-24_urbanbeach_2_1_barcode09.fastq.tar.gz`       | barcode09 |
| UrbanBeach_Rep2_Sub2        | Urban beach      | 20/10/2023    | `SQK-RBK114-24_urbanbeach_2_2_barcode11.fastq.tar.gz`       | barcode11 |
| UrbanBeach_Rep3_Sub1        | Urban beach      | 31/10/2023    | `SQK-RBK114-24_urbanbeach_3_1_barcode12.fastq.tar.gz`       | barcode12 |

*(Note: There seems to be a missing UrbanBeach_Rep2_Sub2 file/barcode10 in the original table; barcode11 is listed instead. This inconsistency is preserved from the original text but might need verification)*.
*(Note 2: Barcodes 01-06 are used for both Greenbelt and Outerbelt samples according to the table. This is possible if they were run on separate flowcells or runs, but worth double-checking if they were run together)*.

### Blank Samples

These blank samples correspond to the different locations and replicate collection times. The sampling date listed is the date associated with the corresponding replicate number for that location.

For the ENA submission, all blank sample FASTQ files are distributed together within a single compressed TAR archive named `blanks_urbanstudy_fastq.tar.gz`. Users downloading this archive will need to unpack it (`tar -xzvf blanks_urbanstudy_fastq.tar.gz`) to access the individual blank sample FASTQ files listed below. If generating FASTQ from POD5, these blanks will be demultiplexed alongside the environmental samples using the specified barcodes.

| Sample name                 | Location            | Sampling Date | File Name (within ENA archive or from Demux)              | Barcode   |
| :-------------------------- | :------------------ | :------------ | :-------------------------------------------------------- | :-------- |
| CityCenter_Blank_Rep1       | City center Blank   | 19/10/2023    | `2024 SQK-RBK114-24_citycenter_blank_1_barcode04.fastq`     | barcode04 |
| CityCenter_Blank_Rep2       | City center Blank   | 31/10/2023    | `2024 SQK-RBK114-24_citycenter_blank_2_barcode05.fastq`     | barcode05 |
| CityCenter_Blank_Rep3       | City center Blank   | 03/11/2023    | `2024 SQK-RBK114-24_citycenter_blank_3_barcode06.fastq`     | barcode06 |
| Greenbelt_Blank_Rep1        | Green belt Blank    | 18/10/2023    | `2024 SQK-RBK114-24_greenbelt_blank_1_barcode07.fastq`      | barcode07 |
| Greenbelt_Blank_Rep2        | Green belt Blank    | 30/10/2023    | `2024 SQK-RBK114-24_greenbelt_blank_2_barcode08.fastq`      | barcode08 |
| Greenbelt_Blank_Rep3        | Green belt Blank    | 03/11/2023    | `2024 SQK-RBK114-24_greenbelt_blank_3_barcode09.fastq`      | barcode09 |
| Outerbelt_Blank_Rep1        | Outer belt Blank    | 17/10/2023    | `2024 SQK-RBK114-24_outerbelt_blank_1_barcode13.fastq`      | barcode13 |
| Outerbelt_Blank_Rep2        | Outer belt Blank    | 20/10/2023    | `2024 SQK-RBK114-24_outerbelt_blank_2_barcode14.fastq`      | barcode14 |
| Outerbelt_Blank_Rep3        | Outer belt Blank    | 02/11/2023    | `2024 SQK-RBK114-24_outerbelt_blank_3_barcode15.fastq`      | barcode15 |
| ResidentialArea_Blank_Rep1  | Residential Blank   | 17/10/2023    | `2024 SQK-RBK114-24_residentialarea_blank_1_barcode01.fastq`| barcode01 |
| ResidentialArea_Blank_Rep2  | Residential Blank   | 19/10/2023    | `2024 SQK-RBK114-24_residentialarea_blank_2_barcode02.fastq`| barcode02 |
| ResidentialArea_Blank_Rep3  | Residential Blank   | 30/10/2023    | `2024 SQK-RBK114-24_residentialarea_blank_3_barcode03.fastq`| barcode03 |
| UrbanBeach_Blank_Rep1       | Urban beach Blank   | 18/10/2023    | `2024 SQK-RBK114-24_urbanbeach_blank_1_barcode10.fastq`     | barcode10 |
| UrbanBeach_Blank_Rep2       | Urban beach Blank   | 20/10/2023    | `2024 SQK-RBK114-24_urbanbeach_blank_2_barcode11.fastq`     | barcode11 |
| UrbanBeach_Blank_Rep3       | Urban beach Blank   | 31/10/2023    | `2024 SQK-RBK114-24_urbanbeach_blank_3_barcode12.fastq`     | barcode12 |

*(Note: Barcodes 04-06 are used for CityCenter blanks, but also listed for Greenbelt_Rep2/3 and Outerbelt_Rep2/3 samples. Barcodes 07-09 are used for Greenbelt blanks but also listed for UrbanBeach_Rep1/2 samples. Barcodes 13-15 are used for Outerbelt blanks but also for ResidentialArea_Rep1/2 samples. Barcodes 01-03 are used for ResidentialArea blanks but also for Greenbelt_Rep1/2 and Outerbelt_Rep1/2 samples. Barcodes 10-12 are used for UrbanBeach blanks but also listed for UrbanBeach_Rep2/3 samples. This extensive reuse of barcodes across sample types and blanks within the same kit (assuming RBK114-24 with 24 barcodes) indicates potential issues or requires careful interpretation, possibly relying on separate sequencing runs to avoid barcode clashes. Please verify the run structure and barcode assignments.)*
