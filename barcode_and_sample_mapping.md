# Sample Information and Barcode Mapping (ENA: PRJEB76446)

This project analyzed air samples collected from various locations using Oxford Nanopore Technologies (ONT) sequencing. Samples were prepared and barcoded using the ONT Rapid Barcoding Kit (likely SQK-RBK114-24 based on filenames), allowing for multiplexed sequencing on a single flow cell. Blank samples were also collected and processed alongside the environmental samples.

The raw sequencing data for ENA submission is provided as TAR archives containing gzipped FASTQ files (`.fastq.tar.gz`). The data is presented below in two sections: the main environmental samples and the associated blank samples.

## Urban Study Samples

These are the main environmental samples collected for the study. The original filenames (prior to archival for ENA) encode information about the sample origin:

* **Location**: e.g., `citycenter`, `greenbelt`, `outerbelt`, `residentialarea`, `urbanbeach`
* **Replicate Number**: The first number after the location (e.g., `_1_`, `_2_`, `_3_`) likely indicates the biological or site replicate visit.
* **Sub-sample/Technical Replicate**: The second number (e.g., `_1`, `_2`) might indicate a technical replicate or fraction, denoted here as Sub1 or Sub2.
* **Barcode**: The assigned ONT barcode (e.g., `barcode01`, `barcode19`).

Due to file size considerations, each sample's FASTQ data was archived individually into its own `.fastq.tar.gz` file for submission.

| Sample name                | Location         | Sampling Date | File Name                                                   | Barcode   |
| :------------------------- | :--------------- | :------------ | :---------------------------------------------------------- | :-------- |
| CityCenter_Rep1_Sub1       | City center      | 19/10/2023    | `SQK-RBK114-24_citycenter_1_1_barcode19.fastq.tar.gz`       | barcode19 |
| CityCenter_Rep1_Sub2       | City center      | 19/10/2023    | `SQK-RBK114-24_citycenter_1_2_barcode20.fastq.tar.gz`       | barcode20 |
| CityCenter_Rep2_Sub1       | City center      | 31/10/2023    | `SQK-RBK114-24_citycenter_2_1_barcode21.fastq.tar.gz`       | barcode21 |
| CityCenter_Rep2_Sub2       | City center      | 31/10/2023    | `SQK-RBK114-24_citycenter_2_2_barcode22.fastq.tar.gz`       | barcode22 |
| CityCenter_Rep3_Sub1       | City center      | 03/11/2023    | `SQK-RBK114-24_citycenter_3_1_barcode23.fastq.tar.gz`       | barcode23 |
| CityCenter_Rep3_Sub2       | City center      | 03/11/2023    | `SQK-RBK114-24_citycenter_3_2_barcode24.fastq.tar.gz`       | barcode24 |
| Greenbelt_Rep1_Sub1        | Green belt       | 18/10/2023    | `SQK-RBK114-24_greenbelt_1_1_barcode01.fastq.tar.gz`        | barcode01 |
| Greenbelt_Rep1_Sub2        | Green belt       | 18/10/2023    | `SQK-RBK114-24_greenbelt_1_2_barcode02.fastq.tar.gz`        | barcode02 |
| Greenbelt_Rep2_Sub1        | Green belt       | 30/10/2023    | `SQK-RBK114-24_greenbelt_2_1_barcode03.fastq.tar.gz`        | barcode03 |
| Greenbelt_Rep2_Sub2        | Green belt       | 30/10/2023    | `SQK-RBK114-24_greenbelt_2_2_barcode04.fastq.tar.gz`        | barcode04 |
| Greenbelt_Rep3_Sub1        | Green belt       | 03/11/2023    | `SQK-RBK114-24_greenbelt_3_1_barcode05.fastq.tar.gz`        | barcode05 |
| Greenbelt_Rep3_Sub2        | Green belt       | 03/11/2023    | `SQK-RBK114-24_greenbelt_3_2_barcode06.fastq.tar.gz`        | barcode06 |
| Outerbelt_Rep1_Sub1        | Outer belt       | 17/10/2023    | `SQK-RBK114-24_outerbelt_1_1_barcode01.fastq.tar.gz`        | barcode01 |
| Outerbelt_Rep1_Sub2        | Outer belt       | 17/10/2023    | `SQK-RBK114-24_outerbelt_1_2_barcode02.fastq.tar.gz`        | barcode02 |
| Outerbelt_Rep2_Sub1        | Outer belt       | 20/10/2023    | `SQK-RBK114-24_outerbelt_2_1_barcode03.fastq.tar.gz`        | barcode03 |
| Outerbelt_Rep2_Sub2        | Outer belt       | 20/10/2023    | `SQK-RBK114-24_outerbelt_2_2_barcode04.fastq.tar.gz`        | barcode04 |
| Outerbelt_Rep3_Sub1        | Outer belt       | 02/11/2023    | `SQK-RBK114-24_outerbelt_3_1_barcode05.fastq.tar.gz`        | barcode05 |
| Outerbelt_Rep3_Sub2        | Outer belt       | 02/11/2023    | `SQK-RBK114-24_outerbelt_3_2_barcode06.fastq.tar.gz`        | barcode06 |
| ResidentialArea_Rep1_Sub1  | Residential area | 17/10/2023    | `SQK-RBK114-24_residentialarea_1_1_barcode13.fastq.tar.gz`  | barcode13 |
| ResidentialArea_Rep1_Sub2  | Residential area | 17/10/2023    | `SQK-RBK114-24_residentialarea_1_2_barcode14.fastq.tar.gz`  | barcode14 |
| ResidentialArea_Rep2_Sub1  | Residential area | 19/10/2023    | `SQK-RBK114-24_residentialarea_2_1_barcode15.fastq.tar.gz`  | barcode15 |
| ResidentialArea_Rep2_Sub2  | Residential area | 19/10/2023    | `SQK-RBK114-24_residentialarea_2_2_barcode16.fastq.tar.gz`  | barcode16 |
| ResidentialArea_Rep3_Sub1  | Residential area | 30/10/2023    | `SQK-RBK114-24_residentialarea_3_1_barcode17.fastq.tar.gz`  | barcode17 |
| ResidentialArea_Rep3_Sub2  | Residential area | 30/10/2023    | `SQK-RBK114-24_residentialarea_3_2_barcode18.fastq.tar.gz`  | barcode18 |
| UrbanBeach_Rep1_Sub1       | Urban beach      | 18/10/2023    | `SQK-RBK114-24_urbanbeach_1_1_barcode07.fastq.tar.gz`       | barcode07 |
| UrbanBeach_Rep1_Sub2       | Urban beach      | 18/10/2023    | `SQK-RBK114-24_urbanbeach_1_2_barcode08.fastq.tar.gz`       | barcode08 |
| UrbanBeach_Rep2_Sub1       | Urban beach      | 20/10/2023    | `SQK-RBK114-24_urbanbeach_2_1_barcode09.fastq.tar.gz`       | barcode09 |
| UrbanBeach_Rep2_Sub2       | Urban beach      | 20/10/2023    | `SQK-RBK114-24_urbanbeach_2_2_barcode11.fastq.tar.gz`       | barcode11 |
| UrbanBeach_Rep3_Sub1       | Urban beach      | 31/10/2023    | `SQK-RBK114-24_urbanbeach_3_1_barcode12.fastq.tar.gz`       | barcode12 |

## Blank Samples

These blank samples correspond to the different locations and replicate collection times. The sampling date listed is the date associated with the corresponding replicate number for that location (e.g., Blank Rep 1 date matches Sample Rep 1 date for that location).

All blank sample files are distributed together within a single compressed TAR archive named `blanks_urbanstudy_fastq.tar.gz`. Users downloading this archive will need to unpack it (`tar -xzvf blanks_urbanstudy_fastq.tar.gz`) to access the individual blank sample files listed below.

| Sample name                  | Location             | Sampling Date | File Name (within archive)                                | Barcode   |
| :--------------------------- | :------------------- | :------------ | :-------------------------------------------------------- | :-------- |
| CityCenter_Blank_Rep1        | City center Blank    | 19/10/2023    | `2024 SQK-RBK114-24_citycenter_blank_1_barcode04.fastq`     | barcode04 |
| CityCenter_Blank_Rep2        | City center Blank    | 31/10/2023    | `2024 SQK-RBK114-24_citycenter_blank_2_barcode05.fastq`     | barcode05 |
| CityCenter_Blank_Rep3        | City center Blank    | 03/11/2023    | `2024 SQK-RBK114-24_citycenter_blank_3_barcode06.fastq`     | barcode06 |
| Greenbelt_Blank_Rep1         | Green belt Blank     | 18/10/2023    | `2024 SQK-RBK114-24_greenbelt_blank_1_barcode07.fastq`      | barcode07 |
| Greenbelt_Blank_Rep2         | Green belt Blank     | 30/10/2023    | `2024 SQK-RBK114-24_greenbelt_blank_2_barcode08.fastq`      | barcode08 |
| Greenbelt_Blank_Rep3         | Green belt Blank     | 03/11/2023    | `2024 SQK-RBK114-24_greenbelt_blank_3_barcode09.fastq`      | barcode09 |
| Outerbelt_Blank_Rep1         | Outer belt Blank     | 17/10/2023    | `2024 SQK-RBK114-24_outerbelt_blank_1_barcode13.fastq`      | barcode13 |
| Outerbelt_Blank_Rep2         | Outer belt Blank     | 20/10/2023    | `2024 SQK-RBK114-24_outerbelt_blank_2_barcode14.fastq`      | barcode14 |
| Outerbelt_Blank_Rep3         | Outer belt Blank     | 02/11/2023    | `2024 SQK-RBK114-24_outerbelt_blank_3_barcode15.fastq`      | barcode15 |
| ResidentialArea_Blank_Rep1   | Residential Blank    | 17/10/2023    | `2024 SQK-RBK114-24_residentialarea_blank_1_barcode01.fastq`| barcode01 |
| ResidentialArea_Blank_Rep2   | Residential Blank    | 19/10/2023    | `2024 SQK-RBK114-24_residentialarea_blank_2_barcode02.fastq`| barcode02 |
| ResidentialArea_Blank_Rep3   | Residential Blank    | 30/10/2023    | `2024 SQK-RBK114-24_residentialarea_blank_3_barcode03.fastq`| barcode03 |
| UrbanBeach_Blank_Rep1        | Urban beach Blank    | 18/10/2023    | `2024 SQK-RBK114-24_urbanbeach_blank_1_barcode10.fastq`     | barcode10 |
| UrbanBeach_Blank_Rep2        | Urban beach Blank    | 20/10/2023    | `2024 SQK-RBK114-24_urbanbeach_blank_2_barcode11.fastq`     | barcode11 |
| UrbanBeach_Blank_Rep3        | Urban beach Blank    | 31/10/2023    | `2024 SQK-RBK114-24_urbanbeach_blank_3_barcode12.fastq`     | barcode12 |
