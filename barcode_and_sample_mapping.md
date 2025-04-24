## Sample Information and Barcode Mapping

This project analyzed air samples collected from various locations using Oxford Nanopore Technologies (ONT) sequencing. Samples were prepared and barcoded using the ONT Rapid Barcoding Kit (likely SQK-RBK114-24 based on filenames), allowing for multiplexed sequencing on a single flow cell.

The raw sequencing data is provided as FASTQ files, compressed using GZIP (`.fastq.gz`). The filenames encode information about the sample origin:

* **Location:** e.g., `citycenter`, `greenbelt`, `outerbelt`, `residentialarea`, `urbanbeach`
* **Replicate Number:** The first number after the location (e.g., `_1_`, `_2_`, `_3_`) likely indicates the biological or site replicate visit.
* **Sub-sample/Technical Replicate:** The second number (e.g., `_1`, `_2`) might indicate a technical replicate or fraction, denoted here as `Sub1` or `Sub2`.
* **Barcode:** The assigned ONT barcode (e.g., `barcode01`, `barcode19`).

The following table provides a mapping between the derived sample names, the sampling dates, the original (compressed) FASTQ filenames, and the assigned barcodes.

| Sample name                 | Sampling Date   | File Name                                                     | Barcode   |
| :-------------------------- | :-------------- | :------------------------------------------------------------ | :-------- |
| CityCenter_Rep1_Sub1        | 2023-10-19      | SQK-RBK114-24_citycenter_1_1_barcode19.fastq.gz             | barcode19 |
| CityCenter_Rep1_Sub2        | 2023-10-19      | SQK-RBK114-24_citycenter_1_2_barcode20.fastq.gz             | barcode20 |
| CityCenter_Rep2_Sub1        | 2023-10-31      | SQK-RBK114-24_citycenter_2_1_barcode21.fastq.gz             | barcode21 |
| CityCenter_Rep2_Sub2        | 2023-10-31      | SQK-RBK114-24_citycenter_2_2_barcode22.fastq.gz             | barcode22 |
| CityCenter_Rep3_Sub1        | 2023-11-03      | SQK-RBK114-24_citycenter_3_1_barcode23.fastq.gz             | barcode23 |
| CityCenter_Rep3_Sub2        | 2023-11-03      | SQK-RBK114-24_citycenter_3_2_barcode24.fastq.gz             | barcode24 |
| Greenbelt_Rep1_Sub1         | 2023-10-18      | SQK-RBK114-24_greenbelt_1_1_barcode01.fastq.gz              | barcode01 |
| Greenbelt_Rep1_Sub2         | 2023-10-18      | SQK-RBK114-24_greenbelt_1_2_barcode02.fastq.gz              | barcode02 |
| Greenbelt_Rep2_Sub1         | 2023-10-30      | SQK-RBK114-24_greenbelt_2_1_barcode03.fastq.gz              | barcode03 |
| Greenbelt_Rep2_Sub2         | 2023-10-30      | SQK-RBK114-24_greenbelt_2_2_barcode04.fastq.gz              | barcode04 |
| Greenbelt_Rep3_Sub1         | 2023-11-03      | SQK-RBK114-24_greenbelt_3_1_barcode05.fastq.gz              | barcode05 |
| Greenbelt_Rep3_Sub2         | 2023-11-03      | SQK-RBK114-24_greenbelt_3_2_barcode06.fastq.gz              | barcode06 |
| Outerbelt_Rep1_Sub1         | 2023-10-17      | SQK-RBK114-24_outerbelt_1_1_barcode01.fastq.gz              | barcode01 |
| Outerbelt_Rep1_Sub2         | 2023-10-17      | SQK-RBK114-24_outerbelt_1_2_barcode02.fastq.gz              | barcode02 |
| Outerbelt_Rep2_Sub1         | 2023-10-20      | SQK-RBK114-24_outerbelt_2_1_barcode03.fastq.gz              | barcode03 |
| Outerbelt_Rep2_Sub2         | 2023-10-20      | SQK-RBK114-24_outerbelt_2_2_barcode04.fastq.gz              | barcode04 |
| Outerbelt_Rep3_Sub1         | 2023-11-02      | SQK-RBK114-24_outerbelt_3_1_barcode05.fastq.gz              | barcode05 |
| Outerbelt_Rep3_Sub2         | 2023-11-02      | SQK-RBK114-24_outerbelt_3_2_barcode06.fastq.gz              | barcode06 |
| ResidentialArea_Rep1_Sub1 | 2023-10-17      | SQK-RBK114-24_residentialarea_1_1_barcode13.fastq.gz      | barcode13 |
| ResidentialArea_Rep1_Sub2 | 2023-10-17      | SQK-RBK114-24_residentialarea_1_2_barcode14.fastq.gz      | barcode14 |
| ResidentialArea_Rep2_Sub1 | 2023-10-19      | SQK-RBK114-24_residentialarea_2_1_barcode15.fastq.gz      | barcode15 |
| ResidentialArea_Rep2_Sub2 | 2023-10-19      | SQK-RBK114-24_residentialarea_2_2_barcode16.fastq.gz      | barcode16 |
| ResidentialArea_Rep3_Sub1 | 2023-10-30      | SQK-RBK114-24_residentialarea_3_1_barcode17.fastq.gz      | barcode17 |
| ResidentialArea_Rep3_Sub2 | 2023-10-30      | SQK-RBK114-24_residentialarea_3_2_barcode18.fastq.gz      | barcode18 |
| UrbanBeach_Rep1_Sub1      | 2023-10-18      | SQK-RBK114-24_urbanbeach_1_1_barcode07.fastq.gz           | barcode07 |
| UrbanBeach_Rep1_Sub2      | 2023-10-18      | SQK-RBK114-24_urbanbeach_1_2_barcode08.fastq.gz           | barcode08 |
| UrbanBeach_Rep2_Sub1      | 2023-10-20      | SQK-RBK114-24_urbanbeach_2_1_barcode09.fastq.gz           | barcode09 |
| UrbanBeach_Rep2_Sub2      | 2023-10-20      | SQK-RBK114-24_urbanbeach_2_2_barcode11.fastq.gz           | barcode11 |
| UrbanBeach_Rep3_Sub1      | 2023-10-31      | SQK-RBK114-24_urbanbeach_3_1_barcode12.fastq.gz           | barcode12 |
