# 1. Guppy Basecaller

```plaintext
Guppy Basecaller is a software tool to convert the raw electrical signal data from nanopore sequencing into DNA sequences. 
Here's an example of a commonly used command with Guppy Basecaller:

guppy_basecaller -i /input/directory -s /output/directory --config configuration.cfg

In this command:

- -i /input/directory specifies the input directory where the raw data files (.fast5 files) are located.
- -s /output/directory specifies the output directory where the basecalled reads will be written.
- --config configuration.cfg is used to specify a configuration file that contains the basecalling model and other settings.

This command will basecall the raw data in the specified input directory using the settings from the configuration file, and write the output to the specified output directory.
```


# 2. Porechop

```plaintext
Porechop is a tool developed for Oxford Nanopore sequencing data. It is used for finding and removing adapters from Oxford Nanopore reads. Adapters on the ends of reads are trimmed off, and when a read has an adapter in its middle, it is treated as chimeric and chopped into separate reads.

Here's an example of a commonly used command with Porechop:

porechop -i input.fastq -o output.fastq

In this command:

- -i input.fastq specifies the input file (in FASTQ format) that you want to trim adapters from.
- -o output.fastq specifies the output file where the trimmed reads will be written.

This command reads the input FASTQ file, trims adapters from the reads, and writes the resulting trimmed reads to the output file.
```


# NanoFilt

```plaintext
NanoFilt is a simple tool to filter Oxford Nanopore sequencing data. It reads in a FASTQ file (or stdin), filters reads based on a minimum quality and/or a minimum length, and writes out the filtered reads to stdout.

Here's an example of a commonly used command with NanoFilt:

gunzip -c input.fastq.gz | NanoFilt -q 10 -l 500 | gzip > output.fastq.gz

In this command:

- gunzip -c input.fastq.gz is used to decompress the input FASTQ file.
- NanoFilt -q 9 -l 500 filters reads based on a minimum quality of 9 and a minimum length of 100.
- gzip > output.fastq.gz compresses the filtered reads and writes them to the output file.

This command decompresses the input FASTQ file, filters the reads based on the specified minimum quality and length, and writes the filtered reads to the output file in compressed format.
```


# Flye

```plaintext
Flye is a de novo assembler for single-molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies. It is designed for a wide range of datasets, from small bacterial projects to large mammalian genomes. The Flye assembler provides accurate, fast, and scalable solutions to assembly problems.

In the context of metagenomic projects, Flye can be run in 'meta' mode. Here's an example of a commonly used command with Flye in 'meta' mode:

flye --meta --nano-raw input.fastq --out-dir output_directory 

In this command:

- --meta indicates that Flye should be run in 'meta' mode for metagenomic data.
- --nano-raw input.fastq specifies the input file containing the reads. Depending on the source of the sequencing data change to nano-hq or nano-corr
- --out-dir output_directory specifies the directory where the assembly result will be written.

This command runs Flye in 'meta' mode on the input reads, assembles the reads, and writes the assembly result to the specified output directory.
```
