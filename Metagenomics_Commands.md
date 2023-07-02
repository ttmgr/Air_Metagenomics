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
