```
## Assembly
-Flye

Flye is a de novo assembler for long-reads. Flye also produces a polished consensus sequence for the assembly which significantly reduces the error rate.

```shell
$ flye --meta --nano-hq <input.fastq> -o /output_dir
```
The above command uses the high-quality polished nanopore reads and assembles the genome.
