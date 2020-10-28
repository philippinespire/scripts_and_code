# scripts_and_code

Bash and R scripts to analyze Next-Gen sequence data. 

Contact: Eric Garcia, e1garcia@odu.edu

---

### Multi_FASTQC.sh

Simple script to create sequence quality reports ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/)) in parallel with a single command on a remote HPC

Usage - Place script in same directory as files to be processed and execute specifying the extension of files between quotations. For example:
```
sbatch Multi_FASTQC.sh "fq.gz"
```

`Multi_FASTQC.sh` has been tested in "fq", "fq.gz" and "bam" files.

Open script for more details

---

### base_calculator.sh

This script counts the number of each base in DNA fragments from single-digest Sbf1 RADseq, pair-end fq.gz files according to the window size and region determined by user. 

Usage:
1.- Set slurm options according to your system
2.- Set "User Variables" 
3.- Check that the ls statement in line 69 will list your input files
4.- Execute in command with: sbatch <script name> <"readDir">
```
sbatch base_calculator.sh "F"
```

Open script for more details


### base_proportions.R

R script to calculate, and plot, base pair proportions and mean base pair proportion of DNA fragments position by position

Uses the output of "base_calculator.sh", or tsv files with base pair counts from single digest RADseq, paired-End sequencing data, as input.

Open script for more details

---
