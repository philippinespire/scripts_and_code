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

---

### base_proportions.R

R script to calculate, and plot, base pair proportions and mean base pair proportion of DNA fragments position by position

Uses the output of "base_calculator.sh", or tsv files with base pair counts from single digest RADseq, paired-End sequencing data, as input

Open script for more details

---

### read_calculator.sh

read_caltulator.sh counts the number of reads in compressed (default) or uncompressed FQ files (open script for details).

Usage:  
1.- Place script in the same directory with FQ files to be processed. 

Open script and

2.- Set slurm options according to your system

3.- Set "User Variables"

4.- Execute
```
sbatch read_calculator.sh
```

---

### concat_fqFiles_diffLanes.sh

Simple script to concatenate files of the same individuals but from multiple sequencing lanes

Usage: Open the script can enter your files' read and lane info 

---

### subsetting_VCF_files.dat

User friendly list of steps and code to succesfully subset VCF files while maintaining functionality

---

### ssh_config_stay_connected.txt

Use this code to maintain a stable ssh connection if your sessions are becoming idle and/or terminated very quickly, after 1-2 min of inactivity.

---
