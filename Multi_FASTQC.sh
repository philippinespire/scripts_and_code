#!/bin/bash -l

#SBATCH --job-name=Multi_fastqc
#SBATCH -o Multi_fastqc-%j.out
#SBATCH --time=144:00:00
#SBATCH -p main
#SBATCH -c 4
#SBATCH --mail-user=<your_email>
#SBATCH --mail-type=begin
#SBATCH --mail-type=END


############# Multi_FASTQC.sh ###################
## runs FASTQC and MultiQC reports in parallel ##
##   contact: Eric Garcia, e1garcia@odu.edu    ##
#################################################

## Requirements: parallel, fastqc, and multiqc in current session 
## To execute, place script in the same directory as the files to be processed and type in the command line:
# sbatch Multi_FASTQC.sh "<extension of files to be processed in quotations>"


#### Details

# Multi_FASTQC.sh is a simple sbatch script that runs FASTQC and MultiQC reports in parallel with a single command 
# Results will be directed to a newly created sub-directory called Multi_FASTQC 
# For FASQC options use <fasqc --help> or visit  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# For MultiQC options use <multiqc --help> or visit  https://multiqc.info/

## Script Usage:
# 1.- Set the above slurm settings (#SBATCH) according to your system 
# 2.- Load parallel, fastqc and multiqc according to your system. Example:
```
#enable_lmod
#module load parallel
#module load container_env multiqc
#module load container_env fastqc
```
# 3.- Execute the script
# in the command line, type "sbatch", the name of the script <Multi_FASTQC.sh>, and the suffix identifying the files to be analyzed in quotations. The last can be file extensions or any other shared file identifier at the end of the files' names
# example: <sbatch Multi_FASTQC.sh ".fq.gz">

# Alternately, the suffix can be replaced by any regex expression that correctly identifies the files to be processed.
# If such regex does not occur at the end of file names, you'll need to remove the wild card " * " in the first fastqc statement in line 50

# Multi_FASTQC.sh has been tested in "fq", "fq.gz" and "bam" files.


#### Script

#run fastqc in parallel 
ls *$1 | parallel "crun fastqc {}" &&

# run multiqc with specific report and subdirectory names
crun multiqc . -n multiqc_report_$1.html -o Multi_FASTQC &&

# move fastqc files to new subdirectory
ls *fastqc.html | parallel "mv {} Multi_FASTQC" &&
ls *fastqc.zip | parallel "mv {} Multi_FASTQC"

