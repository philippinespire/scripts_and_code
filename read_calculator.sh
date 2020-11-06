#!/bin/bash -l

#SBATCH --job-name=read_calculator
#SBATCH -o read_calculator-%j.out
#SBATCH --time=144:00:00
#SBATCH -p main
#SBATCH -c 4
#SBATCH --mail-user=<youremail>
#SBATCH --mail-type=begin
#SBATCH --mail-type=END

# set above slurm options according to your system

###############################################
###            read_calculator.sh           ###
### counts reads in parallel from fq files  ###
###############################################

# read_caltulator.sh counts the number of reads in compressed or uncompressed FQ files.
# It is set up to process compressed FQ files by default, i.e. filetype=fq.gz and greptype=zgrep
# If you have uncompressed FQ files instead, set filetype=fq and greptype=grep  
# Adjust the number of cores/nthreads according to your system

# User Variables
filetype=fq.gz
greptype=zgrep
nthreads=20

# Load parallel and record variables 
enable_lmod
module load parallel
parallel --record-env

# Code
cat <(echo file,reads) <(\
      ls *$filetype | parallel --no-notice -kj$nthreads "echo -n {}',' && $greptype '^@' {} | wc -l") > numReads_$filetype.csv
