#  mapDamage tutorial for Wahab, ODU

Read about the program, its output, and follow instructions for the setup in the [mapDamaga's main page](http://ginolhac.github.io/mapDamage/)

mapDamage and dependencies have already been installed in Wahab. 

Previous examples can be found in /home/e1garcia/PIRE_data/Ssp/mkBAM_0LeadTrim/deamination/mapDamage

## Start here

Download mapDamage from git into new desired directory 
```sh
git clone https://github.com/ginolhac/mapDamage.git
```

Get a working node if you haven't done this yet
```sh
salloc
```

Load the container with mapDamage
```sh
enable_lmod
module load container_env mapdamage2
```

mapDamage creates plots using R, which is already installed in the same container

Load required R libraries in your current session
```sh
crun R		# Access R inside the container

# Load required libraries
library("inline")
library("gam")
library("Rcpp")
library("ggplot2")
library("RcppGSL")
```

You can check that they were loaded by
```sh
sessionInfo()
````

Quit R and save the session
```sh
quit()
y
```

Run using `crun` while specifying a BAM (or SAM) file and the corresponding reference in fasta format using the `-i` and `-r` flags, respectively. Example:
```sh
crun mapDamage -i ../../PIRE2019-Ssp-A-Atu_010-Plate1Pool5Seq1-1E-L4.5.5-RAW -r ../../reference.5.5.fasta
```

This creates a directory with results named: results<name_of_bam_file>

A successful run delivers a message similar to this:
```sh
Started with the command: /opt/mapdamage2/bin/mapDamage -i ../../PIRE2019-Ssp-A-Atu_010-Plate1Pool5Seq1-1E-L4.5.5-RAW.bam -r ../../reference.5.5.fasta
WARNING: Alignment contains a large number of reference sequences (28525)!
  This may lead to excessive memory/disk usage.
  Consider using --merge-reference-sequences

	Reading from '../../PIRE2019-Ssp-A-Atu_010-Plate1Pool5Seq1-1E-L4.5.5-RAW.bam'
	Writing results to 'results_PIRE2019-Ssp-A-Atu_010-Plate1Pool5Seq1-1E-L4.5.5-RAW/'
pdf results_PIRE2019-Ssp-A-Atu_010-Plate1Pool5Seq1-1E-L4.5.5-RAW/Fragmisincorporation_plot.pdf generated
No length distributions are available, plotting length distribution only works for single-end reads
Warning: DNA damage levels are too low, the Bayesian computation should not be performed (0.001559 < 0.01)

Performing Bayesian estimates
Warning, To few substitutions to assess the nick frequency, using constant nick frequency instead
Starting grid search, starting from random values
Adjusting the proposal variance iteration  1  
Adjusting the proposal variance iteration  2  
Adjusting the proposal variance iteration  3  
Adjusting the proposal variance iteration  4  
Adjusting the proposal variance iteration  5  
Adjusting the proposal variance iteration  6  
Adjusting the proposal variance iteration  7  
Adjusting the proposal variance iteration  8  
Adjusting the proposal variance iteration  9  
Adjusting the proposal variance iteration  10  
Done burning, starting the iterations
Done with the iterations, finishing up
Writing and plotting to files
Warning messages:
1: package ‘inline’ was built under R version 4.0.2 
2: package ‘gam’ was built under R version 4.0.2 
Successful run
```

I then created a sbatch script, mapDamage_loop.sh, to run all albatross files in parallel:
```sh
#!/bin/bash -l

#SBATCH --job-name=mapDamage
#SBATCH -o mapDamage-%j.out
#SBATCH -p main
#SBATCH -c 40
#SBATCH --mail-user=youremail
#SBATCH --mail-type=begin
#SBATCH --mail-type=END

enable_lmod
module load container_env mapdamage2
module load paralle

ls ../../PIRE2019-Ssp-A*bam | parallel --no-notice -kj40 "crun mapDamage -i{} -r ../../reference.5.5.fasta"

# Running a few jobs with contemporary files for comparisons as well

crun mapDamage -i../../PIRE2019-Ssp-C-Gub_002-Plate1Pool3Seq1-2D-L4.5.5-RAW.bam -r ../../reference.5.5.fasta
crun mapDamage -i../../PIRE2019-Ssp-C-Gub_005-Plate1Pool12Seq1-3E-L4.5.5-RAW.bam -r ../../reference.5.5.fasta
```

**Don't forget you load R and required libraries before running a script like this**

*message*
```sh
crun mapDamage -i../../PIRE2019-Ssp-C-Gub_005-Plate1Pool12Seq1-3E-L4.5.5-RAW.bam -r ../../reference.5.5.fasta
Started with the command: /opt/mapdamage2/bin/mapDamage -i../../PIRE2019-Ssp-C-Gub_005-Plate1Pool12Seq1-3E-L4.5.5-RAW.bam -r ../../reference.5.5.fasta
WARNING: Alignment contains a large number of reference sequences (28525)!
  This may lead to excessive memory/disk usage.
  Consider using --merge-reference-sequences

	Reading from '../../PIRE2019-Ssp-C-Gub_005-Plate1Pool12Seq1-3E-L4.5.5-RAW.bam'
	Writing results to 'results_PIRE2019-Ssp-C-Gub_005-Plate1Pool12Seq1-3E-L4.5.5-RAW/'
pdf results_PIRE2019-Ssp-C-Gub_005-Plate1Pool12Seq1-3E-L4.5.5-RAW/Fragmisincorporation_plot.pdf generated
No length distributions are available, plotting length distribution only works for single-end reads
Warning: DNA damage levels are too low, the Bayesian computation should not be performed (0.005677 < 0.01)

Performing Bayesian estimates
Warning, To few substitutions to assess the nick frequency, using constant nick frequency instead
Starting grid search, starting from random values
Adjusting the proposal variance iteration  1  
Adjusting the proposal variance iteration  2  
Adjusting the proposal variance iteration  3  
Adjusting the proposal variance iteration  4  
Adjusting the proposal variance iteration  5  
Adjusting the proposal variance iteration  6  
Adjusting the proposal variance iteration  7  
Adjusting the proposal variance iteration  8  
Adjusting the proposal variance iteration  9  
Adjusting the proposal variance iteration  10  
Done burning, starting the iterations
Done with the iterations, finishing up
Writing and plotting to files
Warning messages:
1: package ‘inline’ was built under R version 4.0.2 
2: package ‘gam’ was built under R version 4.0.2 
Successful run

```

## Updated SLURM Script by CBIRD

name: runMAPDMG.sbatch

example: `sbatch runMAPDMG.sbatch "*NoWGA_???-A*RG.bam" ref.fasta`

```bash
#!/bin/bash -l

#SBATCH --job-name=mpdmg
#SBATCH -o mapDamage-%j.out
#SBATCH -p main
#SBATCH -c 32

enable_lmod
module load container_env mapdamage2
module load parallel

PATTERN=$1
REF=$2

echo FILES=$PATTERN
echo REF=$REF

ls $PATTERN | parallel --no-notice -j 32 "crun mapDamage -i {} -r $REF"
```
