#!/bin/bash

#SBATCH --job-name=motif_calculator
#SBATCH -o motif_calculator-%j.out
#SBATCH --time=144:00:00
#SBATCH -p main
#SBATCH -c 40
#SBATCH --mail-user=<youremail>
#SBATCH --mail-type=begin
#SBATCH --mail-type=END

# Above sbatch and module loading statements are optimized for ODU's Wahab super computer. Modify these to fit your system if using sbatch.
# Contact: Eric Garcia, e1garcia@odu.edu

###############################################
###   	      motif_calculator.sh	        ###
### identifies and counts repeated motifs   ###
###############################################

### Requirements: parallel
### Input files: fq or fq.gz files from single-digest Sbf1 RADseq, pair-end sequencing data.
### Output files: files listing the frequency of repeated and never-repeated motifs of specified length. Summary file of repeated motifs by length and the number of reads where they occur. 


# DETAILS:
# motif_calculator.sh identifies motifs that are present at the beginning of multiple reads and counts their frequencies.
# This is useful to characterize sequence patterns and remove sequences with repeated motifs
# The script expects all forward reads to begin with the Sbf1 restriction site motif (TGCAGG) and all reverse reads not to begin with it.
# The code will automatically remove any exceptions to the last sentence as these will shift the position of motifs and might represent sequencing mistakes. 
# The Sbf1 restriction site motif (TGCAGG) is also automatically removed from forwards reads so that the reported "position 1" is at the beginning all reads, forward and reverse.
# motif_calculator can be used in combination with "fq_repeat_cleaner.sh" to remove reads with repeated motifs.
# Keep in mind that the reads with repeated large motifs will also be counted by all smaller motifs searches.
# For example, if your search for 20bp repeated motifs returns X-reads, the same reads will also be counted if you search for 10bp motifs since the first 10bp of the 20bp motifs are also repeated.
# In previous analysis searching for signs of DNA deamination, we have previously use 10bp as the threshold to remove reads with repeated motifs

### USAGE
# Users set file and read info as well as the maximum length (bp) of motifs to be counted
#and the calculator reports frequencies of all motifs within the size range from "position 1", or the first bp, to the specified maximum length (from beginning of reads only).

# User variables options include:
# FILE_DIRECTION ("forward" or "reverse") 
# DIRECTION_SUFFIX ("F","R","R1","R2", etc)  
# FILE_EXTENSION ("fq" or "fq.gz" for uncompressed and compressed files, respectively) 
# MAX_motif_length (digit)(maximum motif size (in bp) to search for repeats)
# THREADS (number of threads according to your system)


# Set user variables here:
FILE_DIRECTION=forward
DIRECTION_SUFFIX=F
FILE_EXTENSION=fq.gz
MAX_MOTIF_length=20
THREADS=40


### Script

# Concatenate all input files (uncompressed or compressed) and transpose rows to columns into a single temporary file where to search motif frequencies at once
zcat -f $(ls *$DIRECTION_SUFFIX.$FILE_EXTENSION) | paste - - - - > INPUT.temp &&

# Removed reads with the restriction site motif from forward files and those without it from reverse files 
if [[ $FILE_DIRECTION == 'forward' ]]
then
  awk '{print $2"\t"$1"\t"$3"\t"$4}' INPUT.temp | grep '^TGCAGG' | cut -c 7- > $DIRECTION_SUFFIX.INPUT.clean  && rm INPUT.temp 
elif [[ $FILE_DIRECTION == 'reverse' ]]
then
  awk '{print $2"\t"$1"\t"$3"\t"$4}' INPUT.temp | grep -v '^TGCAGG' > $DIRECTION_SUFFIX.INPUT.clean  && rm INPUT.temp 
else
  echo "Please provide file direction as 'forward' or 'reverse'"
fi && 

# Copy the bash variables to parallel
parallel --record-env

# Create files listing the repeated motifs and their frequencies. Files are created according to the motif size, from 1bp motif to the specify "MAX_MOTIF_length" 
seq 1 $MAX_MOTIF_length | \
	 parallel --no-notice -k "cat $DIRECTION_SUFFIX.INPUT.clean | cut -c 1-{} | sort | uniq -c | tr -s ' ' '\t' | cut -f2-3 | sort -nk1 | sed 's/\t/ /' | grep -v '^1 ' | sed 's/ /\t/' \
		> {}_bp.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads"
seq 1 $MAX_MOTIF_length | \
	parallel --no-notice -k "cat $DIRECTION_SUFFIX.INPUT.clean | cut -c 1-{} | sort | uniq -c | tr -s ' ' '\t' | cut -f2-3 | sort -nk1 | sed 's/\t/ /' | grep '^1 ' | sed 's/ /\t/' \
		> {}_bp.$DIRECTION_SUFFIX.Motifs_never_repeated"

# Create summary for repeated motifs
echo -e "motif_file\tnReads_with_motifs\tmofit_length(bp)" > header.motif 
ls *bp.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads | sort -g | parallel --no-notice -k "echo -n {}',' && cut -f1 {} | paste -sd+ | bc" | sed 's/,/\t/g' > body.motif1  &&
ls *bp.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads | sort -g | sed 's/_.*//' > body.motif2 &&
cat header.motif <(\
	paste body.motif1 body.motif2) | \
 awk '{print $3"\t"$2"\t"$1}' > summary.$DIRECTION_SUFFIX.repeated_Motifs.tsv

# Remove extra files
rm $DIRECTION_SUFFIX.INPUT.clean header.motif body.motif1 body.motif2

# Report summary and send files to subdirectory
cat summary.$DIRECTION_SUFFIX.repeated_Motifs.tsv
mkdir calculator.$DIRECTION_SUFFIX.Motifs && 
mv *bp.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads *bp.$DIRECTION_SUFFIX.Motifs_never_repeated summary.$DIRECTION_SUFFIX*  calculator.$DIRECTION_SUFFIX.Motifs


