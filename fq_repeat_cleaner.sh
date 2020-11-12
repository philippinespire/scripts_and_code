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
###   	      fq_repeat_cleaner.sh	        ###
###        concatenates FQ files and 		###
### removes sequences with repeated motifs  ###
###############################################

### Requirements: parallel
### Input files: fq or fq.gz files from single-digest Sbf1 RADseq, pair-end sequencing data.


# DETAILS:
# fq_repeat_cleaner.sh removes sequences with repeated motifs at the beginning of the read (likely duplicated reads).
# This code was initially created to explore sequence patterns in FQ files and to search for signs of DNA deamination.
# All input FQ files are first concatenated into temporary single file where the scrips identifies all unique motifs from the first base to the specify maximum motif lenght and creates files listing these motifs and their frequencies
# A summary of motif frequencies is then printed on the terminal to help the user to decide on the motif length to use in sequence removal (typically, a short motif with high repetion among sequences). We used 10bp motifs in the past. 
# Then, the script prompts  the user to enter the desired motif length and procedes to remove all reads with common motifs of the determined size 

# This script expects all forward reads to begin with the Sbf1 restriction site motif (TGCAGG) and all reverse reads not to begin with it.
# The code will automatically remove any exceptions to the last sentence as these will shift the position of motifs and might represent sequencing mistakes. 
# The Sbf1 restriction site motif (TGCAGG) is temporarly removed from forwards reads so that the reported "position 1" is at the beginning all reads, forward and reverse.
# When choosing a motif cutoff, keep in mind that the reads with repeated large motifs will also be counted by all smaller motifs searches.
# For example, if your search for 20bp repeated motifs returns n-reads, the same reads will also be counted if you search for 10bp motifs since the first 10bp of the 20bp motifs are also repeated.
# In previous analysis searching for signs of DNA deamination, we have previously use 10bp as the threshold to remove reads with repeated motifs

### USAGE
# Set file and read info as well as the maximum length (bp) of motifs to be counted in the user variables below
# Save and execute the script
# The script will report frequencies of all motifs within the size range from "position 1", or the first bp, to the specified maximum length (from beginning of reads only) at the terminal
# and the user enters the desire motif length


# User variables options include:
# FILE_DIRECTION ("forward" or "reverse") 
# DIRECTION_SUFFIX ("F","R","R1","R2", etc)  
# FILE_EXTENSION ("fq" or "fq.gz" for uncompressed and compressed files, respectively) 
# MAX_motif_length (maximum motif size (in bp) to search for repeats). Enter only the number as a digit
# THREADS (number of threads according to your system)


# Set user variables here:
FILE_DIRECTION=forward
DIRECTION_SUFFIX=F
FILE_EXTENSION=fq.gz
MAX_MOTIF_length=20
THREADS=40
OUTPUT_BASENAME=Ssp-AAtu

### Script

### Calculating motif frequencies
echo "Calculating motif frequencies..."

# Concatenate all input files (uncompressed or compressed) and transpose rows to columns into a single temporary file where to search motif frequencies at once
zcat -f $(ls *$DIRECTION_SUFFIX.$FILE_EXTENSION) | paste - - - - > INPUT.temp &&

# Removed reads with the restriction site motif from forward files and those without it from reverse files 
if [[ $FILE_DIRECTION == 'forward' ]]
then
  awk '{print $2"\t"$1"\t"$3"\t"$4}' INPUT.temp | grep '^TGCAGG' | cut -c 7- > $DIRECTION_SUFFIX.INPUT.clean && rm INPUT.temp 
elif [[ $FILE_DIRECTION == 'reverse' ]]
then
  awk '{print $2"\t"$1"\t"$3"\t"$4}' INPUT.temp | grep -v '^TGCAGG' > $DIRECTION_SUFFIX.INPUT.clean && rm INPUT.temp 
else
  echo "Please provide file direction as 'forward' or 'reverse'"
fi && 

# Copy the bash variables to parallel
parallel --record-env

# Create files listing the repeated motifs and their frequencies. Files are created according to the motif size, from 1bp motif to the specify "MAX_MOTIF_length" 
seq 1 $MAX_MOTIF_length | \
	 parallel --no-notice -k "cat $DIRECTION_SUFFIX.INPUT.clean | cut -c 1-{} | sort | uniq -c | tr -s ' ' '\t' | cut -f2-3 | sort -nk1 | sed 's/\t/ /' | grep -v '^1 ' | sed 's/ /\t/' \
		> {}-bp.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads"
seq 1 $MAX_MOTIF_length | \
	parallel --no-notice -k "cat $DIRECTION_SUFFIX.INPUT.clean | cut -c 1-{} | sort | uniq -c | tr -s ' ' '\t' | cut -f2-3 | sort -nk1 | sed 's/\t/ /' | grep '^1 ' | sed 's/ /\t/' | cut -f2 \
		> {}-bp.$DIRECTION_SUFFIX.Motifs_never_repeated"

# Create summary for repeated motifs
echo -e "motif_file\tnReads_with_motifs\tmofit_length(bp)" > header.motif 
ls *bp.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads | sort -g | parallel --no-notice -k "echo -n {}',' && cut -f1 {} | paste -sd+ | bc" | sed 's/,/\t/g' > body.motif1  &&
ls *bp.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads | sort -g | sed 's/-.*//' > body.motif2 &&
cat header.motif <(\
	paste body.motif1 body.motif2) | \
 awk '{print $3"\t"$2"\t"$1}' > summary.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads.tsv

# Remove extra files
rm header.motif body.motif1 body.motif2

# Report summary and ask for the desired motif length to be used in sequence removal
echo -e "\nSummary of motifs in multiple sequences:\n"
cat summary.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads.tsv &&
echo -e "\n\n" &&
read -p "What repeated motif length (from above) would you like to use in sequence removal? Please enter digit : `echo $'\n> '`" selected_motif 

### Starting read removal process
echo -e "\nStarting read removal process..."

if [[ $FILE_DIRECTION == 'forward' ]]
then
  adjusted_motif=$((selected_motif + 6)) 
elif [[ $FILE_DIRECTION == 'reverse' ]]
then
  adjusted_motif=$selected_motif 
else
  echo "Please provide file direction as 'forward' or 'reverse'"
fi 

# Create a temporary file with motifs of selected size as the first column
awk -v var=$selected_motif -vFS="" -vOFS="" '{$var=$var"\t"}1' $DIRECTION_SUFFIX.INPUT.clean > INPUT.temp2 &&

# Extract lines for which the first columns match in the file listing the motifs and the read fil
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1]' $selected_motif-bp.$DIRECTION_SUFFIX.Motifs_never_repeated INPUT.temp2 > out_noRepeats.temp &&

# Restore the restriction site
paste <(yes "TGCAGG" | head -n $(cat out_noRepeats.temp | wc -l)) out_noRepeats.temp > out_noRepeats.rest.temp &&

# Restore FQ format and write final clean file
awk '{print $4"\t"$1$2$3"\t"$5"\t"$6}' out_noRepeats.rest.temp | sed 's/\t/\n/g' > $OUTPUT_BASENAME.$DIRECTION_SUFFIX.$selected_motif-bp_noRepeats.fq

# Remove extra files
rm $DIRECTION_SUFFIX.INPUT.clean INPUT.temp2 out_noRepeats.temp out_noRepeats.rest.temp

# Send motif files to a subdirectory
mkdir 1to$MAX_MOTIF_length-bp.$DIRECTION_SUFFIX.Motifs && 
mv *bp.$DIRECTION_SUFFIX.Motifs_inmultiple_Reads *bp.$DIRECTION_SUFFIX.Motifs_never_repeated summary.$DIRECTION_SUFFIX*  1to$MAX_MOTIF_length-bp.$DIRECTION_SUFFIX.Motifs

echo "Read removal completed"


