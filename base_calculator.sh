#!/bin/bash

#SBATCH --job-name=base_calculator
#SBATCH -o base_calculator-%j.out
#SBATCH --time=144:00:00
#SBATCH -p main
#SBATCH -c 40
#SBATCH --mail-user=<youremail>
#SBATCH --mail-type=begin
#SBATCH --mail-type=END

# Above sbatch and module loading statements are optimized for ODU's Wahab super computer. Modify these to fit your system.
# Contact: Eric Garcia, e1garcia@odu.edu

############################################
###   		  base_calculator.sh	     ###
############################################

# Requirements: parallel
# Input files: fq.gz files from single-digest Sbf1 RADseq, pair-end sequencing data.

# This script counts the number of each base in DNA fragments from fq.gz files according to the window size and region determined by user. 
# Returns a tab delimitated file (.tsv) listing file names, base counts and read info.
# The script has been tested in single-digest RADseq, pair-end data, and it is meant to be ran separately for each read direction, forward (F) and reverse (R).
# The remanning Sbf1 restriction site motif 'TGCAGG' should in theory be present at the beginning of all forward reads ("restrct_end") and absent from the reverse reads ("sheared_end"). 
# However, exceptions do occur and cause incorrect position comparisons.
# base_calculator.sh reports base counts for all, restriction and sheared ends in forward and reverse reads. 
# The user should eliminate "sheared end" records from forward data and "restrct_end" from reverse data for further analysis.

#### To run:

# 1. Make sure you have modified the slurm options above according to your system specifications.
# 2.-Set the "User Variable" bellow to the correct species (3-letter code), protocol (baited or unbaited), read direction (F or R), desired window size (1 if you want every position), and region (starting position, ending position). 
# In the example below, the script analyzes every position (i.e. win_size=1) in the entire length of the reads, which in this case was 140 bp (i.e. st_pos=1, en_pos=140) for foward reads and 150 bp for reverse.

# User Variables. Modify these according to your analysis. For readDir=$1, the "$1" means you will provide the direction directly in the command line. Thus, do not modify here.
sp=Ssp
protocol=unbaited
readDir=$1
win_size=1
st_pos=1
en_pos=150

# 3.- Check that the first "ls" statement in the script below (line 69) will list your input files. As default, it is set to "$1.fq.gz" files, the $1 represents the read direction, which will be entered the command line as "F", "R", "R1", or "R2".  
# if you have uncompressed files, delete the ".gz" from the ls statement and the change "zcat" to "cat" inside the parallel statement later in the same line.

# 4.- Save your changes to this script and execute in the command line as:
# sbatch <script name> <"readDir">

# example: sbatch base_calculator.sh "F"


### Script ### Do not modify unless you need to modify line 69 (see above) or you want to change the script :) #######################################

st_pos=$((st_pos - 1))
en_pos=$((en_pos - 1))
st_pos_win=$((st_pos + 1))
en_pos_win=$((st_pos_win + win_size - 1))

output="allbase_prop_${sp}_${protocol}_${readDir}_win${win_size}_st$((st_pos + 1))_en$((en_pos + 1))_baseCalculator.tsv"

export NPROC="${SLURM_CPUS_PER_TASK:=1}"
export HEADER="$(echo 'file A C G T N positions frag_end readDir' | sed 's: :\t:g')"
export CACHE_DIR="$(mktemp -d)"

cache_create() {
	mkdir -p "$CACHE_DIR"

	ls *$1.fq.gz | parallel --no-notice -k -j"$NPROC" "zcat {} |  paste - - - - | cut -f2  > $CACHE_DIR/{}.processed"
	ls $CACHE_DIR/*.processed | parallel --no-notice -k -j"$NPROC"  "grep '^TGCAGG' {}    > {}.TGCAGG"
	ls $CACHE_DIR/*.processed | parallel --no-notice -k -j"$NPROC"  "grep -v '^TGCAGG' {} > {}.none_TGCAGG"

	rm $CACHE_DIR/*.processed 
}

cache_query() {
	file="$1"
	TGCAGG="$2"
	cut_beg="$3"
	cut_end="$4"
	search=($5)

	f="$(mktemp)"
	cut -c "${cut_beg}-${cut_end}" "$CACHE_DIR/$file.processed.$TGCAGG" > "$f"

	for pattern in "${search[@]}"; do
		grep -o "$pattern" $f | wc -l 
	done | tr '\n' "\t"

	rm "$f"
}
export -f cache_query

cache_cleanup() {
	rm -rf "$CACHE_DIR"
}
trap cache_cleanup EXIT

base_prop_func(){
	start_pos_win="$(($1 + $3))"
	end_pos_win="$(($2 + $3))"
	readDir="$4"

	#echo -e "$HEADER" 

	for f in `ls *$readDir.fq`; do
		echo -e -n "$f\t"
		cache_query "$f" 'TGCAGG' "${start_pos_win}" "${end_pos_win}" 'A C G T N'
		echo -e "${start_pos_win}-${end_pos_win}\trestrct_end\t${readDir}"
	done

	for f in `ls *$readDir.fq`; do
		echo -e -n "$f\t"
		cache_query "$f" 'none_TGCAGG' "${start_pos_win}" "${end_pos_win}" 'A C G T N'
		echo -e "${start_pos_win}-${end_pos_win}\tsheared_end\t${readDir}"
	done
}
export -f base_prop_func

parallel --record-env

cache_create $readDir

echo -e "$HEADER" > "$output"
seq ${st_pos} ${en_pos} | parallel --env _ --no-notice -k -j$NPROC \
	"base_prop_func $st_pos_win $en_pos_win {} $readDir"  >> "$output"
