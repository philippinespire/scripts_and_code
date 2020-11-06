# Example file names of a single individual (052) with data from different lanes (L7 and L8):
#`20200331_PIRE-Ssp-CGub_052-Capture_CKDL200150313-1a-D712-AK1545_H7T2LBBXX_L7.R1.fq.gz`
#`20200331_PIRE-Ssp-CGub_052-Capture_CKDL200150313-1a-D712-AK1545_H7T2LBBXX_L8.R1.fq.gz`

# User variables: Match the FILE_DIRECTION and set the number of threads according to your system

FILE_DIRECTION=R1
THREADS=40

# Also match the Lane info in the script below to your files (L7 and L8 in this example)

#Script:

# step 1: list all forward fq.gz files
# ls *$FILE_DIRECTION.fq.gz | 

# step 2: remove file extension from Lane to end, and remove duplicate file basenames
# sed 's/L[78].$FILE_DIRECTION.fq.gz//' | uniq | 

# step 3: parallel cat files together
# parallel --no-notice -j $THREADS "cat {}L[78].$FILE_DIRECTION.fq.gz > {}L78.$FILE_DIRECTION.fq.gz"



ls *$FILE_DIRECTION.fq.gz | sed 's/L[78].$FILE_DIRECTION.fq.gz//' | uniq | parallel --no-notice -j $THREADS "cat {}L[78].$FILE_DIRECTION.fq.gz > {}L78.$FILE_DIRECTION.fq.gz"

ls *R1.fq.gz | sed 's/L[78].R1.fq.gz//' | less -S

