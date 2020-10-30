# example file name `20200331_PIRE-Ssp-CGub_052-Capture_CKDL200150313-1a-D712-AK1545_H7T2LBBXX_L7.R1.fq.gz`

FILE_DIRECTION=R1
THREADS=40

# step 1: list all foward fq.gz files
# ls *$FILE_DIRECTION.fq.gz | 

# step 2: remove file extension from Lane to end, and remove duplicate file basenames
# sed 's/L[78].$FILE_DIRECTION.fq.gz//' | uniq | 

# step 3: parallel cat files together
# parallel --no-notice -j $THREADS "cat {}L[78].$FILE_DIRECTION.fq.gz > {}L78.$FILE_DIRECTION.fq.gz" 


ls *$FILE_DIRECTION.fq.gz | sed 's/L[78].$FILE_DIRECTION.fq.gz//' | uniq | parallel --no-notice -j $THREADS "cat {}L[78].$FILE_DIRECTION.fq.gz > {}L78.$FILE_DIRECTION.fq.gz" 

ls *R1.fq.gz | sed 's/L[78].R1.fq.gz//' | less -S

