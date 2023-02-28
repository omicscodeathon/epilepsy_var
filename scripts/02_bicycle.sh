#!/usr/bin/bash

cd ~/my_shared_data_folder/epilepsy/data/Bisulfite-seq/Bicycle

# Declaring an alias
alias bicycle="docker run --rm -v `pwd`/data:/data -u `id -u
\`whoami\`` -it singgroup/bicycle bicycle"

# Declare the docker data directories
REFERENCE_GENOME="data/referenceGenome"
TARGET_REGIONS="data/target_regions.bed"
SAMPLES_DIR="data/raw_data"
PROJECT_DIR="data/myproject"

# Create project
bicycle create-project -p $PROJECT_DIR -r $REFERENCE_GENOME -f
$SAMPLES_DIR --paired-mate1-regexp _1.fastq

# Align reads to both references
bicycle align -p $PROJECT_DIR -t 4 --bowtie2-quals phred33 --bowtie2-I
0 --bowtie2-X 500

# Perform methylation analysis and methylcytosine calling
bicycle analyze-methylation -p $PROJECT_DIR -t 4 –n 1 --removeambiguous --only-with-one-alignment -b $TARGET_REGIONS

# Perform the differential methylation analysis:
$bicycle analyze-differential-methylation -p $PROJECT_DIR -c
SRR2052487,SRR2052488,SRR2052489,SRR2052490,SRR2052491 -t
SRR2052492,SRR2052493,SRR2052494,SRR2052495,SRR2052496 -b
$TARGET_REGIONS
