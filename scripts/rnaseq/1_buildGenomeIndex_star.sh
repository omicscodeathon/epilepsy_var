#!/usr/bin/bash


mkdir -p ~/my_shared_data_folder/epilepsy/data/rnaseq/ref/star

REF_DIR="~/my_shared_data_folder/epilepsy/data/rnaseq/ref"

# Indexing reference genome 
samtools faidx ${REF_DIR}/GRCh38.fna

NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${REF_DIR}/GRCh38.fai`


# Run STAR genome generate (index)
STAR \
    --runThreadN 10 \
    --runMode genomeGenerate \
    --genomeDir ${REF_DIR}/star \
    --genomeFastaFiles ${REF_DIR}/GRCh38.fna \
    --sjdbGTFfile ${REF_DIR}/GRCh38.gtf \
    --genomeSAindexNbases $NUM_BASES
