#!/usr/bin/bash

FASTQ_DIR="~/my_shared_data_folder/epilepsy/data/bisulfiteseq/"

SAMPLES="SRR10493735 SRR10493736 SRR10493737 SRR10493738 SRR10493739 SRR10493740"

for SAMPLE in $SAMPLES; do

awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' ${SAMPLE}_1_val_1.fq > ${SAMPLE}_1.fa

awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' ${SAMPLE}_2_val_2.fq > ${SAMPLE}_2.fa

done
