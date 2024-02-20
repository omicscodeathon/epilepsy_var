#!/usr/bin/bash

cd ~/my_shared_data_folder/epilepsy/data/bisulfiteseq

mkdir -p ~/my_shared_data_folder/epilepsy/data/bisulfiteseq/results/trimmed-mRNA

FASTQ_DIR="~/my_shared_data_folder/epilepsy/data/bisulfiteseq"

OUT_DIR="~/my_shared_data_folder/epilepsy/data/bisulfiteseq//results/trimmed-mRNA"

SAMPLES="SRR10493735 SRR10493736 SRR10493737 SRR10493738 SRR10493739 SRR10493740"

for SAMPLE in $SAMPLES; do

trim_galore --paired --cores 8 ${FASTQ_DIR}/${SAMPLE}_1.fastq.gz ${FASTQ_DIR}/${SAMPLE}_2.fastq.gz

done
