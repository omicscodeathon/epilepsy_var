#!/usr/bin/bash

cd ~/my_shared_data_folder/epilepsy/data/Bisulfite-seq

SAMPLES="SRR10493735 SRR10493736 SRR10493737 SRR10493738 SRR10493739 SRR10493740"
for SAMPLE in $SAMPLES; do

prefetch ${SAMPLE} --max-size u && fastq-dump --split-files ${SAMPLE}

done
