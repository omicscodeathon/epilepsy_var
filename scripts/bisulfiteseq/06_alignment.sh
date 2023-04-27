#!/usr/bin/bash

SAMPLES="SRR10493738 SRR10493739 SRR10493740"

for SAMPLE in $SAMPLES; do

bismark --non_directional -f ~/my_shared_data_folder/epilepsy/data/bisulfiteseq/ref -1 ${SAMPLE}_1.fa -2 ${SAMPLE}_2.fa report.txt

done
