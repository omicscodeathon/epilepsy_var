#!/usr/bin/bash

SAMPLES="SRR10493738 SRR10493739 SRR10493740"

for SAMPLE in $SAMPLES; do

bismark --non_directional --genome ~/my_shared_data_folder/epilepsy/data/bisulfiteseq/ref -1 ${SAMPLE}_1_val_1.fq.gz -2 ${SAMPLE}_2_val_2.fq.gz -o output_directory


done
