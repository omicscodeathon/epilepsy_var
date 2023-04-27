#!/usr/bin/bash

cd ~/my_shared_data_folder/epilepsy/data/bisulfiteseq

SAMPLES="SRR10493737 SRR10493738 SRR10493739 SRR10493740"

for SAMPLE in $SAMPLES; do

bismark_methylation_extractor -p --merge_non_CpG --bedGraph --counts --buffer_size 10G ${SAMPLE}_1.fa_bismark_bt2_pe.bam

done