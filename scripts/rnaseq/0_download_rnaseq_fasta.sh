#!/usr/bin/bash

# Extract fastqs from SRR from a text file containing SRR IDs

while read line; do
    
    prefetch $line --max-size u
    fastq-dump --gzip --split-files $line

done < GSE134697_SRR_Acc_List.txt