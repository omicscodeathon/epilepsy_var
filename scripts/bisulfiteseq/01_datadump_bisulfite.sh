#!/usr/bin/env bash

# This script downloads fastq.qz files (bisulfite sequencing data)
# Dependencies: fastq-dump

SAMPLES=$(cat ../accessions/Bisulfite-data-accessions.txt)

for SAMPLE in $SAMPLES; do

 prefetch ${SAMPLE} --max-size u && fastq-dump --gzip --split-files ${SAMPLE} && rm -Rv ${SAMPLE}

done
