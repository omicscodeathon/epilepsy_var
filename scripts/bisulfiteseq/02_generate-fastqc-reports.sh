#!/usr/bin/bash

cd ~/my_shared_data_folder/epilepsy/data/bisulfiteseq
#Script to generate FastQC reports using the FastQC tool

# fastq files directory
FASTQ_DIR="~/my_shared_data_folder/epilepsy/data/bisulfiteseq"

for file in $FASTQ_DIR/*.fastq; do
        fastqc ${file}
done

multiqc .