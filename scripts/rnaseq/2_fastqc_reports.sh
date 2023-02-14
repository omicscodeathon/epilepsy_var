#!/usr/bin/bash


# create QC report output directory
mkdir -p ~/my_shared_data_folder/epilepsy/rnaseq_results/QC_report
mkdir ~/my_shared_data_folder/epilepsy/rnaseq_results/MQC_report


RES_DIR="~/my_shared_data_folder/epilepsy/rnaseq_results"
FASTQ_DIR="~/my_shared_data_folder/epilepsy/data/rnaseq/fastq"

# Sample listes
SAMPLE="SRR9733947 SRR9733948 SRR9733949 SRR9733950 SRR9733951 SRR9733952 SRR9733953 
SRR9733954 SRR9733955 SRR9733956 SRR9733957 SRR9733958 SRR9733959 SRR9733960 SRR9733961 
SRR9733962 SRR9733963 SRR9733964 SRR9733965 SRR9733966 SRR9733967 SRR9733968 SRR9733969 
SRR9733970 SRR9733971 SRR9733972 SRR9733973 SRR9733974 SRR9733975 SRR9733976 SRR9733977 
SRR9733978 SRR9733979 SRR9733980 SRR9733981 SRR9733982"


# Run quality control on each fastq files

for reads in $FASTQ_DIR/*.fastq; 
    do
        fastqc ${reads} -o ${RES_DIR}/QC_report
done

# Run MultiQc on each forward (_1.fastq) QC files, this will output a html file
# combining all samples

multiqc ${RES_DIR}/QC_report/*_1.fastq -o ${RES_DIR}/MQC_report