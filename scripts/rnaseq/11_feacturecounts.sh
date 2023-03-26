#!/usr/bin/bash

mkdir /srv/data/my_shared_data_folder/epilepsy/rnaseq_results/features


ALGN_DIR="/srv/data/my_shared_data_folder/epilepsy/rnaseq_results/alignment"
REF_DIR="/srv/data/my_shared_data_folder/epilepsy/data/rnaseq/ref"
COUNT_DIR="/srv/data/my_shared_data_folder/epilepsy/rnaseq_results/features"



SAMPLE_ID="SRR9733947 SRR9733948 SRR9733949 SRR9733950 SRR9733951 SRR9733952 SRR9733953 
SRR9733954 SRR9733955 SRR9733956 SRR9733957 SRR9733958 SRR9733959 SRR9733960 SRR9733961 
SRR9733962 SRR9733963 SRR9733964 SRR9733965 SRR9733966 SRR9733967 SRR9733968 SRR9733969 
SRR9733970 SRR9733971 SRR9733972 SRR9733973 SRR9733974 SRR9733975 SRR9733976 SRR9733977 
SRR9733978 SRR9733979 SRR9733980 SRR9733981 SRR9733982"


#=======READS QUANTIFICATION

# Run fragment count on TotalRNA (pair-end) alignment files

#Featurescounts for GSE127871

featureCounts \
    -T 32 \
    -s 2 \
    --minOverlap 10 \
    -B \
    -C \
    -Q 30 \
    -p \
    -g gene_id \
    -a ${REF_DIR}/GRCh38.gtf \
    -o ${COUNT_DIR}/featureCounts_GSE134697.txt \
    ${ALGN_DIR}/SRR8*_str_Aligned.sortedByCoord.out.bam


#Featurescounts for GSE134697

featureCounts \
    -T 32 \
    -s 2 \
    --minOverlap 10 \
    -B \
    -C \
    -Q 30 \
    -p \
    -g gene_id \
    -a ${REF_DIR}/GRCh38.gtf \
    -o ${COUNT_DIR}/featureCounts_GSE134697.txt \
    ${ALGN_DIR}/SRR9*_str_Aligned.sortedByCoord.out.bam

