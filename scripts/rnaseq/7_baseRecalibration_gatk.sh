#!/usr/bin/bash

mkdir /srv/data/my_shared_data_folder/epilepsy/rnaseq_results/alignment


ALGN_DIR="/srv/data/my_shared_data_folder/epilepsy/rnaseq_results/alignment"
REF_DIR="/srv/data/my_shared_data_folder/epilepsy/data/rnaseq/ref"


SAMPLE_ID="SRR9733947"

#========BASE RECAL======
echo "Running Base recalibration"

for SAMPLE in $SAMPLE_ID; do
        gatk --java-options -Xmx8g BaseRecalibrator \
                -I ${ALGN_DIR}/${SAMPLE}_spl.bam \
                -R ${REF_DIR}/GRCh38.fasta \
                --known-sites ${REF_DIR}/GRCh38.known_sites.vcf.gz \
                -O ${ALGN_DIR}/${SAMPLE}_recal.table 
done

echo "Running APLYBQSR"

for SAMPLE in $SAMPLE_ID; do
        gatk --java-options -Xmx8g ApplyBQSR \
                -R ${REF_DIR}/GRCh38.fasta \
                -I ${ALGN_DIR}/${SAMPLE}_spl.bam \
                --bqsr-recal-file ${ALGN_DIR}/${SAMPLE}_recal.table \
                -O ${ALGN_DIR}/${SAMPLE}_bqsr.bam
        
        #rm -r ${ALGN_DIR}/${SAMPLE}_spl.bam

done