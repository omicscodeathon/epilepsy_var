#!/usr/bin/bash

ALGN_DIR="/srv/data/my_shared_data_folder/epilepsy/rnaseq_results/alignment"
REF_DIR="/srv/data/my_shared_data_folder/epilepsy/data/rnaseq/ref"



SAMPLE_ID="SRR9733947"
#SRR9733948 SRR9733949 SRR9733950 SRR9733951 SRR9733952 SRR9733953 
#SRR9733954 SRR9733955 SRR9733956 SRR9733957 SRR9733958 SRR9733959 SRR9733960 SRR9733961 
#SRR9733962 SRR9733963 SRR9733964 SRR9733965 SRR9733966 SRR9733967 SRR9733968 SRR9733969 
#SRR9733970 SRR9733971 SRR9733972 SRR9733973 SRR9733974 SRR9733975 SRR9733976 SRR9733977 
#SRR9733978 SRR9733979 SRR9733980 SRR9733981 SRR9733982

#=======Split N CIGAR============

for SAMPLE in $SAMPLE_ID; do
    gatk --java-options -Xmx8g SplitNCigarReads \
        -R ${REF_DIR}/GRCh38.fna \
        -I ${ALGN_DIR}/${SAMPLE}_str_Aligned.sortedByCoord.out.bam \
        -O ${ALGN_DIR}/${SAMPLE}_spl.bam

    #rm -f ${ALGN_DIR}/${SAMPLE}_dup.bam
done


