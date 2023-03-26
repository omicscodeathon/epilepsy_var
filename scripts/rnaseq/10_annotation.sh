#!/bin/bash


mkdir /srv/data/my_shared_data_folder/epilepsy/rnaseq_results/vcf_anno


VAR_DIR="/srv/data/my_shared_data_folder/epilepsy/rnaseq_results/varcall"
humandb="/srv/data/my_shared_data_folder/epilepsy/data/rnaseq/ref/humandb"
ANNO_DIR="/srv/data/my_shared_data_folder/epilepsy/rnaseq_results/vcf_anno"


SAMPLE_ID="SRR9733947 SRR9733948 SRR9733949 SRR9733950 SRR9733951 SRR9733952 SRR9733953 
SRR9733954 SRR9733955 SRR9733956 SRR9733957 SRR9733958 SRR9733959 SRR9733960 SRR9733961 
SRR9733962 SRR9733963 SRR9733964 SRR9733965 SRR9733966 SRR9733967 SRR9733968 SRR9733969 
SRR9733970 SRR9733971 SRR9733972 SRR9733973 SRR9733974 SRR9733975 SRR9733976 SRR9733977 
SRR9733978 SRR9733979 SRR9733980 SRR9733981 SRR9733982"

#=====ANNOTATION===============

for SAMPLE in $SAMPLE_ID; do
    perl table_annovar.pl \
    ${VAR_DIR}/${SAMPLE}.vcf.gz \
    ${humandb}/ \
    -buildver hg38 \
    -out ${ANNO_DIR}/${SAMPLE} \
    -remove \
    -protocol refGene,knownGene,cytoBand,exac03,avsnp150,dbnsfp42a,,gnomad211_exome,clinvar_20221231 \
    -operation g,g,r,f,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish
done
