# Integrated Analysis of Genetic Variation, Gene Expression and Methylation Changes in Epilepsy

## Table of Contents
- [Background](#Background)
- [Objectives](#Objectives)
- [Data](#Data)
- [Methods](#Methods)
- [Results](#Results)
- [Reproducibility](#Reproducibility)
- [Discussion](#Discussion)
- [Team](#Team)
- [Acknowledgment](#Acknowledgment)
- [References](#References)

## Background
Epilepsy is a neurological disorder characterized by recurrent seizures caused by abnormal brain activity. While seizures are the hallmark of epilepsy, not all seizures are due to epilepsy. Idiopathic epilepsy presents only with seizures, while in symptomatic epilepsy, seizures are a sign of an underlying brain condition. Studies have shown that epileptic activity can result in significant behavioral and cognitive impairments. These cognitive and personality deficiencies are strongly linked to the duration and frequency of seizures, as well as earlier age of onset of epilepsy.

Genetic factors play a significant role in epilepsy with over 500 genes linked to the disorder. Several studies have identified specific gene mutations associated with different types of epilepsy. Epigenetic changes, such as DNA methylation, have also been found to play a role in epilepsy and several differentially methylated genes have been linked to specific types of epilepsy. Despite these findings, over 30% of epilepsy patients cannot effectively manage their seizures with medication.

## Objectives
To analyse genetic variation, gene expression and methylation changes in transcriptomic and bisulphite whole genome sequences of patients with epilepsy

## Data

## Methods
Variant calling and analysis of differentially expressed genes (DEGs) were carried out on two sets of temporal lobe epilepsy (TLE) transcriptomic data from. Bisulfite mapping, methylation calling, and differential methylation (DMR) analysis were also conducted on two sets of TLE whole genome bisulphite sequence (WGBS) data. For RNA-seq variant calling, the Haplotypecaller tools from GATK4 were utilized, and variant annotation was performed using ANNOVAR. Differential expression analysis was executed using DESeq2 in R, with annotation facilitated by the AnnotationDbi and org.Hs.eg.db packages. DMR analysis of bisulfite sequences was accomplished using the Bismark tool. Bisulphite mapping and methylation calling of the TLE WGBS data were carried out using Bismark tool, while DMRichR was used to perform differential methylation analysis.

## Analysis workflow

<p align="center">
<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/dd96bd28-5ee8-4364-b6d8-c8a224858cfc" width="800" height="700">
</p>

## Results
<p align="center">
  
<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/9ae4d3b9-a1e1-4c62-bfc0-c10de063dfca">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/f34d7175-73c9-41bd-b975-816209f5d232">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/b792c233-9720-4076-b1b0-b3037849b51a">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/950e65d9-711c-4c89-a7c9-3f9cbbf410b9">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/d4397e92-6e5b-46c7-bef6-819164b11182">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/1b18929d-ed99-4b11-89ad-cb8cfedc5861">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/bf4d0e56-2d5b-4fce-a25b-84dcece2c532">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/43ecb76b-e5f6-4531-afe4-736427f767dd">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/f47865e0-74b1-4f2c-91f0-def94e83b94d">

<img src="https://github.com/omicscodeathon/epilepsy_var/assets/116915872/cfb2e920-8e72-42c1-bbe7-19baccf56bcd">

</p>


## References

## Acknowledgment


## Team

Marion Nyaboke - MSc. Bioinformatics

Shamim Osata -MSc. Bioinformatics

Modibo K. Goita - MSc. Bioinformatics

Awe Olaitan - African Society for Bioinformatics and Computational Biology, South Africa
