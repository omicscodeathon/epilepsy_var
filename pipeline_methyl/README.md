# Methylomic Variations in Temporal Lobe Epilepsy Subtypes: Focus on Hippocampal Sclerosis analysis pipeline

This repository contains the reproducibility steps for the differential methylation analysis of WGBS data of TLE patients with(3) and without(3) hippocampal sclerosis. The pipeline downloads WGBS 
data from the sra explorer, performs quality checks using fastqc. The downloaded dataset is then trimmed using trimgalore and aligned to the bismark prepared reference human genome hg38 followed 
by methylation calling with bismark. 
tool.
 
## SETTING UP THE ENVIRONMENT

### CONDA

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

### MAMBA

```
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh
bash Mambaforge-$(uname)-$(uname -m).sh -b -p $HOME/mambaforge
```

### INSTALLATION OF TOOLS

```
mamba install -y -c conda-forge -c bioconda \
    fastqc=0.11.9 \
    multiqc=1.13 \
    bismark=0.23.1 \
    trim-galore=0.6.7 \
    bedtools=2.30.0 \
    samtools=1.15.1 \
    metilene=0.2.8
```

## DATA DOWNLOADING AND PREPROCESSING

. Download Reference genome

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```  

. 
