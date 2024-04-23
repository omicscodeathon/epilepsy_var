#!/bin/bash
#
#SBATCH --job-name=DMR_Epi_Var_trial
#SBATCH --ntasks=20 # Number of cores/threads
#SBATCH --mem=128000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --time=5-12:00:00
#SBATCH --output=job_output_%j.txt
#SBATCH --error=job_error_%j.txt

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

################
# Load Modules #
################

module load R/4.0
module load homer

########
# DM.R #
########

call="Rscript \
--vanilla \
/share/lasallelab/Ensi/project/dmrichr/DM.main.R \
--genome hg38 \
--coverage 1 \
--perGroup '0.1' \
--minCpGs 5 \
--maxPerms 6 \
--maxBlockPerms 6 \
--cutoff '0.1' \
--testCovariate 'Diagnosis' \
--sexCheck FALSE \
--GOfuncR TRUE \
--EnsDb FALSE \
--cores 1"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
