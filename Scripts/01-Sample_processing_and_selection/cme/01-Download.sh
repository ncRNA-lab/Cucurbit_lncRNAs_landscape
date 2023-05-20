#!/bin/bash

#SBATCH --job-name=cmeS1							# Job name.
#SBATCH --output=cme_STEP1.log						# Standard output and error log.
#SBATCH --qos=medium								# Partition (queue)
#SBATCH --ntasks=90								# Run on one mode. 
#SBATCH --cpus-per-task=1							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=15gb							# Job memory request.


########################################################################################
#
# This script is used to download all the libraries of a 
# particular specie present in the txt file provided.
#
########################################################################################


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/01-Sample_processing_and_selection"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/01-Sample_processing_and_selection/Functions.sh"

####### NEW VARIABLES
WD1=$WD"/"$specie
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List.txt"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD1
mkdir -p $WD1/01-Raw_data
mkdir -p $WD1/01-Raw_data/Fastqc
mkdir -p $WD1/01-Raw_data/Multiqc
mkdir -p $WD1/01-Raw_data/Outputs


####### PIPELINE: STEP 1

cd $WD1/01-Raw_data

### DOWNLOAD
## Download each library using prefetch.
echo -e "\nDOWNLOAD LIBRARIES FROM SRA DATABASE..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/01-Raw_data/Outputs/stdout_Prefetch_$SRR.log --quiet --exclusive $F task_Prefetch $SRR &
done
wait

### FASTQC
## Execute Fastqc.
echo -e "\nFASTQC..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/01-Raw_data/Outputs/stdout_Fastqc_$SRR.log --quiet --exclusive $F task_Fastqc1 $SRR $SLURM_CPUS_PER_TASK &
done
wait

### MULTIQC
## Execute Multiqc.
echo -e "\nMULTIQC..."

multiqc ./Fastqc -o ./Multiqc -f >> ./Outputs/stdout_Multiqc.log 2>&1


