#!/bin/bash

#SBATCH --job-name=cmeS1							# Job name.
#SBATCH --output=cme_STEP1.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=90								# Run on one mode. 
#SBATCH --cpus-per-task=1							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=14gb							# Job memory request.


####### VARIABLES
specie="cme"
update="Update_2023_03_29"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Scripts/"$update"/01-Samples_processing/Functions.sh"

####### NEW VARIABLES
WD1=$WD"/"$update"/01-Samples_processing"
WD1_spe=$WD1"/"$specie
Acc_list=$AI"/sra-info/"$update"/accession_list/"$specie"-SRR_Acc_List.txt"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD"/"$update
mkdir -p $WD1
mkdir -p $WD1_spe
mkdir -p $WD1_spe/01-Raw_data
mkdir -p $WD1_spe/01-Raw_data/fastqc
mkdir -p $WD1_spe/01-Raw_data/multiqc
mkdir -p $WD1_spe/01-Raw_data/outputs


####### PIPELINE: STEP 1

echo -e "\n\n\n####################"
echo -e "###### STEP 1 ######"
echo -e "####################\n"

cd $WD1_spe/01-Raw_data

### DOWNLOAD
## Download each library using prefetch.
echo -e "\nDOWNLOAD LIBRARIES FROM SRA DATABASE..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/01-Raw_data/outputs/stdout_Prefetch_$SRR.log --quiet --exclusive $F task_Prefetch $SRR &
done
wait

### FASTQC
## Execute Fastqc.
echo -e "\nFASTQC..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/01-Raw_data/outputs/stdout_Fastqc_$SRR.log --quiet --exclusive $F task_Fastqc1 $SRR $SLURM_CPUS_PER_TASK &
done
wait

### MULTIQC
## Execute Multiqc.
echo -e "\nMULTIQC..."

>./outputs/stdout_Multiqc.log
multiqc ./fastqc -o ./multiqc -f >> ./outputs/stdout_Multiqc.log 2>&1


