#!/bin/bash

#SBATCH --job-name=carS1							# Job name.
#SBATCH --output=car_STEP1.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=10								# Run on one mode. 
#SBATCH --cpus-per-task=10							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb							# Job memory request.


####### VARIABLES
specie="car"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/01-sample_processing_and_selection"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/01-sample_processing_and_selection/Functions.sh"

####### NEW VARIABLES
WD1=$WD"/"$specie
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List.txt"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD1
mkdir -p $WD1/01-Raw_data
mkdir -p $WD1/01-Raw_data/fastqc
mkdir -p $WD1/01-Raw_data/multiqc
mkdir -p $WD1/01-Raw_data/outputs


####### PIPELINE: STEP 1

cd $WD1/01-Raw_data

### DOWNLOAD
## Download each library using prefetch.
echo -e "\nDOWNLOAD LIBRARIES FROM SRA DATABASE..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/01-Raw_data/outputs/stdout_Prefetch_$SRR.log --quiet --exclusive $F task_Prefetch $SRR &
done
wait

### FASTQC
## Execute Fastqc.
echo -e "\nFASTQC..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/01-Raw_data/outputs/stdout_Fastqc_$SRR.log --quiet --exclusive $F task_Fastqc1 $SRR $SLURM_CPUS_PER_TASK &
done
wait

### MULTIQC
## Execute Multiqc.
echo -e "\nMULTIQC..."

multiqc ./fastqc -o ./multiqc -f >> ./outputs/stdout_Multiqc.log 2>&1


