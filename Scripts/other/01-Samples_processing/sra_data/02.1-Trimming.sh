#!/bin/bash

#SBATCH --job-name=cmeS2.1							# Job name.
#SBATCH --output=cme_STEP2.1.log						# Standard output and error log.
#SBATCH --qos=medium								# Partition (queue)
#SBATCH --ntasks=16								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=5-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2500mb							# Job memory request.


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/01-Samples_processing/sra_data/Functions.sh"
trim_adapters="/storage/ncRNA/Softwares/BBMap_38.90/bbmap/resources"

####### NEW VARIABLES
WD1=$WD"/01-Samples_processing/sra_data"
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List.txt"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD/01-Samples_processing
mkdir -p $WD1
mkdir -p $WD1/02-Trimmed_data
mkdir -p $WD1/02-Trimmed_data/outputs


####### PIPELINE: STEP 2.1

echo -e "\n\n\n######################"
echo -e "###### STEP 2.1 ######"
echo -e "######################\n"

cd $WD1/02-Trimmed_data

### TRIMMING
## Trim the libraries.
echo -e "\nTRIMMING..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/02-Trimmed_data/outputs/stdout_Trimming_$SRR.log --quiet --exclusive $F task_Trimming $SRR $trim_adapters $SLURM_CPUS_PER_TASK &
done
wait

