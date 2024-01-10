#!/bin/bash

#SBATCH --job-name=RMask					# Job name.
#SBATCH --output=Repeat_Masker.log				# Standard output and error log.
#SBATCH --qos=medium						# Partition (queue)
#SBATCH --ntasks=1						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=64					# Number of tasks = cpus.
#SBATCH --time=4-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb					# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Scripts/08-TEs_and_genomic_repeats/Functions.sh"
Specie="vvi"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD/08-TEs_and_genomic_repeats/$Specie

####### DIRECTORY
mkdir -p $WD/08-TEs_and_genomic_repeats
mkdir -p $WD1_spe
mkdir -p $WD1_spe/01-Repeat_calling
mkdir -p $WD1_spe/01-Repeat_calling/02-RepeatMasker
mkdir -p $WD1_spe/01-Repeat_calling/02-RepeatMasker/Logs

####### PIPELINE

### REPEAT CALLING
## RepeatMasker.
echo -e "\n\nRepeat calling: RepeatMasker...\n"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_RepeatMasker $Specie $AI $WD1_spe $SLURM_CPUS_PER_TASK


