#!/bin/bash

#SBATCH --job-name=cmeS10RMask					# Job name.
#SBATCH --output=cme_STEP10_RepeatMasker.log				# Standard output and error log.
#SBATCH --qos=medium							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=64						# Number of tasks = cpus.
#SBATCH --time=4-00:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb						# Job memory request.


####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-TEs_and_genomic_repeats"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/08-TEs_and_genomic_repeats/Functions.sh"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1_spe
mkdir -p $WD1_spe/01-Repeat_calling
mkdir -p $WD1_spe/01-Repeat_calling/02-RepeatMasker
mkdir -p $WD1_spe/01-Repeat_calling/02-RepeatMasker/Outputs

####### PIPELINE

### REPEAT CALLING: REPEATMASKER
echo -e "\nREPEAT CALLING: REPEATMASKER..."

cd $WD1_spe/01-Repeat_calling/02-RepeatMasker
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/01-Repeat_calling/02-RepeatMasker/Outputs/stdout_RepeatMasker.log --quiet --exclusive $F task_RepeatMasker $specie $AI $WD1_spe $SLURM_CPUS_PER_TASK


