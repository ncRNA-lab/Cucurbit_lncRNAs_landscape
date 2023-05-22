#!/bin/bash

#SBATCH --job-name=Maskcma						# Job name.
#SBATCH --output=Repeat_Masker_RMfiltered_cma.log			# Standard output and error log.
#SBATCH --partition=medium						# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=16						# Number of tasks = cpus.
#SBATCH --time=4-00:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb						# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/11-TEs_and_genomic_repeats/Functions.sh"
Species_list="cma"

####### NEW AND OTHER VARIABLES
WD1=$WD/11-TEs_and_genomic_repeats

####### DIRECTORY
mkdir -p $WD/11-TEs_and_genomic_repeats
mkdir -p $WD1/01-Repeat_calling
mkdir -p $WD1/01-Repeat_calling/06-RepeatMasker_RMfiltered
mkdir -p $WD1/01-Repeat_calling/06-RepeatMasker_RMfiltered/Logs

####### PIPELINE

### REPEAT CALLING
## RepeatMasker.
echo -e "\n\nRepeat calling: RepeatMasker...\n"
for spe in $Species_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_RepeatMasker_RMfiltered $spe $AI $WD1 $SLURM_CPUS_PER_TASK &
done
wait


