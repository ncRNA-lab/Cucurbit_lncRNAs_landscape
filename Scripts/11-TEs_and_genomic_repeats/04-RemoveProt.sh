#!/bin/bash

#SBATCH --job-name=RemProt					# Job name.
#SBATCH --output=RemoveProt.log				# Standard output and error log.
#SBATCH --partition=medium					# Partition (queue)
#SBATCH --ntasks=3						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=15					# Number of tasks = cpus.
#SBATCH --time=4-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb					# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/11-TEs_and_genomic_repeats/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"

####### NEW AND OTHER VARIABLES
WD1=$WD/11-TEs_and_genomic_repeats

####### DIRECTORY
mkdir -p $WD/11-TEs_and_genomic_repeats
mkdir -p $WD1/01-Repeat_calling
mkdir -p $WD1/01-Repeat_calling/04-RemoveProt_unknown_repeats

####### PIPELINE

### REPEAT CALLING
## Unknown repeats filtered are used to search against a plant protein database where proteins from transposons are excluded..
echo -e "\n\nExclusion of gene fragments...\n"
cd $WD1/01-Repeat_calling/04-RemoveProt_unknown_repeats
for spe in $Species_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_RemoveProt $spe $WD1 $SLURM_CPUS_PER_TASK &
done
wait

