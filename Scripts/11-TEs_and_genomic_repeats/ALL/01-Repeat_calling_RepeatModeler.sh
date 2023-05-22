#!/bin/bash

#SBATCH --job-name=RepMod				# Job name.
#SBATCH --output=Repeat_Modeler.log			# Standard output and error log.
#SBATCH --partition=medium				# Partition (queue)
#SBATCH --ntasks=2					# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=16				# Number of tasks = cpus.
#SBATCH --time=7-00:00:00				# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb				# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/11-TEs_and_genomic_repeats/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"

####### NEW AND OTHER VARIABLES
WD1=$WD/11-TEs_and_genomic_repeats

####### DIRECTORY
mkdir -p $WD/11-TEs_and_genomic_repeats
mkdir -p $WD1/01-Repeat_calling
mkdir -p $WD1/01-Repeat_calling/01-RepeatModeler
mkdir -p $WD1/01-Repeat_calling/01-RepeatModeler/Logs

####### PIPELINE

### REPEAT CALLING
## RepeatModeler.
echo -e "\n\nRepeat calling: RepeatModeler...\n"
for spe in $Species_list; do
	srun -N1 -n1 -c16 --quiet --exclusive $F task_RepeatModeler $spe $AI $WD1 4 &
done
wait


