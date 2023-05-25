#!/bin/bash

#SBATCH --job-name=dendist					# Job name.
#SBATCH --output=dendist.log					# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=9						# Run on one mode.
#SBATCH --cpus-per-task=2					# Number of tasks = cpus. 
#SBATCH --time=0-01:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb					# Job memory request.


####### MODULES
module load R/4.1.2
module load biotools

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/10-Distribution/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/10-Distribution/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"

####### NEW AND OTHER VARIABLES
WD1=$WD/05-predict_lncRNAs
WD2=$WD/10-Distribution

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD/10-Distribution
mkdir -p $WD2/Transcript_density
mkdir -p $WD2/Transcript_density/r
mkdir -p $WD2/Transcript_density/nr


####### PIPELINE

## NON-REDUNDANT: Visualize the lncRNAs and genes distribution.
echo -e "\nNON-REDUNDANT:"
for spe in $Species_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_transcript_density_non_redundant $spe $WD1 $WD2 $AI $AS &
done
wait

## REDUNDANT: Visualize the lncRNAs and genes distribution.
echo -e "\nREDUNDANT:"
for spe in $Species_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_transcript_density_redundant $spe $WD1 $WD2 $AI $AS &
done
wait


