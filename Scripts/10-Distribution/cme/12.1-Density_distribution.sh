#!/bin/bash

#SBATCH --job-name=distvvi					# Job name.
#SBATCH --output=distribution_vvi.log				# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=2						# Run on one mode.
#SBATCH --cpus-per-task=2					# Number of tasks = cpus. 
#SBATCH --time=0-01:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb					# Job memory request.


####### MODULES
module load R/4.1.2
module load biotools

####### VARIABLES
specie="vvi"
WD1="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results/10-Distribution"
AI="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Scripts/10-Distribution/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Scripts/10-Distribution/Functions.sh"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2_spe
mkdir -p $WD2_spe/nr
mkdir -p $WD2_spe/r
mkdir -p $WD2_spe/LOGS


####### PIPELINE

## NON-REDUNDANT: Visualize the lncRNAs and genes distribution.
echo -e "\nNON-REDUNDANT:"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2_spe/LOGS/stdout_Dist_$specie\_NR.log --quiet --exclusive $F task_transcript_density_non_redundant $specie $WD1_spe $WD2_spe $AI $AS &

## REDUNDANT: Visualize the lncRNAs and genes distribution.
echo -e "\nREDUNDANT:"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2_spe/LOGS/stdout_Dist_$specie\_R.log --quiet --exclusive $F task_transcript_density_redundant $specie $WD1_spe $WD2_spe $AI $AS &
wait


