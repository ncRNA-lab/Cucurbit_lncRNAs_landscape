#!/bin/bash

#SBATCH --job-name=cmeS12dist						# Job name.
#SBATCH --output=cme_STEP12_dist.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=2							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus. 
#SBATCH --time=0-01:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb						# Job memory request.


####### MODULES
module load R/4.1.2
module load biotools

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/10-Genomic_distribution"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/10-Genomic_distribution/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/10-Genomic_distribution/Functions.sh"
flag_list="nr"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2_spe


####### PIPELINE
for flag in $flag_list; do
	
	## Directory.
	mkdir -p $WD2_spe/$flag
	mkdir -p $WD2_spe/$flag/Outputs

	## Visualize the lncRNAs and genes distribution.
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2_spe/$flag/Outputs/stdout_Dist_$specie.log --quiet --exclusive $F task_transcript_density $specie $WD1_spe $WD2_spe $AI $AS $flag &
	
done
wait


