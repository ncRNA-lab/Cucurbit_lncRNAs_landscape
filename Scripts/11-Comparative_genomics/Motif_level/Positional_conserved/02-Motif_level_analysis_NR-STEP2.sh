#!/bin/bash

#SBATCH --job-name=PosNR2							# Job name.
#SBATCH --output=Positional_NR_2.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=2								# Run on one mode. 
#SBATCH --cpus-per-task=20							# Number of tasks = cpus.
#SBATCH --time=00-02:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=200mb							# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Motif_level/Positional_conserved/Functions_NR.sh"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High Medium Low"
#strictness_list="ORIGINAL RELAXED STRICT MORE-STRICT"
#nonmatch_list="no yes"
strictness_list="ORIGINAL"
nonmatch_list="no"
n_iter=50

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Positional_level/Approach_2/nr/04-Families
WD2=$WD/08-comparative_genomics/Motif_level/nr/Positional_conserved

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Motif_level
mkdir -p $WD/08-comparative_genomics/Motif_level/nr
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/02-Preparation


####### PIPELINE

### PREPARE FAMILY FILES FOR MOTIF LEVEL ANALYSIS
## For each nonmatch level, strictness level, confidence level and class code, create a fasta file by conserved family.
cd $WD2/02-Preparation

echo -e "\n\nPREPARATION: Create a fasta file by conserved family...\n"

for strictness in $strictness_list; do
	mkdir -p $strictness
	echo -e $strictness"..."
	for nonmatch in $nonmatch_list; do
		mkdir -p $strictness/$nonmatch
		echo -e "\t"$nonmatch"..."
		for confidence in $Confidence_levels_list; do
			mkdir -p $strictness/$nonmatch/$confidence
			echo -e "\t\t"$confidence"..."
			for class in $Classes_list; do
				srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Prepare_motif_level_analysis $strictness $nonmatch $confidence $class $WD1 $WD2 $n_iter $SLURM_CPUS_PER_TASK &
			done
			wait
		done
	done
done


