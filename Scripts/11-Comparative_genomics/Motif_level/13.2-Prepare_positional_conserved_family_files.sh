#!/bin/bash

#SBATCH --job-name=S13mot2							# Job name.
#SBATCH --output=STEP13_motif_2.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=2								# Run on one mode. 
#SBATCH --cpus-per-task=20							# Number of tasks = cpus.
#SBATCH --time=00-02:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=200mb							# Job memory request.


####### MODULES
module load anaconda

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Motif_level/Functions.sh"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High"
strictness_list="ORIGINAL"
nonmatch_list="no"
n_iter=50

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1/Motif_level


####### PIPELINE: STEP 13.2

### ANALYSIS OF CONSERVATION AT MOTIF LEVEL
echo -e "\nANALYSIS OF CONSERVATION AT MOTIF LEVEL..."

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	## Variable.
	I=$WD1/Positional_level/$flag/04-Families
	O=$WD1/Motif_level/$flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/02-Preparation
	
	### PREPARE FAMILY FILES FOR MOTIF LEVEL ANALYSIS
	## For each nonmatch level, strictness level, confidence level and class code, create a fasta file by conserved family.
	cd $O/02-Preparation

	echo -e "\n\nPREPARATION: Create a fasta file by conserved family...\n"

	for strictness in $strictness_list; do
		mkdir -p $strictness
		echo -e $strictness"..."
		for nonmatch in $nonmatch_list; do
			mkdir -p $strictness/$nonmatch
			echo -e "\t"$nonmatch"..."
			for confidence in $confidence_list; do
				mkdir -p $strictness/$nonmatch/$confidence
				echo -e "\t\t"$confidence"..."
				for class in $class_list; do
					srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Prepare_motif_level_analysis $strictness $nonmatch $confidence $class $I $O $n_iter $SLURM_CPUS_PER_TASK &
				done
				wait
			done
		done
	done
done

