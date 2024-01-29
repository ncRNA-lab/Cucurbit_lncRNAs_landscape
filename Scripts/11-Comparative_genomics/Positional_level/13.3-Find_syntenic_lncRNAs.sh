#!/bin/bash

#SBATCH --job-name=S13pos3					# Job name.
#SBATCH --output=STEP13_positional_3.log			# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=5						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb					# Job memory request.


####### MODULES
module load anaconda

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Positional_level/Functions.sh"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High Medium Low"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1/Positional_level


####### PIPELINE

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	## Variable.
	O=$WD1/Positional_level/$flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/03-Synteny_analysis

	## SYNTENY ANALYSIS
	## For each confidence level and class code, execute synteny analysis.
	cd $O/03-Synteny_analysis



	echo -e "\n\nORIGINAL synteny analysis: genesNearby=3, minOverlap=3, minSideOverlap=1, nonMatch='no'...\n"

	param='ORIGINAL'
	genesNearby=3
	minOverlap=3
	minSideOverlap=1
	nonMatch='no'

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $O $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
		done
		wait
	done

	echo -e "\n\nORIGINAL synteny analysis: genesNearby=3, minOverlap=3, minSideOverlap=1, nonMatch='yes'...\n"

	param='ORIGINAL'
	genesNearby=3
	minOverlap=3
	minSideOverlap=1
	nonMatch='yes'

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $O $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
		done
		wait
	done



	echo -e "\n\nRELAXED synteny analysis: genesNearby=4, minOverlap=2, minSideOverlap=1, nonMatch='no'...\n"

	param='RELAXED'
	genesNearby=4
	minOverlap=2
	minSideOverlap=1
	nonMatch='no'

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $O $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
		done
		wait
	done

	echo -e "\n\nRELAXED synteny analysis: genesNearby=4, minOverlap=2, minSideOverlap=1, nonMatch='yes'...\n"

	param='RELAXED'
	genesNearby=4
	minOverlap=2
	minSideOverlap=1
	nonMatch='yes'

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $O $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
		done
		wait
	done



	echo -e "\n\nSTRICT synteny analysis: genesNearby=3, minOverlap=4, minSideOverlap=2, nonMatch='no'...\n"

	param='STRICT'
	genesNearby=3
	minOverlap=4
	minSideOverlap=2
	nonMatch='no'

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $O $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
		done
		wait
	done

	echo -e "\n\nSTRICT synteny analysis: genesNearby=3, minOverlap=4, minSideOverlap=2, nonMatch='yes'...\n"

	param='STRICT'
	genesNearby=3
	minOverlap=4
	minSideOverlap=2
	nonMatch='yes'

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $O $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
		done
		wait
	done



	echo -e "\n\nMORE-STRICT synteny analysis: genesNearby=4, minOverlap=6, minSideOverlap=3, nonMatch='no'...\n"

	param='MORE-STRICT'
	genesNearby=4
	minOverlap=6
	minSideOverlap=3
	nonMatch='no'

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $O $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
		done
		wait
	done

	echo -e "\n\nMORE-STRICT synteny analysis: genesNearby=4, minOverlap=6, minSideOverlap=3, nonMatch='yes'...\n"

	param='MORE-STRICT'
	genesNearby=4
	minOverlap=6
	minSideOverlap=3
	nonMatch='yes'

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $O $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
		done
		wait
	done
done


