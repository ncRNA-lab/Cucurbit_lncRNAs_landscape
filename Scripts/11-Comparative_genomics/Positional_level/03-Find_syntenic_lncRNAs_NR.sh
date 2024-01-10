#!/bin/bash

#SBATCH --job-name=NRP3					# Job name.
#SBATCH --output=NR_P3.log					# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=5						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb					# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Positional_level/Approach_2/Functions.sh"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High Medium Low"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Positional_level/Approach_2/nr

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Positional_level
mkdir -p $WD/08-comparative_genomics/Positional_level/Approach_2
mkdir -p $WD/08-comparative_genomics/Positional_level/Approach_2/nr
mkdir -p $WD1/03-Synteny_analysis


####### PIPELINE

## SYNTENY ANALYSIS
## For each confidence level and class code, execute synteny analysis.
cd $WD1/03-Synteny_analysis



echo -e "\n\nORIGINAL synteny analysis: genesNearby=3, minOverlap=3, minSideOverlap=1, nonMatch='no'...\n"

param='ORIGINAL'
genesNearby=3
minOverlap=3
minSideOverlap=1
nonMatch='no'

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $WD1 $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
	done
	wait
done

echo -e "\n\nORIGINAL synteny analysis: genesNearby=3, minOverlap=3, minSideOverlap=1, nonMatch='yes'...\n"

param='ORIGINAL'
genesNearby=3
minOverlap=3
minSideOverlap=1
nonMatch='yes'

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $WD1 $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
	done
	wait
done



echo -e "\n\nRELAXED synteny analysis: genesNearby=4, minOverlap=2, minSideOverlap=1, nonMatch='no'...\n"

param='RELAXED'
genesNearby=4
minOverlap=2
minSideOverlap=1
nonMatch='no'

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $WD1 $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
	done
	wait
done

echo -e "\n\nRELAXED synteny analysis: genesNearby=4, minOverlap=2, minSideOverlap=1, nonMatch='yes'...\n"

param='RELAXED'
genesNearby=4
minOverlap=2
minSideOverlap=1
nonMatch='yes'

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $WD1 $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
	done
	wait
done



echo -e "\n\nSTRICT synteny analysis: genesNearby=3, minOverlap=4, minSideOverlap=2, nonMatch='no'...\n"

param='STRICT'
genesNearby=3
minOverlap=4
minSideOverlap=2
nonMatch='no'

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $WD1 $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
	done
	wait
done

echo -e "\n\nSTRICT synteny analysis: genesNearby=3, minOverlap=4, minSideOverlap=2, nonMatch='yes'...\n"

param='STRICT'
genesNearby=3
minOverlap=4
minSideOverlap=2
nonMatch='yes'

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $WD1 $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
	done
	wait
done



echo -e "\n\nMORE-STRICT synteny analysis: genesNearby=4, minOverlap=6, minSideOverlap=3, nonMatch='no'...\n"

param='MORE-STRICT'
genesNearby=4
minOverlap=6
minSideOverlap=3
nonMatch='no'

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $WD1 $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
	done
	wait
done

echo -e "\n\nMORE-STRICT synteny analysis: genesNearby=4, minOverlap=6, minSideOverlap=3, nonMatch='yes'...\n"

param='MORE-STRICT'
genesNearby=4
minOverlap=6
minSideOverlap=3
nonMatch='yes'

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_synteny_analysis $confidence $class $WD1 $genesNearby $minOverlap $minSideOverlap $nonMatch $param &
	done
	wait
done

