#!/bin/bash

#SBATCH --job-name=NRP4					# Job name.
#SBATCH --output=NR_P4.log					# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=5						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=10					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=500mb					# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Positional_level/Approach_2/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
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
mkdir -p $WD1/04-Families


####### PIPELINE

## CLASSIFY LNCRNAS INTO FAMILIES
## For each confidence level and class code, classify syntenic lncRNAs into conservated families using home-made script Classify_into_families_mod.py.
cd $WD1/04-Families



for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	for class in $Classes_list; do
		mkdir -p $confidence/$class
		>$confidence/$class/ids_by_specie.tsv
		for spe in $Species_list; do
			awk -v a="$spe" '{print a"\t"$0"-"a}' $WD1/01-LncRNAs_and_Genes/$confidence/$class/$spe\_lncRNA_ids.txt >> $confidence/$class/ids_by_specie.tsv
		done
	done
done



echo -e "\n\nORIGINAL (families): genesNearby=3, minOverlap=3, minSideOverlap=1, nonMatch='no'...\n"

for confidence in $Confidence_levels_list; do
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $WD1 $SLURM_CPUS_PER_TASK 'ORIGINAL_no' &
	done
	wait
done

echo -e "\n\nORIGINAL (families): genesNearby=3, minOverlap=3, minSideOverlap=1, nonMatch='yes'...\n"

for confidence in $Confidence_levels_list; do
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $WD1 $SLURM_CPUS_PER_TASK 'ORIGINAL_yes' &
	done
	wait
done



echo -e "\n\nRELAXED (families): genesNearby=4, minOverlap=2, minSideOverlap=1, nonMatch='no'...\n"

for confidence in $Confidence_levels_list; do
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $WD1 $SLURM_CPUS_PER_TASK 'RELAXED_no' &
	done
	wait
done

echo -e "\n\nRELAXED (families): genesNearby=4, minOverlap=2, minSideOverlap=1, nonMatch='yes'...\n"

for confidence in $Confidence_levels_list; do
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $WD1 $SLURM_CPUS_PER_TASK 'RELAXED_yes' &
	done
	wait
done



echo -e "\n\nSTRICT (families): genesNearby=3, minOverlap=4, minSideOverlap=2, nonMatch='no'...\n"

for confidence in $Confidence_levels_list; do
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $WD1 $SLURM_CPUS_PER_TASK 'STRICT_no' &
	done
	wait
done

echo -e "\n\nSTRICT (families): genesNearby=3, minOverlap=4, minSideOverlap=2, nonMatch='yes'...\n"

for confidence in $Confidence_levels_list; do
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $WD1 $SLURM_CPUS_PER_TASK 'STRICT_yes' &
	done
	wait
done



echo -e "\n\nMORE-STRICT (families): genesNearby=4, minOverlap=6, minSideOverlap=3, nonMatch='no'...\n"

for confidence in $Confidence_levels_list; do
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $WD1 $SLURM_CPUS_PER_TASK 'MORE-STRICT_no' &
	done
	wait
done

echo -e "\n\nMORE-STRICT (families): genesNearby=4, minOverlap=6, minSideOverlap=3, nonMatch='yes'...\n"

for confidence in $Confidence_levels_list; do
	echo -e $confidence"..."
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $WD1 $SLURM_CPUS_PER_TASK 'MORE-STRICT_yes' &
	done
	wait
done


