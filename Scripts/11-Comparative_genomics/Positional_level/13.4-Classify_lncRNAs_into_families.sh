#!/bin/bash

#SBATCH --job-name=S13pos4					# Job name.
#SBATCH --output=STEP13_positional_4.log			# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=5						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=10					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=500mb					# Job memory request.


####### MODULES
module load anaconda

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Positional_level/Functions.sh"
specie_list="car cla cma cme cmo cpe csa lsi mch"
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
	mkdir -p $O/04-Families

	## CLASSIFY LNCRNAS INTO FAMILIES
	## For each confidence level and class code, classify syntenic lncRNAs into conservated families using home-made script Classify_into_families_mod.py.
	cd $O/04-Families



	for confidence in $confidence_list; do
		mkdir -p $confidence
		for class in $class_list; do
			mkdir -p $confidence/$class
			>$confidence/$class/ids_by_specie.tsv
			for spe in $specie_list; do
				awk -v a="$spe" '{print a"\t"$0"-"a}' $O/01-LncRNAs_and_Genes/$confidence/$class/$spe\_lncRNA_ids.txt >> $confidence/$class/ids_by_specie.tsv
			done
		done
	done



	echo -e "\n\nORIGINAL (families): genesNearby=3, minOverlap=3, minSideOverlap=1, nonMatch='no'...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $O $SLURM_CPUS_PER_TASK 'ORIGINAL_no' &
		done
		wait
	done

	echo -e "\n\nORIGINAL (families): genesNearby=3, minOverlap=3, minSideOverlap=1, nonMatch='yes'...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $O $SLURM_CPUS_PER_TASK 'ORIGINAL_yes' &
		done
		wait
	done



	echo -e "\n\nRELAXED (families): genesNearby=4, minOverlap=2, minSideOverlap=1, nonMatch='no'...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $O $SLURM_CPUS_PER_TASK 'RELAXED_no' &
		done
		wait
	done

	echo -e "\n\nRELAXED (families): genesNearby=4, minOverlap=2, minSideOverlap=1, nonMatch='yes'...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $O $SLURM_CPUS_PER_TASK 'RELAXED_yes' &
		done
		wait
	done



	echo -e "\n\nSTRICT (families): genesNearby=3, minOverlap=4, minSideOverlap=2, nonMatch='no'...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $O $SLURM_CPUS_PER_TASK 'STRICT_no' &
		done
		wait
	done

	echo -e "\n\nSTRICT (families): genesNearby=3, minOverlap=4, minSideOverlap=2, nonMatch='yes'...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $O $SLURM_CPUS_PER_TASK 'STRICT_yes' &
		done
		wait
	done



	echo -e "\n\nMORE-STRICT (families): genesNearby=4, minOverlap=6, minSideOverlap=3, nonMatch='no'...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $O $SLURM_CPUS_PER_TASK 'MORE-STRICT_no' &
		done
		wait
	done

	echo -e "\n\nMORE-STRICT (families): genesNearby=4, minOverlap=6, minSideOverlap=3, nonMatch='yes'...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_classify_into_families $confidence $class $O $SLURM_CPUS_PER_TASK 'MORE-STRICT_yes' &
		done
		wait
	done
done


