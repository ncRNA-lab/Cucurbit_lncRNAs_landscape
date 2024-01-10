#!/bin/bash

#SBATCH --job-name=NRP1					# Job name.
#SBATCH --output=NR_P1.log					# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=9						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=1					# Number of tasks = cpus.
#SBATCH --time=0-01:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb					# Job memory request.


####### MODULES
module load R/4.1.2

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
mkdir -p $WD1/01-LncRNAs_and_Genes


####### PIPELINE

### LISTS OF LNCRNAS AND GENES 
## For each confidence level, class code and specie, select the lncRNAs. Then, create a TXT file with the ids of the selected lncRNAs and the genes.
cd $WD1/01-LncRNAs_and_Genes

echo -e "\n\nCreate list of lncRNAs (selected) and genes...\n"

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		mkdir -p $confidence/$class
		echo -e "\t"$class"..."
		for spe in $Species_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_lncRNAS_NR_and_genes_list $spe $confidence $class $WD $WD1 $AS &
		done
		wait
	done
done


