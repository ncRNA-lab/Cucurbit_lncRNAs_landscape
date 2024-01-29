#!/bin/bash

#SBATCH --job-name=S13pos1					# Job name.
#SBATCH --output=STEP13_positional_1.log			# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=9						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=1					# Number of tasks = cpus.
#SBATCH --time=0-01:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb					# Job memory request.


####### MODULES
module load anaconda
module load R/4.1.2

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/11-Comparative_genomics/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/11-Comparative_genomics/Positional_level/Functions.sh"
specie_list="car cla cma cme cmo cpe csa lsi mch"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High Medium Low"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2/Positional_level


####### PIPELINE: STEP 13.1

### ANALYSIS OF CONSERVATION AT POSITIONAL LEVEL
echo -e "\nANALYSIS OF CONSERVATION AT POSITIONAL LEVEL..."

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	## Variable.
	O=$WD2/Positional_level/$flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/01-LncRNAs_and_Genes

	### LISTS OF LNCRNAS AND GENES 
	## For each confidence level, class code and specie, select the lncRNAs. Then, create a TXT file with the ids of the selected lncRNAs and the genes.
	cd $O/01-LncRNAs_and_Genes

	echo -e "\n\nCreate list of lncRNAs (selected) and PCGs...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			for spe in $specie_list; do
				srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_lncRNA_and_PCG_list $spe $confidence $class $WD1 $O $flag $AS &
			done
			wait
		done
	done
done


