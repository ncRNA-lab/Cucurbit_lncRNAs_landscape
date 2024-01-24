#!/bin/bash

#SBATCH --job-name=S14							# Job name.
#SBATCH --output=STEP14_tissue_specificity.log			# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus. It depends on the number of process of your parallelization.
#SBATCH --time=1-00:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=5gb						# Job memory request.


####### MODULES
module load R/4.2.1
module load biotools

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/12-Tissue-specificity"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/12-Tissue-specificity/additional_scripts"
flag_list="nr"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie
WD3_spe=$WD3/$specie

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD3
mkdir -p $WD3_spe
mkdir -p $WD3_spe/ALL


####### PIPELINE
for flag in $flag_list; do
	
	# Directory.
	mkdir -p $WD3_spe/ALL/$flag
	mkdir -p $WD3_spe/ALL/$flag/STEP1
	
	echo -e "\n\n#############################"
	echo -e "########### STEP 1 ##########"
	echo -e "#############################\n\n"

	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/01-ALL-Tissue-specificity.R $WD1_spe $WD2_spe $WD3_spe $AI $flag $specie


	echo -e "\n\n#############################"
	echo -e "########### STEP 2 ##########"
	echo -e "#############################\n\n"

	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/02-ALL-Tissue-specificity.R $WD3_spe $flag $specie


	echo -e "\n\n#############################"
	echo -e "########### STEP 3 ##########"
	echo -e "#############################\n\n"

	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/03-ALL-Tissue-specificity.R $WD3_spe $flag $specie


	echo -e "\n\n#############################"
	echo -e "########### STEP 4 ##########"
	echo -e "#############################\n\n"




	echo -e "\n\n#############################"
	echo -e "############ END ############"
	echo -e "#############################\n\n"

done



