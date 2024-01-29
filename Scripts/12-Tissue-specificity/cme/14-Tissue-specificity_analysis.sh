#!/bin/bash

#SBATCH --job-name=cmeS14						# Job name.
#SBATCH --output=cme_STEP14.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
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
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/12-Tissue-specificity/Additional_scripts"
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


####### PIPELINE: STEP 14

### TISSUE-SPECIFICITY ANALYSIS
echo -e "\nTISSUE-SPECIFICITY ANALYSIS..."

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	O=$WD3_spe/ALL/$flag
	
	mkdir -p $O
	mkdir -p $O/STEP1
	mkdir -p $O/STEP2
	mkdir -p $O/STEP3
	mkdir -p $O/STEP4
	mkdir -p $O/Outputs
	
	echo -e "\n\t-STEP 1"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $O/Outputs/stdout_TS_STEP1.log --quiet --exclusive Rscript $AS/01-ALL-Tissue-specificity.R $WD1_spe $WD2_spe $WD3_spe $AI $flag $specie

	echo -e "\n\t-STEP 2"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $O/Outputs/stdout_TS_STEP2.log --quiet --exclusive 02-ALL-Tissue-specificity.py --path $WD3_spe --flag $flag --specie $specie

	echo -e "\n\t-STEP 3"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $O/Outputs/stdout_TS_STEP3.log --quiet --exclusive Rscript $AS/03-ALL-Tissue-specificity.R $WD2_spe $WD3_spe $flag $specie

	echo -e "\n\t-STEP 4"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $O/Outputs/stdout_TS_STEP4.log --quiet --exclusive Rscript $AS/04-ALL-Tissue-specificity.R $WD2_spe $WD3_spe $flag $specie

done


