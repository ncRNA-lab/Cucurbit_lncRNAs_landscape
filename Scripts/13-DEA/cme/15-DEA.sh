#!/bin/bash

#SBATCH --job-name=cmeS15				# Job name.
#SBATCH --output=cme_STEP15.log			# Standard output and error log.
#SBATCH --qos=short					# Partition (queue)
#SBATCH --ntasks=1					# Run on one mode.
#SBATCH --cpus-per-task=2				# Number of tasks = cpus.
#SBATCH --time=1-00:00:00				# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=5gb				# Job memory request.


####### MODULES
module load R/4.2.1
module load biotools

####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
WD4="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/12-Tissue-specificity"
WD5="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/13-DEA"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/13-DEA/Additional_scripts"
flag_list="nr"
alpha_value=0.05

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie
WD5_spe=$WD5/$specie
Metadata=$AI/sra-info/metadata/Tables_DEA/$specie

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD5
mkdir -p $WD5_spe


####### PIPELINE
for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	mkdir -p $WD5_spe/$flag/01-Metadata_EA
	mkdir -p $WD5_spe/$flag/02-EA
	mkdir -p $WD5_spe/$flag/03-Metadata_DEA
	mkdir -p $WD5_spe/$flag/04-DEA
	mkdir -p $WD5_spe/$flag/05-Tables_and_Figures
	mkdir -p $WD5_spe/$flag/Outputs

	echo -e "\n\t-STEP 1: FILTER AND PREPARE THE METADATA"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD5_spe/$flag/Outputs/stdout_STEP1.log --quiet --exclusive Rscript $AS/Filter_and_prepare_metadata.R $Metadata $WD5_spe $WD2_spe/ALL/$flag/03-Quant

	echo -e "\n\t-STEP 2 and 3: EXPLORATORY ANALYSIS"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD5_spe/$flag/Outputs/stdout_STEP2-3.log --quiet --exclusive Rscript $AS/Exploratory_analysis.R $WD5_spe $WD1_spe/STEP-FINAL/Files/Joined/ALL/$flag $WD2_spe/ALL/$flag/03-Quant

	echo -e "\n\t-STEP 4: DEA"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD5_spe/$flag/Outputs/stdout_STEP4.log --quiet --exclusive Rscript $AS/DEA.R $WD5_spe $WD1_spe/STEP-FINAL/Files/Joined/ALL/$flag $WD2_spe/ALL/$flag/03-Quant $alpha_value

	echo -e "\n\t-STEP 5: TABLES AND FIGURES"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD5_spe/$flag/Outputs/stdout_STEP5.log --quiet --exclusive Rscript $AS/Tables_and_figures.R $WD5_spe $WD1_spe/STEP-FINAL/Database/Database_LncRNAs_${flag^^}.tsv $WD1_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv $WD3/Positional_level/Approach_2/$flag/05-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv $WD4/approach_1/ALL/$flag/STEP3/mean-TAU.tsv $alpha_value $specie

done

