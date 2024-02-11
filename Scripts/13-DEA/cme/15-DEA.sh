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

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/13-DEA"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/13-DEA/Additional_scripts"
flag_list="nr"
alpha_value=0.05

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie
WD3_spe=$WD3/$specie
Metadata=$AI/DEA_metadata/$specie

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD3
mkdir -p $WD3_spe


####### PIPELINE: STEP 15

### DIFFERENTIAL EXPRESSION ANALYSIS
echo -e "\nDIFFERENTIAL EXPRESSION ANALYSIS..."

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	ID=$WD1_spe/STEP-FINAL/Database
	IAF=$WD1_spe/STEP-FINAL/Files/ALL/$flag
	IGF=$WD1_spe/STEP-FINAL/Files/Genes
	IAQ=$WD2_spe/ALL/$flag/03-Quant
	O=$WD3_spe/$flag
	
	mkdir -p $O
	mkdir -p $O/01-Metadata_EA
	mkdir -p $O/02-EA
	mkdir -p $O/03-Metadata_DEA
	mkdir -p $O/04-DEA
	mkdir -p $O/05-Tables_and_Figures
	mkdir -p $O/Outputs

	echo -e "\n\t-STEP 1: FILTER AND PREPARE THE METADATA"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $O/Outputs/stdout_STEP1.log --quiet --exclusive Rscript $AS/Filter_and_prepare_metadata.R $Metadata $O $IAQ

	echo -e "\n\t-STEP 2 and 3: EXPLORATORY ANALYSIS"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $O/Outputs/stdout_STEP2-3.log --quiet --exclusive Rscript $AS/Exploratory_analysis.R $O $IAF $IAQ

	echo -e "\n\t-STEP 4: DEA"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $O/Outputs/stdout_STEP4.log --quiet --exclusive Rscript $AS/DEA.R $O $IAF $IAQ $alpha_value

	echo -e "\n\t-STEP 5: TABLES AND FIGURES"
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $O/Outputs/stdout_STEP5.log --quiet --exclusive Rscript $AS/Tables_and_figures.R $O $ID/Database_LncRNAs_${flag^^}.tsv $IGF/ORIGINAL_GENES.tsv $alpha_value $specie

done


