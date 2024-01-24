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
mkdir -p $WD5_spe/01-Metadata_EA
mkdir -p $WD5_spe/02-EA
mkdir -p $WD5_spe/03-Metadata_DEA
mkdir -p $WD5_spe/04-DEA
mkdir -p $WD5_spe/05-Tables_and_Figures


####### PIPELINE

echo -e "\n\n#############################"
echo -e "########### STEP 1 ##########"
echo -e "#############################\n\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/Filter_and_prepare_metadata.R $Metadata $WD5_spe $WD2_spe/Salmon/ALL/nr/03-Quant


echo -e "\n\n#############################"
echo -e "########### STEP 2 ##########"
echo -e "#############################\n\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/Exploratory_analysis.R $WD5_spe $WD1_spe/STEP-FINAL/Files/Joined/ALL/nr $WD2_spe/Salmon/ALL/nr/03-Quant


echo -e "\n\n#############################"
echo -e "########### STEP 3 ##########"
echo -e "#############################\n\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/DEA.R $WD5_spe $WD1_spe/STEP-FINAL/Files/Joined/ALL/nr $WD2_spe/Salmon/ALL/nr/03-Quant $alpha_value


echo -e "\n\n#############################"
echo -e "########### STEP 4 ##########"
echo -e "#############################\n\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/Tables_and_figures.R $WD5_spe $WD1_spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv $WD1_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv $WD3/Positional_level/Approach_2/nr/05-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv $WD4/approach_1/ALL/nr/STEP3/mean-TAU.tsv $alpha_value $specie


echo -e "\n\n#############################"
echo -e "############ END ############"
echo -e "#############################\n\n"



