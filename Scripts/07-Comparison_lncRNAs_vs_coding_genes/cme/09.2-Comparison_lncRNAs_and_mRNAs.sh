#!/bin/bash

#SBATCH --job-name=cmeS9.2						# Job name.
#SBATCH --output=cme_STEP9_2.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-02:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=4gb						# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/07-Comparison_lncRNAs_vs_coding_genes"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/07-Comparison_lncRNAs_vs_coding_genes/Additional_scripts"

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


####### PIPELINE
echo -e "\n\nREDUNDANT SEQUENCES: CREATING FIGURES\n"
Rscript Comparison_lncRNAs_and_PCGs.R $WD1_spe $WD2_spe $WD3_spe "r"
echo -e "\n\nNON-REDUNDANT SEQUENCES: CREATING FIGURES\n"
Rscript Comparison_lncRNAs_and_PCGs.R $WD1_spe $WD2_spe $WD3_spe "nr"

