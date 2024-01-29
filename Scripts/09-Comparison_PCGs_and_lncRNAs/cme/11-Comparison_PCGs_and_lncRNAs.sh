#!/bin/bash

#SBATCH --job-name=cmeS11							# Job name.
#SBATCH --output=cme_STEP11.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=1								# Run on one mode.
#SBATCH --cpus-per-task=2							# Number of tasks = cpus.
#SBATCH --time=0-02:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=4gb							# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
specie="cme"
specie_long="C. melo"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/07-Get_intergenic_regions"
WD4="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-TEs_and_genomic_repeats"
WD5="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/09-Comparison_PCGs_and_lncRNAs"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/09-Comparison_PCGs_and_lncRNAs/Additional_scripts"
flag_list="nr"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie
WD3_spe=$WD3/$specie
WD4_spe=$WD4/$specie
WD5_spe=$WD5/$specie

####### DIRECTORY
mkdir -p $WD5
mkdir -p $WD5_spe
mkdir -p $WD5_spe/Outputs


####### PIPELINE: STEP 11

### COMPARISON: PCG, LNCRNA AND INTERGENIC REGIONS
echo -e "\nCOMPARISON: PCG, LNCRNA AND INTERGENIC REGIONS..."

for flag in $flag_list; do
	echo -e "\nFLAG: "$flag
	Rscript $AS/Comparison_PCGs_and_lncRNAs.R $WD1_spe $WD2_spe $WD3_spe $WD4_spe $WD5_spe $flag "$specie_long" >> $WD5_spe/Outputs/stdout_comparison_PCGs_and_lncRNAs.log 2>&1
done


