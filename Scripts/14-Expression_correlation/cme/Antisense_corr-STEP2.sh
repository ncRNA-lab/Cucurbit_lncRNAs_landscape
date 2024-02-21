#!/bin/bash

#SBATCH --job-name=cmeAcor2						# Job name.
#SBATCH --output=cme_antisense_corr_2.log				# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-02:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=5gb						# Job memory request.


####### MODULES
module load R/4.2.1

####### VARIABLES
specie="cme"
specie_long="C. melo"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/17-Correlation_DEFINITIVE"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/17-Correlation_DEFINITIVE/Additional_scripts"
flag_list="nr"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### NEW VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2_spe
mkdir -p $WD2_spe/antisense


####### PIPELINE

### CORRELATION
echo -e "\nCORRELATION..."

echo -e "\n\n##############################"
echo -e "########### STEP 2 ###########"
echo -e "##############################\n\n"

for flag in $flag_list; do

	echo -e "\nFLAG: "$flag
	
	## Variable.
	O=$WD2_spe/antisense/$flag
	O1=$WD2_spe/antisense/$flag/STEP1
	O2=$WD2_spe/antisense/$flag/STEP2
	
	## Directory.
	mkdir -p $O
	mkdir -p $O1
	mkdir -p $O2
	
	cd $O2

	Rscript $AS/Antisense-STEP2-Create_Cis_Table.R $O1 $O2 $WD1_spe "$specie_long" $flag
	
done


