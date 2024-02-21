#!/bin/bash

#SBATCH --job-name=cmeAcorr6						# Job name.
#SBATCH --output=cme_antisense_corr_6.log				# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-00:30:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=5gb						# Job memory request.


####### MODULES
module load R/4.2.1

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/17-Correlation_DEFINITIVE"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/17-Correlation_DEFINITIVE/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/17-Correlation_DEFINITIVE/Functions.sh"
flag_list="nr"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### NEW VARIABLES
WD1_spe=$WD1/$specie

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1_spe
mkdir -p $WD1_spe/antisense


####### PIPELINE

### CORRELATION
echo -e "\nCORRELATION..."

echo -e "\n\n##############################"
echo -e "########### STEP 6 ###########"
echo -e "##############################\n\n"

for flag in $flag_list; do

	echo -e "\nFLAG: "$flag
	
	## Variable.
	O=$WD1_spe/antisense/$flag
	O5=$WD1_spe/antisense/$flag/STEP5
	O6=$WD1_spe/antisense/$flag/STEP6
	
	## Directory.
	mkdir -p $O
	mkdir -p $O5
	mkdir -p $O6

	cd $O6

	Rscript $AS/Antisense-STEP6-Create_Final_Table-RANDOM.R $O5 $O6

done


