#!/bin/bash

#SBATCH --job-name=cmeS16X4						# Job name.
#SBATCH --output=cme_STEP16_X4.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-02:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=20gb						# Job memory request.


####### MODULES
module load R/4.2.1

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/14-Expression_correlation/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/14-Expression_correlation/Functions.sh"
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
echo -e "########### STEP 4 ###########"
echo -e "##############################\n\n"

for flag in $flag_list; do

	echo -e "\nFLAG: "$flag
	
	## Variable.
	O=$WD1_spe/antisense/$flag
	O3=$WD1_spe/antisense/$flag/STEP3
	O4=$WD1_spe/antisense/$flag/STEP4
	
	## Directory.
	mkdir -p $O
	mkdir -p $O3
	mkdir -p $O4

	cd $O4

	Rscript $AS/Antisense-STEP4-Join_cis_correlation_tables.R $O3 $O4

done


