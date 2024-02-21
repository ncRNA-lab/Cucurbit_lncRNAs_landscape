#!/bin/bash

#SBATCH --job-name=cmeUcor2						# Job name.
#SBATCH --output=cme_intergenic_corr_2.log				# Standard output and error log.
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
dist_list="500 1000 2000 5000 10000 20000 50000 100000"
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
mkdir -p $WD2_spe/intergenic


####### PIPELINE

### CORRELATION
echo -e "\nCORRELATION..."

echo -e "\n\n##############################"
echo -e "########### STEP 2 ###########"
echo -e "##############################\n\n"

for flag in $flag_list; do

	echo -e "\nFLAG: "$flag
	
	## Variable.
	O=$WD2_spe/intergenic/$flag
	O1=$WD2_spe/intergenic/$flag/STEP1
	O2=$WD2_spe/intergenic/$flag/STEP2
	
	## Directory.
	mkdir -p $O
	mkdir -p $O1
	mkdir -p $O2
	
	cd $O2

	Rscript $AS/Intergenic-STEP2-Create_Cis_Table.R $O1 $O2 $WD1_spe "$dist_list" "$specie_long" $flag
	
done


