#!/bin/bash

#SBATCH --job-name=cmeAcorr3						# Job name.
#SBATCH --output=cme_antisense_corr_3.log				# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=5							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-05:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=20gb						# Job memory request.


####### MODULES
module load R/4.2.1

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-quantification"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/16-DEA"
WD4="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/17-Correlation_DEFINITIVE"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/17-Correlation_DEFINITIVE/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/17-Correlation_DEFINITIVE/Functions.sh"
flag_list="nr"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### NEW VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie
WD3_spe=$WD3/$specie
WD4_spe=$WD4/$specie

####### DIRECTORY
mkdir -p $WD4
mkdir -p $WD4_spe
mkdir -p $WD4_spe/antisense


####### PIPELINE

### CORRELATION
echo -e "\nCORRELATION..."

echo -e "\n\n##############################"
echo -e "########### STEP 3 ###########"
echo -e "##############################\n\n"

for flag in $flag_list; do

	echo -e "\nFLAG: "$flag
	
	## Variable.
	O=$WD4_spe/antisense/$flag
	O2=$WD4_spe/antisense/$flag/STEP2
	O3=$WD4_spe/antisense/$flag/STEP3
	
	## Directory.
	mkdir -p $O
	mkdir -p $O2
	mkdir -p $O3
	mkdir -p $O3/Outputs
	
	cd $O3

	if [ -d "$WD3_spe/03-Metadata_DEA/" ] && [ "$(ls -A "$WD3_spe/03-Metadata_DEA/")" ]; then
		experiment_list=$(find "$WD3_spe/03-Metadata_DEA/" -type f -exec basename {} \; | sed 's/.tsv//g' | sort)
		for exp in $experiment_list; do
			echo -e "\t-Processing experiment: $exp"
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive --output $O3/Outputs/$exp.log $F task_STEP3-ANTISENSE $exp $O2 $O3 $WD1_spe $WD2_spe $WD3_spe $flag $AS &
		done
		wait
	else
		echo "\tDEA folder doesn't exist. Cannot proceed."
	fi
	
done


