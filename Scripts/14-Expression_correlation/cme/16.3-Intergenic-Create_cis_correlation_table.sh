#!/bin/bash

#SBATCH --job-name=cmeS16U3						# Job name.
#SBATCH --output=cme_STEP16_U3.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=3							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=1-00:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=80gb						# Job memory request.


####### MODULES
module load R/4.2.1

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/13-DEA"
WD4="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/14-Expression_correlation/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/14-Expression_correlation/Functions.sh"
dist_list="500 1000 2000 5000 10000 20000 50000 100000"
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
mkdir -p $WD4_spe/intergenic


####### PIPELINE

### CORRELATION
echo -e "\nCORRELATION..."

echo -e "\n\n##############################"
echo -e "########### STEP 3 ###########"
echo -e "##############################\n\n"

for flag in $flag_list; do

	echo -e "\nFLAG: "$flag
	
	## Variable.
	O=$WD4_spe/intergenic/$flag
	O2=$WD4_spe/intergenic/$flag/STEP2
	O3=$WD4_spe/intergenic/$flag/STEP3
	
	## Directory.
	mkdir -p $O
	mkdir -p $O2
	mkdir -p $O3
	mkdir -p $O3/Outputs
	
	cd $O3

	echo -e "-Closest..."

	if [ -d "$WD3_spe/03-Metadata_DEA/" ] && [ "$(ls -A "$WD3_spe/03-Metadata_DEA/")" ]; then
		experiment_list=$(find "$WD3_spe/03-Metadata_DEA/" -type f -exec basename {} \; | sed 's/.tsv//g' | sort)
		for exp in $experiment_list; do
			echo -e "\t-Processing experiment: $exp"
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive --output $O3/Outputs/Closest-$exp.log $F task_STEP3_CLOSEST-INTERGENIC $exp $O2 $O3 $WD1_spe $WD2_spe $WD3_spe $flag $AS &
		done
		wait
	else
		echo "\tDEA folder doesn't exist. Cannot proceed."
	fi

	echo -e "-Range..."

	if [ -d "$WD3_spe/03-Metadata_DEA/" ] && [ "$(ls -A "$WD3_spe/03-Metadata_DEA/")" ]; then
		experiment_list=$(find "$WD3_spe/03-Metadata_DEA/" -type f -exec basename {} \; | sed 's/.tsv//g' | sort)
		for exp in $experiment_list; do
			echo -e "\t-Processing experiment: $exp"
			for dist in $dist_list; do
				echo -e "\t\t-Processing distance: $dist"
				srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive --output $O3/Outputs/Range-$exp\-$dist.log $F task_STEP3_RANGE-INTERGENIC $exp $dist $O2 $O3 $WD1_spe $WD2_spe $WD3_spe $flag $AS &
			done
		done
		wait
	else
		echo "\tDEA folder doesn't exist. Cannot proceed."
	fi
	
done


