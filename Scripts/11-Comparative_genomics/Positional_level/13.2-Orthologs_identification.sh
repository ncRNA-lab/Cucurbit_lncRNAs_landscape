#!/bin/bash

#SBATCH --job-name=S13pos2					# Job name.
#SBATCH --output=STEP13_positional_2.log			# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=1						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=40					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=500mb					# Job memory request.


####### MODULES
module load anaconda

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Additional_scripts"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
specie_list="cma cme"
flag_list="nr"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### SOFTWARES
### ORTHOFINDER
export OFPATH=$SP/OrthoFinder
export PATH=$PATH:${OFPATH}
export PATH=$PATH:${OFPATH}/bin
export PATH=$PATH:${OFPATH}/tools

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1/Positional_level


####### PIPELINE

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag"\n"
	
	## Variable.
	O=$WD1/Positional_level/$flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/02-Orthologs
	mkdir -p $O/02-Orthologs/Inference
	mkdir -p $O/02-Orthologs/Inference/Outputs
	mkdir -p $O/02-Orthologs/Tables_orthologs_1_to_1

	### IDENTIFY ORTHOLOGS
	echo -e "\n\nIDENTIFY ORTHOLOGS...\n"

	## ORTHOFINDER
	echo -e "\n\nOrthologs: OrthoFinder...\n"

	cd $O/02-Orthologs/Inference
	if [ -d "Orthofinder" ]; then
		rm -r Orthofinder
	fi
	mkdir Orthofinder
	cd Orthofinder

	mkdir Proteomes_temp
	for spe in $specie_list; do
		cp $AI/Proteomes/$spe.fa ./Proteomes_temp/
	done

	>$O/02-Orthologs/Inference/Outputs/stdout_orthofinder_all.log
	orthofinder -f Proteomes_temp -o All -t $SLURM_CPUS_PER_TASK >> $O/02-Orthologs/Inference/Outputs/stdout_orthofinder_all.log 2>&1

	rm -r Proteomes_temp

	echo -e "\nOrthologs: Select 1:1 orthologs...\n"

	Extract_1_to_1_orthologs_across_all-Approach_2-OrthoFinder.py \
		--path-orthofinder $O/02-Orthologs/Inference/Orthofinder \
		--output $O/02-Orthologs/Tables_orthologs_1_to_1/Orthotable_1_to_1_orthofinder.tsv
done


