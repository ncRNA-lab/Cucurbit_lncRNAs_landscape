#!/bin/bash

#SBATCH --job-name=RP2						# Job name.
#SBATCH --output=R_P2.log					# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=1						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=40					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=500mb					# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Positional_level/Approach_2/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Positional_level/Approach_2/r

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
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Positional_level
mkdir -p $WD/08-comparative_genomics/Positional_level/Approach_2
mkdir -p $WD/08-comparative_genomics/Positional_level/Approach_2/r
mkdir -p $WD1/02-Orthologs
mkdir -p $WD1/02-Orthologs/Inference
mkdir -p $WD1/02-Orthologs/Inference/Logs
mkdir -p $WD1/02-Orthologs/Tables_orthologs_1_to_1


####### PIPELINE

### IDENTIFY ORTHOLOGS

echo -e "\n\nIDENTIFY ORTHOLOGS...\n"

## ORTHOFINDER
echo -e "\n\nOrthologs: OrthoFinder...\n"

cd $WD1/02-Orthologs/Inference
if [ -d "$WD1/02-Orthologs/Inference/Orthofinder" ]; then
	rm -r Orthofinder
fi
mkdir Orthofinder
cd Orthofinder

mkdir Proteomes_temp
for spe in $Species_list; do
	cp $AI/Proteomes/$spe.fa ./Proteomes_temp/
done

>$WD1/02-Orthologs/Inference/Logs/Orthofinder_stdout_all.log
orthofinder -f Proteomes_temp -o All -t $SLURM_CPUS_PER_TASK >> $WD1/02-Orthologs/Inference/Logs/Orthofinder_stdout_all.log 2>&1

rm -r Proteomes_temp

echo -e "\nOrthologs: Select 1:1 orthologs...\n"

Extract_1_to_1_orthologs_across_all-Approach_2-OrthoFinder.py \
	--path-orthofinder $WD1/02-Orthologs/Inference/Orthofinder \
	--output $WD1/02-Orthologs/Tables_orthologs_1_to_1/Orthotable_1_to_1_orthofinder.tsv


