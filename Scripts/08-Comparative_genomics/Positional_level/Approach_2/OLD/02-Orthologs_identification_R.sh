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
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Positional_level/Approach_2/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Positional_level/Approach_2/r

####### SOFTWARES
### ORTHOFINDER
export OFPATH=$SP/OrthoFinder
export PATH=$PATH:${OFPATH}
export PATH=$PATH:${OFPATH}/bin
export PATH=$PATH:${OFPATH}/tools

### INPARANOID
# It's necessary to be in the directory where the executable is found.

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Positional_level
mkdir -p $WD/08-comparative_genomics/Positional_level/Approach_2
mkdir -p $WD/08-comparative_genomics/Positional_level/Approach_2/r
mkdir -p $WD1/02-Orthologs
mkdir -p $WD1/02-Orthologs/Logs


####### PIPELINE

### IDENTIFY ORTHOLOGS

echo -e "\n\nIDENTIFY ORTHOLOGS...\n"

## ORTHOFINDER
echo -e "\n\nOrthologs: OrthoFinder...\n"

cd $WD1/02-Orthologs
if [ -d "$WD1/02-Orthologs/Orthofinder" ]; then
	rm -r Orthofinder
fi
mkdir Orthofinder
cd Orthofinder

mkdir Proteomes_temp
for spe in $Species_list; do
	cp $AI/Proteomes/$spe.fa ./Proteomes_temp/
done

>../Logs/Orthofinder_stdout_all.log
orthofinder -f Proteomes_temp -o All -t $SLURM_CPUS_PER_TASK >> ../Logs/Orthofinder_stdout_all.log 2>&1

rm -r Proteomes_temp


## INPARANOID
echo -e "\n\nOrthologs: Inparanoid...\n"

cd $WD1/02-Orthologs
if [ -d "$WD1/02-Orthologs/Inparanoid" ]; then
	rm -r Inparanoid
fi
mkdir Inparanoid
cd Inparanoid

mkdir Proteomes_temp
for spe in $Species_list; do
	cp $AI/Proteomes/$spe.fa ./Proteomes_temp/
done

>../Logs/Inparanoid_stdout_all.log
cd $SP/inparanoid
perl inparanoid.pl \
	-input-dir $WD1/02-Orthologs/Inparanoid/Proteomes_temp \
	-out-dir $WD1/02-Orthologs/Inparanoid/All \
	-out-stats \
	-out-table \
	-out-sqltable \
	-cores $SLURM_CPUS_PER_TASK \
	-cores-diamond $SLURM_CPUS_PER_TASK >> $WD1/02-Orthologs/Logs/Inparanoid_stdout_all.log 2>&1

cd $WD1/02-Orthologs/Inparanoid
rm -r Proteomes_temp


