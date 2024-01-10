#!/bin/bash

#SBATCH --job-name=OrthoNR				# Job name.
#SBATCH --output=Orthologs_NR.log			# Standard output and error log.
#SBATCH --partition=short				# Partition (queue)
#SBATCH --ntasks=1					# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=25				# Number of tasks = cpus.
#SBATCH --time=1-00:00:00				# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb				# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
Species_list="car cla cma cme cmo cpe csa lsi mch"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Synteny/nr

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
mkdir -p $WD/08-comparative_genomics/Synteny
mkdir -p $WD/08-comparative_genomics/Synteny/nr
mkdir -p $WD1/02-Orthologs
mkdir -p $WD1/02-Orthologs/Logs


####### PIPELINE

### CREATE PAIRWISE COMBINATIONS FILE
Combinations=$WD1/02-Orthologs/Species_combination.txt
set -- $Species_list
for a; do
    shift
    for b; do
        printf "%s %s\n" "$a" "$b"
    done
done > $Combinations


### IDENTIFY ORTHOLOGS BETWEEN SPECIES
## ORTHOFINDER
echo -e "\n\nOrthologs: OrthoFinder...\n"

cd $WD1/02-Orthologs
if [ -d "$WD1/02-Orthologs/Orthofinder" ]; then
	rm -r Orthofinder
fi
mkdir Orthofinder
cd Orthofinder

# Pairwise
n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb | cut -d" " -f1)
	spe2=$(echo $comb | cut -d" " -f2)
	echo -e "\tCombination: "$spe1"-"$spe2"..."
	mkdir Proteomes_temp
	cp $AI/Proteomes/$spe1.fa ./Proteomes_temp/
	cp $AI/Proteomes/$spe2.fa ./Proteomes_temp/
	orthofinder -f Proteomes_temp -o $spe1\-$spe2 -t $SLURM_CPUS_PER_TASK >> ../Logs/Orthofinder_stdout_$spe1\-$spe2.log 2>&1
	rm -r Proteomes_temp
done

# All
echo -e "\n\tCombination: All..."
mkdir Proteomes_temp
for spe in $Species_list; do
	cp $AI/Proteomes/$spe.fa ./Proteomes_temp/
done
orthofinder -f Proteomes_temp -o All -t $SLURM_CPUS_PER_TASK >> ../Logs/Orthofinder_stdout_all.log 2>&1
rm -r Proteomes_temp


## INPARANOID
echo -e "\n\nOrthologs: Inparanoid...\n"

cd $WD1/02-Orthologs
if [ -d "$WD1/02-Orthologs/Inparanoid" ]; then
	rm -r Inparanoid
fi
mkdir Inparanoid

cd $SP/inparanoid

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb | cut -d" " -f1)
	spe2=$(echo $comb | cut -d" " -f2)
	echo -e "\tCombination: "$spe1"-"$spe2"..."
	perl inparanoid.pl \
		-f1 $AI/Proteomes/$spe1.fa \
		-f2 $AI/Proteomes/$spe2.fa \
		-out-dir $WD1/02-Orthologs/Inparanoid/$spe1\-$spe2/ \
		-out-stats \
		-out-table \
		-out-sqltable \
		-cores $SLURM_CPUS_PER_TASK \
		-cores-diamond $SLURM_CPUS_PER_TASK >> $WD1/02-Orthologs/Logs/Inparanoid_stdout_$spe1\-$spe2.log 2>&1
done



