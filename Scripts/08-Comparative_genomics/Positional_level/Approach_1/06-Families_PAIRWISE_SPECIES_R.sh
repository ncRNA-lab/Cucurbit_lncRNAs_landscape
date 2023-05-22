#!/bin/bash

#SBATCH --job-name=FamPR							# Job name.
#SBATCH --output=Family_PAIRWISE_SPECIES_R.log				# Standard output and error log.
#SBATCH --partition=short							# Partition (queue)
#SBATCH --ntasks=1								# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2							# Number of tasks = cpus. It depends on the number of process of your parallelization.
#SBATCH --time=0-10:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem=10gb								# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
Species_list="car cla cma cme cmo cpe csa lsi mch"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High Medium Low"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Synteny/r

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD1/05-Families
mkdir -p $WD1/05-Families/PAIRWISE_SPECIES


####### PIPELINE

### FAMILIES
## Classify lncRNAs into families.
cd $WD1/05-Families/PAIRWISE_SPECIES

echo -e "\n\nClassify lncRNAs into families...\n"

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		mkdir -p $confidence/$class
		echo -e "\t"$class"..."
		>./$confidence/$class/ids_by_specie.tsv
		for spe in $Species_list; do
			awk -v a="$spe" '{print a"\t"$0}' $WD1/01-LncRNAs/$confidence/$class/$spe\_ids.txt >> ./$confidence/$class/ids_by_specie.tsv
		done
		
		# Classify into families.
		$AS/Classify_into_families.py \
			--pred-lncRNAs $WD1/05-Families/PAIRWISE_SPECIES/$confidence/$class/ids_by_specie.tsv \
			--reciprocal-hits $WD1/04-LncRNAs_inside_syntenic_blocks/PAIRWISE_SPECIES/CloudHits-$class\-$confidence.tsv \
			--fam $WD1/05-Families/PAIRWISE_SPECIES/$confidence/$class/fam.tsv \
			--gen $WD1/05-Families/PAIRWISE_SPECIES/$confidence/$class/gen.tsv
	done
done



