#!/bin/bash

#SBATCH --job-name=ILPNR							# Job name.
#SBATCH --output=Identify_LncRNAs_PAIRWISE_SPECIES_NR.log			# Standard output and error log.
#SBATCH --partition=short							# Partition (queue)
#SBATCH --ntasks=36								# Run on one mode.
#SBATCH --cpus-per-task=2							# Number of tasks = cpus.
#SBATCH --time=0-10:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=5gb							# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Synteny/Functions.sh"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High Medium Low"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Synteny/nr

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD1/04-LncRNAs_inside_syntenic_blocks
mkdir -p $WD1/04-LncRNAs_inside_syntenic_blocks/PAIRWISE_SPECIES


####### PIPELINE: PAIRWISE_SPECIES
echo -e "\n\nPAIRWISE_SPECIES..."

## STEP 6: Identify LncRNAs inside syntenic blocks
echo -e "\n\nSTEP 6: Identify LncRNAs inside syntenic blocks...\n"
cd $WD1/04-LncRNAs_inside_syntenic_blocks/PAIRWISE_SPECIES

Combinations=$WD1/03-Adhore/PAIRWISE_SPECIES/Species_combination.txt

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb | cut -d" " -f1)
	spe2=$(echo $comb | cut -d" " -f2)
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_identify_LncRNAs_inside_syntenic_blocks $spe1 $spe2 $WD1/03-Adhore/PAIRWISE_SPECIES $WD1/04-LncRNAs_inside_syntenic_blocks/PAIRWISE_SPECIES $WD/05-predict_lncRNAs 'nr' $AS &
done
wait

for confidence in $Confidence_levels_list; do
	for class in $Classes_list; do
		cat */output_cloud/CloudHits-$class\-$confidence.tsv > CloudHits-$class\-$confidence.tsv
	done
done


