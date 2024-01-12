#!/bin/bash

#SBATCH --job-name=cmeS10Comp						# Job name.
#SBATCH --output=cme_STEP10_Compare.log				# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=6							# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-05:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb						# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/08-TEs_and_genomic_repeats/Functions.sh"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/08-TEs_and_genomic_repeats/Additional_scripts"
Specie="cme"
confidences="High Medium Low"
flags="NR R"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD/08-TEs_and_genomic_repeats/$Specie
WD2_spe=$WD/05-LncRNAs_prediction/$Specie
WD3_spe=$WD/07-Get_intergenic_regions/$Specie

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD/08-TEs_and_genomic_repeats
mkdir -p $WD1_spe
mkdir -p $WD1_spe/02-Comparison_PCGs_LncRNAs


####### PIPELINE

## Convert RepeatMasker output to bed format and intersect genes, lncRNAs and intergenic regions with the repetitive regions found.
echo -e "\n\nIntersect genes, lncRNAs and intergenic regions with the repetitive regions found...\n"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Intersect $Specie $WD1_spe $WD2_spe $WD3_spe

## Get masked genome percentages.
echo -e "\n\nGet masked genome percentages...\n"
Get_percentage_of_masked_genome.py \
	--path $WD1_spe/01-Repeat_calling/02-RepeatMasker \
	--specie $Specie

## Create final tables.
echo -e "\n\nCreate the final tables...\n"
mkdir -p $WD1_spe/02-Comparison_PCGs_LncRNAs/Final_tables
for flag in $flags; do
	for confidence in $confidences; do
		 srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Final_tables $Specie $WD1_spe/02-Comparison_PCGs_LncRNAs $WD1_spe/01-Repeat_calling/02-RepeatMasker $WD2_spe $AS $flag $confidence &
	done
done
wait


