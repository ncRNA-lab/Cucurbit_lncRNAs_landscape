#!/bin/bash

#SBATCH --job-name=cmeS10Int						# Job name.
#SBATCH --output=cme_STEP10_Intersect.log				# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=6							# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-05:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb						# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/07-Get_intergenic_regions"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-TEs_and_genomic_repeats"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/08-TEs_and_genomic_repeats/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/08-TEs_and_genomic_repeats/Functions.sh"
confidence_list="High Medium Low"
flag_list="nr"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie
WD3_spe=$WD3/$specie

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD3
mkdir -p $WD3_spe
mkdir -p $WD3_spe/02-Intersection
mkdir -p $WD3_spe/02-Intersection/Intersection
mkdir -p $WD3_spe/02-Intersection/Final_tables
mkdir -p $WD3_spe/02-Intersection/Outputs


####### PIPELINE

### INTERSECTION DIFFERENT CLASSES (PCG, LNCRNA AND INTERGENIC REGIONS) WITH REPETITIVE REGIONS
echo -e "\nINTERSECTION DIFFERENT CLASSES (PCG, LNCRNA AND INTERGENIC REGIONS) WITH REPETITIVE REGIONS..."

## Convert RepeatMasker output to bed format and intersect genes, lncRNAs and intergenic regions with the repetitive regions found.
echo -e "\nIntersect genes, lncRNAs and intergenic regions with the repetitive regions found..."
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Intersect $specie $WD3_spe $WD1_spe $WD2_spe

## Create final tables.
echo -e "\nCreate the final tables..."
for flag in $flag_list; do
	for confidence in $confidence_list; do
		 srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Final_tables $specie $WD3_spe/02-Intersection $WD3_spe/01-Repeat_calling/02-RepeatMasker $WD1_spe $AS $flag $confidence &
	done
done
wait


