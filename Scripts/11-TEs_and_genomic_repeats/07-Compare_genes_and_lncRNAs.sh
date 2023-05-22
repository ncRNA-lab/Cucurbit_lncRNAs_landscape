#!/bin/bash

#SBATCH --job-name=Comp					# Job name.
#SBATCH --output=Compare.log					# Standard output and error log.
#SBATCH --partition=short					# Partition (queue)
#SBATCH --ntasks=3						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb					# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/11-TEs_and_genomic_repeats/Functions.sh"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/11-TEs_and_genomic_repeats/additional_scripts"
Species_list="car cla cma cme cmo cpe csa lsi mch"
confidences="High Medium Low"
flags="NR R"

####### NEW AND OTHER VARIABLES
WD1=$WD/11-TEs_and_genomic_repeats
WD2=$WD/05-predict_lncRNAs
WD3=$WD/07-comparison_lncRNAs_vs_coding_genes

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD/11-TEs_and_genomic_repeats
mkdir -p $WD1/02-Comparison_Genes_LncRNAs
mkdir -p $WD1/02-Comparison_Genes_LncRNAs/Figures_and_Tables
mkdir -p $WD1/02-Comparison_Genes_LncRNAs/Figures_and_Tables/A
mkdir -p $WD1/02-Comparison_Genes_LncRNAs/Figures_and_Tables/B
mkdir -p $WD1/02-Comparison_Genes_LncRNAs/Figures_and_Tables/C

####### PIPELINE

## Convert RepeatMasker output to bed format.
echo -e "\n\nConvert RepeatMasker output to bed format...\n"
for spe in $Species_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_intersect $spe $WD1 $WD2 $WD3 &
done
wait

## Get masked genome percentages.
echo -e "\n\nGet masked genome percentages...\n"
Get_percentage_of_masked_genome.py \
	--path $WD1/01-Repeat_calling/05-RepeatMasker \
	--species $Species_list

## Create the Figures and Tables.
echo -e "\n\nCreate the Figures and Tables...\n"
for flag in $flags; do
	for confidence in $confidences; do
		 srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Figures_and_tables $WD $WD1 "02-Comparison_Genes_LncRNAs" "05-RepeatMasker" $flag $confidence $AS &
	done
	wait
done



