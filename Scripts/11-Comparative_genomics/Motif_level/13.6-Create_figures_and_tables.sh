#!/bin/bash

#SBATCH --job-name=S13mot6							# Job name.
#SBATCH --output=STEP13_motif_6.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=1								# Run on one mode. 
#SBATCH --cpus-per-task=2							# Number of tasks = cpus.
#SBATCH --time=0-08:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=10gb							# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Motif_level/Functions.sh"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High"
strictness_list="ORIGINAL"
nonmatch_list="no"
width_list="6-15 6-50"
mode_list="oops"
n_sim=50

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1/Motif_level


####### PIPELINE

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	## Variable.
	O=$WD1/Motif_level/$flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/06-Figures_and_tables

	### CREATE FIGURES AND TABLES
	echo -e "\n\nCreate figures and tables...\n"
	Rscript $AS/Create_figures_and_tables_Motif.R $O/05-Summary $WD1/Positional_level/$flag/05-Figures_and_tables $O/06-Figures_and_tables "$class_list" "$confidence_list" "$strictness_list" "$nonmatch_list" "$width_list" "$mode_list" $n_sim
	Rscript $AS/Create_figures_and_tables_Motif_additional.R $O/05-Summary $WD1/Positional_level/$flag/05-Figures_and_tables $O/06-Figures_and_tables "$class_list" "$confidence_list" "$strictness_list" "$nonmatch_list" "$width_list" "$mode_list" $n_sim

done


