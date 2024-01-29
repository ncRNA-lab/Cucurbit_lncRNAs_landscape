#!/bin/bash

#SBATCH --job-name=S13pos5					# Job name.
#SBATCH --output=STEP13_positional_5.log			# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=1						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2					# Number of tasks = cpus.
#SBATCH --time=0-02:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb					# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Positional_level/Functions.sh"
specie_list="car cla cma cme cmo cpe csa lsi mch"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High Medium Low"
strictness_list="ORIGINAL RELAXED STRICT MORE-STRICT"
nonmatch_list="no yes"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1/Positional_level


####### PIPELINE: STEP 13.5

### ANALYSIS OF CONSERVATION AT POSITIONAL LEVEL
echo -e "\nANALYSIS OF CONSERVATION AT POSITIONAL LEVEL..."

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	## Variable.
	O=$WD1/Positional_level/$flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/05-Figures_and_tables

	## CREATE FINAL FIGURES AND TABLES
	cd $O/05-Figures_and_tables
	
	echo -e "\n\nCreate figures and tables..."
	
	Rscript $AS/Create_figures_and_tables_Synteny.R $O/04-Families $O/05-Figures_and_tables "$specie_list" "$class_list" "$confidence_list" "$strictness_list" "$nonmatch_list"
	
done


