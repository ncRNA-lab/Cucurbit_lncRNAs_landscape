#!/bin/bash

#SBATCH --job-name=PosNR3L							# Job name.
#SBATCH --output=Positional_NR_3_lncLOOM.log					# Standard output and error log.
#SBATCH --qos=long								# Partition (queue)
#SBATCH --ntasks=90								# Run on one mode. 
#SBATCH --cpus-per-task=1							# Number of tasks = cpus.
#SBATCH --time=15-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb							# Job memory request.


####### MODULES
source activate LncLOOMv2

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Motif_level/Positional_conserved/Functions_NR.sh"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High Medium Low"
#strictness_list="ORIGINAL RELAXED STRICT MORE-STRICT"
#nonmatch_list="no yes"
strictness_list="ORIGINAL"
nonmatch_list="no"
widths="6-15 6-50"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Positional_level/Approach_2/nr/04-Families
WD2=$WD/08-comparative_genomics/Motif_level/nr/Positional_conserved

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Motif_level
mkdir -p $WD/08-comparative_genomics/Motif_level/nr
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/03-MotifFinder
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/03-MotifFinder/lncLOOM


####### PIPELINE

### MOTIF FINDER (lncLOOM)
## For each nonmatch level, strictness level, confidence level and class code, find the different motif and make an enrichment analysis.
cd $WD2/03-MotifFinder/lncLOOM

echo -e "\n\nMOTIF FINDER (lncLOOM): Find motivs...\n"

for strictness in $strictness_list; do
	mkdir -p $strictness
	echo -e $strictness"..."
	for nonmatch in $nonmatch_list; do
		mkdir -p $strictness/$nonmatch
		echo -e "\t"$nonmatch"..."
		for confidence in $Confidence_levels_list; do
			mkdir -p $strictness/$nonmatch/$confidence
			echo -e "\t\t"$confidence"..."
			for class in $Classes_list; do
				mkdir -p $strictness/$nonmatch/$confidence/$class
				echo -e "\t\t\t"$class"..."
				for width in $widths; do
					mkdir -p $strictness/$nonmatch/$confidence/$class/$width
					echo -e "\t\t\t\t"$width"..."
					
					DIR_A="$WD2/02-Preparation/$strictness/$nonmatch/$confidence/$class"
					DIR_B="$WD2/03-MotifFinder/lncLOOM/$strictness/$nonmatch/$confidence/$class/$width"
					
					cd $DIR_B
					
					if [ -d "$DIR_B/outputs" ]; then
						rm -r outputs
					fi
					mkdir outputs
					
					# REAL.
					# First, clean the directory.
					if [ -d "$DIR_B/real" ]; then
						rm -r real
					fi
					mkdir real
					
					# Second, execute lncLOOM by lncRNA family.
					echo -e "\t\t\t\t\tREAL"
					families=$(ls $DIR_A/real/ | grep ".fasta")
					>$DIR_B/outputs/real_stdout.log
					for fam in $families; do 
						srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_lncLOOM_1 $fam $DIR_A/real $DIR_B/real $SLURM_CPUS_PER_TASK $width $DIR_B/outputs/real_stdout.log &
					done
					wait
					
					# SIMULATIONS.
					# First, clean the directory.
					if [ -d "$DIR_B/simulations" ]; then
						rm -r simulations
					fi
					mkdir simulations
					
					# Second, execute lncLOOM by lncRNA family in each iteration.
					iterations=$(ls -1 $DIR_A/simulations/ | wc -l)
					for i in $(seq $iterations); do
						echo -e "\t\t\t\t\tSIMULATION ($i/$iterations)"
						mkdir simulations/iter_$i
						families=$(ls $DIR_A/simulations/iter_$i/ | grep ".fasta")
						for fam in $families; do 
							srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_lncLOOM_1 $fam $DIR_A/simulations/iter_$i $DIR_B/simulations/iter_$i $SLURM_CPUS_PER_TASK $width $DIR_B/outputs/simulations_stdout.log &
						done
						wait
					done
					
					cd $WD2/03-MotifFinder/lncLOOM
				done
			done
		done
	done
done


