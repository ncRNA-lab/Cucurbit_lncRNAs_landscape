#!/bin/bash

#SBATCH --job-name=S13mot4							# Job name.
#SBATCH --output=STEP13_motif_4.log						# Standard output and error log.
#SBATCH --qos=medium								# Partition (queue)
#SBATCH --ntasks=20								# Run on one mode. 
#SBATCH --cpus-per-task=2							# Number of tasks = cpus.
#SBATCH --time=6-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb							# Job memory request.


####### MODULES
module load meme/5.5.1/serial

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Motif_level/Functions.sh"
MEME_databases="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares/MEME_databases"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High"
strictness_list="ORIGINAL"
nonmatch_list="no"
width_list="6-15 6-50"
mode_list="oops"

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
	mkdir -p $O/04-MotifEnrichment
	mkdir -p $O/04-MotifEnrichment/GOMO

	### MOTIF ENRICHMENT (GOMO)
	## For each nonmatch level, strictness level, confidence level and class code, determine if any motif is significantly associated with genes linked to one or more Genome Ontology (GO) terms.
	cd $O/04-MotifEnrichment/GOMO

	echo -e "\n\nMOTIF ENRICHMENT (GOMO): Determine if any motif is significantly associated with genes linked to one or more Genome Ontology (GO) terms...\n"

	for strictness in $strictness_list; do
		mkdir -p $strictness
		echo -e $strictness"..."
		for nonmatch in $nonmatch_list; do
			mkdir -p $strictness/$nonmatch
			echo -e "\t"$nonmatch"..."
			for confidence in $confidence_list; do
				mkdir -p $strictness/$nonmatch/$confidence
				echo -e "\t\t"$confidence"..."
				for class in $class_list; do
					mkdir -p $strictness/$nonmatch/$confidence/$class
					echo -e "\t\t\t"$class"..."
					for mode in $mode_list; do
						mkdir -p $strictness/$nonmatch/$confidence/$class/$mode
						echo -e "\t\t\t\t"$mode"..."
						for width in $width_list; do
							mkdir -p $strictness/$nonmatch/$confidence/$class/$mode/$width
							echo -e "\t\t\t\t\t"$width"..."
							
							DIR_A="$O/03-MotifFinder/MEME/$strictness/$nonmatch/$confidence/$class/$mode/$width"
							DIR_B="$O/04-MotifEnrichment/GOMO/$strictness/$nonmatch/$confidence/$class/$mode/$width"
							
							cd $DIR_B
							
							# OUTPUTS
							# First, clean the directory.
							if [ -d "Outputs" ]; then
								rm -r Outputs
							fi
							mkdir Outputs
							
							# REAL.
							# First, clean the directory.
							if [ -d "real" ]; then
								rm -r real
							fi
							mkdir real
							
							# Second, execute gomo by lncRNA family.
							echo -e "\t\t\t\t\t\tREAL"
							>$DIR_B/Outputs/stdout_REAL.log
							srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_GOMO_2 $DIR_A/real $DIR_B/real $MEME_databases/gomo_databases $DIR_B/Outputs/stdout_REAL.log &
							
							# SIMULATIONS.
							# First, clean the directory.
							if [ -d "simulations" ]; then
								rm -r simulations
							fi
							mkdir simulations
							
							# Second, execute gomo by lncRNA family in each iteration.
							echo -e "\t\t\t\t\t\tSIMULATIONS"
							iterations=$(ls -1 $DIR_A/simulations/ | wc -l)
							for i in $(seq $iterations); do
								mkdir simulations/iter_$i
								>$DIR_B/Outputs/stdout_SIMULATION_$i.log
								srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_GOMO_2 $DIR_A/simulations/iter_$i $DIR_B/simulations/iter_$i $MEME_databases/gomo_databases $DIR_B/Outputs/stdout_SIMULATION_$i.log &
							done
							
							cd $O/04-MotifEnrichment/GOMO
						done
					done
				done
			done
		done
	done
	wait
done

