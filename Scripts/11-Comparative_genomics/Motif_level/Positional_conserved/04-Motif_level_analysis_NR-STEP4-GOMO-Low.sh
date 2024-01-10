#!/bin/bash

#SBATCH --job-name=PosNR4GL							# Job name.
#SBATCH --output=Positional_NR_4_GOMO_Low.log					# Standard output and error log.
#SBATCH --qos=long-mem								# Partition (queue)
#SBATCH --ntasks=80								# Run on one mode. 
#SBATCH --cpus-per-task=1							# Number of tasks = cpus.
#SBATCH --time=2-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb							# Job memory request.


####### MODULES
module load meme/5.5.1/serial

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Motif_level/Positional_conserved/Functions_NR.sh"
MEME_databases="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares/MEME_databases"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="Low"
#strictness_list="ORIGINAL RELAXED STRICT MORE-STRICT"
#nonmatch_list="no yes"
strictness_list="ORIGINAL"
nonmatch_list="no"
widths="6-15 6-50"
modes="oops"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Motif_level/nr/Positional_conserved

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Motif_level
mkdir -p $WD/08-comparative_genomics/Motif_level/nr
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/04-MotifEnrichment
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/04-MotifEnrichment/GOMO


####### PIPELINE

### MOTIF ENRICHMENT (GOMO)
## For each nonmatch level, strictness level, confidence level and class code, determine if any motif is significantly associated with genes linked to one or more Genome Ontology (GO) terms.
cd $WD1/04-MotifEnrichment/GOMO

echo -e "\n\nMOTIF ENRICHMENT (GOMO): Determine if any motif is significantly associated with genes linked to one or more Genome Ontology (GO) terms...\n"

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
				for mode in $modes; do
					mkdir -p $strictness/$nonmatch/$confidence/$class/$mode
					echo -e "\t\t\t\t"$mode"..."
					for width in $widths; do
						mkdir -p $strictness/$nonmatch/$confidence/$class/$mode/$width
						echo -e "\t\t\t\t\t"$width"..."
						
						DIR_A="$WD1/03-MotifFinder/MEME/$strictness/$nonmatch/$confidence/$class/$mode/$width"
						DIR_B="$WD1/04-MotifEnrichment/GOMO/$strictness/$nonmatch/$confidence/$class/$mode/$width"
						
						cd $DIR_B
						
						# OUTPUTS
						# First, clean the directory.
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
						
						# Second, execute gomo by lncRNA family.
						echo -e "\t\t\t\t\t\tREAL"
						>$DIR_B/outputs/stdout_REAL.log
						srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_GOMO_2 $DIR_A/real $DIR_B/real $MEME_databases/gomo_databases $DIR_B/outputs/stdout_REAL.log &
						
						# SIMULATIONS.
						# First, clean the directory.
						if [ -d "$DIR_B/simulations" ]; then
							rm -r simulations
						fi
						mkdir simulations
						
						# Second, execute gomo by lncRNA family in each iteration.
						echo -e "\t\t\t\t\t\tSIMULATIONS"
						iterations=$(ls -1 $DIR_A/simulations/ | wc -l)
						for i in $(seq $iterations); do
							mkdir simulations/iter_$i
							>$DIR_B/outputs/stdout_SIMULATION_$i.log
							srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_GOMO_2 $DIR_A/simulations/iter_$i $DIR_B/simulations/iter_$i $MEME_databases/gomo_databases $DIR_B/outputs/stdout_SIMULATION_$i.log &
						done
						
						cd $WD1/04-MotifEnrichment/GOMO
					done
				done
			done
		done
	done
done

wait

