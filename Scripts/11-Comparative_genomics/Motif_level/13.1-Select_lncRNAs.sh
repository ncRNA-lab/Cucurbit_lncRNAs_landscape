#!/bin/bash

#SBATCH --job-name=S13mot1							# Job name.
#SBATCH --output=STEP13_motif_1.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=9								# Run on one mode. 
#SBATCH --cpus-per-task=2							# Number of tasks = cpus.
#SBATCH --time=0-03:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb							# Job memory request.


####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Motif_level/Functions.sh"
specie_list="cma cme"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High"

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2/Motif_level


####### PIPELINE: STEP 13.1

### ANALYSIS OF CONSERVATION AT MOTIF LEVEL
echo -e "\nANALYSIS OF CONSERVATION AT MOTIF LEVEL..."

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	## Variable.
	O=$WD2/Motif_level/$flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/01-LncRNAs

	### SELECTION
	## For each confidence level and class code, create a TXT file with the identifiers of each lncRNA and a FASTA file with the sequences of the lncRNAs. We have to modify lncRNA IDs by adding the specie.
	cd $O/01-LncRNAs

	echo -e "\n\nSELECTION: Select lncRNAs and modify lncRNA IDs by adding the specie...\n"

	for confidence in $confidence_list; do
		mkdir -p $confidence
		echo -e $confidence"..."
		for class in $class_list; do
			echo -e "\t"$class"..."
			
			# First, clean the directory.
			cd $confidence
			if [ -d "$class" ]; then
				rm -r $class
			fi
			mkdir $class
			cd ..
			
			# Second, create fasta file (Sequences) and txt file (IDs) by specie and 
			# then join them in final fasta and txt files.
			for spe in $specie_list; do
				srun -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Select $spe $confidence $class $WD1 $flag &
			done
			wait
			cat $confidence/$class/*.fasta > $confidence/$class/LncRNAs.fasta
			cat $confidence/$class/*.txt > $confidence/$class/LncRNAs_ids.txt
		done
	done
done

