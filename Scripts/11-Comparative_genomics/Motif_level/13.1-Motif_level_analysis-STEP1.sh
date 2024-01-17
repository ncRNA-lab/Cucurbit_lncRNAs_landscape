#!/bin/bash

#SBATCH --job-name=PosNR1							# Job name.
#SBATCH --output=Positional_NR_1.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=9								# Run on one mode. 
#SBATCH --cpus-per-task=2							# Number of tasks = cpus.
#SBATCH --time=0-03:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb							# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Motif_level/Positional_conserved/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Positional_level/Approach_2/nr/04-Families
WD2=$WD/08-comparative_genomics/Motif_level/nr/Positional_conserved

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Motif_level
mkdir -p $WD/08-comparative_genomics/Motif_level/nr
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/01-LncRNAs


####### PIPELINE

### SELECTION
## For each confidence level and class code, create a TXT file with the identifiers of each lncRNA and a FASTA file with the sequences of the lncRNAs. We have to modify lncRNA IDs by adding the specie.
cd $WD2/01-LncRNAs

echo -e "\n\nSELECTION: Select lncRNAs and modify lncRNA IDs by adding the specie...\n"

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		echo -e "\t"$class"..."
		
		# First, clean the directory.
		cd $confidence
		if [ -d "$WD2/01-LncRNAs/$confidence/$class" ]; then
			rm -r $class
		fi
		cd ..
		mkdir $confidence/$class
		
		# Second, create fasta file (Sequences) and txt file (IDs) by specie and then join them in final fasta and txt files.
		for spe in $Species_list; do
			srun -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Select $spe $confidence $class $WD &
		done
		wait
		cat $confidence/$class/*.fasta > $confidence/$class/LncRNAs.fasta
		cat $confidence/$class/*.txt > $confidence/$class/LncRNAs_ids.txt
	done
done


