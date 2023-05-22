#!/bin/bash

#SBATCH --job-name=LNR					# Job name.
#SBATCH --output=LncRNAs_NR.log			# Standard output and error log.
#SBATCH --partition=short				# Partition (queue)
#SBATCH --ntasks=1					# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=4				# Number of tasks = cpus.
#SBATCH --time=1-00:00:00				# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb				# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
Species_list="car cla cma cme cmo cpe csa lsi mch"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High Medium Low"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Synteny/nr

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Synteny
mkdir -p $WD/08-comparative_genomics/Synteny/nr
mkdir -p $WD1/01-LncRNAs


####### PIPELINE

### SELECTION
## Select lncRNAs by class code.
cd $WD1/01-LncRNAs

echo -e "\n\nSelect lncRNAs by class code and confidence level...\n"

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		mkdir -p $confidence/$class
		echo -e "\t"$class"..."
		for spe in $Species_list; do
			echo -e "\t\t"$spe"..."
			if [[ $class = "intergenic" ]]
			then
				tail -n +2 $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$confidence '$31 == var && $8 == "u" {print $1}' > ./$confidence/intergenic/$spe\_ids.txt
				seqtk subseq $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta ./$confidence/intergenic/$spe\_ids.txt > ./$confidence/intergenic/$spe.fasta
			elif [[ $class = "antisense" ]]
			then
				tail -n +2 $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$confidence '$31 == var && $8 == "x" {print $1}' > ./$confidence/antisense/$spe\_ids.txt
				seqtk subseq $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta ./$confidence/antisense/$spe\_ids.txt > ./$confidence/antisense/$spe.fasta
			elif [[ $class = "intronic" ]]
			then
				tail -n +2 $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$confidence '$31 == var && $8 == "i" {print $1}' > ./$confidence/intronic/$spe\_ids.txt
				seqtk subseq $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta ./$confidence/intronic/$spe\_ids.txt > ./$confidence/intronic/$spe.fasta
			elif [[ $class = "sense" ]]
			then
				tail -n +2 $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$confidence '$31 == var && ($8 == "o" || $8 == "e") {print $1}' > ./$confidence/sense/$spe\_ids.txt
				seqtk subseq $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta ./$confidence/sense/$spe\_ids.txt > ./$confidence/sense/$spe.fasta
			else
				tail -n +2 $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$confidence '$31 == var {print $1}' > ./$confidence/ALL/$spe\_ids.txt
				seqtk subseq $WD/05-predict_lncRNAs/$spe/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta ./$confidence/ALL/$spe\_ids.txt > ./$confidence/ALL/$spe.fasta
			fi
		done
	done
done

