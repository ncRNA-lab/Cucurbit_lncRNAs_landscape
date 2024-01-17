#!/bin/bash

#SBATCH --job-name=PosNR5MG							# Job name.
#SBATCH --output=Positional_NR_5_MEME-GOMO.log				# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=1								# Run on one mode. 
#SBATCH --cpus-per-task=50							# Number of tasks = cpus.
#SBATCH --time=0-08:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=140mb							# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Motif_level/Positional_conserved/Functions_NR.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High"
strictness_list="ORIGINAL"
nonmatch_list="no"
widths="6-15 6-50"
modes="oops"

####### NEW AND OTHER VARIABLES
WD1=$WD/05-predict_lncRNAs
WD2=$WD/08-comparative_genomics/Motif_level/nr/Positional_conserved

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Motif_level
mkdir -p $WD/08-comparative_genomics/Motif_level/nr
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/05-Summary
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/05-Summary/LNCRNAS_DB
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/05-Summary/MEME
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/05-Summary/GOMO


####### PIPELINE

### GLOBAL LNCRNAS DATABASE
## Generate a global database of lncRNAs.
echo -e "\n\nGenerate a global database of lncRNAs...\n"

cd $WD2/05-Summary/LNCRNAS_DB

>Database_LncRNAs_NR.tsv
n_species=$(echo $Species_list | wc -w)
for i in $(seq $n_species); do
	spe=$(echo $Species_list | cut -d ' ' -f $i)
	if [[ "$i" == 1 ]]
	then
		head -n 1 $WD1/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv >> Database_LncRNAs_NR.tsv
		tail -n +2 $WD1/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$spe '{print $1"-"var"\t"$0}' | awk '{$0=gensub(/\s*\S+/,"",2)}1' >> Database_LncRNAs_NR.tsv
	else
		tail -n +2 $WD1/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$spe '{print $1"-"var"\t"$0}' | awk '{$0=gensub(/\s*\S+/,"",2)}1' >> Database_LncRNAs_NR.tsv
	fi
done

### SUMMARY TABLE MEME
## For each nonmatch level, strictness level, confidence level and class code, generate a summary table using MEME results.
echo -e "\n\nSUMMARY TABLE: Generate summary table using MEME  results...\n"

cd $WD2/05-Summary/MEME

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
				
						DIR_A="$WD2/02-Preparation/$strictness/$nonmatch/$confidence/$class"
						DIR_B="$WD2/03-MotifFinder/MEME/$strictness/$nonmatch/$confidence/$class/$mode/$width"
						DIR_C="$WD2/05-Summary/MEME/$strictness/$nonmatch/$confidence/$class/$mode/$width"
						
						if [[ "$mode" == "oops" ]]
						then
							Generate_summary_oops_MEME_results.py \
								--path-meme $DIR_B \
								--path-fastas $DIR_A \
								--path-out $DIR_C \
								--n-proc $SLURM_CPUS_PER_TASK
						fi
					done
				done
			done
		done
	done
done

### SUMMARY TABLE GOMO
## For each nonmatch level, strictness level, confidence level and class code, generate a summary table using GOMO results.
echo -e "\n\nSUMMARY TABLE: Generate summary table using GOMO results...\n"

cd $WD2/05-Summary/GOMO

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
				
						DIR_A="$WD2/04-MotifEnrichment/GOMO/$strictness/$nonmatch/$confidence/$class/$mode/$width"
						DIR_B="$WD2/05-Summary/GOMO/$strictness/$nonmatch/$confidence/$class/$mode/$width"
						
						if [[ "$mode" == "oops" ]]
						then
							Generate_summary_GOMO_results.py \
								--path-gomo $DIR_A \
								--path-out $DIR_B \
								--n-proc $SLURM_CPUS_PER_TASK
						fi
					done
				done
			done
		done
	done
done


