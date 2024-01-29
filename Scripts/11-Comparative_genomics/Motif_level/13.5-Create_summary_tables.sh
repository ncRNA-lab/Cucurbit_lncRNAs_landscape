#!/bin/bash

#SBATCH --job-name=S13mot5							# Job name.
#SBATCH --output=STEP13_motif_5.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=1								# Run on one mode. 
#SBATCH --cpus-per-task=10							# Number of tasks = cpus.
#SBATCH --time=0-08:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb							# Job memory request.


####### MODULES
module load anaconda

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Motif_level/Functions.sh"
specie_list="cma cme"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High"
strictness_list="ORIGINAL"
nonmatch_list="no"
width_list="6-15 6-50"
mode_list="oops"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2/Motif_level


####### PIPELINE: STEP 13.5

### ANALYSIS OF CONSERVATION AT MOTIF LEVEL
echo -e "\nANALYSIS OF CONSERVATION AT MOTIF LEVEL..."

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	## Variable.
	O=$WD2/Motif_level/$flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/05-Summary
	mkdir -p $O/05-Summary/LNCRNAS_DB
	mkdir -p $O/05-Summary/MEME
	mkdir -p $O/05-Summary/GOMO

	### GLOBAL LNCRNAS DATABASE
	## Generate a global database of lncRNAs.
	echo -e "\n\nGenerate a global database of lncRNAs...\n"

	cd $O/05-Summary/LNCRNAS_DB

	>Database_LncRNAs_${flag^^}.tsv
	n_species=$(echo $specie_list | wc -w)
	for i in $(seq $n_species); do
		spe=$(echo $specie_list | cut -d ' ' -f $i)
		if [[ "$i" == 1 ]]
		then
			head -n 1 $WD1/$spe/STEP-FINAL/Database/Database_LncRNAs_${flag^^}.tsv >> Database_LncRNAs_${flag^^}.tsv
			tail -n +2 $WD1/$spe/STEP-FINAL/Database/Database_LncRNAs_${flag^^}.tsv | awk -v var=$spe '{print $1"-"var"\t"$0}' | awk '{$0=gensub(/\s*\S+/,"",2)}1' >> Database_LncRNAs_${flag^^}.tsv
		else
			tail -n +2 $WD1/$spe/STEP-FINAL/Database/Database_LncRNAs_${flag^^}.tsv | awk -v var=$spe '{print $1"-"var"\t"$0}' | awk '{$0=gensub(/\s*\S+/,"",2)}1' >> Database_LncRNAs_${flag^^}.tsv
		fi
	done

	### SUMMARY TABLE MEME
	## For each nonmatch level, strictness level, confidence level and class code, generate a summary table using MEME results.
	echo -e "\n\nSUMMARY TABLE: Generate summary table using MEME  results...\n"

	cd $O/05-Summary/MEME

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
					
							DIR_A="$O/02-Preparation/$strictness/$nonmatch/$confidence/$class"
							DIR_B="$O/03-MotifFinder/MEME/$strictness/$nonmatch/$confidence/$class/$mode/$width"
							DIR_C="$O/05-Summary/MEME/$strictness/$nonmatch/$confidence/$class/$mode/$width"
							
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

	cd $O/05-Summary/GOMO

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
					
							DIR_A="$O/04-MotifEnrichment/GOMO/$strictness/$nonmatch/$confidence/$class/$mode/$width"
							DIR_B="$O/05-Summary/GOMO/$strictness/$nonmatch/$confidence/$class/$mode/$width"
							
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
done


