#!/bin/bash

####### VARIABLES
WD="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results"
AS="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
SP="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
F="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Motif_level/Positional_conserved/Functions_NR.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
#Classes_list="intergenic antisense intronic sense ALL"
#Confidence_levels_list="High Medium Low"
#strictness_list="ORIGINAL RELAXED STRICT MORE-STRICT"
#nonmatch_list="no yes"
Classes_list="antisense"
Confidence_levels_list="High"
strictness_list="ORIGINAL"
nonmatch_list="no"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Positional_level/Approach_2/nr/04-Families
WD2=$WD/08-comparative_genomics/Motif_level/nr/Positional_conserved

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### SOFTWARES
### MEME
#https://meme-suite.org/meme/doc/install.html?man_type=web#quick_src
export MEMEPATH=$SP/meme-5.5.1/meme
export PATH=$PATH:${MEMEPATH}/bin
export PATH=$PATH:${MEMEPATH}/libexec/meme-5.5.1

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Motif_level
mkdir -p $WD/08-comparative_genomics/Motif_level/nr
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/01-LncRNAs
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/02-Preparation
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/03-MotifFinder
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/03-MotifFinder/MEME
mkdir -p $WD/08-comparative_genomics/Motif_level/nr/Positional_conserved/03-MotifFinder/lncLOOM


####### PIPELINE

### SELECTION
## For each confidence level and class code, create a TXT file with the identifiers of each lncRNA and a FASTA file with the sequences of the lncRNAs. We have to modify lncRNA IDs by adding the specie.
cd $WD2/01-LncRNAs

echo -e "\n\nSELECTION: Select lncRNAs and modify lncRNA IDs by adding the specie...\n"

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	echo -e $confidence"..."
	for class in $Classes_list; do
		mkdir -p $confidence/$class
		echo -e "\t"$class"..."
		
		if [[ -f "$confidence/$class/LncRNAs.fasta" ]]; then
		    rm $confidence/$class/LncRNAs.fasta
		fi
		if [[ -f "$confidence/$class/LncRNAs_ids.txt" ]]; then
		    rm $confidence/$class/LncRNAs_ids.txt
		fi
		
		for spe in $Species_list; do
			echo -e "\t\t"$spe"..."
	
			DB=$WD/05-predict_lncRNAs/$spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv
			FASTA=$WD/05-predict_lncRNAs/$spe/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta
			DIR=./$confidence/$class
			
			if [[ $class = "intergenic" ]]
			then
				tail -n +2 $DB | awk -v var=$confidence '$31 == var && $8 == "u" {print $1}' > $DIR/$spe\_ids_temp.txt
				seqtk subseq $FASTA $DIR/$spe\_ids_temp.txt > $DIR/$spe\_temp.fasta
				awk -v var=$spe '{print $0"-"var}' $DIR/$spe\_ids_temp.txt > $DIR/$spe\_ids.txt
				rm $DIR/$spe\_ids_temp.txt
				awk -F '>' -v var=$spe '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$spe\_temp.fasta > $DIR/$spe.fasta
				rm $DIR/$spe\_temp.fasta
			elif [[ $class = "antisense" ]]
			then
				tail -n +2 $DB | awk -v var=$confidence '$31 == var && $8 == "x" {print $1}' > $DIR/$spe\_ids_temp.txt
				seqtk subseq $FASTA $DIR/$spe\_ids_temp.txt > $DIR/$spe\_temp.fasta
				awk -v var=$spe '{print $0"-"var}' $DIR/$spe\_ids_temp.txt > $DIR/$spe\_ids.txt
				rm $DIR/$spe\_ids_temp.txt
				awk -F '>' -v var=$spe '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$spe\_temp.fasta > $DIR/$spe.fasta
				rm $DIR/$spe\_temp.fasta
			elif [[ $class = "intronic" ]]
			then
				tail -n +2 $DB | awk -v var=$confidence '$31 == var && $8 == "i" {print $1}' > $DIR/$spe\_ids_temp.txt
				seqtk subseq $FASTA $DIR/$spe\_ids_temp.txt > $DIR/$spe\_temp.fasta
				awk -v var=$spe '{print $0"-"var}' $DIR/$spe\_ids_temp.txt > $DIR/$spe\_ids.txt
				rm $DIR/$spe\_ids_temp.txt
				awk -F '>' -v var=$spe '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$spe\_temp.fasta > $DIR/$spe.fasta
				rm $DIR/$spe\_temp.fasta
			elif [[ $class = "sense" ]]
			then
				tail -n +2 $DB | awk -v var=$confidence '$31 == var && ($8 == "o" || $8 == "e") {print $1}' > $DIR/$spe\_ids_temp.txt
				seqtk subseq $FASTA $DIR/$spe\_ids_temp.txt > $DIR/$spe\_temp.fasta
				awk -v var=$spe '{print $0"-"var}' $DIR/$spe\_ids_temp.txt > $DIR/$spe\_ids.txt
				rm $DIR/$spe\_ids_temp.txt
				awk -F '>' -v var=$spe '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$spe\_temp.fasta > $DIR/$spe.fasta
				rm $DIR/$spe\_temp.fasta
			else
				tail -n +2 $DB | awk -v var=$confidence '$31 == var {print $1}' > $DIR/$spe\_ids_temp.txt
				seqtk subseq $FASTA $DIR/$spe\_ids_temp.txt > $DIR/$spe\_temp.fasta
				awk -v var=$spe '{print $0"-"var}' $DIR/$spe\_ids_temp.txt > $DIR/$spe\_ids.txt
				rm $DIR/$spe\_ids_temp.txt
				awk -F '>' -v var=$spe '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$spe\_temp.fasta > $DIR/$spe.fasta
				rm $DIR/$spe\_temp.fasta
			fi
		done
		cat $confidence/$class/*.fasta > $confidence/$class/LncRNAs.fasta
		cat $confidence/$class/*.txt > $confidence/$class/LncRNAs_ids.txt
	done
done

### PREPARE FAMILY FILES FOR MOTIF LEVEL ANALYSIS
## For each nonmatch level, strictness level, confidence level and class code, create a fasta file by conserved family.
cd $WD2/02-Preparation

echo -e "\n\nPREPARATION: Create a fasta file by conserved family...\n"

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
				echo -e "\t\t\t"$class"..."
	
				cd $WD2/02-Preparation/$strictness/$nonmatch/$confidence/
				if [ -d "$WD2/02-Preparation/$strictness/$nonmatch/$confidence/$class" ]; then
					rm -r $class
				fi
				mkdir $class
				mkdir $class/outputs
				cd $WD2/02-Preparation
				
				# Prepare directory to execute motif level analysis.
				>$strictness/$nonmatch/$confidence/$class/outputs/stdout.log
				Prepare_motif_level_analysis.py \
					--input-table $WD1/$confidence/$class/gen_$strictness\_$nonmatch.tsv \
					--input-fasta $WD2/01-LncRNAs/$confidence/$class/LncRNAs.fasta \
					--path-out $WD2/02-Preparation/$strictness/$nonmatch/$confidence/$class \
					--n-iter 50 \
					--n-proc 10 >> $strictness/$nonmatch/$confidence/$class/outputs/stdout.log 2>&1
			done
		done
	done
	
done


### MOTIF FINDER (MEME)
## For each nonmatch level, strictness level, confidence level and class code, find the different motif and make an enrichment analysis.
cd $WD2/03-MotifFinder/MEME

echo -e "\n\nMOTIF FINDER (MEME): Find motivs...\n"

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
				echo -e "\t\t\t"$class"..."
				
				DIR_A="$WD2/02-Preparation/$strictness/$nonmatch/$confidence/$class"
				DIR_B="$WD2/03-MotifFinder/MEME/$strictness/$nonmatch/$confidence/$class"
	
				mkdir -p $DIR_B
				cd $DIR_B
				if [ -d "$DIR_B/outputs" ]; then
					rm -r outputs
				fi
				mkdir outputs
				
				# REAL.
				if [ -d "$DIR_B/real" ]; then
					rm -r real
				fi
				mkdir real
				
				families=$(ls $DIR_A/real/ | grep ".fasta")
				>$DIR_B/outputs/stdout_REAL.log
				for fam in $families; do 
					name=$(echo $fam | cut -f 1 -d '.')
					echo -e $name
					#-revcomp \ In my opinion, it doesn't make sense search for motifs in the reverse strand because lncRNAs are transcripts (RNA) although the alphabet be dna here. It could make sense if you study lncRNA promoters.
					meme \
						$DIR_A/real/$fam \
						-oc $DIR_B/real/$name \
						-nmotifs 1 \
						-mod oops \
						-evt 0.05 \
						-dna \
						-nostatus \
						-p 20		
				done >> $DIR_B/outputs/stdout_REAL.log 2>&1
				
				# SIMULATIONS.
				if [ -d "$DIR_B/simulations" ]; then
					rm -r simulations
				fi
				mkdir simulations
				
				iterations=$(ls -1 $DIR_A/simulations/ | wc -l)
				for i in $(seq $iterations); do
					>$DIR_B/outputs/stdout_SIMULATION_$i.log
					families=$(ls $DIR_A/simulations/iter_$i/ | grep ".fasta")
					for fam in $families; do 
						name=$(echo $fam | cut -f 1 -d '_')
						echo -e $name
						#-revcomp \ In my opinion, it doesn't make sense search for motifs in the reverse strand because lncRNAs are transcripts (RNA) although the alphabet be dna here. It could make sense if you study lncRNA promoters.
						meme \
							$DIR_A/simulations/iter_$i/$fam \
							-oc $DIR_B/simulations/iter_$i/$name \
							-nmotifs 1 \
							-mod oops \
							-evt 0.05 \
							-dna \
							-nostatus \
							-p 20	
					done >> $DIR_B/outputs/stdout_SIMULATION_$i.log 2>&1
				done 
			done
		done
	done
	
done

<<comment
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
				srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_lncLOOM $strictness $nonmatch $confidence $class $WD1 $WD2 &
			done
			wait
		done
	done
	
done
comment

