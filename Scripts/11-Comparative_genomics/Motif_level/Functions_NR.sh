#!/bin/bash

####### FUNCTIONS

task_Select(){
	echo -e "\t\t"$1"..."
	
	DB=$4/05-predict_lncRNAs/$1/STEP-FINAL/Database/Database_LncRNAs_NR.tsv
	FASTA=$4/05-predict_lncRNAs/$1/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta
	DIR=./$2/$3
	
	if [[ $3 = "intergenic" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "u" {print $1}' > $DIR/$1\_ids_temp.txt
		seqtk subseq $FASTA $DIR/$1\_ids_temp.txt > $DIR/$1\_temp.fasta
		awk -v var=$1 '{print $0"-"var}' $DIR/$1\_ids_temp.txt > $DIR/$1\_ids.txt
		rm $DIR/$1\_ids_temp.txt
		awk -F '>' -v var=$1 '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$1\_temp.fasta > $DIR/$1.fasta
		rm $DIR/$1\_temp.fasta
	elif [[ $3 = "antisense" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "x" {print $1}' > $DIR/$1\_ids_temp.txt
		seqtk subseq $FASTA $DIR/$1\_ids_temp.txt > $DIR/$1\_temp.fasta
		awk -v var=$1 '{print $0"-"var}' $DIR/$1\_ids_temp.txt > $DIR/$1\_ids.txt
		rm $DIR/$1\_ids_temp.txt
		awk -F '>' -v var=$1 '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$1\_temp.fasta > $DIR/$1.fasta
		rm $DIR/$1\_temp.fasta
	elif [[ $3 = "intronic" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "i" {print $1}' > $DIR/$1\_ids_temp.txt
		seqtk subseq $FASTA $DIR/$1\_ids_temp.txt > $DIR/$1\_temp.fasta
		awk -v var=$1 '{print $0"-"var}' $DIR/$1\_ids_temp.txt > $DIR/$1\_ids.txt
		rm $DIR/$1\_ids_temp.txt
		awk -F '>' -v var=$1 '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$1\_temp.fasta > $DIR/$1.fasta
		rm $DIR/$1\_temp.fasta
	elif [[ $3 = "sense" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && ($8 == "o" || $8 == "e") {print $1}' > $DIR/$1\_ids_temp.txt
		seqtk subseq $FASTA $DIR/$1\_ids_temp.txt > $DIR/$1\_temp.fasta
		awk -v var=$1 '{print $0"-"var}' $DIR/$1\_ids_temp.txt > $DIR/$1\_ids.txt
		rm $DIR/$1\_ids_temp.txt
		awk -F '>' -v var=$1 '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$1\_temp.fasta > $DIR/$1.fasta
		rm $DIR/$1\_temp.fasta
	else
		tail -n +2 $DB | awk -v var=$2 '$31 == var {print $1}' > $DIR/$1\_ids_temp.txt
		seqtk subseq $FASTA $DIR/$1\_ids_temp.txt > $DIR/$1\_temp.fasta
		awk -v var=$1 '{print $0"-"var}' $DIR/$1\_ids_temp.txt > $DIR/$1\_ids.txt
		rm $DIR/$1\_ids_temp.txt
		awk -F '>' -v var=$1 '$2~/^MSTRG/ { print $0"-"var; next }1' $DIR/$1\_temp.fasta > $DIR/$1.fasta
		rm $DIR/$1\_temp.fasta
	fi
}

task_Prepare_motif_level_analysis(){
	echo -e "\t\t\t"$4"..."
	
	# First, clean the directory.
	cd $6/02-Preparation/$1/$2/$3/
	if [ -d "$6/02-Preparation/$1/$2/$3/$4" ]; then
		rm -r $4
	fi
	mkdir $4
	mkdir $4/outputs
	cd $6/02-Preparation
	
	# Second, prepare directory to execute motif level analysis.
	>$1/$2/$3/$4/outputs/stdout.log
	Prepare_motif_level_analysis.py \
		--input-table $5/$3/$4/gen_$1\_$2.tsv \
		--input-fasta $6/01-LncRNAs/$3/$4/LncRNAs.fasta \
		--path-out $6/02-Preparation/$1/$2/$3/$4 \
		--n-iter $7 \
		--n-proc $8 >> $1/$2/$3/$4/outputs/stdout.log 2>&1
}

task_MEME_1(){
	name=$(echo $1 | cut -f 1 -d '.')
	min=$(echo $5 | cut -f 1 -d '-')
	max=$(echo $5 | cut -f 2 -d '-')
	meme $2/$1 -oc $3/$name -mod $4 -maxw $max -minw $min -evt 0.05 -dna -nostatus
}

task_MEME_2(){
	families=$(ls $1 | grep ".fasta")
	for fam in $families; do
		name=$(echo $fam | cut -f 1 -d '.')
		min=$(echo $4 | cut -f 1 -d '-')
		max=$(echo $4 | cut -f 2 -d '-')
		echo -e "$name" >> $5
		meme $1/$fam -oc $2/$name -mod $3 -maxw $max -minw $min -evt 0.05 -dna -nostatus
	done
}

task_lncLOOM_1(){
	name=$(echo $1 | cut -f 1 -d '.')
	min=$(echo $5 | cut -f 1 -d '-')
	max=$(echo $5 | cut -f 2 -d '-')
	LncLOOM --fasta $2/$1 --startw $max --stopw $min --iterations 100 --outdir $3 --pname $name --multiprocess $4 >> $6 2>&1
}

task_lncLOOM_2(){
	families=$(ls $1 | grep ".fasta")
	for fam in $families; do
		name=$(echo $fam | cut -f 1 -d '.')
		min=$(echo $4 | cut -f 1 -d '-')
		max=$(echo $4 | cut -f 2 -d '-')
		echo -e "\n\n\n$name\n\n" >> $5
		LncLOOM --fasta $1/$fam --startw $max --stopw $min --iterations 100 --outdir $2 --pname $name --multiprocess $3 >> $5 2>&1
	done	
}

task_GOMO_1(){
	if grep -q "SUMMARY OF MOTIFS" $1/meme.txt; then
		mkdir $2
		ama $1/meme.txt $3/plant_arabidopsis_1000_199.na $3/plant_arabidopsis_1000_199.na.bfile --oc $2/ama_ath --pvalues >> $4 2>&1
		ama $1/meme.txt $3/plant_oryza_sativa_1000_199.na $3/plant_oryza_sativa_1000_199.na.bfile --oc $2/ama_osa --pvalues >> $4 2>&1
		ama $1/meme.txt $3/plant_populus_trichocarpa_1000_199.na $3/plant_populus_trichocarpa_1000_199.na.bfile --oc $2/ama_ptr --pvalues >> $4 2>&1
		ama $1/meme.txt $3/plant_sorghum_brachupodium_1000_199.na $3/plant_sorghum_brachupodium_1000_199.na.bfile --oc $2/ama_sbr --pvalues >> $4 2>&1
		ama $1/meme.txt $3/plant_brachypodium_1000_199.na $3/plant_brachypodium_1000_199.na.bfile --oc $2/ama_bra --pvalues >> $4 2>&1
		gomo --oc $2/gomo_out --t 0.05 --shuffle_scores 1000 --nostatus --dag $3/go.dag --motifs $1/meme.txt $3/plant_arabidopsis_1000_199.na.csv $2/ama_ath/ama.xml $2/ama_osa/ama.xml $2/ama_ptr/ama.xml $2/ama_sbr/ama.xml $2/ama_bra/ama.xml >> $4 2>&1
	fi
}

task_GOMO_2(){
	families=$(ls $1)
	for fam in $families; do
		if grep -q "SUMMARY OF MOTIFS" $1/$fam/meme.txt; then
			mkdir $2/$fam
			ama $1/$fam/meme.txt $3/plant_arabidopsis_1000_199.na $3/plant_arabidopsis_1000_199.na.bfile --oc $2/$fam/ama_ath --pvalues >> $4 2>&1
			ama $1/$fam/meme.txt $3/plant_oryza_sativa_1000_199.na $3/plant_oryza_sativa_1000_199.na.bfile --oc $2/$fam/ama_osa --pvalues >> $4 2>&1
			ama $1/$fam/meme.txt $3/plant_populus_trichocarpa_1000_199.na $3/plant_populus_trichocarpa_1000_199.na.bfile --oc $2/$fam/ama_ptr --pvalues >> $4 2>&1
			ama $1/$fam/meme.txt $3/plant_sorghum_brachupodium_1000_199.na $3/plant_sorghum_brachupodium_1000_199.na.bfile --oc $2/$fam/ama_sbr --pvalues >> $4 2>&1
			ama $1/$fam/meme.txt $3/plant_brachypodium_1000_199.na $3/plant_brachypodium_1000_199.na.bfile --oc $2/$fam/ama_bra --pvalues >> $4 2>&1
			gomo --oc $2/$fam/gomo_out --t 0.05 --shuffle_scores 1000 --nostatus --dag $3/go.dag --motifs $1/$fam/meme.txt $3/plant_arabidopsis_1000_199.na.csv $2/$fam/ama_ath/ama.xml $2/$fam/ama_osa/ama.xml $2/$fam/ama_ptr/ama.xml $2/$fam/ama_sbr/ama.xml $2/$fam/ama_bra/ama.xml >> $4 2>&1
		fi
	done
}

"$@"

