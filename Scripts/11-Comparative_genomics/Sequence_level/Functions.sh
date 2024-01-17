#!/bin/bash

####### FUNCTIONS

task_Select(){
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

task_Makeblastdb(){
	mkdir -p $2/$3
	mkdir -p $2/$3/outputs
	
	>$2/$3/outputs/$1\_stdout.log
	makeblastdb \
		-in ../01-LncRNAs/$2/$3/$1.fasta \
		-dbtype nucl \
		-parse_seqids \
		-out ./$2/$3/$1 \
		-title $1 >> $2/$3/outputs/$1\_stdout.log 2>&1
}

task_Blastn(){	
	mkdir -p $5/$6
	mkdir -p $5/$6/outputs
	echo -e "Combination ($9/${10}): $1 - $2; Confidence-level: $5; Class: $6"
	>$5/$6/outputs/$3\_to\_$4\_stdout.log
	blastn \
		-query ../01-LncRNAs/$5/$6/$1.fasta \
		-db ../02-Makeblastdb/$5/$6/$2 \
		-task blastn \
		-num_threads $8 \
		-out $5/$6/Blast$3\_$4\_temp1.txt \
		-strand plus \
		-outfmt 6 \
		-evalue $7 >> $5/$6/outputs/$3\_to\_$4\_stdout.log 2>&1
	# Now, it's necessary to choose the best hit. First, we will choose the best hsp found for each qseqid-sseqid alignment taking into account the e-value,
	# pidet and length of the alignment. Then, we will choose the best hit found for each qseqid taking into account also the e-value, pidet and length of 
	# the alignment. 
	# As OrthoFinder also requires the results of the blastn of a species against itself but we do not want any paralogues, in those cases in which the blastn 
	# execution is done against the same species, those hits in the blastn table that are between different lncRNAs will be removed to avoid the presence of 
	# any paralog in the results.
	Rscript ${12}/Choose_x_best_blast_hits.R $3 $4 ${11}/03-Blastn/$5/$6/Blast$3\_$4\_temp1.txt ${11}/03-Blastn/$5/$6/Blast$3\_$4\_temp2.txt 1
}

task_OrthoFinder(){
	mkdir -p $1/$2
	
	# Directories.
	DIR1=$3/01-LncRNAs/$1/$2
	DIR2=$3/03-Blastn/$1/$2
	DIR3=$3/04-OrthoFinder/$1/$2
	
	# First, run orthofinder with the '-op' option. This option will prepare the files in the format required by OrthoFinder and print the set of BLAST commands that need to be run. This is useful if you want to manage the BLAST searches yourself.
	if [ -d "$DIR3/Prepare_blast" ]; then
		rm -r $DIR3/Prepare_blast
	fi
	>$DIR3/stdout_prepare_blast.log
	orthofinder -f $DIR1 -d -o $DIR3/Prepare_blast -op -t $6 >> $DIR3/stdout_prepare_blast.log
	
	# Second, copy the blastn results coming from the 03-Blastn step into the folder created (Prepare_blast/$folder_p/WorkingDirectory/).
	folder_p=$(ls $DIR3/Prepare_blast) 
	cp $DIR2/*_temp2.txt $DIR3/Prepare_blast/$folder_p/WorkingDirectory/
	
	# Third, modify lncRNA ids again. In this case, we have to convert lncRNA IDs to OrthoFinder format. We use a table created in the first step of 04-OrthoFinder step.
	$4/Convert_lncRNA_ID_to_OrthoFinder_ID.py \
		--path $DIR3/Prepare_blast/$folder_p/WorkingDirectory \
		--comb $5
	rm $DIR3/Prepare_blast/$folder_p/WorkingDirectory/*_temp2.txt
	
	# Fourth, compress the blastn files.
	n_comb=$(wc -l $5 | cut -d" " -f1)
	for j in `seq 1 $n_comb`; do
		comb_cod=$(cat $5 | head -n $j | tail -n 1)
		cod1=$(echo $comb_cod | cut -d" " -f1)
		cod2=$(echo $comb_cod | cut -d" " -f2)
		gzip -c $DIR3/Prepare_blast/$folder_p/WorkingDirectory/Blast$cod1\_$cod2.txt > $DIR3/Prepare_blast/$folder_p/WorkingDirectory/Blast$cod1\_$cod2.txt.gz
		rm $DIR3/Prepare_blast/$folder_p/WorkingDirectory/Blast$cod1\_$cod2.txt
	done
	
	# Fifth, execute OrthoFinder with the '-og' opion to cluster the lncRNAs into families. The MCL clustering is used to cluster the blastn reciprocal results. OrthoFinder automatically converts the OrthoFinder lncRNA IDs to the initial lncRNA IDs.
	if [ -d "$DIR3/Clustering" ]; then
		rm -r $DIR3/Clustering
	fi
	>$DIR3/stdout_clustering.log
	orthofinder -b $DIR3/Prepare_blast/$folder_p/WorkingDirectory/ -d -og -t $6 >> $DIR3/stdout_clustering.log
	cp -R $DIR3/Prepare_blast/$folder_p/WorkingDirectory/OrthoFinder $DIR3/
	mv $DIR3/OrthoFinder $DIR3/Clustering
	
	# Sixth, liberate space.
	cd $DIR3/Prepare_blast/$folder_p/WorkingDirectory/
	rm -r OrthoFinder
	folder_c=$(ls $DIR3/Clustering) 
	cd $DIR3/Clustering/$folder_c/
	rm -r Orthogroup_Sequences
	rm -r Single_Copy_Orthologue_Sequences
}

task_Get_families(){
	mkdir -p $1/$2/outputs
	
	folder_c=$(ls $3/04-OrthoFinder/$1/$2/Clustering)
	# Get families.
	>$1/$2/outputs/stdout.log
	$4/Get_families_from_OrthoFinder.py \
		--pred-lncRNAs $3/05-Families/$1/$2/ids_by_specie.tsv \
		--orthofinder $3/04-OrthoFinder/$1/$2/Clustering/$folder_c/Orthogroups/Orthogroups.txt \
		--fam $3/05-Families/$1/$2/fam.tsv \
		--gen $3/05-Families/$1/$2/gen.tsv >> $1/$2/outputs/stdout.log 2>&1
}

"$@"

