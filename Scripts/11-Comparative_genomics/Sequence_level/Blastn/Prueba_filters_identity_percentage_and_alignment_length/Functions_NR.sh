#!/bin/bash

####### FUNCTIONS

task_Select(){
	DB=$4/05-predict_lncRNAs/$1/STEP-FINAL/Database/Database_LncRNAs_NR.tsv
	FASTA=$4/05-predict_lncRNAs/$1/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta
	DIR=./$2/$3
	
	if [[ $3 = "intergenic" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "u" {print $1}' > $DIR/$1\_ids.txt
		seqtk subseq $FASTA $DIR/$1\_ids.txt > $DIR/$1.fasta
	elif [[ $3 = "antisense" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "x" {print $1}' > $DIR/$1\_ids.txt
		seqtk subseq $FASTA $DIR/$1\_ids.txt > $DIR/$1.fasta
	elif [[ $3 = "intronic" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "i" {print $1}' > $DIR/$1\_ids.txt
		seqtk subseq $FASTA $DIR/$1\_ids.txt > $DIR/$1.fasta
	elif [[ $3 = "sense" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && ($8 == "o" || $8 == "e") {print $1}' > $DIR/$1\_ids.txt
		seqtk subseq $FASTA $DIR/$1\_ids.txt > $DIR/$1.fasta
	else
		tail -n +2 $DB | awk -v var=$2 '$31 == var {print $1}' > $DIR/$1\_ids.txt
		seqtk subseq $FASTA $DIR/$1\_ids.txt > $DIR/$1.fasta
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
	mkdir -p $3/$4
	mkdir -p $3/$4/outputs
	
	# First possibility.
	echo -e "Combination ($7/$8): $1 - $2; Confidence-level: $3; Class: $4"
	>$3/$4/outputs/$1\_to\_$2\_stdout.log
	blastn \
		-query ../01-LncRNAs/$3/$4/$1.fasta \
		-db ../02-Makeblastdb/$3/$4/$2 \
		-task blastn \
		-num_threads $6 \
		-out $3/$4/$1\_to\_$2\_temp1.tsv \
		-strand plus \
		-outfmt 6 \
		-evalue $5 >> $3/$4/outputs/$1\_to\_$2\_stdout.log 2>&1
	# Filter by pident (>= 20%) and length (>= 50).
	awk '$3 >= 20 && $4 >= 50 {print $0}' $3/$4/$1\_to\_$2\_temp1.tsv > $3/$4/$1\_to\_$2\_temp2.tsv
	rm $3/$4/$1\_to\_$2\_temp1.tsv
	# Now, it's necessary to choose the best hit. First, we will choose the best hsp found for each qseqid-sseqid alignment taking into account the e-value,
	# pidet and length of the alignment. Then, we will choose the best hit found for each qseqid taking into account also the e-value, pidet and length of 
	# the alignment.
	Rscript ${10}/Choose_x_best_blast_hits.R $1 $2 $9/03-Blastn/$3/$4/$1\_to\_$2\_temp2.tsv $9/03-Blastn/$3/$4/$1\_to\_$2.tsv 1
	rm $3/$4/$1\_to\_$2\_temp2.tsv
	
	# Second possibility.
	echo -e "Combination ($7/$8): $2 - $1; Confidence-level: $3; Class: $4"
	>$3/$4/outputs/$2\_to\_$1\_stdout.log
	blastn \
		-query ../01-LncRNAs/$3/$4/$2.fasta \
		-db ../02-Makeblastdb/$3/$4/$1 \
		-task blastn \
		-num_threads $6 \
		-out $3/$4/$2\_to\_$1\_temp1.tsv \
		-strand plus \
		-outfmt 6 \
		-evalue $5 >> $3/$4/outputs/$2\_to\_$1\_stdout.log 2>&1
	# Filter by pident (>= 20%) and length (>= 50).
	awk '$3 >= 20 && $4 >= 50 {print $0}' $3/$4/$2\_to\_$1\_temp1.tsv > $3/$4/$2\_to\_$1\_temp2.tsv
	rm $3/$4/$2\_to\_$1\_temp1.tsv
	# Now, it's necessary to choose the best hit. First, we will choose the best hsp found for each qseqid-sseqid alignment taking into account the e-value,
	# pidet and length of the alignment. Then, we will choose the best hit found for each qseqid taking into account also the e-value, pidet and length of 
	# the alignment.
	Rscript ${10}/Choose_x_best_blast_hits.R $1 $2 $9/03-Blastn/$3/$4/$2\_to\_$1\_temp2.tsv $9/03-Blastn/$3/$4/$2\_to\_$1.tsv 1
	rm $3/$4/$2\_to\_$1\_temp2.tsv
}

task_Reciprocal_hits(){
	mkdir -p $3/$4
	
	echo -e "Combination ($7/$8): $1 - $2 / $2 - $1; Confidence-level: $3; Class: $4..."
	$6/Find_reciprocal_hits.py \
		--spe1 $1 \
		--spe2 $2 \
		--spe1vsspe2 $5/03-Blastn/$3/$4/$1\_to\_$2.tsv \
		--spe2vsspe1 $5/03-Blastn/$3/$4/$2\_to\_$1.tsv \
		--output $5/04-Reciprocal_hits/$3/$4/$1\_and\_$2.tsv
}

task_Classify_into_families(){
	mkdir -p $1/$2/outputs
	
	# Classify into families.
	>$1/$2/outputs/stdout.log
	$4/Classify_into_families.py \
		--pred-lncRNAs $3/05-Families/$1/$2/ids_by_specie.tsv \
		--reciprocal-hits $3/04-Reciprocal_hits/$1/$2/Reciprocal_hits.tsv \
		--fam $3/05-Families/$1/$2/fam.tsv \
		--gen $3/05-Families/$1/$2/gen.tsv \
		--threads $5 >> $1/$2/outputs/stdout.log 2>&1
}

"$@"

