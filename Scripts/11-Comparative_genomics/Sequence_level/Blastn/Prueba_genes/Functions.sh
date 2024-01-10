#!/bin/bash

####### FUNCTIONS

task_Select(){
	DIR=$2/05-predict_lncRNAs/$1/STEP-FINAL/Files/Genes
	
	cp $DIR/ORIGINAL_GENES_ids.txt ./$1\_ids.txt
	seqtk subseq $DIR/ORIGINAL_GENES.fasta $1\_ids.txt > $1.fasta
	$4/Create_temporal_gene_ids.py \
		--spe $1 \
		--ids $3/01-Genes/$1\_ids.txt \
		--fasta $3/01-Genes/$1.fasta \
		--new-ids $3/01-Genes/$1\_ids\_new.txt \
		--new-fasta $3/01-Genes/$1\_new.fasta \
		--database $3/01-Genes/$1\_database.tsv \
		--mode 'Blastn'
	rm $1.fasta $1\_ids.txt
}

task_Makeblastdb(){
	mkdir -p outputs
	
	>outputs/$1\_stdout.log
	makeblastdb \
		-in ../01-Genes/$1\_new.fasta \
		-dbtype nucl \
		-parse_seqids \
		-out $1 \
		-title $1 >> outputs/$1\_stdout.log 2>&1
}

task_Blastn(){
	mkdir -p outputs
	
	comb_spe=$(cat $1 | head -n $4 | tail -n 1)
	spe1=$(echo $comb_spe | cut -d" " -f1)
	spe2=$(echo $comb_spe | cut -d" " -f2)
	
	# First possibility.
	>outputs/$spe1\_to\_$spe2\_stdout.log
	blastn \
		-query ../01-Genes/$spe1\_new.fasta \
		-db ../02-Makeblastdb/$spe2 \
		-task blastn \
		-num_threads $3 \
		-out $spe1\_to\_$spe2\_temp.tsv \
		-strand plus \
		-outfmt 6 \
		-evalue $2 >> outputs/$spe1\_to\_$spe2\_stdout.log 2>&1
	# Now, it's necessary to choose the best hit. First, we will choose the best hsp found for each qseqid-sseqid alignment taking into account the e-value,
	# pidet and length of the alignment. Then, we will choose the best hit found for each qseqid taking into account also the e-value, pidet and length of 
	# the alignment.
	Rscript $6/Choose_x_best_blast_hits.R $spe1 $spe2 $5/03-Blastn/$spe1\_to\_$spe2\_temp.tsv $5/03-Blastn/$spe1\_to\_$spe2.tsv 1
	rm $spe1\_to\_$spe2\_temp.tsv
	
	# Second possibility.
	>outputs/$spe2\_to\_$spe1\_stdout.log
	blastn \
		-query ../01-Genes/$spe2\_new.fasta \
		-db ../02-Makeblastdb/$spe1 \
		-task blastn \
		-num_threads $3 \
		-out $spe2\_to\_$spe1\_temp.tsv \
		-strand plus \
		-outfmt 6 \
		-evalue $2 >> outputs/$spe2\_to\_$spe1\_stdout.log 2>&1
	# Now, it's necessary to choose the best hit. First, we will choose the best hsp found for each qseqid-sseqid alignment taking into account the e-value,
	# pidet and length of the alignment. Then, we will choose the best hit found for each qseqid taking into account also the e-value, pidet and length of 
	# the alignment.
	Rscript $6/Choose_x_best_blast_hits.R $spe1 $spe2 $5/03-Blastn/$spe2\_to\_$spe1\_temp.tsv $5/03-Blastn/$spe2\_to\_$spe1.tsv 1
	rm $spe2\_to\_$spe1\_temp.tsv
}

task_Reciprocal_hits(){
	comb_spe=$(cat $1 | head -n $4 | tail -n 1)
	spe1=$(echo $comb_spe | cut -d" " -f1)
	spe2=$(echo $comb_spe | cut -d" " -f2)
	
	$3/Find_reciprocal_hits.py \
		--spe1 $spe1 \
		--spe2 $spe2 \
		--spe1vsspe2 $2/03-Blastn/$spe1\_to\_$spe2.tsv \
		--spe2vsspe1 $2/03-Blastn/$spe2\_to\_$spe1.tsv \
		--output $2/04-Reciprocal_hits/$spe1\_and\_$spe2.tsv
}

"$@"

