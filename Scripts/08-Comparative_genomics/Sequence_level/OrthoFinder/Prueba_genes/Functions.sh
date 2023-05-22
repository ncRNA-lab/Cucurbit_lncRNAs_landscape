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
		--mode 'OrthoFinder'
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
	
	comb_spe=$(cat $1 | head -n $5 | tail -n 1)
	spe1=$(echo $comb_spe | cut -d" " -f1)
	spe2=$(echo $comb_spe | cut -d" " -f2)
	comb_cod=$(cat $2 | head -n $5 | tail -n 1)
	cod1=$(echo $comb_cod | cut -d" " -f1)
	cod2=$(echo $comb_cod | cut -d" " -f2)
	
	>outputs/$cod1\_to\_$cod2\_stdout.log
	blastn \
		-query ../01-Genes/$spe1\_new.fasta \
		-db ../02-Makeblastdb/$spe2 \
		-task blastn \
		-num_threads $4 \
		-out Blast$cod1\_$cod2\_temp1.txt \
		-strand plus \
		-outfmt 6 \
		-evalue $3 >> outputs/$cod1\_to\_$cod2\_stdout.log 2>&1
	# Now, it's necessary to choose the best hit. First, we will choose the best hsp found for each qseqid-sseqid alignment taking into account the e-value,
	# pidet and length of the alignment. Then, we will choose the best hit found for each qseqid taking into account also the e-value, pidet and length of 
	# the alignment. 
	# As OrthoFinder also requires the results of the blastn of a species against itself but we do not want any paralogues, in those cases in which the blastn 
	# execution is done against the same species, those hits in the blastn table that are between different lncRNAs will be removed to avoid the presence of 
	# any paralog in the results.
	Rscript $7/Choose_x_best_blast_hits.R $cod1 $cod2 $6/03-Blastn/Blast$cod1\_$cod2\_temp1.txt $6/03-Blastn/Blast$cod1\_$cod2\_temp2.txt 1
}

"$@"

