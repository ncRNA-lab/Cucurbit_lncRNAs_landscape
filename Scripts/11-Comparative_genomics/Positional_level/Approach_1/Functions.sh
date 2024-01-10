#!/bin/bash

####### FUNCTIONS
# Extract orthologs 1:1 pairwise.
task_extract_orthologs_pairwise(){
	if [ -d "$3/03-Adhore/PAIRWISE_SPECIES/$1-$2" ]; then
		rm -r $1\-$2
	fi	
	mkdir $1\-$2
	$4/Extract_1_to_1_orthologs_pairwise-Approach_1.py \
		--spe1 $1 \
		--spe2 $2 \
		--path-in $3/02-Orthologs \
		--path-out $3/03-Adhore/PAIRWISE_SPECIES/$1\-$2
}

task_extract_info_GFF3_genes(){
	echo -e "\tCombination: "$1"-"$2"..."
	gffread $3/GFF3_genes/$1.gff3 -F --keep-exon-attrs > $4/$1\-$2/$1.gff3
	gffread $3/GFF3_genes/$2.gff3 -F --keep-exon-attrs > $4/$1\-$2/$2.gff3
	Rscript $5/Extract_info_genes_GFF3.R $1 $2 $4/$1\-$2
}

task_identify_LncRNAs_inside_syntenic_blocks(){
	echo -e "\tCombination: "$1"-"$2"..."
	mkdir -p $1\-$2
	if [ -d "./$1-$2/output_cloud" ]; then
		rm -r ./$1-$2/output_cloud
	fi	
	mkdir ./$1-$2/output_cloud
	$7/Identify_LncRNAs_inside_syntenic_blocks.py \
		--spe1 $1 \
		--spe2 $2 \
		--path1 $5 \
		--path2 $3/$1\-$2 \
		--path3 $4/$1\-$2 \
		--flag $6 >> $4/$1\-$2/output_cloud/stdout.log 2>&1
}

"$@"
