#!/bin/bash

####### FUNCTIONS

task_transcript_density(){
	
	## Create directory.
	cd $3
	if [ -d "$6" ]; then
		rm -r "$6"
	fi
	mkdir $6
	cd $6
	
	## Get chromosome sizes.
	cp $4/Genome/$1.fa ./
	samtools faidx $1.fa
	cut -f1,2 $1.fa.fai > $1.sizes_genome_temp.txt
	awk '{print $1"\t"0"\t"$2}' $1.sizes_genome_temp.txt > $1.sizes_genome.txt
	rm $1.sizes_genome_temp.txt $1.fa $1.fa.fai

	## NON-REDUNDANT: Visualize the lncRNAs and genes distribution.
	LncRNAs_tab="$2/STEP-FINAL/Database/Database_LncRNAs_${6^^}.tsv"
	Genes_tab="$2/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv"
	Rscript $5/Genome-wide_distribution.R $1 $LncRNAs_tab $Genes_tab $3/$6 $4
}

"$@"


