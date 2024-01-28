#!/bin/bash

####### FUNCTIONS

task_lncRNA_and_PCG_list(){
	# Paths.
	DB=$4/$1/STEP-FINAL/Database/Database_LncRNAs_${6^^}.tsv
	G=$4/$1/STEP-FINAL/Files/Genes
	L=$4/$1/STEP-FINAL/Files/LncRNAs/nr
	FOLDER01=$5/01-LncRNAs_and_Genes/$2/$3
	
	mkdir -p $2
	mkdir -p $2/$3
	
	echo -e "Specie: $1; Confidence-level: $2; Class: $3"
	
	# LncRNA ids file.
	if [[ $3 = "intergenic" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "u" {print $1}' > $FOLDER01/$1\_lncRNA_ids.txt
	elif [[ $3 = "antisense" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "x" {print $1}' > $FOLDER01/$1\_lncRNA_ids.txt
	elif [[ $3 = "intronic" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && $8 == "i" {print $1}' > $FOLDER01/$1\_lncRNA_ids.txt
	elif [[ $3 = "sense" ]]
	then
		tail -n +2 $DB | awk -v var=$2 '$31 == var && ($8 == "o" || $8 == "e") {print $1}' > $FOLDER01/$1\_lncRNA_ids.txt
	else
		tail -n +2 $DB | awk -v var=$2 '$31 == var {print $1}' > $FOLDER01/$1\_lncRNA_ids.txt
	fi
	
	# Gene ids file.
	cp $G/ORIGINAL_GENES_ids.txt $FOLDER01/$1\_gene_ids.txt
	
	# Join lncRNA and gene ids files.
	cat $FOLDER01/$1\_lncRNA_ids.txt $FOLDER01/$1\_gene_ids.txt > $FOLDER01/$1\_ids_unsorted.txt
	
	# LncRNAs gtf file.
	cp $L/POTENTIAL_LNCRNAS_pred.gtf $FOLDER01/$1\_lncRNAs_temp.gtf
	Filter_GTF.py \
		--gtf-initial $FOLDER01/$1\_lncRNAs_temp.gtf \
		--gtf-final $FOLDER01/$1\_lncRNAs.gtf \
		--ids $FOLDER01/$1\_lncRNA_ids.txt
	rm $FOLDER01/$1\_lncRNAs_temp.gtf
	
	# Genes gtf file.
	cp $G/ORIGINAL_GENES.gtf $FOLDER01/$1\_genes.gtf 
	
	# Join lncRNA and gene gtf files.
	cat $FOLDER01/$1\_lncRNAs.gtf $FOLDER01/$1\_genes.gtf > $FOLDER01/$1\_unsorted.gtf
	
	# Sort gtf file.
	gffread $FOLDER01/$1\_unsorted.gtf -T -F --keep-exon-attrs > $FOLDER01/$1\_sorted.gtf
	
	# Extract sorted IDs.
	Rscript $7/Extract_transcript_ids_from_GTF.R $FOLDER01/$1\_sorted.gtf $FOLDER01/$1\_ids_sorted_temp.txt
	
	# Add specie tag to lncRNAs.
	awk -v var=$1 '$1 ~ /^MSTRG/ {print $0"-"var; next} {print}' $FOLDER01/$1\_ids_sorted_temp.txt > $FOLDER01/$1\_ids_sorted.txt
	rm $FOLDER01/$1\_ids_sorted_temp.txt
}

task_synteny_analysis(){
	# Directory.
	mkdir -p $1
	mkdir -p $1/$2
	mkdir -p $1/$2/Outputs
	
	echo -e "Confidence-level: $1; Class: $2"
	
	# Run synteny analysis.
	>$1/$2/Outputs/stdout_$8\_$7.log
	Synteny_analysis.py \
		--path-lists $3/01-LncRNAs_and_Genes/$1/$2 \
		--orthotable $3/02-Orthologs/Tables_orthologs_1_to_1/Orthotable_1_to_1_orthofinder.tsv \
		--output $3/03-Synteny_analysis/$1/$2/output_synteny_table_$8\_$7.tsv \
		--genesNearby $4 \
		--minOverlap $5 \
		--minSideOverlap $6 \
		--nonMatch $7 >> $1/$2/Outputs/stdout_$8\_$7.log 2>&1
}

task_classify_into_families(){
	# Directory.
	mkdir -p $1
	mkdir -p $1/$2
	mkdir -p $1/$2/Outputs
	
	echo -e "Confidence-level: $1; Class: $2"
	
	# Classify into families.
	>$1/$2/Outputs/stdout_$5.log
	Classify_into_families-synteny.py \
		--pred-lncRNAs $3/04-Families/$1/$2/ids_by_specie.tsv \
		--reciprocal-hits $3/03-Synteny_analysis/$1/$2/output_synteny_table_$5.tsv \
		--fam $3/04-Families/$1/$2/fam_$5.tsv \
		--gen $3/04-Families/$1/$2/gen_$5.tsv \
		--threads $4 >> $1/$2/Outputs/stdout_$5.log 2>&1
}

"$@"


