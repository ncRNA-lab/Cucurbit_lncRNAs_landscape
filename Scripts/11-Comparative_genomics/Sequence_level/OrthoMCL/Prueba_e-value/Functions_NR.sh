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
	mkdir -p $3/$4
	mkdir -p $3/$4/outputs
	echo -e "Combination ($7/$8): $1 - $2; Confidence-level: $3; Class: $4"
	>$3/$4/outputs/$1\_to\_$2\_stdout.log
	blastn \
		-query ../01-LncRNAs/$3/$4/$1.fasta \
		-db ../02-Makeblastdb/$3/$4/$2 \
		-task blastn \
		-num_threads $6 \
		-max_target_seqs 1 \
		-max_hsps 1 \
		-out $3/$4/$1\_to\_$2\_temp1.tsv \
		-strand plus \
		-outfmt 6 \
		-evalue $5 >> $3/$4/outputs/$1\_to\_$2\_stdout.log 2>&1
	# Now, it's necessary to choose the best hit. First, we will choose the best hsp found for each qseqid-sseqid alignment taking into account the e-value,
	# pidet and length of the alignment. Then, we will choose the best hit found for each qseqid taking into account also the e-value, pidet and length of 
	# the alignment. 
	# As OrthoMCL also requires the results of the blastn of a species against itself but we do not want any paralogues, in those cases in which the blastn 
	# execution is done against the same species, those hits in the blastn table that are between different lncRNAs will be removed to avoid the presence of 
	# any paralog in the results.
	Rscript ${10}/Choose_x_best_blast_hits.R $1 $2 $9/03-Blastn/$3/$4/$1\_to\_$2\_temp1.tsv $9/03-Blastn/$3/$4/$1\_to\_$2\_temp2.tsv 1
}

task_OrthoMCL(){
	mkdir -p $2/$3
	
	Species_list="$(echo $1 | sed 's/-/ /g')"
	
	cd $2/$3
	>Repo_spec.txt
	>Repo_spec.txt.all.GFF3
	>stdout.log
	for spe1 in $Species_list; do
		mkdir -p $spe1
		
		# Directories
		DIR1=$5/05-predict_lncRNAs/$spe1/STEP-FINAL/Files/LncRNAs/nr
		DIR2=$4/03-Blastn/$2/$3
		DIR3=$4/04-OrthoMCL/$2/$3/$spe1
		
		# Store in this folder all blastn tables in which the specie has been used as spe1 in 3-blastn and rename the files.
		if [ -d "$DIR3/RBH_blast_CDS" ]; then
			rm -r $DIR3/RBH_blast_CDS
		fi
		mkdir $DIR3/RBH_blast_CDS
		cp $DIR2/$spe1\_to* $DIR3/RBH_blast_CDS/
		cd $spe1/RBH_blast_CDS
		rm *_temp1.tsv
		
		# Rename files to OrthoMCL format name.
		files="$(ls $DIR3/RBH_blast_CDS/)"
		for file in $files; do
			spe2="$(echo $file | cut -d_ -f3)"
			mv $spe1\_to\_$spe2\_temp2.tsv vs_$spe2.blast.m8
		done
		cd ../../
		
		# Create Repo_spec.txt.all.GFF3.
		cp $DIR1/POTENTIAL_LNCRNAS_pred.gtf $DIR3/$spe1.annotation.gtf
		gffread $DIR3/$spe1.annotation.gtf -o $DIR3/$spe1.annotation.gff3
		tail -n +4 $DIR3/$spe1.annotation.gff3 | cut -d";" -f1 | \
		awk '$3 == "transcript" {print $0}' | sed -e 's/ID=//g' | \
		sed -e "s/StringTie/$spe1/g" | sed -e 's/transcript/mRNA/g' | \
		awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"-"$2}' >> Repo_spec.txt.all.GFF3
		rm $DIR3/$spe1.annotation.gff3
		rm $DIR3/$spe1.annotation.gtf
		
		# Create the Repo_spec.txt.
		echo -e "//" >> Repo_spec.txt
		echo -e "Genome $spe1" >> Repo_spec.txt
		echo -e "Annotation $spe1" >> Repo_spec.txt
	done
	echo -e "//" >> Repo_spec.txt
	
	# Execute OrthoMCL using Synima to cluster the lncRNAs into families using the blastn reciprocal results.
	if [ -d "OMCL_outdir" ]; then
		rm -r OMCL_outdir
	fi
	perl $6/Synima-master/util/Blast_all_vs_all_repo_to_OrthoMCL.pl -r Repo_spec.txt -t CDS -o OMCL_outdir >> stdout.log 2>&1
}

task_Get_families(){
	mkdir -p $1/$2/outputs
		
	# Get families.
	>$1/$2/outputs/stdout.log
	$4/Get_families_from_OrthoMCL.py \
		--pred-lncRNAs $3/05-Families/$1/$2/ids_by_specie.tsv \
		--orthomcl $3/04-OrthoMCL/$1/$2/all_orthomcl.out \
		--fam $3/05-Families/$1/$2/fam.tsv \
		--gen $3/05-Families/$1/$2/gen.tsv >> $1/$2/outputs/stdout.log 2>&1
}

"$@"

