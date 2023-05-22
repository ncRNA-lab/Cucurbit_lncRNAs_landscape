#!/bin/bash

####### FUNCTIONS
# https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
# https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.html#gsc.tab=0

task_RepeatModeler(){
	echo -e "\t"$1"..."
	mkdir -p $3/01-Repeat_calling/01-RepeatModeler/$1
	cd $3/01-Repeat_calling/01-RepeatModeler/$1
	## Run Repeatmodeler.
	# Make up a name for your database, choose your search engine, the number of threads, and the genome file.
	BuildDatabase -name $1 -engine rmblast $2/Genome/$1.fa 1> ../Logs/Database_stdout_$1.log
	RepeatModeler -database $1 -engine rmblast -LTRStruct -pa $4 1> ../Logs/RepeatModeler_stdout_$1.log
}

task_transposonPSI(){
	cd temp/folder_$1
	# Run TransposonPSI.
	transposonPSI.pl db_$1.fasta $2 >> TransposonPSI_stdout_$1.log 2>&1
}

task_RemoveProt(){
	echo -e "\t"$1"..."
	if [ -d "$2/01-Repeat_calling/04-RemoveProt_unknown_repeats/$1" ]; then
		rm -r $1
	fi
	mkdir $1
	cd $1
	
	WD02="$2/01-Repeat_calling/02-TransposonPSI_cucurbitaceae"
	WD03="$2/01-Repeat_calling/03-TransposonPSI_unknown_repeats"
	
	cp $WD02/cucurbitaceae.filtered.fa ./
	sed 's/ /|/g' $WD03/$1/RM_unknowns.filtered.fa > RM_unknowns_mod.filtered.fa
	
	# Blastx
	diamond makedb --in cucurbitaceae.filtered.fa -d cucurbitaceae.dmnd >> Makedb_stdout_$1.log 2>&1
	diamond blastx -d cucurbitaceae.dmnd -q RM_unknowns_mod.filtered.fa -o Blastx_res.tsv --threads $3 --top 5 --strand plus --evalue 1e-10 -f 6 qseqid qlen sseqid slen qstart qend sstart send qframe pident length evalue bitscore >> Blastx_stdout_$1.log 2>&1
	
	# Get fasta.
	awk '{print $1}' Blastx_res.tsv | sort -u > ids.txt
	source activate gaas
	gaas_fasta_removeSeqFromIDlist.pl -f RM_unknowns_mod.filtered.fa -l ids.txt -o RM_unknowns_temp.filtered.noProtFinal.fa
	conda deactivate
	rm ids.txt
	grep ">" RM_unknowns_temp.filtered.noProtFinal.fa | sed 's/>//g' > ids.txt
	seqtk subseq RM_unknowns_temp.filtered.noProtFinal.fa ids.txt | sed 's/|/ /g' > RM_unknowns.filtered.noProtFinal.fa
	rm ids.txt
	
	# Create Final Customed repeat library (CRL).
	cat $WD03/$1/RM_identities.fa $WD03/$1/RM_unknowns_matched_on_transposonPSI.fa $WD03/$1/RM_LTRStruct.fa RM_unknowns.filtered.noProtFinal.fa > $1-families.fa
}

task_RepeatMasker(){
	echo -e "\t"$1"..."
	mkdir -p $3/01-Repeat_calling/05-RepeatMasker/$1
	cd $3/01-Repeat_calling/05-RepeatMasker/$1
	## Run RepeatMasker
	# I moved to a different directory, so I softlinked my classified file.  
	# Make sure you use the consensi.fa.classified file, or your repeats will just be masked by repeatmasker, but unannotated.
	# Make sure you softlink the classified file, otherwise you will not get a table of classified elements after the run.
	ln -s $3/01-Repeat_calling/01-RepeatModeler/$1/$1-families.fa
	ln -s $2/Genome/$1.fa
	# This will produce a gff for the repeat mapping, a masked fasta file, and a table summarizing the repeats found in the genome.
	RepeatMasker -pa $4 -gff -lib $1-families.fa $1.fa 1> ../Logs/RepeatMasker_stdout_$1.log
}

task_RepeatMasker_RMfiltered(){
	echo -e "\t"$1"..."
	mkdir -p $3/01-Repeat_calling/06-RepeatMasker_RMfiltered/$1
	cd $3/01-Repeat_calling/06-RepeatMasker_RMfiltered/$1
	## Run RepeatMasker
	# I moved to a different directory, so I softlinked my classified file.  
	# Make sure you use the consensi.fa.classified file, or your repeats will just be masked by repeatmasker, but unannotated.
	# Make sure you softlink the classified file, otherwise you will not get a table of classified elements after the run.
	ln -s $3/01-Repeat_calling/04-RemoveProt_unknown_repeats/$1/$1-families.fa
	ln -s $2/Genome/$1.fa
	# This will produce a gff for the repeat mapping, a masked fasta file, and a table summarizing the repeats found in the genome.
	RepeatMasker -pa $4 -gff -lib $1-families.fa $1.fa 1> ../Logs/RepeatMasker_stdout_$1.log
}

task_intersect(){
	echo -e "\t"$1"..."
	cd $2/02-Comparison_Genes_LncRNAs/
	if [ -d "$2/02-Comparison_Genes_LncRNAs/$1" ]; then
		rm -r $1
	fi
	mkdir $1
	cd $1
	## Convert RepeatMasker to bed format file.
	rmsk2bed < $2/01-Repeat_calling/05-RepeatMasker/$1/$1.fa.out > sorted-$1.fa.out.bed
	awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$5"\t"$6"\t"$4"\t"$11}' sorted-$1.fa.out.bed > sorted-mod-$1.fa.out.bed
	## Convert GTF to BED.
	cat $3/$1/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > ORIGINAL_GENES.bed
	cat $3/$1/STEP-FINAL/Files/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > POTENTIAL_LNCRNAS_pred_R.bed
	cat $3/$1/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > POTENTIAL_LNCRNAS_pred_NR.bed
	cp $4/$1/Random_IR.bed ./
	## Intersect LncRNAs and Genes BED files with Repeat file. Parameters: strandness and 0.5 minimum overlap
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_R.bed -b sorted-mod-$1.fa.out.bed -s -F 0.5 -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$14"\t"$17"\t"$18}' > POTENTIAL_LNCRNAS_intersect_Rep_R.tsv
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_NR.bed -b sorted-mod-$1.fa.out.bed -s -F 0.5 -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$14"\t"$17"\t"$18}' > POTENTIAL_LNCRNAS_intersect_Rep_NR.tsv
	bedtools intersect -a ORIGINAL_GENES.bed -b sorted-mod-$1.fa.out.bed -s -F 0.5 -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$14"\t"$17"\t"$18}' > ORIGINAL_GENES_intersect_Rep.tsv
	bedtools intersect -a Random_IR.bed -b sorted-mod-$1.fa.out.bed -s -F 0.5 -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$13"\t"$14}' > Random_IR_intersect_Rep.tsv
	## Intersect LncRNAs and Genes BED files with Repeat file. Parameters: strandness and without minimum overlap
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_R.bed -b sorted-mod-$1.fa.out.bed -s -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$12"\t"$13"\t"$14"\t"$17"\t"$18"\t"$19}' > POTENTIAL_LNCRNAS_intersect_Rep_R_19.tsv
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_NR.bed -b sorted-mod-$1.fa.out.bed -s -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$12"\t"$13"\t"$14"\t"$17"\t"$18"\t"$19}' > POTENTIAL_LNCRNAS_intersect_Rep_NR_19.tsv
	bedtools intersect -a ORIGINAL_GENES.bed -b sorted-mod-$1.fa.out.bed -s -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$12"\t"$13"\t"$14"\t"$17"\t"$18"\t"$19}' > ORIGINAL_GENES_intersect_Rep_19.tsv
	bedtools intersect -a Random_IR.bed -b sorted-mod-$1.fa.out.bed -s -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8"\t"$9"\t"$10"\t"$13"\t"$14"\t"$15}' > Random_IR_intersect_Rep_15.tsv
}

task_intersect_RMfiltered(){
	echo -e "\t"$1"..."
	cd $2/03-Comparison_Genes_LncRNAs_RMfiltered/
	if [ -d "$2/03-Comparison_Genes_LncRNAs_RMfiltered/$1" ]; then
		rm -r $1
	fi
	mkdir $1
	cd $1
	## Convert RepeatMasker to bed format file.
	rmsk2bed < $2/01-Repeat_calling/06-RepeatMasker_RMfiltered/$1/$1.fa.out > sorted-$1.fa.out.bed
	awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$5"\t"$6"\t"$4"\t"$11}' sorted-$1.fa.out.bed > sorted-mod-$1.fa.out.bed
	## Convert GTF to BED.
	cat $3/$1/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > ORIGINAL_GENES.bed
	cat $3/$1/STEP-FINAL/Files/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > POTENTIAL_LNCRNAS_pred_R.bed
	cat $3/$1/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > POTENTIAL_LNCRNAS_pred_NR.bed
	cp $4/$1/Random_IR.bed ./
	## Intersect LncRNAs and Genes BED files with Repeat file. Parameters: strandness and 0.2 minimum overlap
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_R.bed -b sorted-mod-$1.fa.out.bed -s -F 0.2 -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$14"\t"$17"\t"$18}' > POTENTIAL_LNCRNAS_intersect_Rep_R.tsv
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_NR.bed -b sorted-mod-$1.fa.out.bed -s -F 0.2 -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$14"\t"$17"\t"$18}' > POTENTIAL_LNCRNAS_intersect_Rep_NR.tsv
	bedtools intersect -a ORIGINAL_GENES.bed -b sorted-mod-$1.fa.out.bed -s -F 0.2 -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$14"\t"$17"\t"$18}' > ORIGINAL_GENES_intersect_Rep.tsv
	bedtools intersect -a Random_IR.bed -b sorted-mod-$1.fa.out.bed -s -F 0.2 -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$13"\t"$14}' > Random_IR_intersect_Rep.tsv
}

task_Figures_and_tables(){
	echo -e "\tFlag: "$5", Confidence: "$6"..."
	Rscript $7/Create_Tables_and_Figures_A.R $1 $3 $4 $5 $6 >> $2"/"$3"/Figures_and_Tables/A/stdout-"$5"-"$6".log" 2>&1
	Rscript $7/Create_Tables_and_Figures_B.R $1 $3 $4 $5 $6 >> $2"/"$3"/Figures_and_Tables/B/stdout-"$5"-"$6".log" 2>&1
	Rscript $7/Create_Tables_and_Figures_C-STEP1.R $1 $3 $4 $5 $6 >> $2"/"$3"/Figures_and_Tables/C/stdout_1-"$5"-"$6".log" 2>&1
	Create_Tables_and_Figures_C-STEP2-Repeat.py --path $2"/"$3"/Figures_and_Tables/C" --flag $5 --confidence $6 >> $2"/"$3"/Figures_and_Tables/C/stdout_2_repeat-"$5"-"$6".log" 2>&1
	Create_Tables_and_Figures_C-STEP2-Transposon.py --path $2"/"$3"/Figures_and_Tables/C" --flag $5 --confidence $6 >> $2"/"$3"/Figures_and_Tables/C/stdout_2_transposon-"$5"-"$6".log" 2>&1
	Rscript $7/Create_Tables_and_Figures_C-STEP3.R $1 $3 $4 $5 $6 >> $2"/"$3"/Figures_and_Tables/C/stdout_3-"$5"-"$6".log" 2>&1
}

"$@"

