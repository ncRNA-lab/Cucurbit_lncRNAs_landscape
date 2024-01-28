#!/bin/bash

####### FUNCTIONS
# https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
# https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.html#gsc.tab=0

task_RepeatModeler(){
	## Run Repeatmodeler.
	# Make up a name for your database, choose your search engine, the number of threads, and the genome file.
	BuildDatabase -name $1 -engine rmblast $2/Genome/$1.fa
	RepeatModeler -database $1 -engine rmblast -LTRStruct -pa $3
}

task_RepeatMasker(){
	## Run RepeatMasker
	# I moved to a different directory, so I softlinked my classified file.  
	# Make sure you use the consensi.fa.classified file, or your repeats will just be masked by repeatmasker, but unannotated.
	# Make sure you softlink the classified file, otherwise you will not get a table of classified elements after the run.
	ln -s $3/01-Repeat_calling/01-RepeatModeler/$1-families.fa
	ln -s $2/Genome/$1.fa
	# This will produce a gff for the repeat mapping, a masked fasta file, and a table summarizing the repeats found in the genome.
	RepeatMasker -pa $4 -gff -lib $1-families.fa $1.fa
}

task_Intersect(){
	cd $2/02-Intersection/Intersection
	## Convert RepeatMasker to bed format file.
	rmsk2bed < $2/01-Repeat_calling/02-RepeatMasker/$1.fa.out > sorted-$1.fa.out.bed
	awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$5"\t"$6"\t"$4"\t"$11}' sorted-$1.fa.out.bed > sorted-mod-$1.fa.out.bed
	## Convert GTF to BED.
	cat $3/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > ORIGINAL_GENES.bed
	cat $3/STEP-FINAL/Files/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > POTENTIAL_LNCRNAS_pred_R.bed
	cat $3/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf | \
		convert2bed -i gtf --attribute-key=transcript_id | \
			awk '$8 == "transcript" {print $0}' > POTENTIAL_LNCRNAS_pred_NR.bed
	cp $4/Random_IR.bed ./
	## Intersect LncRNAs and Genes BED files with Repeat file. Parameters: strandness and without minimum overlap
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_R.bed -b sorted-mod-$1.fa.out.bed -s -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$12"\t"$13"\t"$14"\t"$17"\t"$18"\t"$19}' > POTENTIAL_LNCRNAS_intersect_Rep_R.tsv
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_NR.bed -b sorted-mod-$1.fa.out.bed -s -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$12"\t"$13"\t"$14"\t"$17"\t"$18"\t"$19}' > POTENTIAL_LNCRNAS_intersect_Rep_NR.tsv
	bedtools intersect -a ORIGINAL_GENES.bed -b sorted-mod-$1.fa.out.bed -s -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$12"\t"$13"\t"$14"\t"$17"\t"$18"\t"$19}' > ORIGINAL_GENES_intersect_Rep.tsv
	bedtools intersect -a Random_IR.bed -b sorted-mod-$1.fa.out.bed -s -wo -nonamecheck | \
		awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8"\t"$9"\t"$10"\t"$13"\t"$14"\t"$15}' > Random_IR_intersect_Rep.tsv
}

task_Final_tables(){
	echo -e "\tFlag: "$6", Confidence: "$7"..."
	Rscript $5/STEP1.R $1 $2 $3 $4 ${6^^} $7 >> $2"/Outputs/stdout_1-"$6"-"$7".log" 2>&1
	STEP2-Repeat.py --path $2"/Final_tables" --flag ${6^^} --confidence $7 >> $2"/Outputs/stdout_2_repeat-"$6"-"$7".log" 2>&1
	STEP2-Transposon.py --path $2"/Final_tables" --flag ${6^^} --confidence $7 >> $2"/Outputs/stdout_2_transposon-"$6"-"$7".log" 2>&1
}

"$@"


