#!/bin/bash

####### FUNCTIONS

task_Mapping(){
	
	LIBTYPE=$(tail -n +2 $4 | awk -v a="$1" -F'\t' '$1 == a {print $3}')
	if [ -e $2/04-Selected_data/$1"_tr_1P.fastq.gz" ] && [ -e $2/04-Selected_data/$1"_tr_2P.fastq.gz" ]; then
		if [ $LIBTYPE = "ISR" ]; then
			mkdir -p $1
			echo -e $1" is a strand-specific library (ISR, fr-firststrand, RF).\n"
			## Mapping
			hisat2 \
				-x $3/01-Index/$5 \
				-1 $2/04-Selected_data/$1"_tr_1P.fastq.gz" \
				-2 $2/04-Selected_data/$1"_tr_2P.fastq.gz" \
				-S $1/$1.sam \
				--max-intronlen 10000 \
				--dta \
				-p $6 \
				--rna-strandness RF \
				--known-splicesite-infile $3/01-Index/known_splice_sites.txt \
				--novel-splicesite-outfile $1/novel_splice_sites.txt \
				--summary-file $1/align_summary.txt 
		elif [ $LIBTYPE = "ISF" ]; then
			mkdir -p $1
			echo -e $1" is a strand-specific library (ISF, fr-secondstrand, FR).\n"
			## Mapping
			hisat2 \
				-x $3/01-Index/$5 \
				-1 $2/04-Selected_data/$1"_tr_1P.fastq.gz" \
				-2 $2/04-Selected_data/$1"_tr_2P.fastq.gz" \
				-S $1/$1.sam \
				--max-intronlen 10000 \
				--dta \
				-p $6 \
				--rna-strandness FR \
				--known-splicesite-infile $3/01-Index/known_splice_sites.txt \
				--novel-splicesite-outfile $1/novel_splice_sites.txt \
				--summary-file $1/align_summary.txt
		fi
		
		## Convert sam to bam and sort it.
		sambamba view -t $6 -f bam -F "not unmapped" -S -o $1/$1.bam $1/$1.sam
		sambamba sort -t $6 --tmpdir=$1 -o $1/accepted_hits.bam $1/$1.bam
		rm $1/$1.sam
		rm $1/$1.bam
		
	elif [ -e $2/04-Selected_data/$1"_tr.fastq.gz" ]; then
		if [ $LIBTYPE = "SR" ]; then
			mkdir -p $1
			echo -e $1" is a strand-specific library (SR, fr-firststrand, R).\n"
			## Mapping
			hisat2 \
				-x $3/01-Index/$5 \
				-U $2/04-Selected_data/$1"_tr.fastq.gz" \
				-S $1/$1.sam \
				--max-intronlen 10000 \
				--dta \
				-p $6 \
				--rna-strandness R \
				--known-splicesite-infile $3/01-Index/known_splice_sites.txt \
				--novel-splicesite-outfile $1/novel_splice_sites.txt \
				--summary-file $1/align_summary.txt 
		elif [ $LIBTYPE = "SF" ]; then
			mkdir -p $1
			echo -e $1" is a strand-specific library (SF, fr-secondstrand, F).\n"
			## Mapping
			hisat2 \
				-x $3/01-Index/$5 \
				-U $2/04-Selected_data/$1"_tr.fastq.gz" \
				-S $1/$1.sam \
				--max-intronlen 10000 \
				--dta \
				-p $6 \
				--rna-strandness F \
				--known-splicesite-infile $3/01-Index/known_splice_sites.txt \
				--novel-splicesite-outfile $1/novel_splice_sites.txt \
				--summary-file $1/align_summary.txt
		fi
		
		## Convert sam to bam and sort it.
		sambamba view -t $6 -f bam -F "not unmapped" -S -o $1/$1.bam $1/$1.sam
		sambamba sort -t $6 --tmpdir=$1 -o $1/accepted_hits.bam $1/$1.bam
		rm $1/$1.sam
		rm $1/$1.bam
	else
		echo -e $1" is an unstranded library (IU/U).\n"
	fi
}

"$@"

