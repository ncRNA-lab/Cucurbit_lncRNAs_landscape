#!/bin/bash

####### FUNCTIONS

task_Assembly(){
	
	LIBTYPE=$(tail -n +2 $4 | awk -v a="$1" -F'\t' '$2 == a {print $4}')
	if [ $LIBTYPE = "ISR" ] || [ $LIBTYPE = "SR" ]; then
		echo -e $1" is a strand-specific library (ISR/SR).\n"
		## Assembly
		stringtie \
			-l "STRGOWN" \
			-m 200 \
			-p $6 \
			-G $3/GTF_genes/$5.gtf \
			-o $1.gtf \
			--rf \
			-A $1"_expression.tsv" \
			$2/02-Mapping/$1/accepted_hits.bam
	elif [ $LIBTYPE = "ISF" ] || [ $LIBTYPE = "SF" ]; then
		echo -e $1" is a strand-specific library (ISF/SF).\n"
		## Assembly
		stringtie \
			-l "STRGOWN" \
			-m 200 \
			-p $6 \
			-G $3/GTF_genes/$5.gtf \
			-o $1.gtf \
			--fr \
			-A $1"_expression.tsv" \
			$2/02-Mapping/$1/accepted_hits.bam
	else
		echo -e $1" is unstranded library (IU/U).\n"
	fi
}

"$@"

