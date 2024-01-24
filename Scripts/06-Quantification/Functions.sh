#!/bin/bash

####### FUNCTIONS

## Quantify LncRNAs and genes at the same time by salmon.
task_Salmon_All(){
	LIBTYPE=$(tail -n +2 $5 | awk -v a="$1" -F'\t' '$1 == a {print $3}')
	if [ -e $2/$1"_tr_1P.fastq.gz" ] && [ -e $2/$1"_tr_2P.fastq.gz" ]; then
		if [ $LIBTYPE = "ISR" ] || [ $LIBTYPE = "ISF" ]; then
			echo -e $1" is a strand-specific library (ISR or ISF)."
			salmon quant \
				-i $4/02-Index/$7 \
				-l $LIBTYPE \
				-1 $2/$1"_tr_1P.fastq.gz" \
				-2 $2/$1"_tr_2P.fastq.gz" \
				-g $3/ALL.gtf \
				-p $6 \
				-o $4/03-Quant/$1 \
				--gcBias
		else
			echo -e $1" is unstranded library (IU)."
			salmon quant \
				-i $4/02-Index/$7 \
				-l $LIBTYPE \
				-1 $2/$1"_tr_1P.fastq.gz" \
				-2 $2/$1"_tr_2P.fastq.gz" \
				-g $3/ALL.gtf \
				-p $6 \
				-o $4/03-Quant/$1 \
				--gcBias
		fi	
	else
		if [ $LIBTYPE = "SR" ] || [ $LIBTYPE = "SF" ]; then
			echo -e $1" is a strand-specific library (SR or SF)."
			## Alignment
			salmon quant \
				-i $4/02-Index/$7 \
				-l $LIBTYPE \
		    		-r $2/$1"_tr.fastq.gz" \
		    		-g $3/ALL.gtf \
				-p $6 \
				-o $4/03-Quant/$1 \
				--gcBias
		else
			echo -e $1" is unstranded library (U)."
			salmon quant \
				-i $4/02-Index/$7 \
				-l $LIBTYPE \
		    		-r $2/$1"_tr.fastq.gz" \
		    		-g $3/ALL.gtf \
				-p $6 \
				-o $4/03-Quant/$1 \
				--gcBias
		fi		
	fi
}

## Quantify LncRNAs by salmon.
task_Salmon_LncRNAs(){
	LIBTYPE=$(tail -n +2 $5 | awk -v a="$1" -F'\t' '$1 == a {print $3}')
	if [ -e $2/$1"_tr_1P.fastq.gz" ] && [ -e $2/$1"_tr_2P.fastq.gz" ]; then
		if [ $LIBTYPE = "ISR" ] || [ $LIBTYPE = "ISF" ]; then
			echo -e $1" is a strand-specific library (ISR or ISF)."
			salmon quant \
				-i $4/02-Index/$7 \
				-l $LIBTYPE \
				-1 $2/$1"_tr_1P.fastq.gz" \
				-2 $2/$1"_tr_2P.fastq.gz" \
				-g $3/POTENTIAL_LNCRNAS_pred.gtf \
				-p $6 \
				-o $4/03-Quant/$1 \
				--gcBias
		else
			echo -e $1" is unstranded library (IU)."
			salmon quant \
				-i $4/02-Index/$7 \
				-l $LIBTYPE \
				-1 $2/$1"_tr_1P.fastq.gz" \
				-2 $2/$1"_tr_2P.fastq.gz" \
				-g $3/POTENTIAL_LNCRNAS_pred.gtf \
				-p $6 \
				-o $4/03-Quant/$1 \
				--gcBias
		fi	
	else
		if [ $LIBTYPE = "SR" ] || [ $LIBTYPE = "SF" ]; then
			echo -e $1" is a strand-specific library (SR or SF)."
			## Alignment
			salmon quant \
				-i $4/02-Index/$7 \
				-l $LIBTYPE \
		    		-r $2/$1"_tr.fastq.gz" \
		    		-g $3/POTENTIAL_LNCRNAS_pred.gtf \
				-p $6 \
				-o $4/03-Quant/$1 \
				--gcBias
		else
			echo -e $1" is unstranded library (U)."
			salmon quant \
				-i $4/02-Index/$7 \
				-l $LIBTYPE \
		    		-r $2/$1"_tr.fastq.gz" \
		    		-g $3/POTENTIAL_LNCRNAS_pred.gtf \
				-p $6 \
				-o $4/03-Quant/$1 \
				--gcBias
		fi		
	fi
}

## Quantify Genes by salmon.
task_Salmon_Genes(){
	LIBTYPE=$(tail -n +2 $5 | awk -v a="$1" -F'\t' '$1 == a {print $3}')
	if [ -e $2/04-Selected_data/$1"_tr_1P.fastq.gz" ] && [ -e $2/04-Selected_data/$1"_tr_2P.fastq.gz" ]; then
		if [ $LIBTYPE = "ISR" ] || [ $LIBTYPE = "ISF" ]; then
			echo -e $1" is a strand-specific library (ISR or ISF)."
			salmon quant \
				-i $4/GENES/02-Index/$7 \
				-l $LIBTYPE \
				-1 $2/04-Selected_data/$1"_tr_1P.fastq.gz" \
				-2 $2/04-Selected_data/$1"_tr_2P.fastq.gz" \
				-g $3/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf \
				-p $6 \
				-o $4/GENES/03-Quant/$1 \
				--gcBias
		else
			echo -e $1" is unstranded library (IU)."
			salmon quant \
				-i $4/GENES/02-Index/$7 \
				-l $LIBTYPE \
				-1 $2/04-Selected_data/$1"_tr_1P.fastq.gz" \
				-2 $2/04-Selected_data/$1"_tr_2P.fastq.gz" \
				-g $3/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf \
				-p $6 \
				-o $4/GENES/03-Quant/$1 \
				--gcBias
		fi	
	else
		if [ $LIBTYPE = "SR" ] || [ $LIBTYPE = "SF" ]; then
			echo -e $1" is a strand-specific library (SR or SF)."
			## Alignment
			salmon quant \
				-i $4/GENES/02-Index/$7 \
				-l $LIBTYPE \
		    		-r $2/04-Selected_data/$1"_tr.fastq.gz" \
		    		-g $3/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf \
				-p $6 \
				-o $4/GENES/03-Quant/$1 \
				--gcBias
		else
			echo -e $1" is unstranded library (U)."
			salmon quant \
				-i $4/GENES/02-Index/$7 \
				-l $LIBTYPE \
		    		-r $2/04-Selected_data/$1"_tr.fastq.gz" \
		    		-g $3/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf \
				-p $6 \
				-o $4/GENES/03-Quant/$1 \
				--gcBias
		fi		
	fi
}

"$@"


