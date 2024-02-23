#!/bin/bash

####### FUNCTIONS

task_STEP3_CLOSEST-ANTISENSE(){
	Rscript $8/Antisense-STEP3-Individual_pearson_correlation_table-CLOSEST.R $1 $2 $3 $4 $5 $6 $7
}

task_STEP3_CLOSEST-INTERGENIC(){
	Rscript $8/Intergenic-STEP3-Individual_pearson_correlation_table-CLOSEST.R $1 $2 $3 $4 $5 $6 $7
}

task_STEP3_RANGE-INTERGENIC(){
	Rscript $9/Intergenic-STEP3-Individual_pearson_correlation_table-RANGE.R $1 $2 $3 $4 $5 $6 $7 $8
}

task_STEP5_CLOSEST-ANTISENSE(){
	Rscript $9/Antisense-STEP5-Individual_pearson_correlation_table-CLOSEST-RANDOM.R "$1" $2 $3 $4 $5 $6 $7 $8
}

task_STEP5_CLOSEST-INTERGENIC(){
	Rscript $9/Intergenic-STEP5-Individual_pearson_correlation_table-CLOSEST-RANDOM.R "$1" $2 $3 $4 $5 $6 $7 $8
}

"$@"


