#!/bin/bash

####### FUNCTIONS

task_STEP3_CLOSEST(){
	Rscript $8/Intergenic-STEP3-Individual_Pearson_Corr_Table-CLOSEST.R $1 $2 $3 $4 $5 $6 $7
}

task_STEP3_RANGE(){
	Rscript $9/Intergenic-STEP3-Individual_Pearson_Corr_Table-RANGE.R $1 $2 $3 $4 $5 $6 $7 $8
}

task_STEP5_CLOSEST(){
	Rscript $9/Intergenic-STEP5-Individual_Pearson_Corr_Table-CLOSEST-RANDOM.R "$1" $2 $3 $4 $5 $6 $7 $8
}

task_STEP3-ANTISENSE(){
	Rscript $8/Antisense-STEP3-Individual_Pearson_Corr_Table.R $1 $2 $3 $4 $5 $6 $7
}

task_STEP5-ANTISENSE(){
	Rscript $9/Intergenic-STEP5-Individual_Pearson_Corr_Table-CLOSEST-RANDOM.R "$1" $2 $3 $4 $5 $6 $7 $8
}

"$@"

