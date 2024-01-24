#!/bin/bash

#SBATCH --job-name=cmeS8G						# Job name.
#SBATCH --output=cme_STEP8_GENES.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=50							# Run on one mode.
#SBATCH --cpus-per-task=4						# Number of tasks = cpus.
#SBATCH --time=1-00:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb						# Job memory request.


####### MODULES
module load R/4.1.2
module load anaconda

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/01-Sample_processing_and_selection"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/06-Quantification/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/06-Quantification/Functions.sh"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie
WD3_spe=$WD3/$specie
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_2.txt"
strand_info=$WD1_spe"/03-Strand_detection/Filter_table/strand_info.tsv"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD3
mkdir -p $WD3_spe
mkdir -p $WD3_spe/GENES
mkdir -p $WD3_spe/GENES/01-Ref
mkdir -p $WD3_spe/GENES/02-Index
mkdir -p $WD3_spe/GENES/03-Quant
mkdir -p $WD3_spe/GENES/04-Table


####### PIPELINE
cd $WD3_spe/GENES

## Transcriptome reference.
echo -e "\nSTEP 1: GET TRANSCRIPTOME REFERENCE"
cp $WD2_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES.fasta ./01-Ref/

## Index.
echo -e "\nSTEP 2: INDEX THE TRANSCRIPTOME REFERENCE"
mkdir -p ./02-Index/Outputs
>./02-Index/Outputs/stdout_Index.log
salmon index -i ./02-Index/$specie -t ./01-Ref/ORIGINAL_GENES.fasta >> ./02-Index/Outputs/stdout_Index.log 2>&1

## Quant.
echo -e "\nSTEP 3: QUANTIFY EACH LIBRARY BY SALMON"
mkdir -p ./03-Quant/Outputs
SRR_list=$(sed -e 's/\n/ /g' $Acc_list)
for SRR in $SRR_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD3_spe/GENES/03-Quant/Outputs/stdout_Quant_$SRR.log --quiet --exclusive $F task_Salmon_Genes $SRR $WD1_spe $WD2_spe $WD3_spe $strand_info $SLURM_CPUS_PER_TASK $specie &
done
wait

## Create TPM table.
echo -e "\nSTEP 4: CREATE A GLOBAL TABLE"
mkdir -p ./04-Table/Outputs
Rscript $AS/Create_TPM_table.R $WD3_spe/GENES $WD3_spe/GENES/03-Quant $WD2_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf $Acc_list >> ./04-Table/Outputs/stdout_TPMs_table.log 2>&1


