#!/bin/bash

#SBATCH --job-name=cmeS3							# Job name.
#SBATCH --output=cme_STEP3.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=33								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb							# Job memory request.


####### MODULES
module load anaconda

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/01-Sample_processing_and_selection"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/01-Sample_processing_and_selection/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/01-Sample_processing_and_selection/Functions.sh"

####### NEW VARIABLES
WD1_spe=$WD1"/"$specie
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_1.txt"
New_Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_2.txt"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1_spe
mkdir -p $WD1_spe/03-Strand_detection
mkdir -p $WD1_spe/03-Strand_detection/01-Ref_gen
mkdir -p $WD1_spe/03-Strand_detection/02-Index
mkdir -p $WD1_spe/03-Strand_detection/03-Pseudomapping
mkdir -p $WD1_spe/03-Strand_detection/Filter_table
mkdir -p $WD1_spe/04-Selected_data


####### PIPELINE: STEP 3

### EXTRACT THE REFERENCE TRANSCRIPTOME
## Using the genome and the annotation file GFF3 (only genes), extract the reference transcriptome with the RSEM software. We'll get a fasta file with all gene sequences.
echo -e "\nEXTRACT THE REFERENCE TRANSCRIPTOME..."

cd $WD1_spe/03-Strand_detection/01-Ref_gen
mkdir -p Outputs
>./Outputs/stdout_Ref_gen.log
rsem-prepare-reference --gff3 $AI/GFF3_genes/$specie.gff3 $AI/Genome/$specie.fa $specie >> ./Outputs/stdout_Ref_gen.log 2>&1

### INDEX THE REFERENCE TRANSCRIPTOME
## Index the reference transcriptome using Salmon software.
echo -e "\nINDEX THE REFERENCE TRANSCRIPTOME..."

cd $WD1_spe/03-Strand_detection/02-Index
mkdir -p Outputs
>./Outputs/stdout_Index.log
salmon index -t $WD1_spe/03-Strand_detection/01-Ref_gen/$specie.transcripts.fa -i $specie >> ./Outputs/stdout_Index.log 2>&1

### PSEUDOMAPPING
## Pseudomap each library against the indexed reference transcriptome using Salmon software.
echo -e "\nPSEUDOMAPPING..."

cd $WD1_spe/03-Strand_detection/03-Pseudomapping
mkdir -p Outputs
for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/03-Strand_detection/03-Pseudomapping/Outputs/stdout_Pseudomapping_$SRR.log --quiet --exclusive $F task_Pseudomapping $SRR $WD1_spe $specie $AI $SLURM_CPUS_PER_TASK &
done
wait

### GENERATE A SUMMARY STRAND INFO TABLE
## Generate a summary strand info table using an in-house script.
echo -e "\nGENERATE STRAND INFO..."

cd $WD1_spe/03-Strand_detection/Filter_table
Generate_strand_info.py \
	--path $WD1_spe \
	--acc-list $Acc_list

### SELECT STRAND-SPECIFIC LIBRARIES
## Select and save strand-specific libraries.
echo -e "\nSELECT STRAND SPECIFIC LIBRARIES..."

## Create the new list of accessions according to the results in strand_info.tsv.
tail -n +2 $WD1_spe/03-Strand_detection/Filter_table/strand_info.tsv | awk -F'\t' '$4 == "used for transcript reconstruction" {print $1}' > $New_Acc_list

## Move the selected samples to the directory 04-Files_selected.
cd $WD1_spe/04-Selected_data
for SRR in $(cat $New_Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Select_SS_libraries $SRR $WD1_spe &
done
wait


