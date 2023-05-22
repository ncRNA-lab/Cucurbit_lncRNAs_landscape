#!/bin/bash

#SBATCH --job-name=cmeS4							# Job name.
#SBATCH --output=cme_STEP4.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=33								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb							# Job memory request.


####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/01-Sample_processing_and_selection"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/02-Mapping"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/02-Mapping/Functions.sh"

####### NEW VARIABLES
WD1_spe=$WD1"/"$specie
WD2_spe=$WD2"/"$specie
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_2.txt"
strand_info=$WD1_spe"/03-Strand_detection/Filter_table/strand_info.tsv"

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2_spe
mkdir -p $WD2_spe/01-Index
mkdir -p $WD2_spe/02-Mapping
mkdir -p $AI/GTF_genes


####### PIPELINE: STEP 4

### INDEX THE REFERENCE GENOME
## Index the reference genome using Hisat2 software.
echo -e "\nINDEX THE REFERENCE GENOME..."

cd $WD2_spe/01-Index
mkdir -p Outputs
>./Outputs/stdout_Index.log

## Convert gff3 to gtf.
gffread $AI/GFF3_genes/$specie.gff3 -T -o $AI/GTF_genes/$specie.gtf -v 
## Index the genome.
hisat2-build $AI/Genome/$specie.fa $specie >> ./Outputs/stdout_Index.log 2>&1
## Extract the known splite sites.
hisat2_extract_splice_sites.py $AI/GTF_genes/$specie.gtf > known_splice_sites.txt

### MAPPING
## Map each library against the reference genome using Hisat2 software.
echo -e "\nMAPPING..."

cd $WD2_spe/02-Mapping
mkdir -p Outputs
for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2_spe/02-Mapping/Outputs/stdout_Mapping_$SRR.log --quiet --exclusive $F task_Mapping $SRR $WD1_spe $WD2_spe $strand_info $specie $SLURM_CPUS_PER_TASK &
done
wait


