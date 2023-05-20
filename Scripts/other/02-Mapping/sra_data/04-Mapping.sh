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
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/02-Mapping/sra_data/Functions.sh"

####### NEW VARIABLES
WD1=$WD"/01-Samples_processing/sra_data"
WD2=$WD"/02-Mapping/sra_data"
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_2.txt"
strand_info=$WD1"/03-Strand_detection/Filter_table/Strand_info.tsv"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD/02-Mapping
mkdir -p $WD2
mkdir -p $WD2/01-Index
mkdir -p $WD2/02-Mapping
mkdir -p $AI/GTF_genes


####### PIPELINE: STEP 4

echo -e "\n\n\n####################"
echo -e "###### STEP 4 ######"
echo -e "####################\n"

### INDEX THE REFERENCE GENOME
## Index the reference genome using Hisat2 software.
echo -e "\nINDEX THE REFERENCE GENOME..."

cd $WD2/01-Index
mkdir -p outputs
>./outputs/stdout_Index.log

## Convert gff3 to gtf.
gffread $AI/GFF3_genes/$specie.gff3 -T -o $AI/GTF_genes/$specie.gtf -v >> ./outputs/stdout_Index.log 2>&1
## Index the genome.
hisat2-build $AI/Genome/$specie.fa $specie >> ./outputs/stdout_Index.log 2>&1
## Extract the known splite sites.
hisat2_extract_splice_sites.py $AI/GTF_genes/$specie.gtf > known_splice_sites.txt

### MAPPING
## Map each library against the reference genome using Hisat2 software.
echo -e "\nMAPPING..."

cd $WD2/02-Mapping
mkdir -p outputs
for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2/02-Mapping/outputs/stdout_Mapping_$SRR.log --quiet --exclusive $F task_Mapping $SRR $WD1 $WD2 $strand_info $specie $SLURM_CPUS_PER_TASK &
done
wait


