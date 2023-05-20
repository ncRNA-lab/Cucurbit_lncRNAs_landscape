#!/bin/bash

#SBATCH --job-name=Map								# Job name.
#SBATCH --output=Mapping.log							# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=33								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb							# Job memory request.


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/02-Mapping/own_data/Functions.sh"

####### NEW VARIABLES
WD1=$WD"/01-Samples_processing/own_data"
WD2=$WD"/02-Mapping/own_data"
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

tail -n +2 $strand_info | awk -F'\t' '$5 == "used for transcript reconstruction" {print $2}' > $AI/Summary_own_samples/Samples_SS_filt_conversion.txt
Samples_list=$(cat $AI/Summary_own_samples/Samples_SS_filt_conversion.txt)

for name in $Samples_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2/02-Mapping/outputs/stdout_Mapping_$name.log --quiet --exclusive $F task_Mapping $name $WD1 $WD2 $strand_info $specie $SLURM_CPUS_PER_TASK &
done
wait


