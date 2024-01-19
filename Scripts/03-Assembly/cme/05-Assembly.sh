#!/bin/bash

#SBATCH --job-name=cmeS5							# Job name.
#SBATCH --output=cme_STEP5.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=33								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb							# Job memory request.


####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/01-Sample_processing_and_selection"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/02-Mapping"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/03-Assembly"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/03-Assembly/Functions.sh"

####### NEW VARIABLES
WD1_spe=$WD1"/"$specie
WD2_spe=$WD2"/"$specie
WD3_spe=$WD3"/"$specie
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_2.txt"
strand_info=$WD1_spe"/03-Strand_detection/Filter_table/strand_info.tsv"

####### DIRECTORY
mkdir -p $WD3
mkdir -p $WD3_spe
mkdir -p $WD3_spe/01-Individual_assembly
mkdir -p $WD3_spe/02-Merged_assembly


####### PIPELINE: STEP 5

### ASSEMBLY
## Assembly each trimmed library using Stringtie software.
echo -e "\nASSEMBLY..."

cd $WD3_spe/01-Individual_assembly
mkdir -p Outputs
for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD3_spe/01-Individual_assembly/Outputs/stdout_Assembly_$SRR.log --quiet --exclusive $F task_Assembly $SRR $WD2_spe $AI $strand_info $specie $SLURM_CPUS_PER_TASK &
done
wait

### ASSEMBLY MERGING
## Merge each individual assembly using Stringtie software.
echo -e "\nASSEMBLY MERGING..."

cd $WD3_spe/02-Merged_assembly
mkdir -p Outputs
>./Outputs/stdout_Merging.log

## Create a TXT file with the assembly (GTF file) route for each sample.
awk '{print "../01-Individual_assembly/"$1".gtf"}' $Acc_list > assemblies.txt

## Merge the GTF files.
stringtie \
	--merge \
	-o $specie\_merged.gtf \
	-m 200 \
	-g 50 \
	-F 0.3 \
	-T 0 \
	-p $SLURM_CPUS_PER_TASK \
	-v \
	-G $AI/GTF_genes/$specie.gtf \
	assemblies.txt >> ./Outputs/stdout_Merging.log 2>&1


