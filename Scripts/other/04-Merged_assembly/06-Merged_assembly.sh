#!/bin/bash

#SBATCH --job-name=cmeS6							# Job name.
#SBATCH --output=cme_STEP6.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=1								# Run on one mode. 
#SBATCH --cpus-per-task=20							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb							# Job memory request.


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"

####### NEW VARIABLES
WD1_own=$WD"/03-Individual_assembly/own_data"
WD1_sra=$WD"/03-Individual_assembly/sra_data"
WD2=$WD"/04-Merged_assembly"
Samples_own=$AI"/Summary_own_samples/Samples_SS_filt_conversion.txt"
Samples_sra=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_2.txt"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD2


####### PIPELINE: STEP 6

echo -e "\n\n\n####################"
echo -e "###### STEP 6 ######"
echo -e "####################\n"

### ASSEMBLY MERGING
## Merge each individual assembly using Stringtie software.
echo -e "\nASSEMBLY MERGING..."

cd $WD2
mkdir -p outputs

## Create a TXT file with the assembly (GTF file) route for each sample.
awk -v a="$WD1_own" '{print a"/"$1".gtf"}' $Samples_own > assemblies.txt
awk -v a="$WD1_sra" '{print a"/"$1".gtf"}' $Samples_sra >> assemblies.txt

## Merge the GTF files (Unify the transcriptome -g 250).
echo -e "\n-g 250..."
>./outputs/stdout_Merging_250.log
stringtie \
	--merge \
	-o $specie\_merged_250.gtf \
	-l "MSTRGOWN" \
	-m 200 \
	-g 250 \
	-F 0.3 \
	-T 0 \
	-p $SLURM_CPUS_PER_TASK \
	-v \
	-G $AI/GTF_genes/$specie.gtf \
	assemblies.txt >> ./outputs/stdout_Merging_250.log 2>&1

## Merge the GTF files (Unify the transcriptome -g 50).
echo -e "\n-g 50..."
>./outputs/stdout_Merging_50.log
stringtie \
	--merge \
	-o $specie\_merged_50.gtf \
	-l "MSTRGOWN" \
	-m 200 \
	-g 50 \
	-F 0.3 \
	-T 0 \
	-p $SLURM_CPUS_PER_TASK \
	-v \
	-G $AI/GTF_genes/$specie.gtf \
	assemblies.txt >> ./outputs/stdout_Merging_50.log 2>&1


