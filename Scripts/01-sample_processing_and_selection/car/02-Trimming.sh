#!/bin/bash

#SBATCH --job-name=carS2							# Job name.
#SBATCH --output=car_STEP2.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=10								# Run on one mode. 
#SBATCH --cpus-per-task=10							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb							# Job memory request.


####### VARIABLES
specie="car"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/01-sample_processing_and_selection"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/01-sample_processing_and_selection/Functions.sh"
trim_adapters="/storage/ncRNA/Softwares/Trimmomatic-0.39/adapters"

####### NEW VARIABLES
WD1=$WD"/"$specie
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List.txt"
New_Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_1.txt"
table=$WD1"/02-Trimmed_data/Filter_table/Trim_info.tsv"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD1
mkdir -p $WD1/02-Trimmed_data
mkdir -p $WD1/02-Trimmed_data/fastqc
mkdir -p $WD1/02-Trimmed_data/multiqc
mkdir -p $WD1/02-Trimmed_data/outputs
mkdir -p $WD1/02-Trimmed_data/Filter_table


####### PIPELINE: STEP 2

cd $WD1/02-Trimmed_data

### TRIMMING
## Trim the libraries.
echo -e "\nTRIMMING..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/02-Trimmed_data/outputs/stdout_Trimming_$SRR.log --quiet --exclusive $F task_Trimming $SRR $trim_adapters $SLURM_CPUS_PER_TASK &
done
wait

### FASTQC
## Execute Fastqc.
echo -e "\nFASTQC..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/02-Trimmed_data/outputs/stdout_Fastqc_$SRR.log --quiet --exclusive $F task_Fastqc2 $SRR &
done
wait

### MULTIQC
## Execute Multiqc.
echo -e "\nMULTIQC..."

multiqc ./fastqc -o ./multiqc -f >> ./outputs/stdout_Multiqc.log 2>&1

### SUMMARY TABLE
## Create a summary table about trimming results.
echo -e "\nCREATE A SUMMARY TABLE..."

## Create Trim_info.tsv and SRR_Acc_List-Filter_1.txt
echo -e "sample\tnote\tdecision\tdepth" > $table
>$New_Acc_list

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/02-Trimmed_data/outputs/stdout_Summary_table_trimming_$SRR.log --quiet --exclusive $F task_Summary_table_trimming $SRR $table $New_Acc_list &
done
wait

