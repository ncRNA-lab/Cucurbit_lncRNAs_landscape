#!/bin/bash

#SBATCH --job-name=cmeS2							# Job name.
#SBATCH --output=cme_STEP2.log						# Standard output and error log.
#SBATCH --qos=medium								# Partition (queue)
#SBATCH --ntasks=16								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=5-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2500mb							# Job memory request.


####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/01-Sample_processing_and_selection"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/01-Sample_processing_and_selection/Functions.sh"
trim_adapters="/storage/ncRNA/Softwares/Trimmomatic-0.39/adapters"

####### NEW VARIABLES
WD1_spe=$WD1"/"$specie
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List.txt"
New_Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_1.txt"
table=$WD1_spe"/02-Trimmed_data/Filter_table/trim_info.tsv"

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1_spe
mkdir -p $WD1_spe/02-Trimmed_data
mkdir -p $WD1_spe/02-Trimmed_data/Fastqc
mkdir -p $WD1_spe/02-Trimmed_data/Multiqc
mkdir -p $WD1_spe/02-Trimmed_data/Outputs
mkdir -p $WD1_spe/02-Trimmed_data/Filter_table


####### PIPELINE: STEP 2

cd $WD1_spe/02-Trimmed_data

### TRIMMING
## Trim the libraries.
echo -e "\nTRIMMING..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/02-Trimmed_data/Outputs/stdout_Trimming_$SRR.log --quiet --exclusive $F task_Trimming $SRR $trim_adapters $SLURM_CPUS_PER_TASK &
done
wait

### FASTQC
## Execute Fastqc.
echo -e "\nFASTQC..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/02-Trimmed_data/Outputs/stdout_Fastqc_$SRR.log --quiet --exclusive $F task_Fastqc2 $SRR &
done
wait

### MULTIQC
## Execute Multiqc.
echo -e "\nMULTIQC..."

multiqc ./Fastqc -o ./Multiqc -f >> ./Outputs/stdout_Multiqc.log 2>&1

### SUMMARY TABLE
## Create a summary table about trimming results.
echo -e "\nCREATE A SUMMARY TABLE..."

# Remove temporal directory of Filter_table directory.
if [ -d "$WD1_spe/02-Trimmed_data/Filter_table/Temp" ]; then
	cd $WD1_spe/02-Trimmed_data/Filter_table
	rm -r Temp
	mkdir Temp
	cd $WD1_spe/02-Trimmed_data
else
	mkdir $WD1_spe/02-Trimmed_data/Filter_table/Temp
fi

# Remove temporal directory of accession_list directory.
if [ -d "$AI/sra-info/accession_list/$specie\_Temp" ]; then
	cd $AI/sra-info/accession_list
	rm -r $specie\_Temp
	mkdir $specie\_Temp
	cd $WD1_spe/02-Trimmed_data
else
	mkdir $AI/sra-info/accession_list/$specie\_Temp
fi

# Calculate library depth and determine if samples pass the filter.
for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/02-Trimmed_data/Outputs/stdout_Summary_table_trimming_$SRR.log --quiet --exclusive $F task_Summary_table_trimming $SRR $WD1_spe/02-Trimmed_data/Filter_table/Temp $AI/sra-info/accession_list/$specie\_Temp &
done
wait

# Create a table with the library depth filtering results. 
cd Filter_table
echo -e "sample\tdepth\tnote" > $table
find Temp -type f -name '*.tsv' -print0 | xargs -0 cat >> $table
rm -r Temp

# Create a list with the sample identifiers which pass the library depth filter.
cd $AI/sra-info/accession_list
find $specie\_Temp -type f -name '*.txt' -print0 | xargs -0 cat > $New_Acc_list
rm -r $specie\_Temp


