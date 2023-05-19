#!/bin/bash

#SBATCH --job-name=cmeS2.3							# Job name.
#SBATCH --output=cme_STEP2.3.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=50								# Run on one mode. 
#SBATCH --cpus-per-task=4							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=4gb							# Job memory request.


####### VARIABLES
specie="cme"
update="Update_2023_03_29"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Scripts/"$update"/01-Samples_processing/Functions.sh"

####### NEW VARIABLES
WD1=$WD"/"$update"/01-Samples_processing"
WD1_spe=$WD1"/"$specie
Acc_list=$AI"/sra-info/"$update"/accession_list/"$specie"-SRR_Acc_List.txt"
New_Acc_list=$AI"/sra-info/"$update"/accession_list/"$specie"-SRR_Acc_List-Filter_1.txt"
table=$WD1_spe"/02-Trimmed_data/Filter_table/Trim_info.tsv"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD"/"$update
mkdir -p $WD1
mkdir -p $WD1_spe
mkdir -p $WD1_spe/02-Trimmed_data
mkdir -p $WD1_spe/02-Trimmed_data/fastqc
mkdir -p $WD1_spe/02-Trimmed_data/multiqc
mkdir -p $WD1_spe/02-Trimmed_data/outputs
mkdir -p $WD1_spe/02-Trimmed_data/Filter_table


####### PIPELINE: STEP 2

echo -e "\n\n\n######################"
echo -e "###### STEP 2.3 ######"
echo -e "######################\n"

cd $WD1_spe/02-Trimmed_data

### FASTQC
## Execute Fastqc.
echo -e "\nFASTQC..."

for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/02-Trimmed_data/outputs/stdout_Fastqc_$SRR.log --quiet --exclusive $F task_Fastqc2 $SRR &
done
wait

### MULTIQC
## Execute Multiqc.
echo -e "\nMULTIQC..."

multiqc ./fastqc -o ./multiqc -f >> ./outputs/stdout_Multiqc.log 2>&1

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
if [ -d "$AI/sra-info/$update/accession_list/$specie.temp_1" ]; then
	cd $AI/sra-info/$update/accession_list
	rm -r $specie.temp_1
	mkdir $specie.temp_1
	cd $WD1_spe/02-Trimmed_data
else
	mkdir $AI/sra-info/$update/accession_list/$specie.temp_1
fi

# Calculate library depth and determine if samples pass the filter.
for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/02-Trimmed_data/outputs/stdout_Summary_table_trimming_$SRR.log --quiet --exclusive $F task_Summary_table_trimming $SRR $WD1_spe/02-Trimmed_data/Filter_table/Temp $AI/sra-info/$update/accession_list/$specie.temp_1 &
done
wait

# Create a table with the library depth filtering results. 
cd Filter_table
echo -e "sample\tnote\tdecision\tdepth" > $table
find Temp -type f -name '*[0-9].tsv' -print0 | xargs -0 cat >> $table
rm -r Temp

# Create a list with the sample identifiers which pass the library depth filter.
cd $AI/sra-info/$update/accession_list
find $specie.temp_1 -type f -name '*[0-9].txt' -print0 | xargs -0 cat > $New_Acc_list
rm -r $specie.temp_1

echo -e "\nResults:\n"
cat $table

