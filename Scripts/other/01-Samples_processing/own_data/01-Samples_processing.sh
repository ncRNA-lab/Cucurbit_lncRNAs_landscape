#!/bin/bash

#SBATCH --job-name=SampProc							# Job name.
#SBATCH --output=Sample_Processing.log					# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=33								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2500mb							# Job memory request.


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/01-Samples_processing/own_data/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/01-Samples_processing/own_data/Functions.sh"
trim_adapters="/storage/ncRNA/Softwares/BBMap_38.90/bbmap/resources"

####### NEW VARIABLES
WD1=$WD"/01-Samples_processing/own_data/01-Raw_data"
WD2=$WD"/01-Samples_processing/own_data/02-Trimmed_data"
WD3=$WD"/01-Samples_processing/own_data/03-Strand_detection"
WD4=$WD"/01-Samples_processing/own_data/04-Selected_data"
table_samples_info=$AI"/Summary_own_samples/Summary_samples.tsv"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD1
mkdir -p $WD1/fastqc
mkdir -p $WD1/multiqc
mkdir -p $WD1/outputs
mkdir -p $WD2
mkdir -p $WD2/fastqc
mkdir -p $WD2/multiqc
mkdir -p $WD2/outputs
mkdir -p $WD2/Filter_table
mkdir -p $WD3
mkdir -p $WD3/01-Ref_gen
mkdir -p $WD3/02-Index
mkdir -p $WD3/03-Pseudomapping
mkdir -p $WD3/Filter_table
mkdir -p $WD4



####### PIPELINE STEP 1

echo -e "\n\n\n####################"
echo -e "###### STEP 1 ######"
echo -e "####################\n"

cd $WD1

Samples_list_1=$(cat $table_samples_info | awk -F"\t" '{print $1}')

### PREPARE SAMPLES NAME
## All files should be named as *.fastq and be compressed (*.fastq.gz).
echo -e "\nPREPARE SAMPLES NAME..."

# Convert files named as ".fq" to ".fastq"
for file in *.fq; do
    if [ -f "$file" ]; then
        mv "$file" "${file%.fq}.fastq"
    fi
done

# Compress all fastq files.
for file in *.fastq; do
    if [ -f "$file" ]; then
        pigz -p $SLURM_CPUS_PER_TASK "$file"
    fi
done

# Convert files named as ".fq.gz" to ".fastq.gz"
for file in *.fq.gz; do
    if [ -f "$file" ]; then
        mv "$file" "${file%.fq.gz}.fastq.gz"
    fi
done

### RAW DATA FASTQC
## Execute Fastqc.
echo -e "\nRAW DATA FASTQC..."

for name in $Samples_list_1; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1/outputs/stdout_Fastqc_$name.log --quiet --exclusive $F task_Fastqc1 $name &
done
wait

### RAW DATA MULTIQC
## Execute Multiqc.
echo -e "\nRAW DATA MULTIQC..."

multiqc ./fastqc -o ./multiqc -f >> ./outputs/stdout_Multiqc.log 2>&1



####### PIPELINE STEP 2

echo -e "\n\n\n####################"
echo -e "###### STEP 2 ######"
echo -e "####################\n"

cd $WD2

### TRIMMING
## Trim the libraries.
echo -e "\nRAW DATA TRIMMING..."

for name in $Samples_list_1; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2/outputs/stdout_Trimming_$name.log --quiet --exclusive $F task_Trimming $name $WD1 $WD2 $trim_adapters $SLURM_CPUS_PER_TASK &
done
wait

### CLEAN DATA FASTQC
## Execute Fastqc.
echo -e "\nCLEAN DATA FASTQC..."

for name in $Samples_list_1; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2/outputs/stdout_Fastqc_$name.log --quiet --exclusive $F task_Fastqc2 $name &
done
wait

### MULTIQC
## Execute Multiqc.
echo -e "\nCLEAN DATA MULTIQC..."

multiqc ./fastqc -o ./multiqc -f >> ./outputs/stdout_Multiqc.log 2>&1

### SUMMARY TABLE
## Create a summary table about trimming results.
echo -e "\nCREATE A SUMMARY TABLE..."

# Remove temporal directory of Filter_table directory.
if [ -d "$WD2/Filter_table/Temp" ]; then
	cd $WD2/Filter_table
	rm -r Temp
	mkdir Temp
	cd $WD2
else
	mkdir $WD2/Filter_table/Temp
fi

# Remove temporal directory of accession_list directory.
if [ -d "$AI/Summary_own_samples/Temp_1" ]; then
	cd $AI/Summary_own_samples
	rm -r Temp_1
	mkdir Temp_1
	cd $WD2
else
	mkdir $AI/Summary_own_samples/Temp_1
fi

# Calculate library depth and determine if samples pass the filter.
for name in $Samples_list_1; do
	old_name=$name
	new_name=$(awk -F"\t" -v var="$old_name" '$1 == var {print $2}' $table_samples_info)
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD2/outputs/stdout_Summary_table_trimming_$old_name.log --quiet --exclusive $F task_Summary_table_trimming $old_name $new_name $WD2/Filter_table/Temp $AI/Summary_own_samples/Temp_1 &
done
wait

# Create a table with the library depth filtering results. 
cd Filter_table
echo -e "old_sample\tnew_sample\tnote\tdecision\tdepth" > $WD2/Filter_table/Trim_info.tsv
find Temp -type f -name '*.tsv' -print0 | xargs -0 cat >> $WD2/Filter_table/Trim_info.tsv
rm -r Temp

# Create a list with the sample identifiers which pass the library depth filter.
cd $AI/Summary_own_samples
find Temp_1 -type f -name '*.txt' -print0 | xargs -0 cat > $AI/Summary_own_samples/Samples_trim_filt.txt
rm -r Temp_1

echo -e "\nResults:\n"
cat $WD2/Filter_table/Trim_info.tsv



####### PIPELINE: STEP 3

echo -e "\n\n\n####################"
echo -e "###### STEP 3 ######"
echo -e "####################\n"

Samples_list_2=$(cat $AI/Summary_own_samples/Samples_trim_filt.txt)

### EXTRACT THE REFERENCE TRANSCRIPTOME
## Using the genome and the annotation file GFF3 (only genes), extract the reference transcriptome with the RSEM software. We'll get a fasta file with all gene sequences.
echo -e "\nEXTRACT THE REFERENCE TRANSCRIPTOME..."

cd $WD3/01-Ref_gen
mkdir -p outputs
>./outputs/stdout_Ref_gen.log
rsem-prepare-reference --gff3 $AI/GFF3_genes/$specie.gff3 $AI/Genome/$specie.fa $specie >> ./outputs/stdout_Ref_gen.log 2>&1

### INDEX THE REFERENCE TRANSCRIPTOME
## Index the reference transcriptome using Salmon software.
echo -e "\nINDEX THE REFERENCE TRANSCRIPTOME..."

cd $WD3/02-Index
mkdir -p outputs
>./outputs/stdout_Index.log
salmon index -t $WD3/01-Ref_gen/$specie.transcripts.fa -i $specie >> ./outputs/stdout_Index.log 2>&1

### PSEUDOMAPPING
## Pseudomap each library against the indexed reference transcriptome using Salmon software.
echo -e "\nPSEUDOMAPPING..."

cd $WD3/03-Pseudomapping
mkdir -p outputs
for name in $Samples_list_2; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD3/03-Pseudomapping/outputs/stdout_Pseudomapping_$name.log --quiet --exclusive $F task_Pseudomapping $name $WD2 $WD3 $specie $AI $SLURM_CPUS_PER_TASK &
done
wait

### GENERATE A SUMMARY STRAND INFO TABLE
## Generate a summary strand info table using an in-house script.
echo -e "\nGENERATE STRAND INFO..."

cd $WD3/Filter_table
Generate_strand_info.py \
	--path $WD/01-Samples_processing/own_data \
	--sample-list $AI/Summary_own_samples/Samples_trim_filt.txt

mv $WD3/Filter_table/Strand_info.tsv $WD3/Filter_table/Strand_info_temp.tsv
echo -e "old_sample\tnew_sample\tmapping_rate\tstrandedness\tnote" > $WD3/Filter_table/Strand_info.tsv
for name in $Samples_list_2; do
	old_name=$name
	new_name=$(awk -F"\t" -v var="$old_name" '$1 == var {print $2}' $table_samples_info)
	awk -F"\t" -v var1="$old_name" -v var2="$new_name" '$1 == var1 {print $1"\t"var2"\t"$2"\t"$3"\t"$4}' $WD3/Filter_table/Strand_info_temp.tsv >> $WD3/Filter_table/Strand_info.tsv
done
rm $WD3/Filter_table/Strand_info_temp.tsv

echo -e "\nResults:\n"
cat $WD3/Filter_table/Strand_info.tsv

### SELECT STRAND-SPECIFIC LIBRARIES
## Select and save strand-specific libraries.
echo -e "\nSELECT STRAND SPECIFIC LIBRARIES...\n"

## Create the new list of accessions according to the results in strand_info.tsv.
cd $WD3/Filter_table
tail -n +2 Strand_info.tsv | awk -F'\t' '$5 == "used for transcript reconstruction" {print $1}' > $AI/Summary_own_samples/Samples_SS_filt.txt

## Move the selected samples to the directory 04-Selected_data.
cd $WD4
Samples_list_3=$(cat $AI/Summary_own_samples/Samples_SS_filt.txt)
for name in $Samples_list_3; do
	old_name=$name
	new_name=$(awk -F"\t" -v var="$old_name" '$1 == var {print $2}' $table_samples_info)
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Select_SS_libraries $old_name $new_name $WD2 &
done
wait


