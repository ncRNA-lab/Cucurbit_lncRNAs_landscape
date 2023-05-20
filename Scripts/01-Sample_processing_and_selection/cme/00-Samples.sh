#!/bin/bash

#SBATCH --job-name=cmeS0							# Job name.
#SBATCH --output=cme_STEP0.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=1								# Run on one mode. 
#SBATCH --cpus-per-task=2							# Number of tasks = cpus.
#SBATCH --time=0-02:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=4gb							# Job memory request.


####### VARIABLES
specie="cme"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"

####### NEW VARIABLES
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List.txt"
Metadata=$AI"/sra-info/metadata/"$specie"-SraRunTable.csv"

####### DIRECTORY
mkdir -p $AI"/sra-info/"
mkdir -p $AI"/sra-info/accession_list"
mkdir -p $AI"/sra-info/metadata"


####### PIPELINE: STEP 0

## Download metadata info.
echo -e "\nDownload metadata info from SRA database..."
cd $AI"/sra-info/metadata"
esearch -db sra -query 'Cucumis melo [Organism] AND Illumina [Platform] AND rna seq [Strategy] AND transcriptomic [Source] AND fastq' | efetch -format runinfo > $Metadata

## Extract accession list.
echo -e "\nGenerate accession list..."
cd $AI"/sra-info/accession_list"
tail -n +2 $Metadata | awk -F "," '{print $1}' > $Acc_list


