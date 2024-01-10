#!/bin/bash

#SBATCH --job-name=covvvi					# Job name.
#SBATCH --output=coverage_vvi.log				# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=1						# Run on one mode.
#SBATCH --cpus-per-task=2					# Number of tasks = cpus. 
#SBATCH --time=0-04:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb					# Job memory request.


####### VARIABLES
specie="vvi"
WD1="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results/03-Assembly"
WD2="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results/05-LncRNAs_prediction"
WD3="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results/10-Distribution"
AI="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Scripts/10-Distribution/Additional_scripts"
Confidence_levels_list="High Medium Low"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie
WD3_spe=$WD3/$specie

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD3
mkdir -p $WD3_spe
mkdir -p $WD3_spe/nr


####### PIPELINE

cd $WD3_spe/nr

# Chromosome sizes.
awk -F"\t" '{print $1"\t"$3}' $AI/Chromosomes/$specie\_sizes_genome.txt > $specie\_sizes_genome.txt	

### TRANSCRIPTOME

echo -e "\n\nTRANSCRIPTOME...\n"

# Create bed sorted file.
awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' $WD1_spe/02-Merged_assembly/$specie\_merged_50.gtf > TRANS-$specie.bed
sort -k 1,1 TRANS-$specie.bed > TRANS-$specie\_sorted.bed
# Calculate coverage
bedtools genomecov -i TRANS-$specie\_sorted.bed -g $specie\_sizes_genome.txt > TRANS-$specie\_coverage.tsv	
rm TRANS-$specie.bed TRANS-$specie\_sorted.bed


### GENES

echo -e "\n\nGENES...\n"

# Create bed sorted file.
awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' $WD2_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf > GENES-$specie.bed
sort -k 1,1 GENES-$specie.bed > GENES-$specie\_sorted.bed
# Calculate coverage
bedtools genomecov -i GENES-$specie\_sorted.bed -g $specie\_sizes_genome.txt > GENES-$specie\_coverage.tsv	
rm GENES-$specie.bed GENES-$specie\_sorted.bed


### LNCRNAS

echo -e "\n\nLNCRNAS BY CONFIDENCE LEVEL...\n"

for co in $Confidence_levels_list; do
	echo -e $co"..."
	# Generate ids lists.
	tail -n +2 $WD2_spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$co '$31 == var {print $1}' > LNCRNAS-$co\-$specie\_ids.txt
	# Generate GTF file.
	Filter_GTF.py \
		--gtf-initial $WD2_spe/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf \
		--gtf-final $WD3_spe/nr/LNCRNAS-$co\-$specie.gtf \
		--ids $WD3_spe/nr/LNCRNAS-$co\-$specie\_ids.txt
	# Create bed sorted file.
	awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' LNCRNAS-$co\-$specie.gtf > LNCRNAS-$co\-$specie.bed
	sort -k 1,1 LNCRNAS-$co\-$specie.bed > LNCRNAS-$co\-$specie\_sorted.bed
	# Calculate coverage
	bedtools genomecov -i LNCRNAS-$co\-$specie\_sorted.bed -g $specie\_sizes_genome.txt > LNCRNAS-$co\-$specie\_coverage.tsv
	rm LNCRNAS-$co\-$specie\_ids.txt LNCRNAS-$co\-$specie.gtf LNCRNAS-$co\-$specie.bed LNCRNAS-$co\-$specie\_sorted.bed
done

echo -e "\n\nLNCRNAS...\n"

# Create bed sorted file.
awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' $WD2_spe/STEP-FINAL/Files/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf > LNCRNAS-$specie.bed
sort -k 1,1 LNCRNAS-$specie.bed > LNCRNAS-$specie\_sorted.bed
# Calculate coverage
bedtools genomecov -i LNCRNAS-$specie\_sorted.bed -g $specie\_sizes_genome.txt > LNCRNAS-$specie\_coverage.tsv
rm LNCRNAS-$specie.bed LNCRNAS-$specie\_sorted.bed


### LNCRNAS AND GENES

echo -e "\n\nLNCRNAS AND GENES BY CONFIDENCE LEVEL...\n"

for co in $Confidence_levels_list; do
	echo -e $co"..."
	# Generate ids lists.
	tail -n +2 $WD2_spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk -v var=$co '$31 == var {print $1}' > LNCRNAS-$co\-$specie\_ids.txt
	cp $WD2_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES_ids.txt GENES-$specie\_ids.txt
	cat LNCRNAS-$co\-$specie\_ids.txt GENES-$specie\_ids.txt > BOTH-$co\-$specie\_ids.txt
	# Generate GTF file.
	Filter_GTF.py \
		--gtf-initial $WD2_spe/STEP-FINAL/Files/Joined/ALL/nr/ALL.gtf \
		--gtf-final $WD3_spe/nr/BOTH-$co\-$specie.gtf \
		--ids $WD3_spe/nr/BOTH-$co\-$specie\_ids.txt
	# Create bed sorted file.
	awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' BOTH-$co\-$specie.gtf > BOTH-$co\-$specie.bed
	sort -k 1,1 BOTH-$co\-$specie.bed > BOTH-$co\-$specie\_sorted.bed
	# Calculate coverage
	bedtools genomecov -i BOTH-$co\-$specie\_sorted.bed -g $specie\_sizes_genome.txt > BOTH-$co\-$specie\_coverage.tsv
	rm LNCRNAS-$co\-$specie\_ids.txt GENES-$specie\_ids.txt BOTH-$co\-$specie\_ids.txt BOTH-$co\-$specie.gtf BOTH-$co\-$specie.bed BOTH-$co\-$specie\_sorted.bed
done

echo -e "\n\nLNCRNAS AND GENES...\n"

# Generate ids lists.
tail -n +2 $WD2_spe/STEP-FINAL/Database/Database_LncRNAs_NR.tsv | awk '{print $1}' > LNCRNAS-$specie\_ids.txt
cp $WD2_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES_ids.txt GENES-$specie\_ids.txt
cat LNCRNAS-$specie\_ids.txt GENES-$specie\_ids.txt > BOTH-$specie\_ids.txt
# Generate GTF file.
Filter_GTF.py \
	--gtf-initial $WD2_spe/STEP-FINAL/Files/Joined/ALL/nr/ALL.gtf \
	--gtf-final $WD3_spe/nr/BOTH-$specie.gtf \
	--ids $WD3_spe/nr/BOTH-$specie\_ids.txt
# Create bed sorted file.
awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' BOTH-$specie.gtf > BOTH-$specie.bed
sort -k 1,1 BOTH-$specie.bed > BOTH-$specie\_sorted.bed
# Calculate coverage
bedtools genomecov -i BOTH-$specie\_sorted.bed -g $specie\_sizes_genome.txt > BOTH-$specie\_coverage.tsv
rm LNCRNAS-$specie\_ids.txt GENES-$specie\_ids.txt BOTH-$specie\_ids.txt BOTH-$specie.gtf BOTH-$specie.bed BOTH-$specie\_sorted.bed

