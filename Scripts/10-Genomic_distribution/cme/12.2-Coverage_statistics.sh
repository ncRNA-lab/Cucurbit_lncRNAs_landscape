#!/bin/bash

#SBATCH --job-name=cmeS12cov						# Job name.
#SBATCH --output=cme_STEP12_coverage.log				# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus. 
#SBATCH --time=0-04:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb						# Job memory request.


####### MODULES
module load R/4.1.2
module load biotools

####### VARIABLES
specie="cme"
specie_long="C. melo"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/03-Assembly"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD3="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/10-Genomic_distribution"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/10-Genomic_distribution/Additional_scripts"
flag_list="nr"
confidence_list="High Medium Low"

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
mkdir -p $WD3_spe/Cov


####### PIPELINE

### COVERAGE STATISTICS
echo -e "\nCOVERAGE STATISTICS..."

for flag in $flag_list; do

	echo -e "\nFLAG: "$flag
	
	mkdir -p $WD3_spe/Cov/$flag
	cd $WD3_spe/Cov/$flag

	# Chromosome sizes.
	cp $AI/Genome/$specie.fa ./ 
	samtools faidx $specie.fa
	cut -f1,2 $specie.fa.fai | awk '{print $1"\t"$2}' > $specie.sizes_genome.txt
	rm $specie.fa $specie.fa.fai

	
	## TRANSCRIPTOME
	echo -e "\n\t- TRANSCRIPTOME..."

	# Create bed sorted file.
	awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' $WD1_spe/02-Merged_assembly/$specie\_merged.gtf > TRANS.bed
	sort -k 1,1 TRANS.bed > TRANS-sorted.bed
	# Calculate coverage
	bedtools genomecov -i TRANS-sorted.bed -g $specie.sizes_genome.txt > TRANS-coverage.tsv	
	rm TRANS.bed TRANS-sorted.bed


	## GENES
	echo -e "\n\t- GENES..."

	# Create bed sorted file.
	awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' $WD2_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf > GENES.bed
	sort -k 1,1 GENES.bed > GENES-sorted.bed
	# Calculate coverage
	bedtools genomecov -i GENES-sorted.bed -g $specie.sizes_genome.txt > GENES-coverage.tsv	
	rm GENES.bed GENES-sorted.bed


	## LNCRNAS
	echo -e "\n\t- LNCRNAS..."

	# Create bed sorted file.
	awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' $WD2_spe/STEP-FINAL/Files/LncRNAs/$flag/POTENTIAL_LNCRNAS_pred.gtf > LNCRNAS.bed
	sort -k 1,1 LNCRNAS.bed > LNCRNAS-sorted.bed
	# Calculate coverage
	bedtools genomecov -i LNCRNAS-sorted.bed -g $specie.sizes_genome.txt > LNCRNAS-coverage.tsv
	rm LNCRNAS.bed LNCRNAS-sorted.bed
	
	
	## LNCRNAS BY CONFIDENCE LEVEL
	echo -e "\n\t- LNCRNAS BY CONFIDENCE LEVEL..."

	for co in $confidence_list; do
		echo -e "\t\t- "$co"..."
		# Generate ids lists.
		tail -n +2 $WD2_spe/STEP-FINAL/Database/Database_LncRNAs_${flag^^}.tsv | awk -v var=$co '$31 == var {print $1}' > LNCRNAS-$co\-ids.txt
		# Generate GTF file.
		Filter_GTF.py \
			--gtf-initial $WD2_spe/STEP-FINAL/Files/LncRNAs/$flag/POTENTIAL_LNCRNAS_pred.gtf \
			--gtf-final $WD3_spe/Cov/$flag/LNCRNAS-$co.gtf \
			--ids $WD3_spe/Cov/$flag/LNCRNAS-$co\-ids.txt
		# Create bed sorted file.
		awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' LNCRNAS-$co.gtf > LNCRNAS-$co.bed
		sort -k 1,1 LNCRNAS-$co.bed > LNCRNAS-$co\-sorted.bed
		# Calculate coverage
		bedtools genomecov -i LNCRNAS-$co\-sorted.bed -g $specie.sizes_genome.txt > LNCRNAS-$co\-coverage.tsv
		rm LNCRNAS-$co\-ids.txt LNCRNAS-$co.gtf LNCRNAS-$co.bed LNCRNAS-$co\-sorted.bed
	done

	
	## LNCRNAS AND GENES
	echo -e "\n\t- LNCRNAS AND GENES..."

	# Generate ids lists.
	tail -n +2 $WD2_spe/STEP-FINAL/Database/Database_LncRNAs_${flag^^}.tsv | awk '{print $1}' > LNCRNAS-ids.txt
	cp $WD2_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES_ids.txt GENES-ids.txt
	cat LNCRNAS-ids.txt GENES-ids.txt > BOTH-ids.txt
	# Generate GTF file.
	Filter_GTF.py \
		--gtf-initial $WD2_spe/STEP-FINAL/Files/ALL/$flag/ALL.gtf \
		--gtf-final $WD3_spe/Cov/$flag/BOTH.gtf \
		--ids $WD3_spe/Cov/$flag/BOTH-ids.txt
	# Create bed sorted file.
	awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' BOTH.gtf > BOTH.bed
	sort -k 1,1 BOTH.bed > BOTH-sorted.bed
	# Calculate coverage
	bedtools genomecov -i BOTH-sorted.bed -g $specie.sizes_genome.txt > BOTH-coverage.tsv
	rm LNCRNAS-ids.txt GENES-ids.txt BOTH-ids.txt BOTH.gtf BOTH.bed BOTH-sorted.bed
	
	
	## LNCRNAS AND GENES BY CONFIDENCE LEVEL
	echo -e "\n\t- LNCRNAS AND GENES BY CONFIDENCE LEVEL..."

	for co in $confidence_list; do
		echo -e "\t\t- "$co"..."
		# Generate ids lists.
		tail -n +2 $WD2_spe/STEP-FINAL/Database/Database_LncRNAs_${flag^^}.tsv | awk -v var=$co '$31 == var {print $1}' > LNCRNAS-$co\-ids.txt
		cp $WD2_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES_ids.txt GENES-ids.txt
		cat LNCRNAS-$co\-ids.txt GENES-ids.txt > BOTH-$co\-ids.txt
		# Generate GTF file.
		Filter_GTF.py \
			--gtf-initial $WD2_spe/STEP-FINAL/Files/ALL/$flag/ALL.gtf \
			--gtf-final $WD3_spe/Cov/$flag/BOTH-$co.gtf \
			--ids $WD3_spe/Cov/$flag/BOTH-$co\-ids.txt
		# Create bed sorted file.
		awk -F"\t" '$3 == "transcript" {print $1"\t"$4"\t"$5}' BOTH-$co.gtf > BOTH-$co.bed
		sort -k 1,1 BOTH-$co.bed > BOTH-$co\-sorted.bed
		# Calculate coverage
		bedtools genomecov -i BOTH-$co\-sorted.bed -g $specie.sizes_genome.txt > BOTH-$co\-coverage.tsv
		rm LNCRNAS-$co\-ids.txt GENES-ids.txt BOTH-$co\-ids.txt BOTH-$co.gtf BOTH-$co.bed BOTH-$co\-sorted.bed
	done
	
	## FIGURES
	echo -e "\n\t- CREATE FINAL TABLES AND FIGURES..."
	Rscript $AS/Coverage-Statistics.R $specie "$specie_long" $WD3_spe/Cov/$flag $AI "$confidence_list"

done

