#!/bin/bash

#SBATCH --job-name=cmeS16U1						# Job name.
#SBATCH --output=cme_STEP16_U1.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-00:45:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=3gb						# Job memory request.


####### MODULES
module load biotools

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/14-Expression_correlation/Additional_scripts"
dist_list="500 1000 2000 5000 10000 20000 50000 100000"
flag_list="nr"

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### NEW VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2_spe
mkdir -p $WD2_spe/intergenic

####### PIPELINE

### CORRELATION
echo -e "\nCORRELATION..."

echo -e "\n\n##############################"
echo -e "########### STEP 1 ###########"
echo -e "##############################\n\n"

for flag in $flag_list; do

	echo -e "\nFLAG: "$flag
	
	## Variable.
	O=$WD2_spe/intergenic/$flag
	O1=$WD2_spe/intergenic/$flag/STEP1
	
	## Directory.
	mkdir -p $O
	mkdir -p $O1
	
	cd $O1
	
	## Select intergenic lncRNAs.
	tail -n +2 $WD1_spe/STEP-FINAL/Database/Database_LncRNAs_${flag^^}.tsv | awk -F"\t" '$8 == "u" {print $1}' > LncRNAs_ids.txt
	Filter_GTF.py \
	--gtf-initial $WD1_spe/STEP-FINAL/Files/LncRNAs/$flag/POTENTIAL_LNCRNAS_pred.gtf \
	--gtf-final $O1/LncRNAs.gtf \
	--ids $O1/LncRNAs_ids.txt
	## Convert the LncRNAs GTF file to BED file.
	cat LncRNAs.gtf | convert2bed -i gtf --attribute-key=transcript_id | awk -F"\t" '$8 == "transcript" {print $0}' > LncRNAs_temp.bed
	sort -k1,1 -k2,2n LncRNAs_temp.bed > LncRNAs.bed
	rm LncRNAs_temp.bed
	## Convert the Genes GTF file to BED file.
	cat $WD1_spe/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf | convert2bed -i gtf --attribute-key=transcript_id | awk -F"\t" '$8 == "transcript" {print $0}' > Genes_temp.bed
	sort -k1,1 -k2,2n Genes_temp.bed > Genes.bed
	rm Genes_temp.bed
	## Get the chromosomes sizes.
	cp $AI/Genome/$specie.fa $specie.fa
	samtools faidx $specie.fa
	cut -f1-2 $specie.fa.fai > Genome_chr_length.txt

	
	###### LncRNAs-Genes.

	echo -e "-LncRNAs-Genes..."

	echo -e "\t-Closest..."

	## Closest
	bedtools closest -a LncRNAs.bed -b Genes.bed -D a -t all -k 1 -nonamecheck | awk -F"\t" '$14 != "." {print $1"\t"$4"\t"$2"\t"$3"\t"$6"\t"$14"\t"$12"\t"$13"\t"$16"\t"$21}' > LncRNA_Gene_closest.tsv

	echo -e "\t-Range..."

	## Range: Intersect LncRNAs BED file with ranges and Genes BED file.
	for dist in $dist_list; do
		echo -e "\t\t-"$dist"..."
		## Create a LncRNAs BED file with X kb downstream/upstream ranges.
		bedtools slop -i LncRNAs.bed -g Genome_chr_length.txt -b $(($dist)) > LncRNAs_range_$dist.bed
		## Intersect LncRNAs BED file with ranges and Genes BED file.
		bedtools intersect -a LncRNAs_range_$dist.bed -b Genes.bed -wao -nonamecheck | awk -F"\t" 'int($21) != 0 {print $1"\t"$4"\t"$2"\t"$3"\t"$6"\t"$14"\t"$12"\t"$13"\t"$16"\t"$21}' > LncRNA_Gene_cis_interactions_range_$dist.tsv
		rm LncRNAs_range_$dist.bed
	done

	
	###### Genes-Genes.

	echo -e "-Genes-Genes..."

	echo -e "\t-Closest..."

	## Closest.
	bedtools closest -a Genes.bed -b Genes.bed -io -D a -t all -k 1 -nonamecheck | awk -F"\t" '$14 != "." {print $1"\t"$4"\t"$2"\t"$3"\t"$6"\t"$14"\t"$12"\t"$13"\t"$16"\t"$21}' > Gene_Gene_closest.tsv

	echo -e "\t-Range..."

	## Range: Intersect Genes BED file with ranges and Genes BED file.
	for dist in $dist_list; do
		echo -e "\t\t-"$dist"..."
		## Create a Genes BED file with X kb downstream/upstream ranges.
		bedtools slop -i Genes.bed -g Genome_chr_length.txt -b $(($dist)) > Genes_range_$dist.bed
		## Intersect Genes BED file with ranges and Genes BED file.
		bedtools intersect -a Genes_range_$dist.bed -b Genes.bed -wao -nonamecheck | awk -F"\t" 'int($21) != 0 && $4 != $14 {print $1"\t"$4"\t"$2"\t"$3"\t"$6"\t"$14"\t"$12"\t"$13"\t"$16"\t"$21}' > Gene_Gene_cis_interactions_range_$dist.tsv
		rm Genes_range_$dist.bed
	done

	## Remove.
	rm Genome_chr_length.txt $specie.fa $specie.fa.fai
	
done


