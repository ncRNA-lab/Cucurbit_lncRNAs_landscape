#!/bin/bash

#SBATCH --job-name=cme_IR							# Job name.
#SBATCH --output=cme_Intergenic_regions.log					# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=1								# Run on one mode.
#SBATCH --cpus-per-task=2							# Number of tasks = cpus.
#SBATCH --time=0-01:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=5gb							# Job memory request.


####### MODULES
module load biotools

####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/07-Comparison_lncRNAs_vs_coding_genes"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2_spe


####### PIPELINE
cd $WD2_spe

## Convert the LncRNAs GTF file to BED file.
cat $WD1_spe/STEP1/Potential_lncRNAs/POTENTIAL_LNCRNAS.gtf | convert2bed -i gtf --attribute-key=transcript_id | awk '$8 == "transcript" {print $0}' > LncRNAs.bed
## Convert the Genes GTF file to BED file.
cat $WD1_spe/STEP1/Original_genes/ORIGINAL_GENES.gtf | convert2bed -i gtf --attribute-key=transcript_id | awk '$8 == "transcript" {print $0}' > Genes.bed
## Get the chromosomes sizes.
cp $AI/Genome/$specie.fa $specie.fa
samtools faidx $specie.fa
cut -f1-2 $specie.fa.fai > Genome_chr_length.txt
## Create a random intergenic regions.
bedtools random -l 500 -n 25000 -seed 53458 -g Genome_chr_length.txt > Random_IR.bed
## Filter the random intergenic regions which oberlap with potential lncRNAs and genes by bedtools intersect.
bedtools intersect -a Random_IR.bed -b LncRNAs.bed -s -wao -nonamecheck | awk -F"\t" 'int($17) == 0 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > Random_IR_filt_1.bed
bedtools intersect -a Random_IR_filt_1.bed -b Genes.bed -s -wao -nonamecheck | awk -F"\t" 'int($17) == 0 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > Random_IR_filt_2.bed
rm Random_IR.bed Random_IR_filt_1.bed
## Filter repeated IR.
awk -F"\t" '!seen[$1,$2,$3,$6]++' Random_IR_filt_2.bed > Random_IR.bed
rm Random_IR_filt_2.bed
## Convert the bed file to fasta
bedtools getfasta -fi $specie.fa -bed Random_IR.bed -s -fo Random_IR.fasta
rm $specie.fa
rm $specie.fa.fai
