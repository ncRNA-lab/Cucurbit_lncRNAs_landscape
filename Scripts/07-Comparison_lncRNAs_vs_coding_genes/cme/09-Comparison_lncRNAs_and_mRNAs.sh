#!/bin/bash

#SBATCH --job-name=cmecomp			# Job name.
#SBATCH --output=cme_comparison.log		# Standard output and error log.
#SBATCH --qos=short				# Partition (queue)
#SBATCH --ntasks=1				# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2			# Number of tasks = cpus. It depends on the number of process of your parallelization.
#SBATCH --time=0-02:00:00			# Time limit days-hrs:min:sec.
#SBATCH --mem=15gb				# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"

####### NEW AND OTHER VARIABLES
WD1=$WD/05-predict_lncRNAs/$specie
WD2=$WD/06-quantification/$specie
WD3=$WD/07-comparison_lncRNAs_vs_coding_genes/$specie

####### DIRECTORY
mkdir -p $WD/07-comparison_lncRNAs_vs_coding_genes
mkdir -p $WD/07-comparison_lncRNAs_vs_coding_genes/$specie

####### PIPELINE
echo -e "\n\nREDUNDANT SEQUENCES: CREATING FIGURES\n"
Rscript 02-Comparison_lncRNAs_and_mRNAs_R.R $WD1 $WD2 $WD3
echo -e "\n\nNON-REDUNDANT SEQUENCES: CREATING FIGURES\n"
Rscript 02-Comparison_lncRNAs_and_mRNAs_NR.R $WD1 $WD2 $WD3

