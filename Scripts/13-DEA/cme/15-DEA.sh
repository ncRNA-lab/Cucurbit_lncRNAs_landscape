#!/bin/bash

#SBATCH --job-name=cmedea			# Job name.
#SBATCH --output=cme_dea.log			# Standard output and error log.
#SBATCH --qos=short				# Partition (queue)
#SBATCH --ntasks=1				# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2			# Number of tasks = cpus. It depends on the number of process of your parallelization.
#SBATCH --time=1-00:00:00			# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=5gb			# Job memory request.


####### MODULES
module load R/4.2.1
module load biotools

####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/16-DEA/additional_scripts"
alpha_value=0.05

####### NEW AND OTHER VARIABLES
WD1=$WD/05-predict_lncRNAs/$specie
WD2=$WD/06-quantification/$specie
WD3=$WD/08-comparative_genomics
WD4=$WD/09-Tissue-specificity
WD5=$WD/16-DEA/$specie
Metadata=$AI/sra-info/metadata/Search_studies/$specie/Studies

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD/16-DEA
mkdir -p $WD5
mkdir -p $WD5/01-Metadata_EA
mkdir -p $WD5/02-EA
mkdir -p $WD5/03-Metadata_DEA
mkdir -p $WD5/04-DEA
mkdir -p $WD5/05-Tables_and_Figures


####### PIPELINE

echo -e "\n\n#############################"
echo -e "########### STEP 1 ##########"
echo -e "#############################\n\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/Filter_and_prepare_metadata.R $Metadata $WD5 $WD2/Salmon/ALL/nr/03-Quant


echo -e "\n\n#############################"
echo -e "########### STEP 2 ##########"
echo -e "#############################\n\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/Exploratory_analysis.R $WD5 $WD1/STEP-FINAL/Files/Joined/ALL/nr $WD2/Salmon/ALL/nr/03-Quant


echo -e "\n\n#############################"
echo -e "########### STEP 3 ##########"
echo -e "#############################\n\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/DEA.R $WD5 $WD1/STEP-FINAL/Files/Joined/ALL/nr $WD2/Salmon/ALL/nr/03-Quant 0.05


echo -e "\n\n#############################"
echo -e "########### STEP 4 ##########"
echo -e "#############################\n\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript $AS/Tables_and_figures.R $WD5 $WD1/STEP-FINAL/Database/Database_LncRNAs_NR.tsv $WD1/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv $WD3/Positional_level/Approach_2/nr/05-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv $WD4/approach_1/ALL/nr/STEP3/mean-TAU.tsv 0.05 $specie


echo -e "\n\n#############################"
echo -e "############ END ############"
echo -e "#############################\n\n"



