#!/bin/bash

#SBATCH --job-name=cmeS9				# Job name.
#SBATCH --output=cme_STEP9.log			# Standard output and error log.
#SBATCH --qos=short					# Partition (queue)
#SBATCH --ntasks=1					# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=6				# Number of tasks = cpus.
#SBATCH --time=0-02:00:00				# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb				# Job memory request.


####### VARIABLES
specie="cme"
WD_first_assembly="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
WD_second_assembly="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"

####### NEW AND OTHER VARIABLES
WD1_first_assembly=$WD_first_assembly"/05-predict_lncRNAs/"$specie"/STEP-FINAL/Files/LncRNAs"
WD1_second_assembly=$WD_second_assembly"/06-LncRNAs_prediction/STEP-FINAL/Files/LncRNAs"
WD2=$WD_second_assembly"/06-LncRNAs_prediction/STEP-FINAL/Comparison_2022_01_24_landscape"

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2/r
mkdir -p $WD2/nr

####### PIPELINE: STEP 9

echo -e "\n\n\n####################"
echo -e "###### STEP 9 ######"
echo -e "####################\n"

## REDUNDANT
echo -e "\nREDUNDANT..."

cd $WD2/r

gffcompare $WD1_second_assembly/r/POTENTIAL_LNCRNAS_pred.gtf -o Compared -r $WD1_first_assembly/r/POTENTIAL_LNCRNAS_pred.gtf -V
mv $WD1_second_assembly/r/Compared.POTENTIAL_LNCRNAS_pred.gtf.tmap Compared.POTENTIAL_LNCRNAS_pred.gtf.tmap
mv $WD1_second_assembly/r/Compared.POTENTIAL_LNCRNAS_pred.gtf.refmap Compared.POTENTIAL_LNCRNAS_pred.gtf.refmap
tail -n +2 Compared.POTENTIAL_LNCRNAS_pred.gtf.tmap | awk '{print $3}' | sort | uniq -c

## NON-REDUNDANT
echo -e "\nNON-REDUNDANT..."

cd $WD2/nr

gffcompare $WD1_second_assembly/nr/POTENTIAL_LNCRNAS_pred.gtf -o Compared -r $WD1_first_assembly/nr/POTENTIAL_LNCRNAS_pred.gtf -V
mv $WD1_second_assembly/nr/Compared.POTENTIAL_LNCRNAS_pred.gtf.tmap Compared.POTENTIAL_LNCRNAS_pred.gtf.tmap
mv $WD1_second_assembly/nr/Compared.POTENTIAL_LNCRNAS_pred.gtf.refmap Compared.POTENTIAL_LNCRNAS_pred.gtf.refmap
tail -n +2 Compared.POTENTIAL_LNCRNAS_pred.gtf.tmap | awk '{print $3}' | sort | uniq -c


