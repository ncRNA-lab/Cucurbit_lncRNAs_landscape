#!/bin/bash

#SBATCH --job-name=cmeS6						# Job name.
#SBATCH --output=cme_STEP6.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=1							# Run on one mode.
#SBATCH --cpus-per-task=2						# Number of tasks = cpus.
#SBATCH --time=0-02:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem=5gb							# Job memory request.


####### VARIABLES
specie="cme"
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/03-Assembly"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/04-Transcript_annotation"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"

####### NEW AND OTHER VARIABLES
WD1_spe=$WD1/$specie
WD2_spe=$WD2/$specie

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2_spe


####### PIPELINE: STEP 6

### TRANSCRIPT ANNOTATION
## Annotate each transcript according to the gene's position using gffcompare.
echo -e "\nTRANSCRIPT ANNOTATION..."

cd $WD2_spe
mkdir -p Outputs
>./Outputs/stdout_Annotation.log

## Compare the merged GTF file with the original genome annotation using gffcompare to identify novel transcripts.
gffcompare $WD1_spe/02-Merged_assembly/$specie\_merged.gtf -o $specie\_merged_compared -r $AI/GTF_genes/$specie.gtf -V >> ./Outputs/stdout_Annotation.log 2>&1
mv $WD1_spe/02-Merged_assembly/$specie\_merged_compared.$specie\_merged.gtf.tmap ./
mv $WD1_spe/02-Merged_assembly/$specie\_merged_compared.$specie\_merged.gtf.refmap ./

echo -e "\nResults:\n"
tail -n +2 $specie\_merged_compared.$specie\_merged.gtf.tmap | awk '{print $3}' | sort | uniq -c


