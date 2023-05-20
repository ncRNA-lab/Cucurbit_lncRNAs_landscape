#!/bin/bash

#SBATCH --job-name=cmeS7			# Job name.
#SBATCH --output=cme_STEP7.log		# Standard output and error log.
#SBATCH --qos=short				# Partition (queue)
#SBATCH --ntasks=1				# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2			# Number of tasks = cpus. It depends on the number of process of your parallelization.
#SBATCH --time=0-02:00:00			# Time limit days-hrs:min:sec.
#SBATCH --mem=5gb				# Job memory request.


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"

####### NEW AND OTHER VARIABLES
WD1=$WD"/04-Merged_assembly"
WD2=$WD"/05-Transcript_annotation"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD2


####### PIPELINE: STEP 7

echo -e "\n\n\n####################"
echo -e "###### STEP 7 ######"
echo -e "####################\n"

### TRANSCRIPT ANNOTATION
## Annotate each transcript according to the gene's position using gffcompare.
echo -e "\nTRANSCRIPT ANNOTATION..."

cd $WD2
mkdir -p outputs

## -g 250
echo -e "\n-g 250..."
>./outputs/stdout_Annotation_250.log
## Compare the merged GTF file with the original genome annotation using gffcompare to identify novel transcripts.
gffcompare $WD1/$specie\_merged_250.gtf -o $specie\_merged_compared_250 -r $AI/GTF_genes/$specie.gtf -V >> ./outputs/stdout_Annotation_250.log 2>&1
mv $WD1/$specie\_merged_compared_250.$specie\_merged_250.gtf.tmap ./
mv $WD1/$specie\_merged_compared_250.$specie\_merged_250.gtf.refmap ./
echo -e "\nResults:\n"
tail -n +2 $specie\_merged_compared_250.$specie\_merged_250.gtf.tmap | awk '{print $3}' | sort | uniq -c

## -g 50
echo -e "\n-g 50..."
>./outputs/stdout_Annotation_50.log
## Compare the merged GTF file with the original genome annotation using gffcompare to identify novel transcripts.
gffcompare $WD1/$specie\_merged_50.gtf -o $specie\_merged_compared_50 -r $AI/GTF_genes/$specie.gtf -V >> ./outputs/stdout_Annotation_50.log 2>&1
mv $WD1/$specie\_merged_compared_50.$specie\_merged_50.gtf.tmap ./
mv $WD1/$specie\_merged_compared_50.$specie\_merged_50.gtf.refmap ./
echo -e "\nResults:\n"
tail -n +2 $specie\_merged_compared_50.$specie\_merged_50.gtf.tmap | awk '{print $3}' | sort | uniq -c


