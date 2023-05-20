#!/bin/bash

#SBATCH --job-name=cmecomp			# Job name.
#SBATCH --output=cme_gffcompare.log		# Standard output and error log.
#SBATCH --qos=short				# Partition (queue)
#SBATCH --ntasks=1				# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2			# Number of tasks = cpus. It depends on the number of process of your parallelization.
#SBATCH --time=0-02:00:00			# Time limit days-hrs:min:sec.
#SBATCH --mem=5gb				# Job memory request.


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"

####### NEW AND OTHER VARIABLES
WD1=$WD/03-assembly/$specie
WD2=$WD/04-annotate_transcripts/$specie

####### DIRECTORY
mkdir -p $WD/04-annotate_transcripts
mkdir -p $WD/04-annotate_transcripts/$specie

####### PIPELINE
cd $WD2

### -g 50
echo -e "\n-G 50..."

echo -e "\nAnnotate the transcripts using gffcompare...\n"
## Compare the merged GTF file with the original genome annotation using gffcompare to identify novel transcripts.
gffcompare $WD1/02-Merged_assembly/$specie\_merged_50.gtf -o $specie\_merged_compared_50 -r $AI/GTF_genes/$specie.gtf -V
mv $WD1/02-Merged_assembly/$specie\_merged_compared_50.$specie\_merged_50.gtf.tmap $specie\_merged_compared_50.$specie\_merged_50.gtf.tmap
mv $WD1/02-Merged_assembly/$specie\_merged_compared_50.$specie\_merged_50.gtf.refmap $specie\_merged_compared_50.$specie\_merged_50.gtf.refmap
tail -n +2 $specie\_merged_compared_50.$specie\_merged_50.gtf.tmap | awk '{print $3}' | sort | uniq -c


### -g 250
echo -e "\n-G 250..."

echo -e "\nAnnotate the transcripts using gffcompare...\n"
## Compare the merged GTF file with the original genome annotation using gffcompare to identify novel transcripts.
gffcompare $WD1/02-Merged_assembly/$specie\_merged_250.gtf -o $specie\_merged_compared_250 -r $AI/GTF_genes/$specie.gtf -V
mv $WD1/02-Merged_assembly/$specie\_merged_compared_250.$specie\_merged_250.gtf.tmap $specie\_merged_compared_250.$specie\_merged_250.gtf.tmap
mv $WD1/02-Merged_assembly/$specie\_merged_compared_250.$specie\_merged_250.gtf.refmap $specie\_merged_compared_250.$specie\_merged_250.gtf.refmap
tail -n +2 $specie\_merged_compared_250.$specie\_merged_250.gtf.tmap | awk '{print $3}' | sort | uniq -c

