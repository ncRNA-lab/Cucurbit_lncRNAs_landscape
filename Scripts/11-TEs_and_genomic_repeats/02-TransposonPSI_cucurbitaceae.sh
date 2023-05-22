#!/bin/bash

#SBATCH --job-name=TEPSIa						# Job name.
#SBATCH --output=TransposonPSI_cucurbitaceae.log			# Standard output and error log.
#SBATCH --partition=medium						# Partition (queue)
#SBATCH --ntasks=50							# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=1						# Number of tasks = cpus.
#SBATCH --time=7-00:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=10gb						# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/11-TEs_and_genomic_repeats/Functions.sh"
lines=5000

####### NEW AND OTHER VARIABLES
WD1=$WD/11-TEs_and_genomic_repeats

####### DIRECTORY
mkdir -p $WD/11-TEs_and_genomic_repeats
mkdir -p $WD1/01-Repeat_calling
mkdir -p $WD1/01-Repeat_calling/02-TransposonPSI_cucurbitaceae
mkdir -p $WD1/01-Repeat_calling/02-TransposonPSI_cucurbitaceae/temp

####### PIPELINE

### REPEAT CALLING
## Subdivide cucurbitaceae proteome database.
echo -e "\n\nAccording to the lines, subdivide the cucurbitaceae proteome..."
cd $WD1/01-Repeat_calling/02-TransposonPSI_cucurbitaceae
if [ -d "$WD1/01-Repeat_calling/02-TransposonPSI_cucurbitaceae/temp" ]; then
	rm -r temp
fi
mkdir temp
cd temp
split -a 7 -d -l $lines $AI/cucurbitaceae_proteome/cucurbitaceae_protein_SwissProt_and_RefSeq_mod.fasta
n=$(ls -1 x* | wc -l)
for i in $(seq -s " " -w 0000000 $((n-1))); do
	mkdir folder_$i
	mv x$i ./folder_$i/db_$i.fasta
done

## TransposonPSI.
echo -e "\nIdentifing protein sequences from the cucurbitaceae proteome with homology to proteins encoded by diverse families of transposable elements..."
cd $WD1/01-Repeat_calling/02-TransposonPSI_cucurbitaceae
source activate PSI
for i in $(seq -s " " -w 0000000 $((n-1))); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_transposonPSI $i "prot" &
done
wait
conda deactivate

cat temp/folder_*/*.TPSI.allHits > cucurbitaceae.TPSI.allHits
cat temp/folder_*/*.TPSI.topHits > cucurbitaceae.TPSI.topHits

## Filter the cucurbitaceae proteome file removing the proteins with homology to proteins encoded by diverse families of transposable elements.
echo -e "\nFiltering protein sequences from the cucurbitaceae proteome with homology to proteins encoded by diverse families of transposable elements..."
cd $WD1/01-Repeat_calling/02-TransposonPSI_cucurbitaceae
source activate gaas
awk '{if($0 ~ /^[^\/\/.*]/) print $5}' cucurbitaceae.TPSI.topHits | sort -u > accessions.list
gaas_fasta_removeSeqFromIDlist.pl -f $AI/cucurbitaceae_proteome/cucurbitaceae_protein_SwissProt_and_RefSeq_mod.fasta -l accessions.list -o cucurbitaceae_temp.filtered.fa
conda deactivate

grep ">" cucurbitaceae_temp.filtered.fa | sed 's/>//g' > ids.txt
seqtk subseq cucurbitaceae_temp.filtered.fa ids.txt > cucurbitaceae.filtered.fa
rm cucurbitaceae_temp.filtered.fa
rm ids.txt


