#!/bin/bash

#SBATCH --job-name=TEPSIb						# Job name.
#SBATCH --output=TransposonPSI_unknown_repeats.log			# Standard output and error log.
#SBATCH --partition=medium						# Partition (queue)
#SBATCH --ntasks=40							# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=1						# Number of tasks = cpus.
#SBATCH --time=4-00:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=8gb						# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/11-TEs_and_genomic_repeats/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
lines=40

####### NEW AND OTHER VARIABLES
WD1=$WD/11-TEs_and_genomic_repeats

####### DIRECTORY
mkdir -p $WD/11-TEs_and_genomic_repeats
mkdir -p $WD1/01-Repeat_calling
mkdir -p $WD1/01-Repeat_calling/03-TransposonPSI_unknown_repeats

####### PIPELINE

### REPEAT CALLING
## Subdivide the customed repeat library of each specie.
echo -e "\n\nDivide the customed repeat library file in unknown and known repeat. Then, according to the lines, subdivide the unknown repeat library of each specie in fragments to be faster..."
cd $WD1/01-Repeat_calling/03-TransposonPSI_unknown_repeats
for spe in $Species_list; do
	echo -e "\t-"$spe
	if [ -d "$WD1/01-Repeat_calling/03-TransposonPSI_unknown_repeats/$spe" ]; then
		rm -r $spe
	fi
	mkdir $spe
	cd $spe
	
	# Separate the repeat libraries: Known and Unknown.
	repeatmodeler_parse.pl \
		--fastafile $WD1/01-Repeat_calling/01-RepeatModeler/$spe/$spe-families.fa \
		--identities RM_identities_temp.fa \
		--unknowns RM_unknowns_temp.fa >> repeatmodeler_parse_$spe.log 2>&1
	
	# Modify repeat libraries to work.
	grep ">" RM_identities_temp.fa | sed 's/>//g' > ids.txt
	seqtk subseq RM_identities_temp.fa ids.txt > RM_identities.fa
	rm RM_identities_temp.fa
	rm ids.txt
	
	grep ">" RM_unknowns_temp.fa | sed 's/>//g' > ids.txt
	seqtk subseq RM_unknowns_temp.fa ids.txt > RM_unknowns.fa
	rm RM_unknowns_temp.fa
	rm ids.txt
	
	sed 's/ /|/g' RM_unknowns.fa > RM_unknowns_mod.fa
	
	grep ">ltr-" $WD1/01-Repeat_calling/01-RepeatModeler/$spe/$spe-families.fa | sed 's/>//g' > ids.txt
	seqtk subseq $WD1/01-Repeat_calling/01-RepeatModeler/$spe/$spe-families.fa ids.txt > RM_LTRStruct.fa
	rm ids.txt
	
	# Split the unknown repeat library.
	mkdir temp
	cd temp
	split -a 5 -d -l $lines ../RM_unknowns_mod.fa
	n=$(ls -1 x* | wc -l)
	for i in $(seq -s " " -w 00000 $((n-1))); do
		mkdir folder_$i
		mv x$i ./folder_$i/db_$i.fasta
	done
	
	cd ../../
done

## TransposonPSI.
echo -e "\nIdentifing sequences from the unknown repeat library (RM_unknowns.fasta) of each specie with homology to proteins encoded by diverse families of transposable elements..."
cd $WD1/01-Repeat_calling/03-TransposonPSI_unknown_repeats
source activate PSI
for spe in $Species_list; do
	echo -e "\t-"$spe
	cd $spe
	
	# Execute TransposonPSI.
	n=$(ls -1 temp/ | grep "folder_*" | wc -l)
	for i in $(seq -s " " -w 00000 $((n-1))); do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_transposonPSI $i "nuc" &
	done
	wait
	
	# Join the results.
	cat temp/folder_*/*.TPSI.allHits > $spe.TPSI.allHits
	
	cd ..
done
conda deactivate

## Filter.
echo -e "\nFiltering sequences from the unknown repeat library (RM_unknowns.fasta) of each specie with homology to proteins encoded by diverse families of transposable elements..."
cd $WD1/01-Repeat_calling/03-TransposonPSI_unknown_repeats
for spe in $Species_list; do
	echo -e "\t-"$spe
	cd $spe
	
	# Remove the sequences from the unknown repeat library with match on TransposonPSI database.
	awk '$16 == "Plus" {print $5$14}' $spe.TPSI.allHits | sort -u > accessions.list
	source activate gaas
	gaas_fasta_removeSeqFromIDlist.pl -f RM_unknowns_mod.fa -l accessions.list -o RM_unknowns_temp.filtered.fa
	conda deactivate
	grep ">" RM_unknowns_temp.filtered.fa | sed 's/>//g' > ids.txt
	seqtk subseq RM_unknowns_temp.filtered.fa ids.txt | sed 's/|/ /g' > RM_unknowns.filtered.fa
	rm RM_unknowns_temp.filtered.fa
	rm ids.txt
	
	# Retrieve the sequences from the unknown repeat library with match.
	seqtk subseq RM_unknowns_mod.fa accessions.list | sed 's/|/ /g' > RM_unknowns_matched_on_transposonPSI.fa
	
	# Remove the initial unknown repeat library with modified headers.
	rm RM_unknowns_mod.fa
	
	cd ..
done



