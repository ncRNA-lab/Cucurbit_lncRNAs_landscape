#!/bin/bash

#SBATCH --job-name=BRdef					# Job name.
#SBATCH --output=Blastn_R_definitive.log			# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=25						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=4					# Number of tasks = cpus.
#SBATCH --time=0-02:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb					# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Sequence_level/Blastn/Definitive/Functions_R.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High Medium Low"
evalue=1e-3
Scripts=$(pwd)

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Sequence_level/Blastn/r/Definitive
if [ -d "$WD1" ]; then
	rm -r $WD1
fi

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### SOFTWARES PREDICTION
### BLAST
#https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#http://nebc.nerc.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf
export NCBIBLASTPATH=$SP/ncbi-blast-2.13.0+
export PATH=$PATH:${NCBIBLASTPATH}/bin

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Sequence_level
mkdir -p $WD/08-comparative_genomics/Sequence_level/Blastn
mkdir -p $WD/08-comparative_genomics/Sequence_level/Blastn/r
mkdir -p $WD/08-comparative_genomics/Sequence_level/Blastn/r/Definitive
mkdir -p $WD1/01-LncRNAs
mkdir -p $WD1/02-Makeblastdb
mkdir -p $WD1/03-Blastn
mkdir -p $WD1/04-Reciprocal_hits
mkdir -p $WD1/05-Families



####### PIPELINE

### SELECTION
## For each confidence level, class code and specie, create a TXT file with the identifiers of each lncRNA and a FASTA file with the sequences of the lncRNAs.
cd $WD1/01-LncRNAs

echo -e "\n\nSelect lncRNAs...\n"

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	for class in $Classes_list; do
		mkdir -p $confidence/$class
		for spe in $Species_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Select $spe $confidence $class $WD &
		done
	done
done
wait

### MAKEBLASTDB
## For each confidence level, class code and specie, create a blast database from the FASTA file created in the previous step using the function maskeblastdb.
cd $WD1/02-Makeblastdb

echo -e "\n\nMakeblastdb...\n"

for spe in $Species_list; do
	for confidence in $Confidence_levels_list; do
		mkdir -p $confidence
		for class in $Classes_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Makeblastdb $spe $confidence $class &
		done
	done
done
wait

### BLASTN
## For each confidence level and class code, run blastn for all possible pairwise combinations of the list of species except for combinations with itself.
cd $WD1/03-Blastn

Combinations=$WD1/03-Blastn/Species_combination.txt
set -- $Species_list
for a; do
    shift
    for b; do
        printf "%s %s\n" "$a" "$b"
    done
done > $Combinations

echo -e "\n\nBlastn...\n"

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb_spe=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb_spe | cut -d" " -f1)
	spe2=$(echo $comb_spe | cut -d" " -f2)
	for confidence in $Confidence_levels_list; do
		mkdir -p $confidence
		for class in $Classes_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Blastn $spe1 $spe2 $confidence $class $evalue $SLURM_CPUS_PER_TASK $j $n_comb $WD1 $AS &
		done
	done
done
wait

### RECIPROCAL HITS
## For each confidence level and class code, find the reciprocal blast hits between species, i.e. those hits that can be found in both the spe1-spe2 and spe2-spe1 combination. To do this, we use the home-made script Find_reciprocal_hits.py.
cd $WD1/04-Reciprocal_hits

echo -e "\n\nReciprocal hits...\n"

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb_spe=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb_spe | cut -d" " -f1)
	spe2=$(echo $comb_spe | cut -d" " -f2)
	for confidence in $Confidence_levels_list; do
		mkdir -p $confidence
		for class in $Classes_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Reciprocal_hits $spe1 $spe2 $confidence $class $WD1 $AS $j $n_comb &
		done
	done
done
wait

## For each confidence level and class code, pool the results of all paired species combinations.
for confidence in $Confidence_levels_list; do
	for class in $Classes_list; do
		if [ -f $confidence/$class/Reciprocal_hits.tsv ] ; then
			rm $confidence/$class/Reciprocal_hits.tsv
		fi
		cat $confidence/$class/* > $confidence/$class/Reciprocal_hits.tsv
	done
done

### FAMILIES
## For each confidence level and class code, classify the conserved lncRNAs into conservation families using home-made script Classify_into_families.py.
cd $WD1/05-Families

echo -e "\n\nClassify lncRNAs into families...\n"

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	for class in $Classes_list; do
		mkdir -p $confidence/$class
		>$confidence/$class/ids_by_specie.tsv
		for spe in $Species_list; do
			awk -v a="$spe" '{print a"\t"$0}' ../01-LncRNAs/$confidence/$class/$spe\_ids.txt >> $confidence/$class/ids_by_specie.tsv
		done
	done
done

for confidence in $Confidence_levels_list; do
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Classify_into_families $confidence $class $WD1 $AS $SLURM_CPUS_PER_TASK &
	done
done
wait

### FIGURES
echo -e "\n\nFigures...\n"
Rscript $Scripts/Figures_R.R


