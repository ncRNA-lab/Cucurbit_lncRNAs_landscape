#!/bin/bash

#SBATCH --job-name=Bgen					# Job name.
#SBATCH --output=Blastn_genes.log				# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=25						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=4					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb					# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Sequence_level/Blastn/Prueba_genes/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
evalue=1e-3
Scripts=$(pwd)

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Sequence_level/Blastn/nr/Prueba_genes
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
mkdir -p $WD/08-comparative_genomics/Sequence_level/Blastn/nr
mkdir -p $WD/08-comparative_genomics/Sequence_level/Blastn/nr/Prueba_genes
mkdir -p $WD1/01-Genes
mkdir -p $WD1/02-Makeblastdb
mkdir -p $WD1/03-Blastn
mkdir -p $WD1/04-Reciprocal_hits
mkdir -p $WD1/05-Families



####### PIPELINE

### SELECTION
## After generating temporal gene IDs, create a TXT file with the gene IDs and a FASTA file with the sequences of the genes.
cd $WD1/01-Genes

echo -e "\n\nSelect genes...\n"

for spe in $Species_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Select $spe $WD $WD1 $AS &
done
wait

### MAKEBLASTDB
## Create a blast database from the FASTA file created in the previous step using the function maskeblastdb.
cd $WD1/02-Makeblastdb

echo -e "\n\nMakeblastdb...\n"

for spe in $Species_list; do	
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Makeblastdb $spe &
done
wait

### BLASTN
## Run blastn for all possible pairwise combinations of the list of species except for combinations with itself.
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
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Blastn $Combinations $evalue $SLURM_CPUS_PER_TASK $j $WD1 $AS &
done
wait

### RECIPROCAL HITS
## Find the reciprocal blast hits between species, i.e. those hits that can be found in both the spe1-spe2 and spe2-spe1 combination. To do this, we use the home-made script Find_reciprocal_hits.py.
cd $WD1/04-Reciprocal_hits

echo -e "\n\nReciprocal hits...\n"

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Reciprocal_hits $Combinations $WD1 $AS $j &
done
wait

## Pool the results of all paired species combinations.
if [ -f Reciprocal_hits.tsv ] ; then
	rm Reciprocal_hits.tsv
fi
cat * > Reciprocal_hits.tsv

### FAMILIES
## Classify the conserved genes into conservation families using home-made script Classify_into_families.py.
cd $WD1/05-Families

echo -e "\n\nClassify lncRNAs into families...\n"

>ids_by_specie.tsv
for spe in $Species_list; do
	awk -v a="$spe" '{print a"\t"$0}' ../01-Genes/$spe\_ids\_new.txt >> ids_by_specie.tsv
done

mkdir -p outputs

# Classify into families.
>outputs/stdout.log
$AS/Classify_into_families.py \
	--pred-lncRNAs $WD1/05-Families/ids_by_specie.tsv \
	--reciprocal-hits $WD1/04-Reciprocal_hits/Reciprocal_hits.tsv \
	--fam $WD1/05-Families/fam.tsv \
	--gen $WD1/05-Families/gen.tsv \
	--threads $SLURM_CPUS_PER_TASK >> outputs/stdout.log 2>&1

### FIGURES
echo -e "\n\nFigures...\n"
Rscript $Scripts/Figures_genes.R


