#!/bin/bash

#SBATCH --job-name=OMgen					# Job name.
#SBATCH --output=OrthoMCL_genes.log				# Standard output and error log.
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
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Sequence_level/OrthoMCL/Prueba_genes/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
evalue=1e-3
Scripts=$(pwd)

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Sequence_level/OrthoMCL/nr/Prueba_genes
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
mkdir -p $WD/08-comparative_genomics/Sequence_level/OrthoMCL
mkdir -p $WD/08-comparative_genomics/Sequence_level/OrthoMCL/nr
mkdir -p $WD/08-comparative_genomics/Sequence_level/OrthoMCL/nr/Prueba_genes
mkdir -p $WD1/01-Genes
mkdir -p $WD1/02-Makeblastdb
mkdir -p $WD1/03-Blastn
mkdir -p $WD1/04-OrthoMCL
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
## Run blastn for all possible pairwise combinations of the list of species including combinations with itself. Filenames can be TXT or TSV files.
cd $WD1/03-Blastn

Combinations_species=$WD1/03-Blastn/Species_combination.txt
set -- $Species_list
for a; do
    for b; do
        printf "%s %s\n" "$a" "$b"
    done
done > $Combinations_species

echo -e "\n\nBlastn...\n"

n_comb=$(wc -l $Combinations_species | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Blastn $Combinations_species $evalue $SLURM_CPUS_PER_TASK $j $WD1 $AS &
done
wait

### ORTHOMCL
## Run OrthoMCL to create families of conserved genes using MCL clustering.
cd $WD1/04-OrthoMCL

echo -e "\n\nOrthoMCL...\n"

>Repo_spec.txt
>Repo_spec.txt.all.GFF3
>stdout.log
for spe1 in $Species_list; do
	mkdir -p $spe1
	
	# Directories
	DIR1=$WD1/01-Genes
	DIR2=$WD1/03-Blastn
	DIR3=$WD1/04-OrthoMCL/$spe1
	
	# Store in this folder all blastn tables in which the specie has been used as spe1 in 3-blastn and rename the files.
	if [ -d "$DIR3/RBH_blast_CDS" ]; then
		rm -r $DIR3/RBH_blast_CDS
	fi
	mkdir $DIR3/RBH_blast_CDS
	cp $DIR2/$spe1\_to* $DIR3/RBH_blast_CDS/
	cd $spe1/RBH_blast_CDS
	rm *_temp1.tsv
	
	# Rename files to OrthoMCL format name.
	files="$(ls $DIR3/RBH_blast_CDS/)"
	for file in $files; do
		spe2="$(echo $file | cut -d_ -f3)"
		mv $spe1\_to\_$spe2\_temp2.tsv vs_$spe2.blast.m8
	done
	cd ../../
	
	# Create Repo_spec.txt.all.GFF3.
	cp $DIR1/$spe1\_new.gtf $DIR3/$spe1.annotation.gtf
	gffread $DIR3/$spe1.annotation.gtf -o $DIR3/$spe1.annotation.gff3
	tail -n +4 $DIR3/$spe1.annotation.gff3 | cut -d";" -f1 | \
	awk '$3 == "transcript" {print $0}' | sed -e 's/ID=//g' | \
	sed -e 's/transcript/mRNA/g' | \
	awk -v var="$spe1" '{print $1"\t"var"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' >> Repo_spec.txt.all.GFF3
	rm $DIR3/$spe1.annotation.gff3
	rm $DIR3/$spe1.annotation.gtf
	
	# Create the Repo_spec.txt.
	echo -e "//" >> Repo_spec.txt
	echo -e "Genome $spe1" >> Repo_spec.txt
	echo -e "Annotation $spe1" >> Repo_spec.txt
done
echo -e "//" >> Repo_spec.txt

# Execute OrthoMCL using Synima to cluster the genes into families using the blastn reciprocal results.
if [ -d "OMCL_outdir" ]; then
	rm -r OMCL_outdir
fi
perl $SP/Synima-master/util/Blast_all_vs_all_repo_to_OrthoMCL.pl -r Repo_spec.txt -t CDS -o OMCL_outdir >> stdout.log 2>&1

### FAMILIES
## Get gene families from OrthoMCL.
cd $WD1/05-Families

echo -e "\n\nGet gene families from OrthoMCL...\n"

>ids_by_specie.tsv
for spe in $Species_list; do
	awk -v a="$spe" '{print a"\t"$0}' ../01-Genes/$spe\_ids\_new.txt >> ids_by_specie.tsv
done

mkdir -p outputs

# Get families.
>outputs/stdout.log
$AS/Get_families_from_OrthoMCL.py \
	--pred-lncRNAs $WD1/05-Families/ids_by_specie.tsv \
	--orthomcl $WD1/04-OrthoMCL/all_orthomcl.out \
	--fam $WD1/05-Families/fam.tsv \
	--gen $WD1/05-Families/gen.tsv >> outputs/stdout.log 2>&1

### FIGURES
echo -e "\n\nFigures...\n"
Rscript $Scripts/Figures_genes.R


