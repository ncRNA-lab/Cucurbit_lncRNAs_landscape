#!/bin/bash

#SBATCH --job-name=OFgen					# Job name.
#SBATCH --output=OrthoFinder_genes.log			# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=25						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=4					# Number of tasks = cpus.
#SBATCH --time=1-00:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=2gb					# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Sequence_level/OrthoFinder/Prueba_genes_e-value_5/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
evalue=1e-5
Scripts=$(pwd)

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Sequence_level/OrthoFinder/nr/Prueba_genes_e-value_5
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
### ORTHOFINDER
export OFPATH=$SP/OrthoFinder
export PATH=$PATH:${OFPATH}
export PATH=$PATH:${OFPATH}/bin
export PATH=$PATH:${OFPATH}/tools

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Sequence_level
mkdir -p $WD/08-comparative_genomics/Sequence_level/OrthoFinder
mkdir -p $WD/08-comparative_genomics/Sequence_level/OrthoFinder/nr
mkdir -p $WD/08-comparative_genomics/Sequence_level/OrthoFinder/nr/Prueba_genes_e-value_5
mkdir -p $WD1/01-Genes
mkdir -p $WD1/02-Makeblastdb
mkdir -p $WD1/03-Blastn
mkdir -p $WD1/04-OrthoFinder
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
## Run blastn for all possible pairwise combinations of the list of species including combinations with itself. In addition, OrthoFinder works with filenames containing numbers instead of specie's name. Therefore, we create a code number for each specie. Filenames have to be TXT because if you use TSV files, you will get an error.
cd $WD1/03-Blastn

Combinations_species=$WD1/03-Blastn/Species_combination.txt
set -- $Species_list
for a; do
    for b; do
        printf "%s %s\n" "$a" "$b"
    done
done > $Combinations_species

Combinations_codes=$WD1/03-Blastn/Codes_combination.txt
Codes_list=$(seq 0 $(($(wc -w <<< $Species_list)-1)))
set -- $Codes_list
for a; do
    for b; do
        printf "%s %s\n" "$a" "$b"
    done
done > $Combinations_codes

echo -e "\n\nBlastn...\n"

n_comb=$(wc -l $Combinations_species | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Blastn $Combinations_species $Combinations_codes $evalue $SLURM_CPUS_PER_TASK $j $WD1 $AS &
done
wait

### ORTHOFINDER
## Run OrthoFinder to create families of conserved genes using MCL clustering.
cd $WD1/04-OrthoFinder

echo -e "\n\nOrthoFinder...\n"

# Directories
DIR1=$WD1/01-Genes
DIR2=$WD1/03-Blastn
DIR3=$WD1/04-OrthoFinder

# First, run orthofinder with the '-op' option. This option will prepare the files in the format required by OrthoFinder and print the set of BLAST commands that need to be run. This is useful if you want to manage the BLAST searches yourself.
if [ -d "$DIR3/Prepare_blast" ]; then
	rm -r $DIR3/Prepare_blast
fi
>$DIR3/stdout_prepare_blast.log
orthofinder -f $DIR1 -d -o $DIR3/Prepare_blast -op -t $SLURM_CPUS_PER_TASK >> $DIR3/stdout_prepare_blast.log

# Second, copy the blastn results coming from the 03-Blastn step into the folder created (Prepare_blast/$folder_p/WorkingDirectory/).
folder_p=$(ls $DIR3/Prepare_blast) 
cp $DIR2/*_temp2.txt $DIR3/Prepare_blast/$folder_p/WorkingDirectory/

# Third, modify gene ids again. In this case, we have to convert gene IDs to OrthoFinder format. We use a table created in the first step of 04-OrthoFinder step.
$AS/Convert_lncRNA_ID_to_OrthoFinder_ID.py \
	--path $DIR3/Prepare_blast/$folder_p/WorkingDirectory \
	--comb $Combinations_codes
rm $DIR3/Prepare_blast/$folder_p/WorkingDirectory/*_temp2.txt

# Fourth, compress the blastn files.
n_comb=$(wc -l $Combinations_codes | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb_cod=$(cat $Combinations_codes | head -n $j | tail -n 1)
	cod1=$(echo $comb_cod | cut -d" " -f1)
	cod2=$(echo $comb_cod | cut -d" " -f2)
	gzip -c $DIR3/Prepare_blast/$folder_p/WorkingDirectory/Blast$cod1\_$cod2.txt > $DIR3/Prepare_blast/$folder_p/WorkingDirectory/Blast$cod1\_$cod2.txt.gz
	rm $DIR3/Prepare_blast/$folder_p/WorkingDirectory/Blast$cod1\_$cod2.txt
done

# Fifth, execute OrthoFinder with the '-og' opion to cluster the genes into families. The MCL clustering is used to cluster the blastn reciprocal results. OrthoFinder automatically converts the OrthoFinder gene IDs to the initial gene IDs.
if [ -d "$DIR3/Clustering" ]; then
	rm -r $DIR3/Clustering
fi
>$DIR3/stdout_clustering.log
orthofinder -b $DIR3/Prepare_blast/$folder_p/WorkingDirectory/ -d -og -t $SLURM_CPUS_PER_TASK >> $DIR3/stdout_clustering.log
cp -R $DIR3/Prepare_blast/$folder_p/WorkingDirectory/OrthoFinder $DIR3/
mv $DIR3/OrthoFinder $DIR3/Clustering

# Sixth, liberate space.
cd $DIR3/Prepare_blast/$folder_p/WorkingDirectory/
rm -r OrthoFinder
folder_c=$(ls $DIR3/Clustering) 
cd $DIR3/Clustering/$folder_c/
rm -r Orthogroup_Sequences
rm -r Single_Copy_Orthologue_Sequences

### FAMILIES
## Get gene families from OrthoFinder results.
cd $WD1/05-Families

echo -e "\n\nGet gene families from OrthoFinder...\n"

>ids_by_specie.tsv
for spe in $Species_list; do
	awk -v a="$spe" '{print a"\t"$0}' ../01-Genes/$spe\_ids\_new.txt >> ids_by_specie.tsv
done

mkdir -p outputs

folder_c=$(ls $WD1/04-OrthoFinder/Clustering)
# Get families.
>outputs/stdout.log
$AS/Get_families_from_OrthoFinder.py \
	--pred-lncRNAs $WD1/05-Families/ids_by_specie.tsv \
	--orthofinder $WD1/04-OrthoFinder/Clustering/$folder_c/Orthogroups/Orthogroups.txt \
	--fam $WD1/05-Families/fam.tsv \
	--gen $WD1/05-Families/gen.tsv >> outputs/stdout.log 2>&1

### FIGURES
echo -e "\n\nFigures...\n"
Rscript $Scripts/Figures_genes.R


