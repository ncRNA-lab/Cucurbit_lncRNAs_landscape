#!/bin/bash

#SBATCH --job-name=OFNRP1					# Job name.
#SBATCH --output=OrthoFinder_NR_P1.log			# Standard output and error log.
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
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Sequence_level/OrthoFinder/Prueba_e-value_5/Functions_NR.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"
Classes_list="intergenic antisense intronic sense ALL"
Confidence_levels_list="High Medium Low"
evalue=1e-5
Scripts=$(pwd)

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Sequence_level/OrthoFinder/nr/Prueba_e-value_5
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
mkdir -p $WD/08-comparative_genomics/Sequence_level/OrthoFinder/nr/Prueba_e-value_5
mkdir -p $WD1/01-LncRNAs
mkdir -p $WD1/02-Makeblastdb
mkdir -p $WD1/03-Blastn
mkdir -p $WD1/04-OrthoFinder
mkdir -p $WD1/05-Families



####### PIPELINE

### SELECTION
## For each confidence level, class code and specie, create a TXT file with the identifiers of each lncRNA and a FASTA file with the sequences of the lncRNAs. We have to modify lncRNA IDs by adding the specie. This is important so that after running OrthoFinder we can distinguish which lncRNAs come from which species since all transcripts assembled with stringtie are called MSTRG regardless of the species.
cd $WD1/01-LncRNAs

echo -e "\n\nSelect lncRNAs and modify lncRNA IDs by adding the specie...\n"

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
## For each confidence level and class code, run blastn for all possible pairwise combinations of the list of species including combinations with itself. In addition, OrthoFinder works with filenames containing numbers instead of specie's name. Therefore, we create a code number for each specie. Filenames have to be TXT because if you use TSV files, you will get an error.
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
	comb_spe=$(cat $Combinations_species | head -n $j | tail -n 1)
	spe1=$(echo $comb_spe | cut -d" " -f1)
	spe2=$(echo $comb_spe | cut -d" " -f2)
	comb_cod=$(cat $Combinations_codes | head -n $j | tail -n 1)
	cod1=$(echo $comb_cod | cut -d" " -f1)
	cod2=$(echo $comb_cod | cut -d" " -f2)
	for confidence in $Confidence_levels_list; do
		mkdir -p $confidence
		for class in $Classes_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Blastn $spe1 $spe2 $cod1 $cod2 $confidence $class $evalue $SLURM_CPUS_PER_TASK $j $n_comb $WD1 $AS &
		done
	done
done
wait

### ORTHOFINDER
## For each confidence level and class code, run OrthoFinder to create families of conserved lncRNAs using MCL clustering.
cd $WD1/04-OrthoFinder

echo -e "\n\nOrthoFinder...\n"

for confidence in $Confidence_levels_list; do
	mkdir -p $confidence
	for class in $Classes_list; do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_OrthoFinder $confidence $class $WD1 $AS $Combinations_codes $SLURM_CPUS_PER_TASK &
	done
done
wait

### FAMILIES
## For each confidence level and class code, get lncRNA families from OrthoFinder results.
cd $WD1/05-Families

echo -e "\n\nGet lncRNA families from OrthoFinder...\n"

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
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Get_families $confidence $class $WD1 $AS &
	done
done
wait

### FIGURES
echo -e "\n\nFigures...\n"
Rscript $Scripts/Figures_NR.R


