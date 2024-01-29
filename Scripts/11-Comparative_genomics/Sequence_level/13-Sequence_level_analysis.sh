#!/bin/bash

#SBATCH --job-name=S13seq					# Job name.
#SBATCH --output=STEP13_sequence.log				# Standard output and error log.
#SBATCH --qos=short						# Partition (queue)
#SBATCH --ntasks=5						# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=4					# Number of tasks = cpus.
#SBATCH --time=0-02:00:00					# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb					# Job memory request.


####### MODULES
module load anaconda
module load R/4.1.2

####### VARIABLES
WD1="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction"
WD2="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/11-Comparative_genomics/Sequence_level/Functions.sh"
specie_list="car cla cma cme cmo cpe csa lsi mch"
class_list="intergenic antisense intronic sense"
flag_list="nr"
confidence_list="High Medium Low"
evalue=1e-5

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### SOFTWARES PREDICTION
### BLAST
export NCBIBLASTPATH=$SP/ncbi-blast-2.13.0+
export PATH=$PATH:${NCBIBLASTPATH}/bin
### ORTHOFINDER
export OFPATH=$SP/OrthoFinder
export PATH=$PATH:${OFPATH}
export PATH=$PATH:${OFPATH}/bin
export PATH=$PATH:${OFPATH}/tools

####### DIRECTORY
mkdir -p $WD2
mkdir -p $WD2/Sequence_level


####### PIPELINE: STEP 13

### ANALYSIS OF CONSERVATION AT SEQUENCE LEVEL
echo -e "\nANALYSIS OF CONSERVATION AT SEQUENCE LEVEL..."

for flag in $flag_list; do

	echo -e "\n\nFLAG: "$flag
	
	## Variable.
	O=$WD2/Sequence_level/$flag
	
	## Clean
	cd $WD2/Sequence_level
	if [ -d "$flag" ]; then
		rm -r $flag
	fi
	mkdir $flag
	cd $flag
	
	## Directory.
	mkdir -p $O
	mkdir -p $O/01-LncRNAs
	mkdir -p $O/02-Makeblastdb
	mkdir -p $O/03-Blastn
	mkdir -p $O/04-OrthoFinder
	mkdir -p $O/05-Families
	mkdir -p $O/06-Figures_and_tables
	
	### SELECTION
	## For each confidence level, class code and specie, create a TXT file with the identifiers of each lncRNA and a FASTA file with the sequences of the lncRNAs. We have to modify lncRNA IDs by adding the specie. This is important so that after running OrthoFinder we can distinguish which lncRNAs come from which species since all transcripts assembled with stringtie are called MSTRG regardless of the species.
	cd $O/01-LncRNAs

	echo -e "\n\nSelect lncRNAs and modify lncRNA IDs by adding the specie...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			for spe in $specie_list; do
				srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Select $spe $confidence $class $WD1 $flag &
			done
		done
	done
	wait

	### MAKEBLASTDB
	## For each confidence level, class code and specie, create a blast database from the FASTA file created in the previous step using the function maskeblastdb.
	cd $O/02-Makeblastdb

	echo -e "\n\nMakeblastdb...\n"

	for spe in $specie_list; do
		for confidence in $confidence_list; do
			for class in $class_list; do
				srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Makeblastdb $spe $confidence $class &
			done
		done
	done
	wait

	### BLASTN
	## For each confidence level and class code, run blastn for all possible pairwise combinations of the list of species including combinations with itself. In addition, OrthoFinder works with filenames containing numbers instead of specie's name. Therefore, we create a code number for each specie. Filenames have to be TXT because if you use TSV files, you will get an error.
	cd $O/03-Blastn

	Combinations_species=$O/03-Blastn/Species_combination.txt
	set -- $specie_list
	for a; do
	    for b; do
		printf "%s %s\n" "$a" "$b"
	    done
	done > $Combinations_species

	Combinations_codes=$O/03-Blastn/Codes_combination.txt
	Codes_list=$(seq 0 $(($(wc -w <<< $specie_list)-1)))
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
		for confidence in $confidence_list; do
			for class in $class_list; do
				srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Blastn $spe1 $spe2 $cod1 $cod2 $confidence $class $evalue $SLURM_CPUS_PER_TASK $j $n_comb $O $AS &
			done
		done
	done
	wait

	### ORTHOFINDER
	## For each confidence level and class code, run OrthoFinder to create families of conserved lncRNAs using MCL clustering.
	cd $O/04-OrthoFinder

	echo -e "\n\nOrthoFinder...\n"

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_OrthoFinder $confidence $class $O $AS $Combinations_codes $SLURM_CPUS_PER_TASK &
		done
	done
	wait

	### FAMILIES
	## For each confidence level and class code, get lncRNA families from OrthoFinder results.
	cd $O/05-Families

	echo -e "\n\nGet lncRNA families from OrthoFinder...\n"

	for confidence in $confidence_list; do
		mkdir -p $confidence
		for class in $class_list; do
			mkdir -p $confidence/$class
			>$confidence/$class/ids_by_specie.tsv
			for spe in $specie_list; do
				awk -v a="$spe" '{print a"\t"$0}' ../01-LncRNAs/$confidence/$class/$spe\_ids.txt >> $confidence/$class/ids_by_specie.tsv
			done
		done
	done

	for confidence in $confidence_list; do
		for class in $class_list; do
			srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Get_families $confidence $class $O $AS &
		done
	done
	wait

	### FIGURES
	echo -e "\n\nCreate figures and tables..."
	Rscript $AS/Create_figures_and_tables_OrthoFinder.R $O/05-Families $O/06-Figures_and_tables "$specie_list" "$class_list" "$confidence_list"
	
done


