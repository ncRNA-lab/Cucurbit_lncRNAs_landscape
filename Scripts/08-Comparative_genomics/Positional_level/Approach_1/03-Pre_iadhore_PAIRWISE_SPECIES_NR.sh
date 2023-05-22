#!/bin/bash

#SBATCH --job-name=PreAdNR								# Job name.
#SBATCH --output=Prerequisities_adhore_PAIRWISE_SPECIES_NR.log			# Standard output and error log.
#SBATCH --partition=short								# Partition (queue)
#SBATCH --ntasks=10									# Run on one mode. 
#SBATCH --cpus-per-task=2								# Number of tasks = cpus.
#SBATCH --time=1-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=4gb								# Job memory request.

####### MODULES
module load R/4.1.2

####### VARIABLES
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Synteny/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Synteny/nr

####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}

####### DIRECTORY
mkdir -p $WD/08-comparative_genomics
mkdir -p $WD/08-comparative_genomics/Synteny
mkdir -p $WD/08-comparative_genomics/Synteny/nr
mkdir -p $WD1/03-Adhore
mkdir -p $WD1/03-Adhore/PAIRWISE_SPECIES


####### PIPELINE: PAIRWISE_SPECIES
echo -e "\n\nPAIRWISE_SPECIES..."
# We use all the 1:1 orthologs that we can predict between two species with Orthofinder and we validate it with Inparanoid (Other software to predict orthologs). Then we extract gene, strand and chromosomes info from the GFF3 annotation file, and we create a list by chromosome (lists/chr1.lst) with all gene names that we can find in the chromosome and its strand. Finally, we build the file that we will use to execute i-adhore. This file will contain all the parameter, genome names... (More info in adhore manual)


## STEP 1: Pairwise combination species
echo -e "\n\nSTEP 1: Create combinations pairwise table..."
Combinations=$WD1/03-Adhore/PAIRWISE_SPECIES/Species_combination.txt
set -- $Species_list
for a; do
    shift
    for b; do
        printf "%s %s\n" "$a" "$b"
    done
done > $Combinations

## STEP 2: Select 1:1 orthologs
echo -e "\n\nSTEP 2: Select 1:1 orthologs...\n"
cd $WD1/03-Adhore/PAIRWISE_SPECIES

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb | cut -d" " -f1)
	spe2=$(echo $comb | cut -d" " -f2)
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_extract_orthologs_pairwise $spe1 $spe2 $WD1 $AS &
done
wait

## STEP 3: Extract info from GFF3 genes
echo -e "\n\nSTEP 3: Extract info from GFF3 genes...\n"
cd $WD1/03-Adhore/PAIRWISE_SPECIES

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb | cut -d" " -f1)
	spe2=$(echo $comb | cut -d" " -f2)
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_extract_info_GFF3_genes $spe1 $spe2 $AI $WD1/03-Adhore/PAIRWISE_SPECIES $AS &
done
wait

## STEP 4: Build adhore_parameters.ini
echo -e "\n\nSTEP 4: Build adhore_parameters.ini...\n"
cd $WD1/03-Adhore/PAIRWISE_SPECIES

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb | cut -d" " -f1)
	spe2=$(echo $comb | cut -d" " -f2)
	
	#----- Collinear mode
	echo -e "\tCombination: "$spe1"-"$spe2"..."
	cd $spe1\-$spe2
	>idhore_parameters_collinear.ini
	echo -e "genome="$spe1 >> idhore_parameters_collinear.ini
	ls -1 Lists_$spe1/ | awk -v a="$spe1" '{split($1, id, "."); print id[1]" ./Lists_"a"/"$0}' >> idhore_parameters_collinear.ini
	echo -e "\n" >> idhore_parameters_collinear.ini
	
	echo -e "genome="$spe2 >> idhore_parameters_collinear.ini
	ls -1 Lists_$spe2/ | awk -v a="$spe2" '{split($1, id, "."); print id[1]" ./Lists_"a"/"$0}' >> idhore_parameters_collinear.ini
	echo -e "\n" >> idhore_parameters_collinear.ini
	
	echo -e "blast_table= Orthogroups.tsv" >> idhore_parameters_collinear.ini
	echo -e "table_type= family" >> idhore_parameters_collinear.ini
	echo -e "output_path= output_collinear" >> idhore_parameters_collinear.ini
	echo -e "\n" >> idhore_parameters_collinear.ini
	
	echo -e "cluster_type=colinear" >> idhore_parameters_collinear.ini
	
	echo -e "alignment_method=gg2" >> idhore_parameters_collinear.ini
	echo -e "gap_size=15" >> idhore_parameters_collinear.ini
	echo -e "tandem_gap=10" >> idhore_parameters_collinear.ini
	echo -e "cluster_gap=20" >> idhore_parameters_collinear.ini
	echo -e "max_gaps_in_alignment=20" >> idhore_parameters_collinear.ini
	echo -e "q_value=0.9" >> idhore_parameters_collinear.ini
	echo -e "prob_cutoff=0.001" >> idhore_parameters_collinear.ini
	echo -e "anchor_points=5" >> idhore_parameters_collinear.ini
	echo -e "multiple_hypothesis_correction=FDR" >> idhore_parameters_collinear.ini
	echo -e "level_2_only=true" >> idhore_parameters_collinear.ini
	echo -e "number_of_threads=16" >> idhore_parameters_collinear.ini
	echo -e "visualizeAlignment=false" >> idhore_parameters_collinear.ini
	echo -e "write_stats=true" >> idhore_parameters_collinear.ini
	
	
	#----- Cloud mode
	>idhore_parameters_cloud.ini
	echo -e "genome="$spe1 >> idhore_parameters_cloud.ini
	ls -1 Lists_$spe1/ | awk -v a="$spe1" '{split($1, id, "."); print id[1]" ./Lists_"a"/"$0}' >> idhore_parameters_cloud.ini
	echo -e "\n" >> idhore_parameters_cloud.ini
	
	echo -e "genome="$spe2 >> idhore_parameters_cloud.ini
	ls -1 Lists_$spe2/ | awk -v a="$spe2" '{split($1, id, "."); print id[1]" ./Lists_"a"/"$0}' >> idhore_parameters_cloud.ini
	echo -e "\n" >> idhore_parameters_cloud.ini
	
	echo -e "blast_table= Orthogroups.tsv" >> idhore_parameters_cloud.ini
	echo -e "table_type= family" >> idhore_parameters_cloud.ini
	echo -e "output_path= output_cloud" >> idhore_parameters_cloud.ini
	echo -e "\n" >> idhore_parameters_cloud.ini
	
	echo -e "cluster_type=cloud" >> idhore_parameters_cloud.ini
	
	echo -e "alignment_method=gg2" >> idhore_parameters_cloud.ini
	echo -e "cloud_gap_size=20" >> idhore_parameters_cloud.ini
	echo -e "cloud_cluster_gap=25" >> idhore_parameters_cloud.ini
	echo -e "cloud_filter_method=binomial" >> idhore_parameters_cloud.ini
	echo -e "q_value=0.9" >> idhore_parameters_cloud.ini
	echo -e "prob_cutoff=0.001" >> idhore_parameters_cloud.ini
	echo -e "anchor_points=5" >> idhore_parameters_cloud.ini
	echo -e "multiple_hypothesis_correction=FDR" >> idhore_parameters_cloud.ini
	echo -e "level_2_only=true" >> idhore_parameters_cloud.ini
	echo -e "number_of_threads=16" >> idhore_parameters_cloud.ini
	echo -e "visualizeAlignment=false" >> idhore_parameters_cloud.ini
	echo -e "write_stats=true" >> idhore_parameters_cloud.ini
	
	#----- Hybrid mode
	>idhore_parameters_hybrid.ini
	echo -e "genome="$spe1 >> idhore_parameters_hybrid.ini
	ls -1 Lists_$spe1/ | awk -v a="$spe1" '{split($1, id, "."); print id[1]" ./Lists_"a"/"$0}' >> idhore_parameters_hybrid.ini
	echo -e "\n" >> idhore_parameters_hybrid.ini
	
	echo -e "genome="$spe2 >> idhore_parameters_hybrid.ini
	ls -1 Lists_$spe2/ | awk -v a="$spe2" '{split($1, id, "."); print id[1]" ./Lists_"a"/"$0}' >> idhore_parameters_hybrid.ini
	echo -e "\n" >> idhore_parameters_hybrid.ini
	
	echo -e "blast_table= Orthogroups.tsv" >> idhore_parameters_hybrid.ini
	echo -e "table_type= family" >> idhore_parameters_hybrid.ini
	echo -e "output_path= output_hybrid" >> idhore_parameters_hybrid.ini
	echo -e "\n" >> idhore_parameters_hybrid.ini
	
	echo -e "cluster_type=hybrid" >> idhore_parameters_hybrid.ini
	
	echo -e "alignment_method=gg2" >> idhore_parameters_hybrid.ini
	echo -e "gap_size=15" >> idhore_parameters_hybrid.ini
	echo -e "tandem_gap=10" >> idhore_parameters_hybrid.ini
	echo -e "cluster_gap=20" >> idhore_parameters_hybrid.ini
	echo -e "max_gaps_in_alignment=20" >> idhore_parameters_hybrid.ini
	echo -e "q_value=0.9" >> idhore_parameters_hybrid.ini
	echo -e "prob_cutoff=0.001" >> idhore_parameters_hybrid.ini
	echo -e "cloud_gap_size=20" >> idhore_parameters_hybrid.ini
	echo -e "cloud_cluster_gap=25" >> idhore_parameters_hybrid.ini
	echo -e "cloud_filter_method=binomial" >> idhore_parameters_hybrid.ini
	echo -e "anchor_points=5" >> idhore_parameters_hybrid.ini
	echo -e "multiple_hypothesis_correction=FDR" >> idhore_parameters_hybrid.ini
	echo -e "level_2_only=true" >> idhore_parameters_hybrid.ini
	echo -e "number_of_threads=16" >> idhore_parameters_hybrid.ini
	echo -e "visualizeAlignment=false" >> idhore_parameters_hybrid.ini
	echo -e "write_stats=true" >> idhore_parameters_hybrid.ini
	
	cd ..
done


