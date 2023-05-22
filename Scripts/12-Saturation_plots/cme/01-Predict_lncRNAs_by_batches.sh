#!/bin/bash

#SBATCH --job-name=cmeCS5				# Job name.
#SBATCH --output=cmeCS5.log				# Standard output and error log.
#SBATCH --qos=long-mem					# Partition (queue)
#SBATCH --ntasks=40					# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=2				# Number of tasks = cpus.
#SBATCH --time=7-00:00:00				# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=8gb				# Job memory request.


# exit when any command fails
set -e

####### MODULES
module load anaconda/anaconda3
module load R/4.1.2

####### VARIABLES
specie="cme"
evalue=1e-5
range=5
permutations=10
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/13-saturation_plots/additional_scripts"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/05-predict_lncRNAs/Softwares_prediction"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/13-saturation_plots/additional_scripts/Functions.sh"


####### ADDITIONAL SCRIPTS
export ASPATH=$AS
export PATH=$PATH:${ASPATH}


####### SOFTWARES PREDICTION
### CPC2
#https://github.com/gao-lab/CPC2_standalone
#pip3 install biopython
# If other error exists is convenient to remove the package and install it again. This works for me.
export CPC2PATH=$SP/CPC2_standalone-1.0.1
export PATH=$PATH:${CPC2PATH}/bin

### FEELnc
#https://github.com/tderrien/FEELnc
#sudo cpanm Bio::DB::GenPept --force
#sudo cpanm Bio::DB::GenBank --force
#sudo cpanm Bio::Perl
#sudo cpan Parallel::ForkManager
# I had to install many perl libraries manually. All of them were saved in lib.
# Moreover, I had to change png by svg in the script utils/codpot_randomforest.r manually because an error was launched.
export FEELNCPATH=$SP/FEELnc
export PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib #order is important to avoid &Bio::DB::IndexedBase::_strip_crnl error with bioperl >=v1.7
export PATH=$PATH:${FEELNCPATH}/scripts
export PATH=$PATH:${FEELNCPATH}/utils
export PATH=$PATH:${FEELNCPATH}/bin/LINUX

### CPAT
#http://rna-cpat.sourceforge.net/
#https://github.com/gbgolding/crema
#pip3 install CPAT

### DIAMOND
#https://github.com/bbuchfink/diamond
export DIAMONDPATH=$SP/diamond-linux64
export PATH=$PATH:${DIAMONDPATH}

### TRANSDECODER
#https://github.com/TransDecoder/TransDecoder/wiki
# It requires URI/Escape.pm perl library.
export TRANSDECODERPATH=$SP/TransDecoder-TransDecoder-v5.5.0
export PATH=$PATH:${TRANSDECODERPATH}
export PATH=$PATH:${TRANSDECODERPATH}/util

### HMMER
#https://github.com/EddyRivasLab/hmmer
export HMMERPATH=$SP/hmmer-3.3.2
export PATH=$PATH:${HMMERPATH}
export PATH=$PATH:${HMMERPATH}/bin

### BLAST
#https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#http://nebc.nerc.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf
export NCBIBLASTPATH=$SP/ncbi-blast-2.13.0+
export PATH=$PATH:${NCBIBLASTPATH}/bin

### VIENNARNA (RNAFOLD)
#https://github.com/ViennaRNA/ViennaRNA
export VIENNAPATH=$SP/ViennaRNA-2.5.1/ViennaRNA
export PATH=$PATH:${VIENNAPATH}/bin
export PATH=$PATH:${VIENNAPATH}/lib

### MIRENA
#http://www.lcqb.upmc.fr/mirena/index.html
export MIRENAPATH=$SP/MIReNA-2.0
export PATH=$PATH:${MIRENAPATH}

### CD-HIT
#https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide
export CDHITPATH=$SP/cdhit-4.8.1
export PATH=$PATH:${CDHITPATH}

### AGAT
#https://github.com/NBISweden/AGAT
export AGATPATH=$SP/AGAT-0.9.2
export PATH=$PATH:${AGATPATH}
export PATH=$PATH:${AGATPATH}/bin


####### DIRECTORY
mkdir -p $WD/13-saturation_plots
mkdir -p $WD/13-saturation_plots/$specie
mkdir -p $WD/13-saturation_plots/$specie/Predictions


####### PIPELINE BY PERMUTATION
for permutation in $(seq -s " " $permutations); do
	
	echo -e "\n\n\n\n################################## PERMUTATION $permutation\n\n"
	
	####### NEW AND OTHER VARIABLES
	WD1=$WD/03-assembly/$specie
	WD2=$WD/13-saturation_plots/$specie/Predictions/range_$range-$permutation
	WD3=$WD/13-saturation_plots/$specie/Tables_CS
	
	####### DIRECTORY
	mkdir -p $WD2
	mkdir -p $WD3

	###########################################################################
	##
	## Create random sample lists to predict lncRNAs
	##
	########################################################################### 	

	echo -e "\n\n#############################"
	echo -e "########### STEP 1 ##########"
	echo -e "#############################\n"

	####### DIRECTORY
	mkdir -p $WD2/01-Batches

	####### PIPELINE

	### Create random sample lists to predict lncRNAs
	echo -e "Creating random sample lists (Batches)...\n"
	Create_random_sample_lists.py \
		--input-list $AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_2.txt" \
		--output-path $WD2/01-Batches \
		--step-size $range




	###########################################################################
	##
	## Merge GTF files by batch
	##
	########################################################################### 	

	echo -e "\n\n#############################"
	echo -e "########### STEP 2 ##########"
	echo -e "#############################\n"

	####### DIRECTORY
	mkdir -p $WD2/02-Merged_assemblies
	mkdir -p $WD2/02-Merged_assemblies/LOGS

	####### PIPELINE

	### Merge GTF files by batch
	echo -e "Merge GTF files by batches...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/02-Merged_assemblies/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Merge_assemblies $i $WD1 $WD2 $AI $specie $SLURM_CPUS_PER_TASK &
	done
	wait



	###########################################################################
	##
	## Annotate the transcripts using gffcompare
	##
	########################################################################### 	

	echo -e "\n\n#############################"
	echo -e "########### STEP 3 ##########"
	echo -e "#############################\n"

	####### DIRECTORY
	mkdir -p $WD2/03-Annotate_transcripts
	mkdir -p $WD2/03-Annotate_transcripts/LOGS

	####### PIPELINE

	## Compare the merged GTF file with the original genome annotation using gffcompare to identify novel transcripts.
	echo -e "Annotate the transcripts using gffcompare...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/03-Annotate_transcripts/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c1 -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Annotate_transcripts $i $WD2 $AI $specie &
	done
	wait




	###########################################################################
	##
	## Predict lncRNAs using own pipeline
	##
	########################################################################### 	

	echo -e "\n\n#############################"
	echo -e "########### STEP 4 ##########"
	echo -e "#############################\n"

	####### DIRECTORY
	mkdir -p $WD2/04-Predict_lncRNAs
	mkdir -p $WD2/04-Predict_lncRNAs/STEP1
	mkdir -p $WD2/04-Predict_lncRNAs/STEP1/LOGS
	mkdir -p $WD2/04-Predict_lncRNAs/STEP2
	mkdir -p $WD2/04-Predict_lncRNAs/STEP2/LOGS
	mkdir -p $WD2/04-Predict_lncRNAs/STEP3
	mkdir -p $WD2/04-Predict_lncRNAs/STEP3/LOGS
	mkdir -p $WD2/04-Predict_lncRNAs/STEP4
	mkdir -p $WD2/04-Predict_lncRNAs/STEP4/LOGS
	mkdir -p $WD2/04-Predict_lncRNAs/STEP5
	mkdir -p $WD2/04-Predict_lncRNAs/STEP5/LOGS
	mkdir -p $WD2/04-Predict_lncRNAs/STEP6
	mkdir -p $WD2/04-Predict_lncRNAs/STEP6/LOGS
	mkdir -p $WD2/04-Predict_lncRNAs/STEP-FINAL
	mkdir -p $WD2/04-Predict_lncRNAs/STEP-FINAL/LOGS

	####### PIPELINE

	### STEP1: Create two folders:
	### 	- Original_genes: Save important files about genes such as GTF 
	### 	(Annotation), FASTA (Nucleotide sequences), TSV (Additional info) 
	### 	and TXT (Sequence IDs).
	### 	- Potential_lncRNAs: Save important files about potential lncRNAs 
	### 	such as GTF (Annotation), FASTA (Nucleotide sequences), TSV 
	### 	(Additional info) and TXT (Sequence IDs).

	echo -e "STEP (1/7): PIPELINE...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/04-Predict_lncRNAs/STEP1/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c1 -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Predict_lncRNAs_STEP1 $i $WD2 $AI $specie &
	done
	wait

	### STEP2: Predict the coding potential of the potential lncRNAs using three 
	### tools:
	### 	- CPC2 (https://github.com/gao-lab/CPC2_standalone; v.1.0.1)
	### 	- FEELnc (https://github.com/tderrien/FEELnc; v.0.2)
	### 	- CPAT (https://github.com/liguowang/cpat; v.3.0.2)

	echo -e "STEP (2/7): PIPELINE...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/04-Predict_lncRNAs/STEP2/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Predict_lncRNAs_STEP2 $i $WD2 $SP $SLURM_CPUS_PER_TASK &
	done
	wait

	### STEP3: Predict the coding potential of the potential lncRNAs using two 
	### protein databases:
	### 	- SwissProt (https://www.uniprot.org/uniprot/?query=reviewed:yes; 
	### 	18/02/2022)
	### 	- Pfam-A (http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/;
	### 	31/05/2022)

	echo -e "STEP (3/7): PIPELINE...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/04-Predict_lncRNAs/STEP3/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Predict_lncRNAs_STEP3 $i $WD2 $SP $evalue $SLURM_CPUS_PER_TASK &
	done
	wait

	### STEP4: Predict if it exits any ORF with length greater than 80, 100 or 
	### 120 aa

	echo -e "STEP (4/7): PIPELINE...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/04-Predict_lncRNAs/STEP4/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c1 -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Predict_lncRNAs_STEP4 $i $WD2 &
	done
	wait

	### STEP5: Annotate the potential lncRNAs using RNAcentral (rRNA, tRNA, snRNA, 
	### snoRNA) and miRBase and PmiREN (miRNAs and miRNA-Precursors).

	echo -e "STEP (5/7): PIPELINE...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/04-Predict_lncRNAs/STEP5/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Predict_lncRNAs_STEP5 $i $WD2 $AI $evalue $SLURM_CPUS_PER_TASK &
	done
	wait

	### STEP6: Annotate the potential lncRNAs using three lncRNAs databases:
	### 	- CANTATAdb (http://cantata.amu.edu.pl/)
	### 	- PLncDB (plncdb.tobaccodb.org)
	### 	- GreeNC (http://greenc.sequentiabiotech.com/wiki2/Main_Page)

	echo -e "STEP (6/7): PIPELINE...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/04-Predict_lncRNAs/STEP6/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Predict_lncRNAs_STEP6 $i $WD2 $AI $evalue $SLURM_CPUS_PER_TASK &
	done
	wait

	### STEP-FINAL: Database creation: LncRNA classification and redundancy remove.

	echo -e "STEP-FINAL (7/7): PIPELINE...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	LOG=$WD2/04-Predict_lncRNAs/STEP-FINAL/LOGS
	for i in $(seq -s " " $n_batches); do
		srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $LOG/Batch_stdout-$i.log --quiet --exclusive $F task_Predict_lncRNAs_STEP-FINAL $i $WD2 $AI $AS $SP $specie $SLURM_CPUS_PER_TASK &
	done
	wait




	###########################################################################
	##
	## Create the batch tables from which we will plot the saturation curves
	##
	########################################################################### 	

	echo -e "\n\n#############################"
	echo -e "########### STEP 5 ##########"
	echo -e "#############################\n"

	####### PIPELINE

	## Create Batch tables
	echo -e "Create Batch tables...\n"
	n_batches=$(ls -1 $WD2/01-Batches/*.txt | wc -l)
	n_samples_total=$(wc -l $AI/sra-info/accession_list/$specie-SRR_Acc_List-Filter_2.txt | cut -d " " -f 1)
	Rscript $AS/Batch_tables.R $WD2 $WD3 $n_batches $range $permutation $n_samples_total
done

## Final Batch table (ALL PERMUTATIONS)
echo -e "\n\n\nCreate Final Batch tables...\n"
Rscript $AS/Batch_tables_mean_permutations.R $WD3 $range $permutations


echo -e "\n\n\n\n#############################"
echo -e "############ END ############"
echo -e "#############################\n\n"



