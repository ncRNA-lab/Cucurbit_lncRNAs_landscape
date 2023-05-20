#!/bin/bash

#SBATCH --job-name=cmeS8				# Job name.
#SBATCH --output=cme_STEP8.log			# Standard output and error log.
#SBATCH --qos=short					# Partition (queue)
#SBATCH --ntasks=5					# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=20				# Number of tasks = cpus.
#SBATCH --time=1-00:00:00				# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb				# Job memory request.


# exit when any command fails
set -e

####### MODULES
module load anaconda/anaconda3_2021.07
module load R/4.1.2
module load biotools

####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"
AS="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/06-LncRNAs_prediction/additional_scripts"
F="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/06-LncRNAs_prediction/Functions.sh"
SP="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Scripts/Pascual/05-predict_lncRNAs/Softwares_prediction"
evalue=1e-5

####### NEW AND OTHER VARIABLES
WD1=$WD"/05-Transcript_annotation"
WD2=$WD"/06-LncRNAs_prediction"


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
mkdir -p $WD
mkdir -p $WD2
mkdir -p $WD2/LOGS
mkdir -p $WD2/STEP1
mkdir -p $WD2/STEP2
mkdir -p $WD2/STEP3
mkdir -p $WD2/STEP4
mkdir -p $WD2/STEP5
mkdir -p $WD2/STEP6
mkdir -p $WD2/STEP-FINAL

####### PIPELINE: STEP 8

echo -e "\n\n\n####################"
echo -e "###### STEP 8 ######"
echo -e "####################\n"

### STEP1: Create two folders:
### 	- Original_genes: Save important files about genes such as GTF 
### 	(Annotation), FASTA (Nucleotide sequences), TSV (Additional info) 
### 	and TXT (Sequence IDs).
### 	- Potential_lncRNAs: Save important files about potential lncRNAs 
### 	such as GTF (Annotation), FASTA (Nucleotide sequences), TSV 
### 	(Additional info) and TXT (Sequence IDs).

echo -e "STEP (1/7): PIPELINE...\n"
g=50
srun -N1 -n1 -c1 -o $WD2/LOGS/stdout_STEP1.log --quiet --exclusive $F task_LncRNAs_prediction_STEP1 $WD1 $WD2 $AI $specie $g

### STEP2: Predict the coding potential of the potential lncRNAs using three 
### tools:
### 	- CPC2 (https://github.com/gao-lab/CPC2_standalone; v.1.0.1)
### 	- FEELnc (https://github.com/tderrien/FEELnc; v.0.2)
### 	- CPAT (https://github.com/liguowang/cpat; v.3.0.2)

echo -e "STEP (2/7): PIPELINE...\n"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $WD2/LOGS/stdout_STEP2.log --quiet --exclusive $F task_LncRNAs_prediction_STEP2 $WD2 $SP $SLURM_CPUS_PER_TASK &

### STEP3: Predict the coding potential of the potential lncRNAs using two 
### protein databases:
### 	- SwissProt (https://www.uniprot.org/uniprot/?query=reviewed:yes; 
### 	18/02/2022)
### 	- Pfam-A (http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/;
### 	31/05/2022)

echo -e "STEP (3/7): PIPELINE...\n"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $WD2/LOGS/stdout_STEP3.log --quiet --exclusive $F task_LncRNAs_prediction_STEP3 $WD2 $SP $evalue $SLURM_CPUS_PER_TASK &

### STEP4: Predict if it exits any ORF with length greater than 80, 100 or 
### 120 aa

echo -e "STEP (4/7): PIPELINE...\n"
srun -N1 -n1 -c1 -o $WD2/LOGS/stdout_STEP4.log --quiet --exclusive $F task_LncRNAs_prediction_STEP4 $WD2 &

### STEP5: Annotate the potential lncRNAs using RNAcentral (rRNA, tRNA, snRNA, 
### snoRNA) and miRBase and PmiREN (miRNAs and miRNA-Precursors).

echo -e "STEP (5/7): PIPELINE...\n"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $WD2/LOGS/stdout_STEP5.log --quiet --exclusive $F task_LncRNAs_prediction_STEP5 $WD2 $AI $evalue $SLURM_CPUS_PER_TASK &

### STEP6: Annotate the potential lncRNAs using three lncRNAs databases:
### 	- CANTATAdb (http://cantata.amu.edu.pl/)
### 	- PLncDB (plncdb.tobaccodb.org)
### 	- GreeNC (http://greenc.sequentiabiotech.com/wiki2/Main_Page)

echo -e "STEP (6/7): PIPELINE...\n"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $WD2/LOGS/stdout_STEP6.log --quiet --exclusive $F task_LncRNAs_prediction_STEP6 $WD2 $AI $evalue $SLURM_CPUS_PER_TASK &
wait

### STEP-FINAL: Database creation: LncRNA classification and redundancy remove.

echo -e "STEP-FINAL (7/7): PIPELINE...\n"
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $WD2/LOGS/stdout_STEP7.log --quiet --exclusive $F task_LncRNAs_prediction_STEP-FINAL $WD2 $AI $AS $SP $specie $SLURM_CPUS_PER_TASK
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK -o $WD2/LOGS/stdout_STEP7_circos.log --quiet --exclusive $F task_LncRNAs_prediction_STEP-FINAL_circos $WD2 $AI $AS $specie


