#!/bin/bash

#SBATCH --job-name=cmeS5							# Job name.
#SBATCH --output=cme_STEP5.log						# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=33								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb							# Job memory request.


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/03-Individual_assembly/sra_data/Functions.sh"

####### NEW VARIABLES
WD1=$WD"/01-Samples_processing/sra_data"
WD2=$WD"/02-Mapping/sra_data"
WD3=$WD"/03-Individual_assembly/sra_data"
Acc_list=$AI"/sra-info/accession_list/"$specie"-SRR_Acc_List-Filter_2.txt"
strand_info=$WD1"/03-Strand_detection/Filter_table/Strand_info.tsv"

####### DIRECTORY
mkdir -p $WD
mkdir -p $WD/03-Individual_assembly
mkdir -p $WD3


####### PIPELINE: STEP 5

echo -e "\n\n\n####################"
echo -e "###### STEP 5 ######"
echo -e "####################\n"

### INDIVIDUAL ASSEMBLY
## Assemble each trimmed library using Stringtie software.
echo -e "\nINDIVIDUAL ASSEMBLY..."

cd $WD3
mkdir -p outputs
for SRR in $(cat $Acc_list); do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD3/outputs/stdout_Assembly_$SRR.log --quiet --exclusive $F task_Assembly $SRR $WD2 $AI $strand_info $specie $SLURM_CPUS_PER_TASK &
done
wait

