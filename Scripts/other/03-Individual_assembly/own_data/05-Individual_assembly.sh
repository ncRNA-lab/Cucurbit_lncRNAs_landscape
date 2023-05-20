#!/bin/bash

#SBATCH --job-name=Indass							# Job name.
#SBATCH --output=Individual_assembly.log					# Standard output and error log.
#SBATCH --qos=short								# Partition (queue)
#SBATCH --ntasks=33								# Run on one mode. 
#SBATCH --cpus-per-task=6							# Number of tasks = cpus.
#SBATCH --time=1-00:00:00							# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=1gb							# Job memory request.


####### VARIABLES
specie="cme"
WD="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Results"
AI="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Additional_info"
F="/storage/ncRNA/Projects/lncRNAs/Cme_stress/Scripts/03-Individual_assembly/own_data/Functions.sh"

####### NEW VARIABLES
WD1=$WD"/01-Samples_processing/own_data"
WD2=$WD"/02-Mapping/own_data"
WD3=$WD"/03-Individual_assembly/own_data"
Samples_list=$(cat $AI/Summary_own_samples/Samples_SS_filt_conversion.txt)
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
for name in $Samples_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD3/outputs/stdout_Assembly_$name.log --quiet --exclusive $F task_Assembly $name $WD2 $AI $strand_info $specie $SLURM_CPUS_PER_TASK &
done
wait

