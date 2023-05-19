#!/bin/bash

#SBATCH --job-name=cmeS2.2						# Job name.
#SBATCH --output=cme_STEP2.2.log					# Standard output and error log.
#SBATCH --qos=short							# Partition (queue)
#SBATCH --ntasks=5							# Run on one mode. Don't change unless you know what you are doing.
#SBATCH --cpus-per-task=4						# Number of tasks = cpus. It depends on the number of process of your parallelization.
#SBATCH --time=0-10:00:00						# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=25gb						# Job memory request.


# Para buscar aquellas muestras con errores como un numero diferente de reads hay que ejecutar: 
# grep "WARNNIG" /storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Results/Update_24_01_2022/01-Samples_processing/cme/02-Trimmed_data/outputs/stdout_Trimming_*
#
# Fastp se queda bloqueado cuando ocurre el error:
#
# WARNNIG: different read numbers of the 18903 pack
# Read1 pack size: 407
# Read2 pack size: 1000
#
# Hay que pararlo, buscar los identificadores de las muestras y reparar las muestras con este script.


####### VARIABLES
specie="cme"
update="Update_2023_03_29"
WD="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Results"
F="/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae_updates/Scripts/"$update"/01-Samples_processing/Functions-REPAIR.sh"
trim_adapters="/storage/ncRNA/Softwares/BBMap_38.90/bbmap/resources"
Samples_to_repair="SRR1342457 SRR6986549 SRR6986551 SRR19977856 SRR19977866"

####### NEW VARIABLES
WD1=$WD"/"$update"/01-Samples_processing"
WD1_spe=$WD1"/"$specie


####### PIPELINE: STEP 2 (OPTIONAL)

echo -e "\n\n\n######################"
echo -e "###### STEP 2.2 ######"
echo -e "######################\n"

####################

cd $WD1_spe/01-Raw_data

### REPAIR
## Repair each library using repair.sh.
echo -e "\nREPAIR..."

for SRR in $Samples_to_repair; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/01-Raw_data/outputs/stdout_Repair_$SRR.log --quiet --exclusive $F task_Repair $SRR &
done
wait

### FASTQC
## Execute Fastqc.
echo -e "\nFASTQC..."

for SRR in $Samples_to_repair; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/01-Raw_data/outputs/stdout_Fastqc_$SRR.log --quiet --exclusive $F task_Fastqc1 $SRR $SLURM_CPUS_PER_TASK &
done
wait

### MULTIQC
## Execute Multiqc.
echo -e "\nMULTIQC..."

rm -r ./multiqc
mkdir multiqc
>./outputs/stdout_Multiqc.log
multiqc ./fastqc -o ./multiqc -f >> ./outputs/stdout_Multiqc.log 2>&1

####################

cd $WD1_spe/02-Trimmed_data

### TRIMMING
## Trim the libraries.
echo -e "\nTRIMMING..."

for SRR in $Samples_to_repair; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --output $WD1_spe/02-Trimmed_data/outputs/stdout_Trimming_$SRR.log --quiet --exclusive $F task_Trimming $SRR $trim_adapters $SLURM_CPUS_PER_TASK &
done
wait


