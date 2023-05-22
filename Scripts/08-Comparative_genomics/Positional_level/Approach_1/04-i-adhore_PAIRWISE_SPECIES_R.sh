#!/bin/bash

####### VARIABLES
WD="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results"
SP="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Softwares"
F="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Scripts/Pascual/08-comparative_genomics/Synteny/Functions.sh"
Species_list="car cla cma cme cmo cpe csa lsi mch"

####### NEW AND OTHER VARIABLES
WD1=$WD/08-comparative_genomics/Synteny/r


####### SOFTWARES
### ORTHOFINDER
export APATH=$SP/3.0/i-ADHoRe-3.0/i-adhore
export PATH=$PATH:${APATH}
export PATH=$PATH:${APATH}/bin


####### PIPELINE: PAIRWISE_SPECIES
echo -e "\n\nPAIRWISE_SPECIES..."

## STEP 5: Execute i-adhore
echo -e "\n\nSTEP 5: Execute i-adhore...\n"
cd $WD1/03-Adhore/PAIRWISE_SPECIES

Combinations=$WD1/03-Adhore/PAIRWISE_SPECIES/Species_combination.txt

n_comb=$(wc -l $Combinations | cut -d" " -f1)
for j in `seq 1 $n_comb`; do
	comb=$(cat $Combinations | head -n $j | tail -n 1)
	spe1=$(echo $comb | cut -d" " -f1)
	spe2=$(echo $comb | cut -d" " -f2)
	echo -e "\tCombination: "$spe1"-"$spe2"..."
	cd $WD1/03-Adhore/PAIRWISE_SPECIES/$spe1\-$spe2
	if [ -f "stdout_collinear.log" ]; then
		rm stdout_collinear.log
	fi
	if [ -f "stdout_cloud.log" ]; then
		rm stdout_cloud.log
	fi
	if [ -f "stdout_hybrid.log" ]; then
		rm stdout_hybrid.log
	fi
	i-adhore idhore_parameters_collinear.ini >> stdout_collinear.log 2>&1
	i-adhore idhore_parameters_cloud.ini >> stdout_cloud.log 2>&1
	i-adhore idhore_parameters_hybrid.ini >> stdout_hybrid.log 2>&1
done

