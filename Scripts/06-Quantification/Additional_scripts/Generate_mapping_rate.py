#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### This script will generate mapping_rate.tsv file. 
### It has to be run when the pseudomapping step is finished.


## MODULES

import os

## VARIABLES

WD="/mnt/doctorado/4-MELON_Cristina/Transcriptoma/Results/03-New_analysis/07-Quantification_salmon/LncRNAs/"
samples="/mnt/doctorado/4-MELON_Cristina/Transcriptoma/Additional_info/Samples_list/samples.txt"

## PIPELINE

os.chdir(WD)

file_list = open(samples, "r+")
output = open("mapping_rate.tsv", "w")
    
output.write("sample\tmapping_rate\n")
    
for line in file_list:
    line=line.rstrip()

    log_file = open("./03-Quant/" + line + "/logs/salmon_quant.log", "r+")
    
    for line2 in log_file.readlines():
        if "Mapping rate =" in line2:
            mapping_rate = line2.rstrip().split(" ")[-1]
            output.write(line + "\t" + mapping_rate + "\n")
            
