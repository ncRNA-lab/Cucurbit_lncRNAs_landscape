#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### This script will generate masked_genome_percentage.tsv file. 
### It has to be run when the repeatmasker step is finished.


## MODULES

import os
import sys
import argparse

## VARIABLES

parser = argparse.ArgumentParser(
        prog='GetMaskedGenomePercentage V1',
        description='''This program gets the masked genome percentage.''', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
parser.add_argument(
        "--path", 
        type=str, 
        nargs=1,
        help="Absolute path of directory."
        )
parser.add_argument(
        "--species", 
        type=str, 
        nargs="+",
        help="Species list."
        )
parser.add_argument(
        '--version', 
        action='version', 
        version='%(prog)s 1.0'
        )

args = parser.parse_args()

try:
    WD = args.path[0]
    species = args.species[:]
     
except:
    print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
    parser.print_help()
    sys.exit()

#WD="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/11-TEs_and_genomic_repeats/01-Repeat_calling/02-RepeatMasker"
#species=["car", "cla", "cma", "cme", "cmo", "cpe", "csa", "lsi", "mch"]

## PIPELINE

os.chdir(WD)

output_ = open("masked_genome_percentage.tsv", "w")
output_.write("spe\tmasked_perc\n")
    
for spe in species:
    input_ = open("./" + spe + "/" + spe + ".fa.tbl", "r")
    for line in input_:
        line = line.strip()
        if "bases masked:" in line:
            per = line.split(" ")[-2]
            output_.write(spe + "\t" + per + "\n")
            break
    input_.close()
output_.close()
            
