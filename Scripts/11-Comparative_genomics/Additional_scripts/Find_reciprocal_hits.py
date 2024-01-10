#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- FIND RECIPROCAL HITS.

This script will create a table with the reciprocal hits of two species.

---------

Created on Wed Jul 20 12:30:00 2022

Last Modified:
    - Wed Jul 20 12:30:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='FindReciprocalBlastHits V1', 
            description='''This program creates a TXT table which contains the blast reciprocal hits of two species.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--spe1", 
            type=str, 
            nargs=1,
            help="Specie 1."
            )
    parser.add_argument(
            "--spe2", 
            type=str, 
            nargs=1,
            help="Specie 2."
            )
    parser.add_argument(
            "--spe1vsspe2", 
            type=str, 
            nargs=1,
            help="Absolute path of TSV table containing the blast results from specie 1 vs specie 2 comparison."
            )
    parser.add_argument(
            "--spe2vsspe1", 
            type=str, 
            nargs=1,
            help="Absolute path of TSV table containing the blast results from specie 2 vs specie 1 comparison."
            )
    parser.add_argument(
            "--output", 
            type=str, 
            nargs=1, 
            help="Absolute path of output table."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        spe1 = args.spe1[0]
        spe2 = args.spe2[0]
        spe1tospe2_file = args.spe1vsspe2[0]
        spe2tospe1_file = args.spe2vsspe1[0]
        output_file = args.output[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Find_reciprocal_hits.py 
    #--spe1 cme
    #--spe2 csa
    #--spe1vsspe2 /mnt/doctorado/.../cme_to_csa.tsv
    #--spe2vsspe1 /mnt/doctorado/.../csa_to_cme.tsv
    #--output /mnt/doctorado/.../cme_and_csa.tsv
    
    """
    spe1 = "cma"
    spe2 = "cme"
    spe1tospe2_file = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn_genes/nr/03-Blastn/cma_to_cme.tsv"
    spe2tospe1_file = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn_genes/nr/03-Blastn/cme_to_cma.tsv"
    output_file = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn_genes/nr/04-Reciprocal_hits/cma_and_cme.tsv"
    """
    
    spe1tospe2_file_open = open(spe1tospe2_file, 'r')
    spe2tospe1_file_open = open(spe2tospe1_file, 'r')
    output_file_open = open(output_file, 'w')
    
    spe1tospe2_list = []
    for line in spe1tospe2_file_open:
        line = line.strip().split("\t")
        spe1tospe2_list.append([line[0], line[1]])
        
    spe2tospe1_list = []
    for line in spe2tospe1_file_open:
        line = line.strip().split("\t")
        spe2tospe1_list.append([line[0], line[1]])
    
    
    for x in spe1tospe2_list:
        for y in spe2tospe1_list:
            if x[0] == y[1] and x[1] == y[0]:
                output_file_open.write(spe1 + "\t" + spe2 + "\t" + x[0] + "\t" + x[1] + "\n")
    
    for y in spe2tospe1_list:
        for x in spe1tospe2_list:
            if y[0] == x[1] and y[1] == x[0]:
                output_file_open.write(spe2 + "\t" + spe1 + "\t" + y[0] + "\t" + y[1] + "\n")
    
    output_file_open.close()
    


# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
