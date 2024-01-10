#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- FILTER BEAGLE RESULTS.

This script will create a table with the hits of two species.

---------

Created on Thu Sep 15 11:49:00 2022

Last Modified:
    - Thu Sep 15 11:49:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import pandas as pd


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='FilterBeagleResults V1', 
            description='''This program creates a table with the hits of two species.''', 
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
            "--input", 
            type=str, 
            nargs=1,
            help="Absolute path of Beagle results."
            )
    parser.add_argument(
            "--output", 
            type=str, 
            nargs=1, 
            help="Absolute path of Beagle results filtered."
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
        input_file = args.input[0]
        output_file = args.output[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Find_reciprocal_hits.py 
    #--spe1 cme
    #--spe2 csa
    #--input /mnt/doctorado/.../cme_and_csa_beagle_outfile.txt
    #--output /mnt/doctorado/.../cme_and_csa_beagle_outfile_filtered.txt
    
    
    input_file_open = open(input_file, 'r')
    output_file_open = open(output_file, 'w')
    
    list_results = []
    for line in input_file_open:
        if line[:6] == ">MSTRG":
            line = line[1:].strip().split("|")
            id_1 = line[0]
            id_2 = line[1]
            pval = float(line[6].split(":")[1].replace(",", "."))
            zscore = float(line[7].split(":")[1].replace(",", "."))
            list_results.append([id_1, id_2, pval, zscore])
    
    df = pd.DataFrame(list_results)   
    df.columns=["id_1", "id_2", "pval", "zscore"]
    
    df['zscore'] = df['zscore'].astype(float)
    df['pval'] = df['pval'].astype(float)
    
    df = df[df["pval"]<0.05]
    df = df[df["zscore"]>3]
    
    inds = df.groupby(['id_1'])['zscore'].transform(max) == df['zscore']
    df = df[inds]
    df.reset_index(drop=True, inplace=True)                   
    
    final_list = df.values.tolist()
    
    for hit in final_list:
        output_file_open.write(spe1 + "\t" + spe2 + "\t" + hit[0] + "\t" + hit[1] + "\n")
        output_file_open.write(spe2 + "\t" + spe1 + "\t" + hit[1] + "\t" + hit[0] + "\n")
    
    output_file_open.close()


# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
