#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- CONVERT INPARANOID RESULTS TO A RECORD FORMAT FILE.

This script will create a record format file using the Inparanoid results.

---------

Created on Mon Sep 26 12:00:00 2022

Last Modified:
    - Mon Sep 26 12:00:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import os


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    """
    parser = argparse.ArgumentParser(
            prog='ConvertInparanoidResultsTORecordFormat V1', 
            description='''This program creates a record format file using the inparanoid results.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--I-path", 
            type=str, 
            nargs=1,
            help="Inparanoid directory with the results."
            )
    parser.add_argument(
            "--output-table", 
            type=str, 
            nargs=1,
            help="File with orthologs in record format."
            )
    parser.add_argument(
            "--spe-ref", 
            type=str, 
            nargs=1,
            help="Specie which will be used as a reference to get all the information."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        I_path = args.I_path[0]
        output_table = args.output_table[0]
        spe_ref = args.spe_ref[0]
        
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Convert_Inparanoid_to_record_format.py 
    #--I-path /mnt/doctorado/.../Inparanoid
    #--output-table /mnt/doctorado/.../Record_format_files/Orthologs_record_format_Inparanoid.tsv
    #--spe-ref cme
    """
    #I_path = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/01-Orthologs/Inparanoid" 
    #output_table = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/01-Orthologs/Record_format_files/Orthologs_record_format_Inparanoid.tsv" 
    #spe_ref = "cme" 
    
    
    ## Getting the species name.
    species = [file.replace("table." + spe_ref + ".fa-", "").replace(".fa", "") for file in os.listdir(I_path) if "table." + spe_ref + ".fa-" in file and "SQLtable" not in file]
    
    ## Create a dictionary
    ortho = {}
    
    for i in range(len(species)):
        file = I_path + "/" + "table." + spe_ref + ".fa-" + species[i] + ".fa"
        file_open = open(file, 'r')
        for line in file_open:
            line = line.strip().split("\t")
            if line[0] != "OrtoID":
                ref = line[2].split(" ")
                ref_remove_numbers = [ref[j] for j in range(0, len(ref), 2)]
                other = line[3].split(" ")
                other_remove_numbers = [other[j] for j in range(0, len(other), 2)]
                
                for ID in ref_remove_numbers:
                    if ID not in list(ortho.keys()):
                        ortho[ID] = [[] for spe in species]
                        ortho[ID][i] = other_remove_numbers
                    else:
                        ortho[ID][i] = other_remove_numbers
                        
    ## Filter the dictionary. Orthologue in all species.
    ortho_filt = {}
    for ID in list(ortho.keys()):
        values = ortho[ID]
        values_removimg_empty_lists = [val for val in values if val]
        if len(values_removimg_empty_lists) == len(species):
            ortho_filt[ID] = values
            
    ## Write the record format file.
    output_file_open = open(output_table, 'w')
    
    for ID in list(ortho_filt.keys()):
        output_file_open.write(ID + "\t" + ID + "|" + spe_ref + "\n")
        for i in range(len(species)):
            spe = species[i]
            values = ortho_filt[ID][i]
            for val in values:
                output_file_open.write(spe + "\t" + val + "\t" + val + "|" + spe + "\tInparanoid\n")
        output_file_open.write("=\n")
    output_file_open.close()
                
     

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
