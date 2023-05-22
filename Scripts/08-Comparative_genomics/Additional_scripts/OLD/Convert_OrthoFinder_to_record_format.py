#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- CONVERT ORTHOFINDER RESULTS TO A RECORD FORMAT FILE.

This script will create a record format file using the OrthoFinder results.

---------

Created on Fri Sep 23 10:00:00 2022

Last Modified:
    - Fri Sep 23 10:00:00 2022 --> Initial code.

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
            prog='ConvertOrthoFinderResultsTORecordFormat V1', 
            description='''This program creates a record format file using the orthofinder results.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--OF-path", 
            type=str, 
            nargs=1,
            help="OrthoFinder directory with the results."
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
        OF_path = args.OF_path[0]
        output_table = args.output_table[0]
        spe_ref = args.spe_ref[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    """
    # python3 ./Convert_OrthoFinder_to_record_format.py 
    #--OF-path /mnt/doctorado/.../OrthoFinder
    #--output-table /mnt/doctorado/.../Record_format_files/Orthologs_record_format_OrthoFinder.tsv
    #--spe-ref cme
    
    OF_path = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/01-Orthologs/Orthofinder" 
    output_table = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/01-Orthologs/Record_format_files/Orthologs_record_format_OrthoFinder.tsv" 
    spe_ref = "cme"
    
    ## Other variables.
    orthologs_path = OF_path + "/" + os.listdir(OF_path)[0] + "/Orthologues/Orthologues_" + spe_ref
    
    ## Getting the species name.
    species = [file.replace(spe_ref + "__v__", "").replace(".tsv", "") for file in os.listdir(orthologs_path) if spe_ref in file]
    
    ## Create a dictionary
    ortho = {}
    
    for i in range(len(species)):
        file = orthologs_path + "/" + spe_ref + "__v__" + species[i] + ".tsv"
        file_open = open(file, 'r')
        for line in file_open:
            line = line.strip().split("\t")
            if line[0] != "Orthogroup":
                ref = line[1].split(", ")
                other = line[2].split(", ")
                
                for ID in ref:
                    if ID not in list(ortho.keys()):
                        ortho[ID] = [[] for spe in species]
                        ortho[ID][i] = other
                    else:
                        ortho[ID][i] = other
                        
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
                output_file_open.write(spe + "\t" + val + "\t" + val + "|" + spe + "\tOrthoFinder\n")
        output_file_open.write("=\n")
    output_file_open.close()
                
     

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
