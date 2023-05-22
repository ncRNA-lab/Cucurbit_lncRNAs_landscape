#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- JOIN RECORD FORMAT FILES.

This script will join two o more record format files.

---------

Created on Mon Sep 26 17:00:00 2022

Last Modified:
    - Fri Mon Sep 26 17:00:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import os

# FUNCTIONS

def CreateDictByFile (File):
    
    ortho = {}
    file_open = open(File, 'r')
    for line in file_open:
        line = line.strip().split("\t")
        if len(line) < 4 and line[0] != "=":
            ID = line[1]
            ortho[ID] = []
        elif len(line) == 4:
            ortho[ID].append(line)
    
    return ortho
    


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
            "--input-folder", 
            type=str, 
            nargs=1,
            help="Folder with all the files in recod format that you want to join."
            )
    parser.add_argument(
            "--output_table", 
            type=str, 
            nargs=1,
            help="Folder with all the files in recod format that you want to join."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        input_folder = args.input_folder[0]
        output_table = args.output_table[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Convert_OrthoFinder_to_record_format.py 
    #--input-folder /mnt/doctorado/.../Record_format_files
    #--output-table /mnt/doctorado/.../Orthologs_record_format_final.tsv
    """
    
    input_folder = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/01-Orthologs/Record_format_files"
    output_table = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/01-Orthologs/Orthologs_record_format_FINAL.tsv"
    ## Getting the species name.
    files = os.listdir(input_folder)
    softwares = [file.replace("Orthologs_record_format_", "").replace(".tsv", "") for file in os.listdir(input_folder)]
    
    ## Create a dictionary with orthologs coming from each software.
    ORTHO = {}
    
    for i in range(len(files)):
        ortho = CreateDictByFile (input_folder + "/" + files[i])
        ORTHO[softwares[i]] = ortho
        print("Number of genes with relatioships across all species according to " + softwares[i] + ": " + str(len(list(ortho.keys()))))
        if i == 0:
            Genes_overlap = list(ortho.keys())
        else:
            Genes_overlap = list(set(Genes_overlap) & set(list(ortho.keys())))
    print("Number of genes with overlap: " + str(len(Genes_overlap)))
                     
    ## Filter the dictionary. Genes overlapping all softwares.
    ORTHO_FINAL = {}
    
    for ID in Genes_overlap:
        for i in range(len(softwares)):
            List_initial = [element[2] for element in ORTHO[softwares[i]][ID]]
            if i == 0:
                List_final = List_initial
            else:
                List_final = list(set(List_final) & set(List_initial))
        
        New_List_Final = [[x.split("-")[1], x.split("-")[0], x, ";".join(softwares)] for x in List_final]
        ORTHO_FINAL[ID] = New_List_Final
    
    ## Write the record format file.
    output_file_open = open(output_table, 'w')
    
    for ID in list(ORTHO_FINAL.keys()):
        output_file_open.write(ID.split("-")[0] + "\t" + ID + "\n")
        values = ORTHO_FINAL[ID]
        for val in values:
            output_file_open.write("\t".join(val) + "\n")
        output_file_open.write("=\n")
    output_file_open.close()
     

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
