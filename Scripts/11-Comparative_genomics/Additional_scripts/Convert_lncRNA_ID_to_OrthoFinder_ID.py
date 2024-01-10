#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- CONVERT LNCRNA IDS TO ORTHOFINDER IDS.

This script converts lncRNA IDs to the OrthoFinder IDs, I mean in the correct
format to execute OrthoFinder again.

---------

Created on Tue Dec 06 23:00:00 2022

Last Modified:
    - Tue Dec 06 23:00:00 2022 --> Initial code.

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
            prog='ConvertLncRNAIDToOFID', 
            description='''This program converts lncRNA IDs to the OrthoFinder IDs.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument(
            "--path", 
            type=str, 
            nargs=1,
            help="Absolute path of WorkingDirectory."
            )
    parser.add_argument(
            "--comb", 
            type=str, 
            nargs=1,
            help="Absolute path of codes table."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        WD = args.path[0]
        comb = args.comb[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Convert_lncRNA_ID_to_OrthoFinder_ID.py 
    #--path /mnt/doctorado/...
    #--comb /mnt/doctorado/.../....txt
    
    """
    WD = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/OrthoFinder/04-OrthoFinder/High/intergenic/Results/Results_Dec06/WorkingDirectory"
    comb = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/OrthoFinder/03-Blastn/Codes_combination.txt"
    """
    
    ## Create the dictionary of lncRNA IDs.
    SeqIDs_file_open = open(WD + "/SequenceIDs.txt", 'r')
    IDs_dic = {}
    for line in SeqIDs_file_open:
        spe = line.strip().split(": ")[0].split("_")[0]
        new_ID = line.strip().split(": ")[0]
        old_ID = line.strip().split(": ")[1]
        if spe in IDs_dic.keys():
            IDs_dic[spe][old_ID] = new_ID
        else:
            IDs_dic[spe] = {}
            IDs_dic[spe][old_ID] = new_ID
    SeqIDs_file_open.close()
    
    ## Convert lncRNA IDs to OrthoFinder IDs in each blast table.
    Comb_file_open = open(comb, 'r')
    for line1 in Comb_file_open:
        cod1 = line1.strip().split(" ")[0]
        cod2 = line1.strip().split(" ")[1]
        file_in = WD + "/Blast" + str(cod1) + "_" + str(cod2) + "_temp2.txt"
        file_out = WD + "/Blast" + str(cod1) + "_" + str(cod2) + ".txt"
        
        file_in_open = open(file_in, 'r')
        file_out_open = open(file_out, 'w')
        
        for line2 in file_in_open:
            info = line2.strip().split("\t")
            ID1 = info[0]
            ID2 = info[1]
            file_out_open.write(IDs_dic[cod1][ID1] + "\t" + IDs_dic[cod2][ID2] + "\t" + "\t".join(info[2:]) + "\n")
        
        file_in_open.close()
        file_out_open.close()
    

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()      
            
