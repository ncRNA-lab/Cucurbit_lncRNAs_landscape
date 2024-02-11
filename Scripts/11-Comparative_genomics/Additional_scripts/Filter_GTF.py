#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- FILTER GTF BY TRANSCRIPT_ID.

This script filters a GTF file by transcript_id using a TXT file
which contains a transcript_id list.

---------

Created on Fri Mar 25 17:45:00 2022

Last Modified:
    - Fri Mar 25 17:45:00 2022 --> Initial code.
    - Sun May 29 13:50:00 2022 --> Add arguments command line.

@author: pvbermell - Pascual Villalba Bermell

"""


## MODULES

import sys
import argparse


## FUNCTIONS

def GTFtoDICT (gtf):
    """
    This function generates a dictionary from the information contained in 
    a GTF file keeping the structure 'transcript' - 'exon'.
    """
    
    f_in = open(gtf, 'r')
    Transcripts = {}
    for line in f_in:
        line = line.strip().split("\t")
        
        if line[0][:2] == "##":
            continue
        
        elif line[2] == "transcript":
            ID = line[8].split(";")[0].split(" ")[1].replace('"', '')
            Transcripts[ID] = [line, []]
            
        elif line[2] == "exon":
            ID = line[8].split(";")[0].split(" ")[1].replace('"', "")
            Transcripts[ID][1].append(line)
    
    return Transcripts

def FilterDictGTF (Dict_gtf, IDs, gtf_final_path):
    """
    This function writes the output GTF file keeping only the transcript_ids of
    the list.
    """
    
    f_in = open(IDs, 'r')
    IDs = [line.strip() for line in f_in]
    IDs_sorted = [ID for ID in list(Dict_gtf.keys()) if ID in IDs]
    
    f_out = open(gtf_final_path, 'w')
    for ID in IDs_sorted:
        trans = Dict_gtf[ID][0]
        exons = Dict_gtf[ID][1]
        f_out.write("\t".join(trans) + "\n")
        
        if exons:
            for exon in exons:
                f_out.write("\t".join(exon) + "\n")
    
    f_out.close()
    

# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
        prog='FilterGTF V2',
        description='''This program filters a GTF file by transcript_id''',
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    
    parser.add_argument(
        "-I", 
        "--gtf-initial", 
        type=str, 
        nargs=1,
        help="Absolute path of gtf file."
        )
    parser.add_argument(
        "-F", 
        "--gtf-final", 
        type=str, 
        nargs=1,   
        help="Absolute path of new gtf file."
        )
    parser.add_argument(
        "-i", 
        "--ids", 
        type=str, 
        nargs=1,
        help="Absolute path of ID-transcriptd that you want to keep."
        )
    parser.add_argument(
        '--version', 
        action='version', 
        version='%(prog)s 2.0'
        )
    
    args = parser.parse_args()
    
    try:
        gtf_i = args.gtf_initial[0]
        gtf_f = args.gtf_final[0]
        IDs = args.ids[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # gtf_i = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/01-LncRNAs_and_Genes/High/intergenic/cme_lncRNAs_temp.gtf"
    # gtf_f = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/01-LncRNAs_and_Genes/High/intergenic/cme_lncRNAs.gtf"
    # IDs = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/01-LncRNAs_and_Genes/High/intergenic/cme_lncRNAs_ids.txt" 
    
    ### PIPELINE
    ## Get the information.
    Dict_gtf = GTFtoDICT (gtf_i)
    
    ## Filter the gtf dict by transcript ID and write the new gtf.
    FilterDictGTF(Dict_gtf, IDs, gtf_f)
    

# CALL THE MAIN PROGRAM

if __name__ == '__main__':
    """
    Call the main program.
    """
    main()

