#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- ADD INFO TO TSV FILE COMING FROM SWISSPROT BLASTX (DIAMOND).

The sseqid column contained in the TSV file indicates the id of the swissprot's 
protein, but this is only an ID. Then, this script adds info to know what is 
this ID using the headers of the swissprot database fasta file which can be 
downloaded from https://www.uniprot.org/uniprot/?query=reviewed:yes. If the fasta 
file has not more info than the ID, this script will be useless.

---------

Created on Fri Jun 08 10:45:00 2022

Last Modified:
    - Fri Jun 08 10:45:00 2022 --> Initial code.

@author: pasviber - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse


# FUNCTIONS

def TSVtoLIST (tsv):
    """
    This function generates a list from the information contained in 
    the TSV file coming from diamond blastx using as a reference (source) 
    the SwissProt database.
    """
    
    List = []
    
    f_in = open(tsv, 'r')
    for line in f_in:
        line = line.strip().split("\t")
        List.append(line)
    
    return List


def FASTAtoDICT (fasta):
    """
    This function creates a dictionary using only the header. The info from the
    begining to the first space will be the ID of the dictionary as well as the
    sseqid in the table. The rest of the header will be the info which will be
    added to the TSV table (sinfo).
    """
    
    Dic = {}
    
    f_in = open(fasta, 'r')
    for line in f_in:
        line = line.strip()
        if line[0] == ">":
            header = line[1:]
            list_ = header.split(" ")
            Dic[list_[0]] = " ".join(list_[1:])
    
    return Dic
    

def Write_table (tsv, records, info):
    """
    This function writes the output table.
    """
    
    f_out = open(tsv, 'w')
    f_out.write("qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tpident\tgaps\tlength\tevalue\tbitscore\tsinfo\n")
    for record in records:
        ID = record[2]
        f_out.write("\t".join(record) + "\t" + info[ID] + "\n")
    
    f_out.close()
    

# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='AddInfoToSwissprotBLASTxTSV V1',
            description='''This program adds info to the Swissprot BLASTx TSV table results''',  
            #formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "-i", 
            "--tsv-initial", 
            type=str, 
            nargs=1,
            help="Absolute path of initial Swissprot BLASTx table results."
            )
    parser.add_argument(
            "-d", 
            "--fasta-db", 
            type=str, 
            nargs=1,
            help="Absolute path of RNAcentral fasta file (database)."
            )
    parser.add_argument(
            "-f", 
            "--tsv-final", 
            type=str, 
            nargs=1,
            help="Absolute path of final Swissprot BLASTx table results."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        tsv1 = args.tsv_initial[0]
        fasta = args.fasta_db[0]
        tsv2 = args.tsv_final[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # tsv1 =  "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP3/SwissProt/diamond_output_temp.tsv"
    # fasta = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info/Swissprot/uniprot_sprot.fasta"
    # tsv2 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP3/SwissProt/diamond_output.tsv"
    
    
    ### PIPELINE
    ## Get the information 1.
    records = TSVtoLIST (tsv1)
    
    ## Get the information 2.
    info = FASTAtoDICT (fasta)
    
    ## Build the table.
    Write_table (tsv2, records, info)
    

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
