#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- CREATE TWO FASTA FILES IN MIRENA FORMAT.

This script will create a fasta file from the blastn output ((1) lncRNAs vs
miRNA-precursors and (2) lncRNAs aligned in the previous step vs mature miRNAs) 
and a fasta file which contains all the lncRNAs sequences. The new fasta file
will be in the MIReNA format.

---------

Created on Fri Jun 30 12:30:00 2022

Last Modified:
    - Fri Jun 30 12:30:00 2022 --> Initial code.
    - Sun Jul 03 17:30:00 2022 --> Length of precursor.

@author: pasviber - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
from Bio import SeqIO


# FUNCTIONS

def TSVtoLIST (tsv):
    """
    This function generates a list from the information contained in the TSV 
    file coming from blastn using as a reference (source) a mature miRNA or 
    miRNA precursor database.
    """
    
    List = []
    f_in = open(tsv, 'r')
    n = 0
    for line in f_in:
        if n != 0:
            line = line.strip().split("\t")
            List.append(line)
        n += 1
    
    return List


def FASTAtoDICT (fasta):
    """
    This function creates a dictionary using the header as ID and the sequence 
    as value.
    """
    
    seq_dict = {rec.id : str(rec.seq) for rec in SeqIO.parse(fasta, "fasta")}
    
    return seq_dict


def CreateNewDictFasta(blastn_info_mat, blastn_info_prec, fasta_info, L, remove_IDs_list):
    """
    This function creates a dictionary with the info which will be used to write 
    the fasta file. The length limit of the miRNA precursor is L, so this function
    also cuts the sequences. Moreover, the mature miRNA region detected must be 
    in the middle of the miRNA-precursor region detected.
    """
    
    Dic_fasta = {}
    Unique_combinations_mat = []
    
    for record in blastn_info_mat:
        mismatch = int(record[13])
        length_mat = int(record[14])
        trans_id = record[0]
        
        if mismatch <= 2 and length_mat >= 18 and trans_id not in remove_IDs_list:
            trans_len = int(record[1])
            trans_start_vs_mat = int(record[4])
            trans_end_vs_mat = int(record[5])
            
            ID_mat = trans_id + "_" + str(trans_start_vs_mat) + "_" + str(trans_end_vs_mat)
            if ID_mat not in Unique_combinations_mat:
                Unique_combinations_mat.append(ID_mat)
                
                Unique_combinations_prec = []
                precs = [x for x in blastn_info_prec if x[0] == trans_id]
                for prec in precs:
                    trans_start_vs_prec = int(prec[4])
                    trans_end_vs_prec = int(prec[5])
                    length_prec = int(prec[13])
                    
                    if length_prec < L:
                        diff = L - length_prec
                        add = diff//2
                        trans_start_vs_prec_add = trans_start_vs_prec - add
                        if trans_start_vs_prec_add <= 0:
                            trans_start_vs_prec_add = 1
                        trans_end_vs_prec_add = trans_end_vs_prec + add
                        if trans_end_vs_prec_add > trans_len:
                            trans_end_vs_prec_add = trans_len
                    else:
                        trans_start_vs_prec_add = trans_start_vs_prec
                        trans_end_vs_prec_add = trans_end_vs_prec
                    
                    ID_prec = trans_id + "_" + str(trans_start_vs_prec_add) + "_" + str(trans_end_vs_prec_add)
                    if trans_start_vs_prec_add <= trans_start_vs_mat and trans_end_vs_mat <= trans_end_vs_prec_add:
                        if ID_prec not in Unique_combinations_prec:
                            Unique_combinations_prec.append(ID_prec)
                            
                            trans_seq = fasta_info[trans_id]
                            i = trans_start_vs_prec_add - 1
                            j = trans_end_vs_prec_add
                            trans_seq = trans_seq[i:j]
                            
                            before = trans_start_vs_mat - trans_start_vs_prec_add
                            after = trans_end_vs_prec_add - trans_end_vs_mat
                            
                            header = ">" + trans_id + " miRNA:" + ID_mat + " precursor:" + ID_prec + " before:" + str(before) + " after:" + str(after)
                            Dic_fasta[header] = trans_seq
                    else:
                        print("In " + trans_id + ", " + ID_mat + " is not contained in " + ID_prec + " alignment region.")
    
    return Dic_fasta
    

def Write_fasta (Dic_fasta, path):
    """
    This function writes the fasta files for MIReNA.
    """
    
    f_out = open(path, 'w')
    
    for header in list(Dic_fasta.keys()):
        f_out.write(header + "\n")
        f_out.write(Dic_fasta[header] + "\n")
        
    f_out.close()
    

# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='FastaMIReNACreator V1', 
            description='''This program creates a fasta file which will be used in MIReNA.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "-m", 
            "--blastn-mat", 
            type=str, 
            nargs=1,
            help="Absolute path of BLASTn table results using the mature miRNAs as a reference."
            )
    parser.add_argument(
            "-p", 
            "--blastn-prec", 
            type=str, 
            nargs=1,
            help="Absolute path of BLASTn table results using the precursor-miRNAs as a reference."
            )
    parser.add_argument(
            "-r", 
            "--remove-IDs", 
            type=str, 
            nargs='?',
            help="Absolute path of file with the transcript IDs to remove."
            )
    parser.add_argument(
            "-i", 
            "--fasta-initial", 
            type=str, 
            nargs=1, 
            help="Absolute path of initial fasta file."
            )
    parser.add_argument(
            "-o", 
            "--fasta-MIReNA", 
            type=str, 
            nargs=1, 
            help="Absolute path of MIReNA fasta file."
            )
    parser.add_argument(
            "-l", 
            "--length", 
            type=int, 
            nargs=1, 
            help="miRNA-precursor length."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        tsv_mat = args.blastn_mat[0]
        tsv_prec = args.blastn_prec[0]
        remove_IDs = args.remove_IDs
        fasta_initial = args.fasta_initial[0]
        fasta_MIReNA = args.fasta_MIReNA[0]
        length_prec = args.length[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # tsv_mat = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP5/PmiREN/output_blastn_mat.tsv"
    # tsv_prec = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP5/PmiREN/output_blastn_prec.tsv"
    # remove_IDs = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP5/PmiREN/prec_ID-100-200.txt"
    # fasta_initial = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta"
    # fasta_MIReNA = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP5/PmiREN/MIReNA-300.fasta"
    # length_prec = 300
    
    
    ### PIPELINE
    ## Get the blastn information.
    blastn_info_mat = TSVtoLIST (tsv_mat)
    blastn_info_prec = TSVtoLIST (tsv_prec)
    
    ## Get the fasta information.
    fasta_info = FASTAtoDICT (fasta_initial)
    
    ## Get the remove IDs information.
    if remove_IDs != None:
        f_in = open(remove_IDs, 'r')
        remove_IDs_list = [line for line in f_in]
    else:
        remove_IDs_list = []
    
    ## Create a dict for each new fasta file.
    New_dic_fasta = CreateNewDictFasta(blastn_info_mat, blastn_info_prec, fasta_info, length_prec, remove_IDs_list)
    
    ## Save the fasta files.
    Write_fasta (New_dic_fasta, fasta_MIReNA)
    

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
