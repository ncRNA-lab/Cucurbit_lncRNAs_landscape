#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- PREPARE MOTIF LEVEL ANALYSIS

This script creates a directory where each family inferred in the conservation 
analysis is stored in a separate fasta file. Families can be conserved at position 
or sequence level.

---------

Created on Fri Feb 10 14:00:00 2023

Last Modified:
    - Fri Feb 10 14:00:00 2023 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import os
import random
import multiprocessing


# FUNCTIONS

def FamiliesDict (file):
    
    file_open = open(file, 'r')
    dic = {}
    for line in file_open:
        line = line.strip().split("\t")
        if line[0] not in dic.keys():
            dic[line[0]] = [line[1]]
        else:
            dic[line[0]].append(line[1])
    file_open.close()
    
    dic_filt = {}
    for fam in list(dic.keys()):
        if len(dic[fam]) > 1:
            dic_filt[fam] = dic[fam]
    
    return dic_filt

def LncRNAsAndSequencesDict (file):
    
    file_open = open(file, 'r')
    dic = {}
    for line in file_open:
        line = line.strip()
        if line[0] == ">":
            ID = line[1:]
            dic[ID] = ""
        else:
            dic[ID] += line
    file_open.close()
    
    return dic

def LncRNAsListsBySpecieDict (file):
    
    file_open = open(file, 'r')
    dic = {}
    for line in file_open:
        line = line.strip()
        if line[0] == ">":
            spe = line.split("-")[1]
            ID = line[1:]
            if spe not in list(dic.keys()):
                dic[spe] = [ID]
            else:
                dic[spe].append(ID)
    file_open.close()
    
    return dic

def CreateFastaFilesREAL (Families_dic, LncRNAs_and_sequences_dic, path_out):
    
    for fam, lncs in Families_dic.items():
        fam_fasta = open(path_out + "/real/" + fam + "_real.fasta", "w")
        for lnc in lncs:
            fam_fasta.write(">" + lnc + "\n" + LncRNAs_and_sequences_dic[lnc] + "\n")
        fam_fasta.close()
        
def CreateFastaFilesSIMULATIONS (iterations, Families_dic, LncRNAs_and_sequences_dic, LncRNAs_lists_by_specie_dic, path_out):

    for i in iterations:
        os.makedirs(path_out + "/simulations/iter_" + str(i))
        for fam, lncs in Families_dic.items():
            random_lncs = []
            for lnc in lncs:
                spe = lnc.split("-")[1]
                LncRNAs_spe = LncRNAs_lists_by_specie_dic[spe]
                Res = False
                while Res == False:
                    random_lnc = random.choice(LncRNAs_spe)
                    if random_lnc not in random_lncs and random_lnc != lnc:
                        random_lncs.append(random_lnc)
                        Res = True
                
            fam_fasta = open(path_out + "/simulations/iter_" + str(i) + "/" + fam + "_random.fasta", "w")
            for random_lnc in random_lncs:
                fam_fasta.write(">" + random_lnc + "\n" + LncRNAs_and_sequences_dic[random_lnc] + "\n")
            fam_fasta.close()


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='PrepareMotifLevelAnalysis V1', 
            description='''This script creates a directory where each family inferred \
                in the conservation analysis is stored in a separate fasta file. Families \
                can be conserved at position or sequence level.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--input-table", 
            type=str, 
            nargs=1,
            help="Absolute path of binary table of gene or lncRNA families."
            )
    parser.add_argument(
            "--input-fasta", 
            type=str, 
            nargs=1,
            help="Absolute path of lncRNAs fasta file."
            )
    parser.add_argument(
            "--path-out", 
            type=str, 
            nargs=1,
            help="Absolute output path."
            )
    parser.add_argument(
            "--n-iter", 
            type=int, 
            nargs=1,
            help="Number of simulations."
            )
    parser.add_argument(
            "--n-proc", 
            type=int, 
            nargs=1,
            help="Number of processes to be launched."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        table_in = args.input_table[0]
        fasta_in = args.input_fasta[0]
        path_out = args.path_out[0]
        n_iter = args.n_iter[0]
        n_proc = args.n_proc[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    if n_proc > n_iter:
        print("ERROR: n_proc cannot be greater than n_iter.")
        sys.exit()
    
    # python3 ./Prepare_motif_level_analysis.py 
    #--input-table /mnt/doctorado/.../....tsv
    #--input-fasta /mnt/doctorado/.../....fasta
    #--path-out /mnt/doctorado/.../...
    #--n-iter 100
    #--n-proc 50
    
    """
    table_in = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/04-Families/High/intergenic/gen_ORIGINAL_no.tsv"
    fasta_in = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/01-LncRNAs/High/intergenic/LncRNAs.fasta"
    path_out = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Motif_level/Positional_conserved/02-Preparation/ORIGINAL/no/High/intergenic"  
    n_iter = 100
    processes = 50
    """
    
    ## Create a dictionary with the family IDs as key and a list with the LncRNAs 
    ## (family members) as value.
    Families_dic = FamiliesDict (table_in)
    
    ## Create a dictionary with the LncRNA IDs as key and its sequence as value.
    LncRNAs_and_sequences_dic = LncRNAsAndSequencesDict (fasta_in)
    
    ## Create a dictionary with the spe as key and lncRNA list as value.
    LncRNAs_lists_by_specie_dic = LncRNAsListsBySpecieDict (fasta_in)
    
    ## REAL: Create a fasta file by conserved family.
    print("\n-REAL")
    os.makedirs(path_out + "/real")
    CreateFastaFilesREAL (Families_dic, LncRNAs_and_sequences_dic, path_out)
    
    ## SIMULATION: Create random fasta files.
    print("\n-SIMULATION:")
    os.makedirs(path_out + "/simulations")
    iterations = [i for i in range(1, n_iter + 1)]
    dist = n_iter//n_proc
    processes = []
    pi = 0
    for i in range(n_proc):
        if i == n_proc-1:
            pf = n_iter
        else:
            pf = pi + dist
        processes.append(multiprocessing.Process(target = CreateFastaFilesSIMULATIONS , args=(iterations[pi:pf], Families_dic, LncRNAs_and_sequences_dic, LncRNAs_lists_by_specie_dic, path_out) )) 
        processes[i].start()
        print('\t-Process', i, 'launched.')
        pi = pf
    for p in processes:
        p.join()
    print('\nDONE.') 


# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()      
            
