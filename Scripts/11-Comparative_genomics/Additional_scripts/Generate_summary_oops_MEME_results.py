#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- GENERATE A SUMMARY OF OOPS MEME RESULTS

This script creates a summary file about oops MEME results.

---------

Created on Sun Feb 12 21:00:00 2023

Last Modified:
    - Sun Feb 12 21:00:00 2023 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""


# MODULES

import sys
import argparse
import os
import multiprocessing
import pandas as pd
import numpy as np


# FUNCTIONS

def GenerateSummaryDict (families_list, path_meme, path_fastas, type_):
    
    Summary_list = []
    for fam in families_list:
        meme_res = open(path_meme + "/" + fam + "/meme.txt", "r")
        lines = meme_res.readlines()
        starts = [i for i in range(len(lines)) if lines[i].strip().startswith("MOTIF")]
        ends = [i for i in range(len(lines)) if lines[i].strip().startswith("Time")]
        ids = [i for i in range(1, len(starts) + 1)]
        if len(ids) > 0:
            MOTIF_SECTIONS = {"MOTIF_" + str(ids[i]): lines[starts[i]:ends[i] + 1] for i in range(len(ids))}
            for mot in list(MOTIF_SECTIONS.keys()):
                section = MOTIF_SECTIONS[mot]
                motif_list = []
                for j in range(len(section)):
                    line = section[j].strip()
                    if line.startswith("MOTIF"):
                        if float(line.split("=")[-1]) <= 0.05:
                            motif = line.split()[1]
                            w = int(line.split()[5])
                            sites = int(line.split()[8])
                            evalue = float(line.split()[14])
                        else:
                            break
                    elif "Multilevel" in line:
                        cons_seq = line.split()[-1]
                    elif "sites sorted by position p-value" in line:
                        ini = j + 4
                        fin = j + 4 + sites
                        for z in range(ini, fin):
                            lncRNA = section[z].split()[0]
                            start = int(section[z].split()[1])
                            pvalue = float(section[z].split()[2])
                            seq = section[z].split()[4]
                            motif_list.append([fam.split("_")[0], lncRNA, motif, cons_seq, evalue, seq, w, start, pvalue, sites])
                    elif "regular expression" in line:
                        motif_list_mod = []
                        reg_exp = section[j + 2].split()[0]
                        for x in motif_list:
                            motif_list_mod.append(x[:3] + [reg_exp] + x[3:])
                Summary_list = Summary_list + motif_list_mod
        else:
            motif_list_mod = []
            fasta = open(path_fastas + "/" + fam + ".fasta", "r")
            lncRNAs = [line[1:].strip() for line in fasta if line[0] == ">"]
            for lnc in lncRNAs:
                motif_list_mod.append([fam.split("_")[0], lnc, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
            Summary_list = Summary_list + motif_list_mod
        
    Summary_tab = pd.DataFrame(Summary_list)
    Summary_tab.columns = ["Family", "LncRNA", "Meme_Motif.Identifier", "Meme_Motif.Regex", "Meme_Motif.Conserved", "Meme_Motif.E.value", "Meme_Motif.LncRNA", "Meme_Motif.Width", "Meme_Motif.Start", "Meme_Motif.P.value", "Meme_Motif.Sites"]
    Summary_tab.insert(0, "Type", [type_]*Summary_tab.shape[0])
    
    return Summary_tab

def CallIterations (iterations, path_meme, path_fastas, simulation_dict):
    
    for iter_ in iterations:
        families_list = os.listdir(path_meme + "/" + iter_)
        families_list_sorted = ["Fam" + str(idx) + "_random" for idx in sorted([int(fam.split("_")[0].replace("Fam", "")) for fam in families_list])]
        Summary_tab = GenerateSummaryDict (families_list_sorted, path_meme + "/" + iter_, path_fastas + "/" + iter_, "SIMULATION_" + str(iter_.split("_")[-1]))
        simulation_dict[iter_] = Summary_tab
        
    return simulation_dict


# MAIN PROGRAM

def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='SummaryOopsMEMEResults V1', 
            description='''This script creates a summary about oops MEME results.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--path-meme", 
            type=str, 
            nargs=1,
            help="Absolute oops meme results path."
            )
    parser.add_argument(
            "--path-fastas", 
            type=str, 
            nargs=1,
            help="Absolute fastas path."
            )
    parser.add_argument(
            "--path-out", 
            type=str, 
            nargs=1,
            help="Absolute output path."
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
        path_meme = args.path_meme[0]
        path_fastas = args.path_fastas[0]
        path_out = args.path_out[0]
        n_proc = args.n_proc[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # path_meme = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Motif_level/nr/03-MotifFinder/MEME/ORIGINAL/no/Low/intergenic/oops/6-15"
    # path_fastas = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Motif_level/nr/02-Preparation/Low/intergenic"
    # path_out = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Motif_level/nr/05-Summary/MEME/ORIGINAL/no/Low/intergenic/oops/6-15"  
    # n_proc = 25
    
    ### PIPELINE
    ## MEME INFO: REAL
    families_list = os.listdir(path_meme + "/real")
    families_list_sorted = ["Fam" + str(idx) + "_real" for idx in sorted([int(fam.split("_")[0].replace("Fam", "")) for fam in families_list])]
    Summary_REAL_tab = GenerateSummaryDict (families_list_sorted, path_meme + "/real", path_fastas + "/real", "REAL")
    
    ## MEME INFO: SIMULATIONS
    iterations = os.listdir(path_meme + "/simulations")
    n_iter = len(iterations)
    dist = n_iter//n_proc
    processes = []
    manager = multiprocessing.Manager()
    simulation_dict = manager.dict()
    pi = 0
    for i in range(n_proc):
        if i == n_proc-1:
            pf = n_iter
        else:
            pf = pi + dist
        processes.append(multiprocessing.Process(target = CallIterations , args=(iterations[pi:pf], path_meme + "/simulations", path_fastas + "/simulations", simulation_dict) )) 
        processes[i].start()
        #print('\t-Process', i, 'launched.')
        pi = pf
    for p in processes:
        p.join()
    
    Summary_SIMULATION_tab = pd.DataFrame()
    simulations_list = list(simulation_dict.keys())
    simulations_list_sorted = ["iter_" + str(idx) for idx in sorted([int(iter_.split("_")[1]) for iter_ in simulations_list])]
    for key in simulations_list_sorted:
        tab = simulation_dict[key]
        Summary_SIMULATION_tab = pd.concat([Summary_SIMULATION_tab, tab])
    
    ## Join REAL and SIMULATIONS tables and save.
    Summary_tab = pd.concat([Summary_REAL_tab, Summary_SIMULATION_tab])
    Summary_tab.to_csv(path_out + "/MEME-SUMMARY.tsv", sep = '\t', index = False, header = True)


# CALL THE MAIN PROGRAM

if __name__ == '__main__':
    """
    Call the main program.
    """
    main()      
            
