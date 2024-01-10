#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- GENERATE A SUMMARY OF GOMO RESULTS

This script creates a summary file about GOMO results.

---------

Created on Wed Feb 22 12:00:00 2023

Last Modified:
    - Wed Feb 22 12:00:00 2023 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import os
import pandas as pd
import multiprocessing


# FUNCTIONS

def GenerateSummaryTable (families_list, path_gomo, type_):
    
    Summary_tab = pd.DataFrame()
    for fam in families_list:
        gomo = [line.strip().split() for line in open(path_gomo + "/" + fam + "/gomo_out/gomo.tsv", 'r') if line[0] != "#" and len(line.strip().split()) == 5]
        if len(gomo) > 1:
            gomo_res = gomo[1:]
            tab = pd.DataFrame(gomo_res)
            tab.columns = ["Gomo_Motif.Identifier", "Gomo_GO.Term.Identifier", "Gomo_Score", "Gomo_P.value", "Gomo_Q.value"]
            tab.insert(0, "Family", [fam.split("_")[0]]*tab.shape[0])
            tab.insert(0, "Type", [type_]*tab.shape[0])
            Summary_tab = pd.concat([Summary_tab, tab])
    
    return Summary_tab

def CallIterations (iterations, path_gomo, simulation_dict):
    
    for iter_ in iterations:
        families_list = os.listdir(path_gomo + "/" + iter_)
        families_list_sorted = ["Fam" + str(idx) + "_random" for idx in sorted([int(fam.split("_")[0].replace("Fam", "")) for fam in families_list])]
        Summary_tab = GenerateSummaryTable (families_list_sorted, path_gomo + "/" + iter_, "SIMULATION_" + str(iter_.split("_")[-1]))
        simulation_dict[iter_] = Summary_tab
        
    return simulation_dict


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='SummaryGOMOResults V1', 
            description='''This script creates a summary about GOMO results.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--path-gomo", 
            type=str, 
            nargs=1,
            help="Absolute gomo results path."
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
        path_gomo = args.path_gomo[0]
        path_out = args.path_out[0]
        n_proc = args.n_proc[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Generate_summary_GOMO_results.py
	#--path-gomo $DIR_A
	#--path-out $DIR_B
    #--n-proc $SLURM_CPUS_PER_TASK
	
    
    """
    path_gomo = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Motif_level/nr/Positional_conserved/04-MotifEnrichment/GOMO/ORIGINAL/no/Low/intergenic/oops/6-15"
    path_out = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Motif_level/nr/Positional_conserved/05-Summary/GOMO/ORIGINAL/no/Low/intergenic/oops/6-15"
    n_proc = 25
    """
    
    ### GOMO INFO: REAL
    families_list = os.listdir(path_gomo + "/real")
    families_list_sorted = ["Fam" + str(idx) + "_real" for idx in sorted([int(fam.split("_")[0].replace("Fam", "")) for fam in families_list])]
    Summary_REAL_tab = GenerateSummaryTable (families_list_sorted, path_gomo + "/real", "REAL")
    
    ### GOMO INFO: SIMULATIONS
    iterations = os.listdir(path_gomo + "/simulations")
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
        processes.append(multiprocessing.Process(target = CallIterations , args=(iterations[pi:pf], path_gomo + "/simulations", simulation_dict) )) 
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
    Summary_tab.to_csv(path_out + "/GOMO-SUMMARY.tsv", sep = '\t', index = False, header = True)


# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()      
            
