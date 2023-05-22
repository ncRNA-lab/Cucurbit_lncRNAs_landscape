#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- EXTRACT THE 1:1 ORTHOLOGS BETWEEN TWO SPECIES (PAIRWISE) FROM THE RESULTS COMING 
FROM THE ORTHOFINDER AND INPARANOID SOFTWARES GENERATING A TXT FILE WITH THE ORTHOLOGS
LIST AND THEIR ORTHOGROUPS.

---------

Created on Thur Oct 20 10:00:00 2022

Last Modified:
    - FThur Oct 20 10:00:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import os
import pandas as pd
import numpy as np


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='ExtractOrthologs V1', 
            description='''This program extract 1:1 orthologs between two species º
            from the results coming from the orthofinder and inparanoid softwares.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--spe1", 
            type=str, 
            nargs=1,
            help="Specie one."
            )
    parser.add_argument(
            "--spe2", 
            type=str, 
            nargs=1,
            help="Specie two."
            )
    parser.add_argument(
            "--path-in", 
            type=str, 
            nargs=1,
            help="Path to the Orthologs directory."
            )
    parser.add_argument(
            "--path-out", 
            type=str, 
            nargs=1,
            help="Path to the results."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        spe1 = args.spe1[0]
        spe2 = args.spe2[0]
        path_in = args.path_in[0]
        path_out = args.path_out[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Extract_1_to_1_orthologs_pairwise.py 
    #--spe1 car
    #--spe2 cla
    #--path-in /mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/01-Orthologs
    #--path-out /mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/02-Adhore/PAIRWISE_SPECIES/Individual/car-cla
    
    """
    spe1 = "car"
    spe2 = "cla"
    path_in = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/01-Orthologs" 
    path_out = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/02-Adhore/PAIRWISE_SPECIES/Individual/car-cla" 
    """
    
    ##### PAIRWISE 1:1 ORTHOLOGS
    
    ## Paths and Files.
    OF_genecounts = path_in + "/Orthofinder/" + spe1 + "-" + spe2 + "/" + os.listdir(path_in + "/Orthofinder/" + spe1 + "-" + spe2)[0] +  "/Orthogroups/Orthogroups.GeneCount.tsv"
    OF_genes_assigned = path_in + "/Orthofinder/" + spe1 + "-" + spe2 + "/" + os.listdir(path_in + "/Orthofinder/"+ spe1 + "-" + spe2)[0] + "/Orthogroups/Orthogroups.tsv"
    OF_genes_unassigned = path_in + "/Orthofinder/" + spe1 + "-" + spe2 + "/" + os.listdir(path_in + "/Orthofinder/"+ spe1 + "-" + spe2)[0] + "/Orthogroups/Orthogroups_UnassignedGenes.tsv"
    OG_counts = pd.read_csv(OF_genecounts, sep='\t', header=0)
    OG_genes_assigned = pd.read_csv(OF_genes_assigned, sep='\t', header=0)
    OG_genes_unassigned = pd.read_csv(OF_genes_unassigned, sep='\t', header=0)
    IP_genes_assigned = open(path_in + "/Inparanoid/" + spe1 + "-" + spe2 + "/table." + spe1 + ".fa-" + spe2 + ".fa", 'r')
    
    ## Create two lists from OrthoFinder with orthogroup ids, according to the following conditions:
    ## List_1: 1:1.
    ## List_2: many:1, many:0, 1:many, 0:many and many:many
    subset_1 = OG_counts[(OG_counts[spe1] == 1) & (OG_counts[spe2] == 1)]
    LIST_OG_1 = list(subset_1["Orthogroup"])
    LIST_OG_2 = [OG for OG in list(OG_counts["Orthogroup"]) if OG not in LIST_OG_1]
    
    ### STEP 1
    ## Create a 1:1 orthologs dictionary with orthogroup id as key and the list of the gene ids as a value.
    DICT_Genes_1 = {}
    for OG in LIST_OG_1:
        list_ = []
        subset = OG_genes_assigned[(OG_genes_assigned["Orthogroup"] == OG)]
        list_.append(list(subset[spe1])[0])
        list_.append(list(subset[spe2])[0])
        DICT_Genes_1[OG] = list_
    
    ## Create a list of lists with the orthologs 1:1 of each specie (PAIRED) according to inparanoid.
    LIST_inparanoid = [] 
    for line in IP_genes_assigned:
        line = line.strip().split("\t")
        if line[0] == "OrtoID":
            continue
        else:
            line = line[2:]
            list_temp = []
            for x in line:
                x_list = x.strip().split(" ")
                if len(x_list) > 2:
                    list_temp = []
                    break
                else:
                    list_temp.append(x_list[0])
            if list_temp:
                LIST_inparanoid.append(list_temp)
    
    ## Validate the results coming from orthofinder (DICT_Genes_1) using the inparanoid results (LIST_inparanoid).
    DICT_Genes_1_checked = {}
    i = 0
    for OG in list(DICT_Genes_1.keys()):
        values_conf = ['Yes' for list_genes in LIST_inparanoid if set(list_genes) == set(DICT_Genes_1[OG])]
        if values_conf:
            DICT_Genes_1_checked[OG] = DICT_Genes_1[OG]
        else:
            i += 1
            LIST_OG_2.append(OG)
    print(spe1 + "-" + spe2 + ": " + str(i) + " OGs removed | " + str(len(DICT_Genes_1_checked.keys())) + " OGs kept")            
    
    
    ### STEP 2
    ## Create a many:1, many:0, 1:many, 0:many and many:many dictionary with orthogroup id as key and 
    ## the list of the gene ids as a value.
    DICT_Genes_2 = {}
    i = 1
    for OG in LIST_OG_2:
        list_ = []
        subset = OG_genes_assigned[(OG_genes_assigned["Orthogroup"] == OG)]
        list_ += list(subset[spe1])
        list_ += list(subset[spe2])
        for x in list_:
            if x is np.nan:
                continue
            elif ", " in x:
                x_list = x.split(", ")
                for y in x_list:
                    new_OG = 'OG' + str(i) + '_new'
                    DICT_Genes_2[new_OG] = [y]
                    i += 1
            elif ", " not in x:
                new_OG = 'OG' + str(i) + '_new'
                DICT_Genes_2[new_OG] = [x]
                i += 1
                
              
    ### STEP 3
    ## Create am unassigned gene dictionary with orthogroup id as key and the list of the gene ids as a value.
    subset_1 = OG_genes_unassigned[~OG_genes_unassigned[spe1].isnull()]
    IDs_1 = list(subset_1["Orthogroup"])
    Genes_1 = list(subset_1[spe1])
    subset_2 = OG_genes_unassigned[~OG_genes_unassigned[spe2].isnull()]
    IDs_2 = list(subset_2["Orthogroup"])
    Genes_2 = list(subset_2[spe2])
    
    DICT_Genes_3 = {}
    for i in range(len(IDs_1)):
        DICT_Genes_3[IDs_1[i]] = [Genes_1[i]]
    for i in range(len(IDs_2)):
        DICT_Genes_3[IDs_2[i]] = [Genes_2[i]]
        
    ### STEP 4
    ## Write the results in the format required by iadhore.
    file_out = open(path_out + "/Orthogroups.tsv", 'w')
    for OG in list(DICT_Genes_1_checked.keys()):
        for gene in DICT_Genes_1_checked[OG]:
            file_out.write(gene + "\t" + OG + "\n")
    for OG in list(DICT_Genes_2.keys()):
        for gene in DICT_Genes_2[OG]:
            file_out.write(gene + "\t" + OG + "\n")
    for OG in list(DICT_Genes_3.keys()):
        for gene in DICT_Genes_3[OG]:
            file_out.write(gene + "\t" + OG + "\n")
    file_out.close()
    

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
