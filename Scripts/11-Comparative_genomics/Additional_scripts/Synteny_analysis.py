#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- LOOKING FOR SYNTENIC LNCRNAS

This script is used to obtain lncRNAs that are positionally conserved. It is an 
adaptation of the code used in Pegueroles et al. (2019).

---------

Created on Fri Jan 27 12:00:00 2023

Last Modified:
    - Fri Jan 27 12:00:00 2023 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""


# MODULES

import sys
import argparse
import pandas as pd
import itertools


# FUNCTIONS

def RecodeList (spe, translist, orthodict, nonMatch):
    
    translist_new = []
    
    lncRNAs = 0
    matched = 0
    nonmatched = 0
    
    for x in translist:
        if x.startswith("MSTRG"):
            translist_new.append(x)
            lncRNAs += 1
        elif x in list(orthodict.keys()):
            ID = orthodict[x]
            translist_new.append(ID)
            matched += 1
        else:
            nonmatched += 1
            if nonMatch == 'yes':
                translist_new.append(x)
    
    print("\nNumber of matched geneIDs for " + spe + ": " + str(matched))
    print("Number of non-matched geneIDs for " + spe + ": " + str(nonmatched))
    print("Number of candidate lncRNAs for " + spe + ": " + str(lncRNAs))
    
    return translist_new


def SearchSyntenicLncRNAs (spe1, spe2, Translistsdict_mod, LncRNAs_idx, genesNearby, minOverlap, minSideOverlap):
    
    dict1 = DictOfClusters(LncRNAs_idx[spe1], Translistsdict_mod[spe1], genesNearby)
    dict2 = DictOfClusters(LncRNAs_idx[spe2], Translistsdict_mod[spe2], genesNearby)
    
    myHomologs = []
    homologFound = 'false'
    
    for key1, val1 in dict1.items():
        for key2, val2 in dict2.items():
            if len(set(dict1[key1]['all']).intersection(set(dict2[key2]['all']))) >= minOverlap:   
                if len(set(dict1[key1]['right']).intersection(set(dict2[key2]['right']))) >= minSideOverlap:
                    if len(set(dict1[key1]['left']).intersection(set(dict2[key2]['left']))) >= minSideOverlap: 
                        homologFound = 'true'
                if len(set(dict1[key1]['right']).intersection(set(dict2[key2]['left']))) >= minSideOverlap:
                    if len(set(dict1[key1]['left']).intersection(set(dict2[key2]['right']))) >= minSideOverlap:
                        homologFound = 'true'
                if homologFound == 'true':             
                        mytup = (key1, key2)
                        myHomologs.append(mytup)
                        homologFound = 'false'
    
    print("Number of tuples: " + str(len(myHomologs)))
    
    return myHomologs


def DictOfClusters (myidx, mylist, genesNearby):
    
    mydict = {}
    
    for idx in myidx:
        key = mylist[idx]
        val = {'left':[], 'right':[], 'all':[]}
        if not key in mydict:
            mydict[key] = val
        i = 1
        while i <= genesNearby:
            try:
                if idx - i >= 0:                 
                    mydict[key]['left'].append(mylist[idx - i])
                    mydict[key]['all'].append(mylist[idx - i])
                mydict[key]['right'].append(mylist[idx + i])                                   
                mydict[key]['all'].append(mylist[idx + i])
            except IndexError:
                pass
            i += 1
            
    return mydict
    

# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='Software V1', 
            description='''This program looks for syntenic lncRNAs.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument(
            "--path-lists", 
            type=str, 
            nargs=1,
            help="Path to the directory of gene and lncRNA lists."
            )
    parser.add_argument(
            "--orthotable", 
            type=str, 
            nargs=1,
            help="Orthotable file."
            )
    parser.add_argument(
            "--output", 
            type=str, 
            nargs=1,
            help="Output file."
            )
    parser.add_argument(
            "--genesNearby", 
            type=int, 
            nargs=1,
            help="When comparing the transcripts around a lncRNA with the transcripts \
                around another lncRNA in another specie, number of transcripts to be \
                evaluated on either side."
            )
    parser.add_argument(
            "--minOverlap", 
            type=int, 
            nargs=1,
            help="When comparing the transcripts around a lncRNA with the transcripts \
                around another lncRNA in another specie, minimum number of orthologs \
                genes that must overlap."
            )
    parser.add_argument(
            "--minSideOverlap", 
            type=int, 
            nargs=1,
            help="When comparing the transcripts around a lncRNA with the transcripts \
                around another lncRNA in another specie, minimum number of orthologs \
                genes that must overlap on each side."
            )
    parser.add_argument(
            "--nonMatch", 
            type=str, 
            nargs=1,
            help="When the lists of genes and lncRNAs are recoded using the orthotable, \
                you can keep the non-orthologs genes in the new lists or you can remove \
                them. (no/yes)."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        path_lists = args.path_lists[0]
        path_orthotable = args.orthotable[0]
        path_output = args.output[0]
        genesNearby = args.genesNearby[0]
        minOverlap = args.minOverlap[0]
        minSideOverlap = args.minSideOverlap[0]
        nonMatch = args.nonMatch[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # path_lists = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/01-LncRNAs_and_Genes/High/intergenic"
    # path_orthotable = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/02-Orthologs/Tables_orthologs_1_to_1/Orthotable_1_to_1_orthofinder.tsv"
    # path_output = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/03-Synteny_analysis/High/intergenic/output_synteny_table_ORIGINAL_no.tsv"
    # genesNearby = 3
    # minOverlap = 3
    # minSideOverlap = 1
    # nonMatch = 'no'
    
    ### PIPELINE
    ## Species.
    print("\nSTEP 1: Capture the species of the folder.")
    orthotable = pd.read_csv(path_orthotable, sep = '\t', header = 0)
    species = list(orthotable.columns)[1:]
    
    ## Dictionary with species as key and a subdictionary as value. This subdictionary 
    ## has genes orthologs 1:1 as key and orthogroup ids as value.
    print("\nSTEP 2: Load the table with the orthologous genes 1:1.")
    orthodict = {}
    for spe in species:
        subset = orthotable[['OG', spe]]
        subdict = dict(zip(subset[spe], subset['OG']))
        orthodict[spe] = subdict
    
    ## Dictionary with species as key and lists of genes and lncRNAs as value.
    print("\nSTEP 3: Load lists of genes and lncRNAs.")
    Translistsdict = {}
    for spe in species:
        Translistsdict[spe] = [line.strip() for line in open(path_lists + "/" + spe + "_ids_sorted.txt", 'r')]
    
    ## Recode lists of genes and lncRNAs. Change genes identified as orthologs 
    ## 1:1 by orthogroup id using orthodict. Dictionary with species as key and 
    ## lists of genes (orthologs 1:1 modified) and lncRNAs as value.
    print("\nSTEP 4: Recode lists of genes and lncRNAs using the orthotable.")
    Translistsdict_mod = {}
    for spe in species:
        Translistsdict_mod[spe] = RecodeList (spe, Translistsdict[spe], orthodict[spe], nonMatch)
    
    ## Dictionary with species as key and lists containing the respective lncRNA 
    ## positions in Translistsdict_mod lists.
    print("\nSTEP 5: Get the index of the lncRNAs in the recoded lists.")
    LncRNAs_idx = {}
    for spe in species:
        LncRNAs_idx[spe] = [i for i, trans in enumerate(Translistsdict_mod[spe]) if trans.startswith('MSTRG')]
    
    ## Compare species to look for syntenic lncRNAs and write the results in a file.
    print("\nSTEP 6: Compare species to look for syntenic lncRNAs.")
    combinations = [list(pair) for pair in itertools.combinations(species,2)]
    output = open(path_output, 'w')
    for comb in combinations:
        spe1 = comb[0]
        spe2 = comb[1]
        print("\nCombination: " + spe1 + "-" + spe2)
        for tup in SearchSyntenicLncRNAs (spe1, spe2, Translistsdict_mod, LncRNAs_idx, genesNearby, minOverlap, minSideOverlap):
            output.write(spe1 + "\t" + spe2 + "\t" + '\t'.join(tup) + "\n")
        print("Combination: " + spe2 + "-" + spe1)
        for tup in SearchSyntenicLncRNAs (spe2, spe1, Translistsdict_mod, LncRNAs_idx, genesNearby, minOverlap, minSideOverlap):
            output.write(spe2 + "\t" + spe1 + "\t" + '\t'.join(tup) + "\n")  
    output.close()
    
    
# CALL THE MAIN PROGRAM

if __name__ == '__main__':
    """
    Call the main program.
    """
    main()

