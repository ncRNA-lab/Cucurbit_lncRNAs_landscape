#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- GET LNCRNA FAMILIES FROM ORTHOFINDER.

This script writes all lncRNA families coming from OrthoFinder in the correct 
format.

---------

Created on Wed Dec 07 14:00:00 2022

Last Modified:
    - Wed Dec 07 14:00:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""


# MODULES

import sys
import argparse
import pandas as pd


# FUNCTIONS

def CreateDictFromOrthoFinder (orthofinder):
    
    dictFam = {}
    file_open = open(orthofinder, 'r')
    i = 1
    for line in file_open:
        info = line.strip().split(": ")
        genes = info[1].split(" ")
        dictFam["Fam" + str(i)] = sorted(genes)
        i += 1
    
    return dictFam

    
def FamiliesBinaryTable (dictFam, species, famfile):
    
    BinTab = []
    for key in dictFam.keys():
        fam = dictFam[key]
        famline = [key]
        fam_members_species = [member.rsplit("-", 1)[1] for member in fam]
        for spe in species:
            if spe in fam_members_species:
                famline.append(1)
            else:
                famline.append(0)
        BinTab.append(famline)
    
    df = pd.DataFrame(BinTab)
    columnname = ["Family"] + species
    df.columns = columnname
    df.to_csv(famfile, sep = '\t', index = False)
    

def GenesBinaryTable (dictFam, species, genfile):
    
    BinTab = []
    for key in dictFam.keys():
        fam = dictFam[key]
        fam_members_species = [member.rsplit("-", 1)[1] for member in fam]
        for member in fam:
            memberline = [key, member]
            for spe in species:
                if spe in fam_members_species:
                    memberline.append(1)
                else:
                    memberline.append(0)
            BinTab.append(memberline)
    
    df = pd.DataFrame(BinTab)
    columnname = ["Family", "Member"] + species
    df.columns = columnname
    df.to_csv(genfile, sep = '\t', index = False)
    

# MAIN PROGRAM

def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='GetLncRNAFamiliesOrthoFinder V1', 
            description='''This program gets lncRNA families coming from OrthoFinder.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument(
            "--pred-lncRNAs", 
            type=str, 
            nargs=1,
            help="Absolute path of TSV table containing the ids of lncRNAs predicted by specie."
            )
    parser.add_argument(
            "--orthofinder", 
            type=str, 
            nargs=1,
            help="Absolute path of table containing the lncRNA families coming from OrthoFinder."
            )
    parser.add_argument(
            "--fam", 
            type=str, 
            nargs=1,
            help="Absolute path of binary table of families."
            )
    parser.add_argument(
            "--gen", 
            type=str, 
            nargs=1,
            help="Absolute path of binary table of lncRNAs."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        predicted_lncRNAs_by_specie = args.pred_lncRNAs[0]
        orthofinder = args.orthofinder[0]
        fam = args.fam[0]
        gen = args.gen[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # predicted_lncRNAs_by_specie = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Sequence_level/nr/05-Families/High/intergenic/ids_by_specie.tsv"
    # orthofinder = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Sequence_level/nr/04-OrthoFinder/High/intergenic/Clustering/Results_Dec07/Orthogroups/Orthogroups.txt"
    # fam = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Sequence_level/nr/05-Families/High/intergenic/fam.tsv"
    # gen = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Sequence_level/nr/05-Families/High/intergenic/gen.tsv"
    
    ### PIPELINE
    ## List of species.
    file_open = open(predicted_lncRNAs_by_specie, 'r')
    species = [line.strip().split("\t")[0] for line in file_open]
    species = sorted(list(set(species)))
    file_open.close()
    
    ## Dictionary families.
    print("STEP 1/2")
    try:
        Families_dict_all = CreateDictFromOrthoFinder (orthofinder)
    except FileNotFoundError:
        print("ERROR: The file Orthogroups.txt doesn't exit. Check it.")
        sys.exit()
    
    ## Generate the output files.
    print("STEP 2/2")
    FamiliesBinaryTable (Families_dict_all, species, fam)
    GenesBinaryTable (Families_dict_all, species, gen)
    

# CALL THE MAIN PROGRAM

if __name__ == '__main__':
    """
    Call the main program.
    """
    main()      
            
