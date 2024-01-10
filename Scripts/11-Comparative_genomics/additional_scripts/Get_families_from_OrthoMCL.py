#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- GET LNCRNA FAMILIES FROM ORTHOMCL.

This script write all lncRNA families coming from OrthoMCL in the correct format.

---------

Created on Mon Dec 05 19:00:00 2022

Last Modified:
    - Mon Dec 05 19:00:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import pandas as pd
import re


# FUNCTIONS

def PredictedLncRNAsDict (file):
    
    file_open = open(file, 'r')
    dic = {}
    for line in file_open:
        line = line.strip().split("\t")
        if line[0] not in dic.keys():
            dic[line[0]] = [line[1]]
        else:
            dic[line[0]].append(line[1])
    file_open.close()
    
    return dic
    

def CreateDictFromOrthoMCL (orthoMCL_out):
    
    dictFam = {}
    file_open = open(orthoMCL_out, 'r')
    i = 1
    for line in file_open:
        info = line.strip().split("\t")
        genes = [re.sub(r'\(G[0-9][0-9][0-9]\)', '', gen.split("|")[1]) for gen in info[1].strip().split(" ")]
        dictFam["Fam" + str(i)] = sorted(genes)
        i += 1
    
    return dictFam


def AddIndividualFamiliesToFamDict (dictFam, species, predicted_lncRNAs):
    
    i = len(dictFam.keys()) + 1
    
    all_lncRNAs_classified = [member for fam in dictFam.values() for member in fam]
    
    for spe in species:
        predicted = predicted_lncRNAs[spe]
        for lnc in predicted:
            if lnc not in all_lncRNAs_classified:
                key = "Fam" + str(i)
                dictFam[key] = [lnc]
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
            prog='GetLncRNAFamiliesOrthoMCL V1', 
            description='''This program gets lncRNA families coming from OrthoMCL.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument(
            "--pred-lncRNAs", 
            type=str, 
            nargs=1,
            help="Absolute path of TSV table containing the ids of lncRNAs predicted by specie."
            )
    parser.add_argument(
            "--orthomcl", 
            type=str, 
            nargs=1,
            help="Absolute path of table containing the lncRNA families coming from OrthoMCL."
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
        orthomcl = args.orthomcl[0]
        fam = args.fam[0]
        gen = args.gen[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Get_families_from_OrthoMCL.py 
    #--pred-lncRNAs /mnt/doctorado/.../....tsv
    #--orthomcl /mnt/doctorado/.../....out
    #--fam /mnt/doctorado/.../....tsv
    #--gen /mnt/doctorado/.../....tsv
    
    """
    predicted_lncRNAs_by_specie = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/OrthoMCL/05-Families/High/intergenic/ids_by_specie.tsv"
    orthomcl = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/OrthoMCL/04-OrthoMCL/High/intergenic/all_orthomcl.out"
    fam = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/OrthoMCL/05-Families/High/intergenic/fam.tsv"
    gen = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/OrthoMCL/05-Families/High/intergenic/gen.tsv"
    """
    ## List of species.
    file_open = open(predicted_lncRNAs_by_specie, 'r')
    species = [line.strip().split("\t")[0] for line in file_open]
    species = sorted(list(set(species)))
    file_open.close()
    
    ## Create a dictionary of lists with the predicted lncRNAs (value) by specie
    ## (key).
    print("STEP 1/4")
    predicted_lncRNAs_by_specie_dict = PredictedLncRNAsDict (predicted_lncRNAs_by_specie)
    
    ## Dictionary families.
    print("STEP 2/4")
    try:
        Families_dict = CreateDictFromOrthoMCL (orthomcl)
    except FileNotFoundError:
        print("ERROR: The file all_orthomcl.out doesn't exit. Maybe there is some blast_table without hits and then the OrthoMCL software failed. Check it.")
        sys.exit()
    
    ## Add individual lncRNAs in the Families_dict as a individual families.
    print("STEP 3/4")
    Families_dict_all = AddIndividualFamiliesToFamDict (Families_dict, species, predicted_lncRNAs_by_specie_dict)
    
    ## Generate the output files.
    print("STEP 4/4")
    FamiliesBinaryTable (Families_dict_all, species, fam)
    GenesBinaryTable (Families_dict_all, species, gen)
    

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()      
            
