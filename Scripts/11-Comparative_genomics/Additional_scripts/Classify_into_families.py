#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- CLASSIFY LNCRNAS INTO FAMILIES.

This script will classify all the lncRNAs into families.

---------

Created on Fri Jul 22 12:10:00 2022

Last Modified:
    - Fri Jul 22 12:10:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import pandas as pd
import timeit
import multiprocessing


# FUNCTIONS

def PredictedLncRNAsDict (file):
    
    file_open = open(file, 'r')
    dic = {}
    for line in file_open:
        line = line.strip().split("\t")
        if line[0] not in dic.keys():
            dic[line[0]] = [line[1] + "-" + line[0]]
        else:
            dic[line[0]].append(line[1] + "-" + line[0])
    file_open.close()
    
    return dic


def ReciprocalHitsDict (file):
    
    file_open = open(file, 'r')
    dic = {}
    for line in file_open:
        line = line.strip().split("\t")
        key = line[2] + "-" + line[0]
        if key not in dic:
            dic[key] = [key, line[3] + "-" + line[1]]
        else:
            dic[key].append(line[3] + "-" + line[1])
    file_open.close()
    
    return dic


def FusionLists (list_of_lists, global_list_of_lists_collapsed_temp):
    
    claims = {}
    for lst in list_of_lists:
        glob = lst[:]
        for e in lst:
            if e not in claims:
                claims[e] = glob  # We are the first to claim e.
                continue
    
            # Ignore self-intersection.
            if claims[e] is glob:
                continue
    
            # We found someone else that we intersect with.
            # Take over ALL their territory, not just e.
            glob.extend(claims[e])
            for e2 in claims[e]:
                claims[e2] = glob  # at some point, e2 = e
    
    deduper = {}
    for v in claims.values():
        deduper[id(v)] = v
    
    list_of_lists_collapsed = [list(set(element)) for element in list(deduper.values())]
    global_list_of_lists_collapsed_temp += list_of_lists_collapsed
    
    return global_list_of_lists_collapsed_temp
    

def FindFamiliesFromReciprocalHitsDict (dic, num_process):
    
    ## Create a list of lists with the value of the previous dictionary, I mean
    ## with the families and then remove the repeated lists.
    ## From...
    ## [[LncRNA_1-cma, LncRNA_32-cla, LncRNA_74-car], [LncRNA_32-cla, LncRNA_1-cma, LncRNA_74-car], [LncRNA_74-car, LncRNA_32-cla, LncRNA_1-cma]]
    ## to...
    ## [[LncRNA_1-cma, LncRNA_32-cla, LncRNA_74-car]]
    list_of_lists = list(dic.values())
    
    number_lists = len(list_of_lists)
    distribution = number_lists//num_process
    
    processes = []
    manager = multiprocessing.Manager()
    global_list_of_lists_collapsed_temp = manager.list()
        
    pi = 0
    for i in range(num_process):
        if i == num_process-1:
            pf = number_lists
        else:
            pf = pi + distribution
        processes.append(multiprocessing.Process(target = FusionLists , args=(list_of_lists[pi:pf], global_list_of_lists_collapsed_temp) )) 
        processes[i].start()
        pi = pf
    for p in processes:
        p.join()
    
    global_list_of_lists_collapsed_sorted_temp = [sorted(val) for val in global_list_of_lists_collapsed_temp]
    global_list_of_lists_collapsed_sorted_uniq_temp = [list(t) for t in set(map(tuple, global_list_of_lists_collapsed_sorted_temp))]
    
    
    claims = {}
    for lst in global_list_of_lists_collapsed_sorted_uniq_temp:
        glob = lst[:]
        for e in lst:
            if e not in claims:
                claims[e] = glob  # We are the first to claim e.
                continue
    
            # Ignore self-intersection.
            if claims[e] is glob:
                continue
    
            # We found someone else that we intersect with.
            # Take over ALL their territory, not just e.
            glob.extend(claims[e])
            for e2 in claims[e]:
                claims[e2] = glob  # at some point, e2 = e
    
    deduper = {}
    for v in claims.values():
        deduper[id(v)] = v
    
    global_list_of_lists_collapsed_sorted = [sorted(list(set(element))) for element in list(deduper.values())]
    global_list_of_lists_collapsed_sorted_uniq = [list(t) for t in set(map(tuple, global_list_of_lists_collapsed_sorted))]
    
    ## Create a dictionary of detected families. For example:
    ## {Fam1: [LncRNA_1-cma, LncRNA_32-cla, LncRNA_74-car]}
    dictFam = {}
    i = 1
    for x in global_list_of_lists_collapsed_sorted_uniq:
        dictFam["Fam" + str(i)] = x
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
        
        fam_members_species = [member.split("-")[1] for member in fam]
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
        fam_members_species = [member.split("-")[1] for member in fam]
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
            prog='ClassifyLncRNAsIntoFamilies V1', 
            description='''This program classifies lncRNAs into families using the reciprocal \
                hits coming from blast.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--pred-lncRNAs", 
            type=str, 
            nargs=1,
            help="Absolute path of TSV table containing the ids of lncRNAs predicted by specie."
            )
    parser.add_argument(
            "--reciprocal-hits", 
            type=str, 
            nargs=1,
            help="Absolute path of TSV table containing the reciprocal blast hits between species."
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
            "--threads", 
            type=str, 
            nargs=1,
            help="Number of process or instances to parallelize."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        predicted_lncRNAs_by_specie = args.pred_lncRNAs[0]
        reciprocal_hits = args.reciprocal_hits[0]
        fam = args.fam[0]
        gen = args.gen[0]
        num_process = int(args.threads[0])
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Classify_into_families.py 
    #--pred-lncRNAs /mnt/doctorado/.../....tsv
    #--reciprocal-hits /mnt/doctorado/.../....tsv
    #--fam /mnt/doctorado/.../....tsv
    #--gen /mnt/doctorado/.../....tsv
    
    """
    predicted_lncRNAs_by_specie = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/05-Families/High/intergenic/ids_by_specie.tsv"
    reciprocal_hits = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/04-Reciprocal_hits/High/intergenic/Reciprocal_hits.tsv"
    fam = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/05-Families/High/intergenic/fam.tsv"
    gen = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Blastn/nr/05-Families/High/intergenic/gen.tsv"
    num_process = 10
    """   
    """
    predicted_lncRNAs_by_specie = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Synteny/nr/05-Families/PAIRWISE_SPECIES/High/intergenic/ids_by_specie.tsv"
    reciprocal_hits = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Synteny/nr/04-LncRNAs_inside_syntenic_blocks/PAIRWISE_SPECIES/CloudHits-intergenic-High.tsv"
    fam = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Synteny/nr/05-Families/PAIRWISE_SPECIES/High/intergenic/fam.tsv"
    gen = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Synteny/nr/05-Families/PAIRWISE_SPECIES/High/intergenic/gen.tsv"
    num_process = 10
    """
    
    start = timeit.default_timer()
    
    ## List of species.
    file_open = open(predicted_lncRNAs_by_specie, 'r')
    species = [line.strip().split("\t")[0] for line in file_open]
    species = sorted(list(set(species)))
    file_open.close()
    
    ## Create a dictionary of lists with the predicted lncRNAs (value) by specie
    ## (key).
    print("STEP 1/5")
    predicted_lncRNAs_by_specie_dict = PredictedLncRNAsDict (predicted_lncRNAs_by_specie)
    
    ## Create a dictionary of lists with the reciprocal blastn hits and rename the
    ## lncRNA ids adding the specie. For example:
    ## LncRNA_1-cma:[LncRNA_1-cma, LncRNA_32-cla]
    print("STEP 2/5")
    reciprocal_hits_dict = ReciprocalHitsDict (reciprocal_hits)
    print("Number of lncRNAs involved in reciprocal hits: " + str(len(reciprocal_hits_dict.keys())))
    
    ## Classify LncRNAs into families and add a family code. For example:
    ## LncRNA_1-cma:[LncRNA_1-cma, LncRNA_32-cla]
    ## LncRNA_32-cla:[LncRNA_32-cla, LncRNA_1-cma]
    ## LncRNA_74-car:[LncRNA_74-car, LncRNA_32-cla]
    ##
    ## LncRNA_1-cma:[LncRNA_1-cma, LncRNA_32-cla, LncRNA_74-car]
    ## LncRNA_32-cla:[LncRNA_32-cla, LncRNA_1-cma, LncRNA_74-car]
    ## LncRNA_74-car:[LncRNA_74-car, LncRNA_32-cla, LncRNA_1-cma]
    print("STEP 3/5")
    Families_dict = FindFamiliesFromReciprocalHitsDict (reciprocal_hits_dict, num_process)
    print("Number of families: " + str(len(Families_dict.keys())))
    
    ## Add individual lncRNAs in the Families_dict as a individual families.
    print("STEP 4/5")
    Families_dict_all = AddIndividualFamiliesToFamDict (Families_dict, species, predicted_lncRNAs_by_specie_dict)
    
    ## Generate the output files.
    print("STEP 5/5")
    FamiliesBinaryTable (Families_dict_all, species, fam)
    GenesBinaryTable (Families_dict_all, species, gen)
    
    stop = timeit.default_timer()

    print('Time: ', stop - start)  


# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()      
            
