#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- EXTRACT 1:1 ORTHOLOGS ACROSS ALL SPECIES FROM ORTHOFINDER RESULTS.

---------

Created on Wed Jan 25 17:00:00 2023

Last Modified:
    - Wed Jan 25 17:00:00 2023 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import pandas as pd
import os


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='Software V1', 
            description='''This program extracts 1:1 orthologs across all species from orthofinder results.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument(
            "--path-orthofinder", 
            type=str, 
            nargs=1,
            help="Path to the OrthoFinder directory."
            )
    parser.add_argument(
            "--output", 
            type=str, 
            nargs=1,
            help="Output table (1:1 orthologs table)."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        path_OF = args.path_orthofinder[0]
        output = args.output[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Extract_1_to_1_orthologs_across_all-Approach_2-Orthofinder.py 
    #--path-orthofinder /storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/02-Orthologs/Inference/Orthofinder
    #--output /storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/02-Orthologs/Tables_orthologs_1_to_1/Orthotable_1_to_1_orthofinder.tsv
    
    """
    path_OF = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/02-Orthologs/Inference/Orthofinder"
    output = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/02-Orthologs/Tables_orthologs_1_to_1/Orthotable_1_to_1_orthofinder.tsv"
    """
    
    ##### ORTHOFINDER
    
    ## Paths.
    path_OF_genecounts = path_OF + "/All/" + os.listdir(path_OF + "/All")[0] +  "/Orthogroups/Orthogroups.GeneCount.tsv"
    path_OF_genes_assigned = path_OF + "/All/" + os.listdir(path_OF + "/All")[0] + "/Orthogroups/Orthogroups.tsv"
    
    ## Tables.
    OF_genecounts = pd.read_csv(path_OF_genecounts, sep='\t', header=0)
    OF_genes_assigned = pd.read_csv(path_OF_genes_assigned, sep='\t', header=0)
    
    ## Species
    cols = list(OF_genecounts.columns)
    species = cols[1:-1]
    
    ## Keep only 1:1 orthogroups.
    condition = " == 1 & ".join(species) + " == 1"
    OF_genecounts_filt = OF_genecounts.query(condition)
    
    ## Create a list with 1:1 orthogroup ids.
    OF_OG_filt_list = list(OF_genecounts_filt["Orthogroup"])
    
    ## Create a 1:1 orthologs dictionary with orthogroup id as key and a subdictionary 
    ## as value. This subdictionary has species as key and the gene (ortholog 1:1) as 
    ## value.
    OF_OG_filt_dict = {}
    for OG in OF_OG_filt_list:
        OF_OG_filt_dict_sub = {}
        subset = OF_genes_assigned[(OF_genes_assigned["Orthogroup"] == OG)]
        for spe in species:
            OF_OG_filt_dict_sub[spe] = list(subset[spe])[0]
        OF_OG_filt_dict[OG] = OF_OG_filt_dict_sub
    
    ## Save the results.
    orthotable = open(output, 'w')
    orthotable.write("OG" + "\t" + "\t".join(species) + "\n")
    for OG in list(OF_OG_filt_dict.keys()):
        orthotable.write(OG + "\t" + "\t".join([OF_OG_filt_dict[OG][spe] for spe in species]) + "\n")
    orthotable.close()
    

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
