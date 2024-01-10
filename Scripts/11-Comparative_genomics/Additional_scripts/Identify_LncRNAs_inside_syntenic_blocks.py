#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- EXTRACT THE LNCRNAS CONTAINED INSIDE THE SYNTENIC BLOCKS PREDICTED BY ADHORE
IN A CLOUD MODE.

---------

Created on Thur Nov 24 10:00:00 2022

Last Modified:
    - Thur Nov 24 10:00:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import pandas as pd
import numpy as np


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='LncRNAsIdentificationFromSyntenicBlocks V1', 
            description='''Extract the LncRNAs contained inside the syntenic blocks \
            predicted by adhore in cloud mode''', 
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
            "--path1", 
            type=str, 
            nargs=1,
            help="Path to the predict lncRNAs."
            )
    parser.add_argument(
            "--path2", 
            type=str, 
            nargs=1,
            help="Path to the syntenic blocks."
            )
    parser.add_argument(
            "--path3", 
            type=str, 
            nargs=1,
            help="Path to the LncRNAs got from the syntenic blocks."
            )
    parser.add_argument(
            "--flag", 
            type=str, 
            nargs=1,
            help="redundant (r) or non-redundant (nr)."
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
        path_1 = args.path1[0]
        path_2 = args.path2[0]
        path_3 = args.path3[0]
        flag = args.flag[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Identify_LncRNAs_inside_syntenic_blocks.py 
    #--spe1 cla
    #--spe2 lsi
    #--path1 /mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs
    #--path2 /mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Synteny/nr/03-Adhore/PAIRWISE_SPECIES/car-csa
    #--path3 /mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Synteny/nr/04-LncRNAs_inside_syntenic_blocks/PAIRWISE_SPECIES/car-csa
    #--flag nr
    
    #spe1 = "cme"
    #spe2 = "csa"
    #path_1 = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs" 
    #path_2 = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Synteny/nr/03-Adhore/PAIRWISE_SPECIES/cme-csa" 
    #path_3 = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Synteny/nr/04-LncRNAs_inside_syntenic_blocks/PAIRWISE_SPECIES/cme-csa" 
    #flag = "nr"
    
    ## Load tables: Genes present in each genome.
    db_genes_spe1 = pd.read_csv(path_1 + "/" + spe1 + "/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv", sep='\t', header=0)
    db_genes_spe1 = db_genes_spe1.loc[:, ['ID_transcript', 'Chr', 'Start', 'End']]
    db_genes_spe1['Spe'] = spe1
    db_genes_spe2 = pd.read_csv(path_1 + "/" + spe2 + "/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv", sep='\t', header=0)
    db_genes_spe2 = db_genes_spe2.loc[:, ['ID_transcript', 'Chr', 'Start', 'End']]
    db_genes_spe2['Spe'] = spe2
    db_genes = pd.concat([db_genes_spe1, db_genes_spe2])
    db_genes = db_genes.astype({"Start":"int", "End":"int"})
    
    ## Load tables: LncRNAs present in each genome.
    if flag == 'nr':
        db_LncRNAs_spe1 = pd.read_csv(path_1 + "/" + spe1 + "/STEP-FINAL/Database/Database_LncRNAs_NR.tsv", sep='\t', header=0)
    else:
        db_LncRNAs_spe1 = pd.read_csv(path_1 + "/" + spe1 + "/STEP-FINAL/Database/Database_LncRNAs.tsv", sep='\t', header=0)
    db_LncRNAs_spe1 = db_LncRNAs_spe1.loc[:, ['ID_transcript', 'Chr', 'Start', 'End', 'Class_code', 'Significance_level']]
    if 'Significance_level' in list(db_LncRNAs_spe1.columns):
        db_LncRNAs_spe1 = db_LncRNAs_spe1.rename(columns={'Significance_level': 'Confidence_level'})
    db_LncRNAs_spe1['Spe'] = spe1
    if flag == 'nr':
        db_LncRNAs_spe2 = pd.read_csv(path_1 + "/" + spe2 + "/STEP-FINAL/Database/Database_LncRNAs_NR.tsv", sep='\t', header=0)
    else:
        db_LncRNAs_spe2 = pd.read_csv(path_1 + "/" + spe2 + "/STEP-FINAL/Database/Database_LncRNAs.tsv", sep='\t', header=0)
    db_LncRNAs_spe2 = db_LncRNAs_spe2.loc[:, ['ID_transcript', 'Chr', 'Start', 'End', 'Class_code', 'Significance_level']]
    if 'Significance_level' in list(db_LncRNAs_spe2.columns):
        db_LncRNAs_spe2 = db_LncRNAs_spe2.rename(columns={'Significance_level': 'Confidence_level'})
    db_LncRNAs_spe2['Spe'] = spe2
    db_LncRNAs = pd.concat([db_LncRNAs_spe1, db_LncRNAs_spe2])
    db_LncRNAs = db_LncRNAs.astype({"Start":"int", "End":"int"})
    
    ## Load tables: Adhore results. I mean, tables with the syntenic regions predicted between both species.
    clouds = pd.read_csv(path_2 + "/output_cloud/clouds.txt", sep='\t', header=0)
    cloudsAP = pd.read_csv(path_2 + "/output_cloud/cloudAP.txt", sep='\t', header=0)
    clouds_id = list(set(clouds["id"]))
    
    ## Additional variables
    studies = ["intergenic", "antisense", "intronic", "sense", "ALL"]
    studies_class_codes = [['u'], ['x'], ['i'], ['o', 'e'], ['u', 'x', 'i', 'o', 'e']]
    studies_factor = pd.Series(["intergenic", "antisense", "intronic", "sense", "ALL"])
    class_codes_factor = pd.Series(['u', 'x', 'i', 'o', 'e'])
    confidences = ['High', 'Medium', 'Low']
    confidences_factor = pd.Series(['High', 'Medium', 'Low'])
    
    ## Identify LncRNAs inside the syntenic blocks by study, confidence level and synteny block.
    dic_studies = {}
    dic_hits_A = {}
    TAB_ALL = pd.DataFrame()
    
    for i in range(len(studies)):
        study = studies[i]
        study_class_codes = studies_class_codes[i]
        dic_studies[study] = study_class_codes
        dic_hits_B = {}
        
        for confidence in confidences:
            print("Analizing: Study --> " + study + ", Confidence --> " + confidence)
            TAB = pd.DataFrame()
            list_hits = []
            
            for ID in clouds_id:
                # Get the syntenic block subset according to the ID.
                subset_clouds = clouds[(clouds["id"] == ID)]
                subset_cloudsAP = cloudsAP[(cloudsAP["CloudID"] == ID)]
                
                # Determine the genomes.
                genome_x = list(subset_clouds["genome_x"])[0]
                genome_y = list(subset_clouds["genome_y"])[0]
                
                # Determine the first and last ortholog gene in the syntnic block of each genome.
                gene_min_x = list(subset_cloudsAP[subset_cloudsAP.coord_x == subset_cloudsAP.coord_x.min()]["gene_x"])[0]
                gene_max_x = list(subset_cloudsAP[subset_cloudsAP.coord_x == subset_cloudsAP.coord_x.max()]["gene_x"])[0]
                gene_min_y = list(subset_cloudsAP[subset_cloudsAP.coord_y == subset_cloudsAP.coord_y.min()]["gene_y"])[0]
                gene_max_y = list(subset_cloudsAP[subset_cloudsAP.coord_y == subset_cloudsAP.coord_y.max()]["gene_y"])[0]
                
                # Determine the first and last nucleotide position and chromosome of the syntenic block in each genome.              
                df_gene_min_x = db_genes[(db_genes["Spe"] == genome_x) & (db_genes["ID_transcript"] == gene_min_x)]
                min_pos_x = list(df_gene_min_x["Start"])[0]
                min_chr_x = list(df_gene_min_x["Chr"])[0]
                df_gene_max_x = db_genes[(db_genes["Spe"] == genome_x) & (db_genes["ID_transcript"] == gene_max_x)]
                max_pos_x = list(df_gene_max_x["End"])[0]
                max_chr_x = list(df_gene_max_x["Chr"])[0]
                df_gene_min_y = db_genes[(db_genes["Spe"] == genome_y) & (db_genes["ID_transcript"] == gene_min_y)]
                min_pos_y = list(df_gene_min_y["Start"])[0]
                min_chr_y = list(df_gene_min_y["Chr"])[0]
                df_gene_max_y = db_genes[(db_genes["Spe"] == genome_y) & (db_genes["ID_transcript"] == gene_max_y)]
                max_pos_y = list(df_gene_max_y["End"])[0]
                max_chr_y = list(df_gene_max_y["Chr"])[0]
                
                # Check if the there is some error. For example, the syntenic block is the fusion of different genomes.
                if min_chr_x != max_chr_x or min_chr_y != max_chr_y:
                    sys.exit()
                
                # Get the LncRNAs subset according to the study and confidence level.
                d = db_LncRNAs[(db_LncRNAs["Class_code"].isin(study_class_codes)) & (db_LncRNAs["Confidence_level"] == confidence)]
                # Determine the LncRNAs contained in the syntenic blocks.
                db_LncRNAs_select_x = d[(d["Spe"] == genome_x) & \
                                        (d["Chr"] == min_chr_x) & \
                                        (((d["Start"] < min_pos_x) & (d["End"] >= min_pos_x)) | \
                                         ((d["Start"] >= min_pos_x) & (d["End"] <= max_pos_x)) | \
                                         ((d["Start"] <= max_pos_x) & (d["End"] > max_pos_x)))]
                db_LncRNAs_select_y = d[(d["Spe"] == genome_y) & \
                                        (d["Chr"] == min_chr_y) & \
                                        (((d["Start"] < min_pos_y) & (d["End"] >= min_pos_y)) | \
                                         ((d["Start"] >= min_pos_y) & (d["End"] <= max_pos_y)) | \
                                         ((d["Start"] <= max_pos_y) & (d["End"] > max_pos_y)))]
                
                # Join both subsets.
                tab = pd.concat([db_LncRNAs_select_x, db_LncRNAs_select_y])
                tab['CloudID'] = ID
                tab['Study'] = study
                
                # Add the table of this syntenic block to the table with the rest of sintenic blocks.
                TAB = pd.concat([TAB, tab])
                
                # Add hits to the list of hits. I mean, with the rest of hits found in other syntenic blocks..
                for a in range(len(list(db_LncRNAs_select_x["ID_transcript"]))):
                    for b in range(len(list(db_LncRNAs_select_y["ID_transcript"]))):
                        list_hits.append([genome_x, genome_y, list(db_LncRNAs_select_x["ID_transcript"])[a], list(db_LncRNAs_select_y["ID_transcript"])[b]])
                        list_hits.append([genome_y, genome_x, list(db_LncRNAs_select_y["ID_transcript"])[b], list(db_LncRNAs_select_x["ID_transcript"])[a]])
            
            # Save the table for a particular study and confidence level.
            #TAB.to_csv(path_3 + "/output_cloud/LncRNAs-" + study + "-" + confidence + ".tsv", sep = '\t', index = False, header = True)
            
            # Add the table about all syntenic blocks contained in a particular confidence level and study to the rest of studies and confidence
            # levels of the table.
            TAB_ALL = pd.concat([TAB_ALL, TAB])
            
            # Remove duplicated records.
            my_set = set(tuple(l) for l in list_hits)
            new_list_hits = [list(tup) for tup in my_set]
            dic_hits_B[confidence] = new_list_hits
            
            if len(list_hits) == len(new_list_hits):
                print("There isn't any duplicated lncRNA hit as a consequence of two syntenic blocks of the same specie sharing a syntenic block in the other specie (" + str(len(list_hits)) + "-" + str(len(list_hits)) + ").")
            else:
                print("There are duplicated lncRNA hits as a consequence of two syntenic blocks of the same specie sharing a syntenic block in the other specie (" + str(len(list_hits)) + "-" + str(len(list_hits)) + ").")
            
        dic_hits_A[study] = dic_hits_B
            
    
    # Convert these columns to factors.
    TAB_ALL['Study'] = pd.Categorical(TAB_ALL['Study'], categories = studies_factor)
    TAB_ALL['Confidence_level'] = pd.Categorical(TAB_ALL['Confidence_level'], categories = confidences_factor)
    TAB_ALL['Class_code'] = pd.Categorical(TAB_ALL['Class_code'], categories = class_codes_factor)
    
    # Save the table with all the studies and confidence levels.
    TAB_ALL.to_csv(path_3 + "/output_cloud/LncRNAs-Global.tsv", sep = '\t', index = False, header = True)
    
    # Total LncRNAs by study and confidence level.
    dic_A = {}
    for i in range(len(studies)):
        study = studies[i]
        study_class_codes = studies_class_codes[i]
        dic_B = {}
        for confidence in confidences:
            dic_C = {}
            for spe in [spe1, spe2]:
                d = db_LncRNAs[(db_LncRNAs["Class_code"].isin(study_class_codes)) & (db_LncRNAs["Confidence_level"] == confidence)]
                d_spe = d[(d["Spe"] == spe)]
                total_lncRNAs = len(d_spe.index)
                dic_C[spe] = total_lncRNAs
            dic_B[confidence] = dic_C
        dic_A[study] = dic_B
    
    print("\nCollapse 1...")
    TAB_ALL_1 = TAB_ALL
    TAB_ALL_1 = TAB_ALL_1.drop('CloudID', axis=1)
    TAB_ALL_1 = TAB_ALL_1.drop_duplicates(keep = "first")
    TAB_ALL_1["Count"] = 1
    TAB_ALL_COLLAPSED_1 = TAB_ALL_1.groupby(['Study', 'Confidence_level', 'Class_code', 'Spe'])['Count'].sum().reset_index()
    TAB_ALL_COLLAPSED_1 = TAB_ALL_COLLAPSED_1.sort_values(['Study', 'Confidence_level', 'Class_code'])
    TAB_ALL_COLLAPSED_1['Count'] = TAB_ALL_COLLAPSED_1['Count'].fillna(0)
    TAB_ALL_COLLAPSED_1['Count'] = TAB_ALL_COLLAPSED_1['Count'].astype(int)
    TAB_ALL_COLLAPSED_1['R'] = TAB_ALL_COLLAPSED_1.apply(lambda x: x.Class_code in dic_studies[x.Study], axis=1)
    TAB_ALL_COLLAPSED_1.drop(TAB_ALL_COLLAPSED_1[TAB_ALL_COLLAPSED_1.R == False].index, inplace=True)
    del TAB_ALL_COLLAPSED_1["R"]
    TAB_ALL_COLLAPSED_1['Total'] = TAB_ALL_COLLAPSED_1.apply(lambda x: dic_A[x.Study][x.Confidence_level][x.Spe], axis=1)
    TAB_ALL_COLLAPSED_1['Percentage'] = np.around((np.array(list(TAB_ALL_COLLAPSED_1['Count']))*100)/np.array(list(TAB_ALL_COLLAPSED_1['Total'])), 2)
    TAB_ALL_COLLAPSED_1.to_csv(path_3 + "/output_cloud/LncRNAs-Global-Collapsed-1.tsv", sep = '\t', index = False, header = True)
    
    print("Collapse 2...")
    TAB_ALL_2 = TAB_ALL
    TAB_ALL_2["Count"] = 1
    TAB_ALL_COLLAPSED_2 = TAB_ALL_2.groupby(['Study', 'Confidence_level', 'Class_code', 'Spe', 'CloudID'])['Count'].sum().reset_index()
    TAB_ALL_COLLAPSED_2 = TAB_ALL_COLLAPSED_2.sort_values(['Study', 'Confidence_level'])
    TAB_ALL_COLLAPSED_2['Count'] = TAB_ALL_COLLAPSED_2['Count'].fillna(0)
    TAB_ALL_COLLAPSED_2['Count'] = TAB_ALL_COLLAPSED_2['Count'].astype(int)
    TAB_ALL_COLLAPSED_2['R'] = TAB_ALL_COLLAPSED_2.apply(lambda x: x.Class_code in dic_studies[x.Study], axis=1)
    TAB_ALL_COLLAPSED_2.drop(TAB_ALL_COLLAPSED_2[TAB_ALL_COLLAPSED_2.R == False].index, inplace=True)
    del TAB_ALL_COLLAPSED_2["R"]
    TAB_ALL_COLLAPSED_2['Total'] = TAB_ALL_COLLAPSED_2.apply(lambda x: dic_A[x.Study][x.Confidence_level][x.Spe], axis=1)
    TAB_ALL_COLLAPSED_2['Percentage'] = np.around((np.array(list(TAB_ALL_COLLAPSED_2['Count']))*100)/np.array(list(TAB_ALL_COLLAPSED_2['Total'])), 2)
    TAB_ALL_COLLAPSED_2.to_csv(path_3 + "/output_cloud/LncRNAs-Global-Collapsed-2.tsv", sep = '\t', index = False, header = True)
    
    # Hits
    print("\nHits...")
    for study in studies:
        for confidence in confidences:
            file = open(path_3 + "/output_cloud/CloudHits-" + study + "-" + confidence + ".tsv", 'w')
            list_hits = dic_hits_A[study][confidence]
            for hit in list_hits:
                file.write("\t".join(hit) + "\n")
            
    
    
# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
