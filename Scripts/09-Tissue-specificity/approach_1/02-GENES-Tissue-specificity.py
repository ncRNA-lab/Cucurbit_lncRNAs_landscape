#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GENES: TISSUE SPECIFICITY STUDY: APPROACH 1 - STEP 2

Calculation of metrics for assessing tissue specificity of transcripts (Genes). 
We get TAU metric (General scoring metrics) and TSI metrics (Individualized 
scoring metrics). Therefore, we generate two files.

---------

Created on Wed Jan 04 10:00:00 2022

Last Modified:
    - Wed Jan 04 10:00:00 2022 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES
import os
import tspex
import pandas as pd


# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    ## Own computer.
    species = ["car", "cla", "cma", "cme", "cmo", "cpe", "csa", "lsi", "mch"]
    flag = "nr"
    path = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity/approach_1/INDIVIDUAL"
    
    ### Garnatxa.
    #species = ["car", "cla", "cma", "cme", "cmo", "cpe", "csa", "lsi", "mch"]
    #flag = "nr"
    #path = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity/approach_1/INDIVIDUAL"
    
    for spe in species:
        print("\tSpecie: " + spe)
        
        ## Paths and directories
        path_spe_in = path + "/" + flag + "/STEP1/" + spe
        if not os.path.exists(path + "/" + flag + "/STEP2/"):
            os.makedirs(path + "/" + flag + "/STEP2/")
        path_spe_out = path + "/" + flag + "/STEP2/" + spe
        if not os.path.exists(path_spe_out):
            os.makedirs(path_spe_out)
        
        ## Files
        Files_median = [file for file in os.listdir(path_spe_in) if "_median.tsv" in file and "GENES-" in file]
        Files_mean = [file for file in os.listdir(path_spe_in) if "_mean.tsv" in file and "GENES-" in file]
        
        ## Median
        print("\t\tMedian values...")
        if len(Files_median) > 0:
            for file in Files_median:
                print("\t\t\tSRA.Study: " + file.replace("_median.tsv", "").replace("GENES-", ""))
                
                ## Load the table built in STEP 1.
                tab = pd.read_csv(path_spe_in + "/" + file, sep = '\t', header = 0)
                tab = tab.set_index("ID_transcript")
                cols = list(tab.columns)
                tab_info = tab.loc[:, cols[:2]]
                tab_exp = tab.loc[:, cols[2:]]
                
                ## Filter the transcripts by expression values.
                tab_exp_filt = tab_exp.loc[(tab_exp >= 1).any(axis=1)]
                tab_info_filt = tab_info[tab_info.index.isin(list(tab_exp_filt.index.values))]
                
                ## TAU: General scoring metric. Describe in a single value how tissue-specific 
                ## or ubiquitous is a gene across all tissues (counts, tau, gini, simpson, 
                ## shannon_specificity, roku_specificity, spm_dpm and js_specificity_dpm).
                tso = tspex.TissueSpecificity(tab_exp_filt, 'tau', log = True)
                # expression_data attribute.
                tab_tau_exp = pd.DataFrame(tso.expression_data)
                tab_tau_exp.columns = ['EXP.' + col for col in cols[2:]]
                # tissue_specificity attribute.
                tab_tau_ts = pd.DataFrame(tso.tissue_specificity)
                tab_tau_ts.columns = ['TAU']
                
                ## TSI: Individualized scoring metric. Quantify how specific is the 
                ## expression of each gene to each tissue (tsi, zscore, spm and js_specificity). 
                tso = tspex.TissueSpecificity(tab_exp_filt, 'tsi', log = True)
                # expression_data attribute.
                tab_tsi_exp = pd.DataFrame(tso.expression_data)
                tab_tsi_exp.columns = ['EXP.' + col for col in cols[2:]]
                # tissue_specificity attribute.
                tab_tsi_ts = pd.DataFrame(tso.tissue_specificity)
                tab_tsi_ts.columns = ['TSI.' + col for col in cols[2:]]
                
                ## Convert index to column and reset index.
                tab_info_filt.insert(loc = len(tab_info_filt.columns), column = "ID_transcript", value = list(tab_info_filt.index))
                tab_info_filt.reset_index(inplace = True, drop = True)
                
                tab_tau_ts.insert(loc = len(tab_tau_ts.columns), column = "ID_transcript", value = list(tab_tau_ts.index))
                tab_tau_ts.reset_index(inplace = True, drop = True)
                tab_tsi_ts.insert(loc = len(tab_tsi_ts.columns), column = "ID_transcript", value = list(tab_tsi_ts.index))
                tab_tsi_ts.reset_index(inplace = True, drop = True)
                
                tab_tau_exp.insert(loc = len(tab_tau_exp.columns), column = "ID_transcript", value = list(tab_tau_exp.index))
                tab_tau_exp.reset_index(inplace = True, drop = True)
                tab_tsi_exp.insert(loc = len(tab_tsi_exp.columns), column = "ID_transcript", value = list(tab_tsi_exp.index))
                tab_tsi_exp.reset_index(inplace = True, drop = True)
                
                ## Create final tables.
                tab_tau_final = pd.merge(tab_info_filt, tab_tau_exp, how = "inner", on = "ID_transcript")
                tab_tau_final = pd.merge(tab_tau_final, tab_tau_ts, how = "inner", on = "ID_transcript")
                tab_tau_final = tab_tau_final.loc[:, ["ID_transcript", "Confidence", "Class_code"] + ['EXP.' + col for col in cols[2:]] + ["TAU"]]
                
                tab_tsi_final = pd.merge(tab_info_filt, tab_tsi_exp, how = "inner", on = "ID_transcript")
                tab_tsi_final = pd.merge(tab_tsi_final, tab_tsi_ts, how = "inner", on = "ID_transcript")
                tab_tsi_final = tab_tsi_final.loc[:, ["ID_transcript", "Confidence", "Class_code"] + ['EXP.' + col for col in cols[2:]] + ['TSI.' + col for col in cols[2:]]]
                
                ## Save final tables.
                tab_tau_final.to_csv(path_spe_out + "/" + file.replace(".tsv", "") + "-TAU.tsv", sep = "\t", index = False, header = True)
                tab_tsi_final.to_csv(path_spe_out + "/" + file.replace(".tsv", "") + "-TSI.tsv", sep = "\t", index = False, header = True)
        else:
            print("\t\t\tThere is no SRA.study")
        
        ## Mean
        print("\t\tMean values...")
        if len(Files_mean) > 0:
            for file in Files_mean:
                print("\t\t\tSRA.Study: " + file.replace("_mean.tsv", "").replace("GENES-", ""))
                
                ## Load the table built in STEP 1.
                tab = pd.read_csv(path_spe_in + "/" + file, sep = '\t', header = 0)
                tab = tab.set_index("ID_transcript")
                cols = list(tab.columns)
                tab_info = tab.loc[:, cols[:2]]
                tab_exp = tab.loc[:, cols[2:]]
                
                ## Filter the transcripts by expression values.
                tab_exp_filt = tab_exp.loc[(tab_exp >= 1).any(axis=1)]
                tab_info_filt = tab_info[tab_info.index.isin(list(tab_exp_filt.index.values))]
                
                ## TAU: General scoring metric. Describe in a single value how tissue-specific 
                ## or ubiquitous is a gene across all tissues (counts, tau, gini, simpson, 
                ## shannon_specificity, roku_specificity, spm_dpm and js_specificity_dpm).
                tso = tspex.TissueSpecificity(tab_exp_filt, 'tau', log = True)
                # expression_data attribute.
                tab_tau_exp = pd.DataFrame(tso.expression_data)
                tab_tau_exp.columns = ['EXP.' + col for col in cols[2:]]
                # tissue_specificity attribute.
                tab_tau_ts = pd.DataFrame(tso.tissue_specificity)
                tab_tau_ts.columns = ['TAU']
                
                ## TSI: Individualized scoring metric. Quantify how specific is the 
                ## expression of each gene to each tissue (tsi, zscore, spm and js_specificity). 
                tso = tspex.TissueSpecificity(tab_exp_filt, 'tsi', log = True)
                # expression_data attribute.
                tab_tsi_exp = pd.DataFrame(tso.expression_data)
                tab_tsi_exp.columns = ['EXP.' + col for col in cols[2:]]
                # tissue_specificity attribute.
                tab_tsi_ts = pd.DataFrame(tso.tissue_specificity)
                tab_tsi_ts.columns = ['TSI.' + col for col in cols[2:]]
                
                ## Convert index to column and reset index.
                tab_info_filt.insert(loc = len(tab_info_filt.columns), column = "ID_transcript", value = list(tab_info_filt.index))
                tab_info_filt.reset_index(inplace = True, drop = True)
                
                tab_tau_ts.insert(loc = len(tab_tau_ts.columns), column = "ID_transcript", value = list(tab_tau_ts.index))
                tab_tau_ts.reset_index(inplace = True, drop = True)
                tab_tsi_ts.insert(loc = len(tab_tsi_ts.columns), column = "ID_transcript", value = list(tab_tsi_ts.index))
                tab_tsi_ts.reset_index(inplace = True, drop = True)
                
                tab_tau_exp.insert(loc = len(tab_tau_exp.columns), column = "ID_transcript", value = list(tab_tau_exp.index))
                tab_tau_exp.reset_index(inplace = True, drop = True)
                tab_tsi_exp.insert(loc = len(tab_tsi_exp.columns), column = "ID_transcript", value = list(tab_tsi_exp.index))
                tab_tsi_exp.reset_index(inplace = True, drop = True)
                
                ## Create final tables.
                tab_tau_final = pd.merge(tab_info_filt, tab_tau_exp, how = "inner", on = "ID_transcript")
                tab_tau_final = pd.merge(tab_tau_final, tab_tau_ts, how = "inner", on = "ID_transcript")
                tab_tau_final = tab_tau_final.loc[:, ["ID_transcript", "Confidence", "Class_code"] + ['EXP.' + col for col in cols[2:]] + ["TAU"]]
                
                tab_tsi_final = pd.merge(tab_info_filt, tab_tsi_exp, how = "inner", on = "ID_transcript")
                tab_tsi_final = pd.merge(tab_tsi_final, tab_tsi_ts, how = "inner", on = "ID_transcript")
                tab_tsi_final = tab_tsi_final.loc[:, ["ID_transcript", "Confidence", "Class_code"] + ['EXP.' + col for col in cols[2:]] + ['TSI.' + col for col in cols[2:]]]
                
                ## Save final tables.
                tab_tau_final.to_csv(path_spe_out + "/" + file.replace(".tsv", "") + "-TAU.tsv", sep = "\t", index = False, header = True)
                tab_tsi_final.to_csv(path_spe_out + "/" + file.replace(".tsv", "") + "-TSI.tsv", sep = "\t", index = False, header = True)
        else:
            print("\t\t\tThere is no SRA.study")

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
