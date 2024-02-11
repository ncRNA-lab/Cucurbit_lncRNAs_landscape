#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- STEP 2: CREATE FINAL TABLES (TRANSPOSON)

This script creates one table showing the percentage of bases covered by only
transposons in each category (LncRNAs, PCGs and IR) encompassing all of types pf
transposons as a single group (Table 1). Sequences with an overlap percentage of 
0 do not appear.

---------

Created on Sun Dec 25 20:00:00 2022

Last Modified:
    - Sun Dec 25 20:00:00 2022 --> Initial code.

@author: pasviber - Pascual Villalba Bermell

"""


# MODULES

import sys
import argparse
import pandas as pd
import math


# VARIABLES

parser = argparse.ArgumentParser(
        prog='CollapseTableTransposon V1',
        description='''This program collapses table by transcript id.''', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
parser.add_argument(
        "--path", 
        type=str, 
        nargs=1,
        help="Absolute path of directory."
        )
parser.add_argument(
        "--confidence", 
        type=str, 
        nargs=1,
        help="Confidence."
        )
parser.add_argument(
        "--flag", 
        type=str, 
        nargs=1,
        help="Flag."
        )
parser.add_argument(
        '--version', 
        action='version', 
        version='%(prog)s 1.0'
        )

args = parser.parse_args()

try:
    WD = args.path[0]
    confidence = args.confidence[0]
    flag = args.flag[0]
     
except:
    print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
    parser.print_help()
    sys.exit()

# WD = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-TEs_and_genomic_repeats/cme/02-Intersection/Final_tables"
# confidence = "High"
# flag = "NR"


# PIPELINE

file = WD + "/Final_tab-Transposon-" + flag + "-" + confidence + ".tsv"
TAB = pd.read_csv(file, sep='\t', header=0)
class_codes = list(set(TAB["Class_code"]))
Repeats = list(set(TAB["Repeat_type_2_mod"]))

print("\n\nBY TRANSPOSON:")    
TAB = TAB[(TAB["Repeat_type_2_mod"] == "Transposon")]
List_results = []
for cl in class_codes:
    print("\t-" + cl)
    subset_1 = TAB[(TAB["Class_code"] == cl)]
    T_ids = list(set(subset_1["ID_transcript"]))
    for T_id in T_ids:
        subset_2 = subset_1[(subset_1["ID_transcript"] == T_id)]
        if len(subset_2.index) == 1:
            row = subset_2.values.tolist()[0]
            length = int(row[6]) - int(row[5])
            overlap = int(row[13])
            overlap_per = (overlap*100)/length
            overlap_per_log = math.log(overlap_per + 1, 2)
            row = row[:7] + [overlap, round(overlap_per, 3), round(overlap_per_log, 3)]
        else:
            ti = int(list(set(subset_2["Start"]))[0]) + 1 # convert 1-based.
            tf = int(list(set(subset_2["End"]))[0])
            tpositions = [x for x in range(ti, tf + 1)]
            length_i = len(tpositions)
            for i in range(len(subset_2.index)):
                ri = int(list(subset_2["Start_rep"])[i]) + 1 # convert 1-based.
                rf = int(list(subset_2["End_rep"])[i])
                rpositions = [x for x in range(ri, rf + 1)]
                tpositions = list(set(tpositions).difference(rpositions))
            length_f = len(tpositions)
            overlap = length_i - length_f
            overlap_per = (overlap*100)/length_i
            overlap_per_log = math.log(overlap_per + 1, 2)
            row = subset_2.values.tolist()[0]
            row = row[:7] + [overlap, round(overlap_per, 3), round(overlap_per_log, 3)] 
        List_results.append(row)

column_names = ["Spe", "Class_code", "Chr", "Strand", "ID_transcript", "Start", "End", "Overlap", "Overlap_per", "Log.overlap_per"]
TAB_filtered = pd.DataFrame(List_results, columns = column_names)
TAB_filtered.to_csv(WD + "/Final_tab-Transposon-" + flag + "-" + confidence + "-COLLAPSED_TRANSPOSON.tsv", sep='\t', index = False, header = True)

