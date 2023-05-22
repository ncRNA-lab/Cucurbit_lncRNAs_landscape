#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- GENERATE A SUMMARY STRAND INFO TABLE.

This script creates a summary strand info table from pseudomapping results.

---------

Created on Thu Feb 10 10:00:00 2022

Last Modified:
    - Thu Feb 10 10:00:00 2022 --> Initial code.
    - Thu Mar 09 12:15:00 2023 --> Add argparse code.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse

# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
            prog='SummaryStrandInfo V1', 
            description='''This script creates a summary strand info table from pseudomapping results.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--path", 
            type=str, 
            nargs=1,
            help="Absolute path to work directory."
            )
    parser.add_argument(
            "--acc-list", 
            type=str, 
            nargs=1,
            help="List of accessions (samples) which passed the trimming filter."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        path = args.path[0]
        acc_list = args.acc_list[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Generate_strand_info.py
	#--path $WD1
	#--acc-list $Acc_list
    
        
    ### PIPELINE
    samples = open(acc_list, "r+")
    output = open(path + "/03-Strand_detection/Filter_table/strand_info.tsv", "w")
    
    output.write("sample\tmapping_rate\tstrandedness\tnote\n")
        
    for sample in samples:
        sample = sample.rstrip()
        log_file = open(path + "/03-Strand_detection/03-Pseudomapping/Outputs/stdout_Pseudomapping_" + sample + ".log", "r+")
        
        Lib_info = 0
        mapping_info = 0
        for line in log_file.readlines():
            if "Mapping rate =" in line:
                mapping_info += 1
                mapping_rate = line.rstrip().split(" ")[-1]
            if  "Automatically detected most" in line:
                Lib_info += 1
                strand = line.rstrip().split(" ")[-1]
                
        if Lib_info == 0 and mapping_info == 0:
            output.write(sample + "\tNA\tNA\tnot used for transcript reconstruction\n")
        if Lib_info == 0 and mapping_info == 1:
            output.write(sample + "\t" + mapping_rate + "\tNA\tnot used for transcript reconstruction\n")
        if Lib_info == 1 and mapping_info == 0:
            output.write(sample + "\tNA\t" +  strand + "\tnot used for transcript reconstruction\n")
            print(line + ": ERROR")
        if Lib_info == 1 and mapping_info == 1:
            if "U" in strand:
                output.write(sample + "\t" + mapping_rate + "\t" + strand + "\tnot used for transcript reconstruction\n")
            else:
                output.write(sample + "\t" + mapping_rate + "\t" + strand + "\tused for transcript reconstruction\n")


# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()      
            
