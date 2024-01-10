#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- MODIFY RNAFOLD FILE.

This script will create a new table with the RNAfold results. The free energy
information will be removed.

---------

Created on Wed Aug 31 11:30:00 2022

Last Modified:
    - Wed Aug 31 11:30:00 2022 --> Initial code.

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
            prog='ModifyRNAFoldFile V1', 
            description='''This script will create a new table with the RNAfold results. The free energy information will be removed.''', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    
    parser.add_argument(
            "--input", 
            type=str, 
            nargs=1,
            help="Absolute path of RNAfold results."
            )
    parser.add_argument(
            "--output", 
            type=str, 
            nargs=1,
            help="Absolute path of modified RNAfold results."
            )
    parser.add_argument(
            '--version', 
            action='version', 
            version='%(prog)s 1.0'
            )
    
    args = parser.parse_args()
    
    try:
        input_file = args.input[0]
        output_file = args.output[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Modify_RNAfold_file.py 
    #--input /mnt/doctorado/.../csa_to_cme.tsv
    #--output /mnt/doctorado/.../cme_and_csa.tsv
    
    
    input_file_open = open(input_file, 'r')
    output_file_open = open(output_file, 'w')
    
    lines = input_file_open.readlines()
    New = []
    i = 0
    while i < len(lines):
        New.append([">" + lines[i][1:], lines[i + 1], lines[i + 2].split(" ")[0]])
        i += 3
        
    for x in New:
        if x[2].startswith("("):
            output_file_open.write(x[0] + "N" + x[1] + "." + x[2] + "\n")
        else:
            output_file_open.write("".join(x) + "\n")


# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()



