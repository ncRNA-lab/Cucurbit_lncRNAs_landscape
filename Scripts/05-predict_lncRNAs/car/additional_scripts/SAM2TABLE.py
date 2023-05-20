#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
- CONVERT SAM TO TABLE (TSV,CSV...)

This script converts a SAM file to a TABLE with the separator chosen. 
In addition, all the lines must measure the length of the maximum length 
line. Then, this script fills the lines which measure less than maximum 
length with empty elements "".

---------

Created on Thu Jan 21 13:22:18 2021

Last Modified:
    - Thu Jan 21 13:22:18 2021 --> Initial code.
    - Sun May 29 13:50:00 2022 --> Add arguments command line.

@author: Pascual Villalba-Bermell
"""

## MODULES

import sys
import argparse


# MAIN PROGRAM

def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
        prog='SAM2TABLE V2',
        description='''This program converts a SAM file to TABLE.''',
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    
    parser.add_argument(
        "-i", 
        "--sam", 
        type=str, 
        nargs=1,
        help="Absolute path of sam file."
        )
    parser.add_argument(
        "-o", 
        "--tab", 
        type=str, 
        nargs=1,
        help="Absolute path of table."
        )
    parser.add_argument(
        "-s", 
        "--sep", 
        type=str, 
        nargs=1,
        help="Field separator."
        )
    parser.add_argument(
        '--version', 
        action='version', 
        version='%(prog)s 2.0'
        )
    
    args = parser.parse_args()
    
    try:
        sam = args.sam[0]
        tab = args.tab[0]
        sep = args.sep[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./SAM2TABLE.py 
    #--sam /mnt/doctorado/.../....sam
    #--tab /mnt/doctorado/.../....tsv
    #--sep $'\t' 
    
    
    ## Determine the maximum length line.
    with open(sam, 'r') as f:
        lines = f.read().splitlines()
        max_ = -1
        for line in lines:
            line = line.split("\t")
            if len(line) > max_:
                max_ = len(line)
    
    ## Open the files in read mode.
    f_in = open(sam, 'r')
    f_out = open(tab, 'w')
    
    ## Create the new file.
    for line in f_in:
        line = line.strip().split("\t")
        L = len(line)
        # If the length line is less than max_, it will be filled with "" 
        # to reach max_.
        if L < max_:
            dif = max_-L
            list_ = [""]*dif
            line = line + list_
        # Now write the new line in the new file.
        new_line = sep.join(line)
        f_out.write(new_line + str("\n"))
    
    f_in.close()
    f_out.close()


# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
