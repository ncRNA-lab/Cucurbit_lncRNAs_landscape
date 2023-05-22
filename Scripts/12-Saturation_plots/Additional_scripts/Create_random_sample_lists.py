#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- CREATE RANDOM SAMPLE LISTS.

This software creates random sample lists from other sample list.

---------

Created on Mon Jan 9 20:00:00 2023

Last Modified:
    - Mon Jan 9 20:00:00 2023 --> Initial code.

@author: pvbermell - Pascual Villalba Bermell

"""

## MODULES

import sys
import argparse
import random


## MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
        prog='RandomSampleLists V1',
        description='''This program filters a GTF file by transcript_id''',
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    
    parser.add_argument(
        "-i", 
        "--input-list", 
        type=str, 
        nargs=1,
        help="Absolute path of initial list used to generate the other lists randomly."
        )
    parser.add_argument(
        "-o", 
        "--output-path", 
        type=str, 
        nargs=1,
        help="Absolute path where the new lists are stored."
        )
    parser.add_argument(
        "-S", 
        "--step-size", 
        type=int, 
        nargs=1,   
        help="The number of new samples to be analyzed in each iteration."
        )
    parser.add_argument(
        '--version', 
        action='version', 
        version='%(prog)s 2.0'
        )
    
    args = parser.parse_args()
    
    try:
        input_list = args.input_list[0]
        output_path = args.output_path[0]
        step_size = args.step_size[0]
         
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    # python3 ./Create_random_sample_lists.py 
    #--input-list /mnt/doctorado/.../....txt
    #--output-path /mnt/doctorado/.../...
    #--number 5 
    
    ## PIPELINE
    
    file_in = open(input_list, 'r')
    Samples = [line.strip() for line in file_in]
    file_in.close
    
    number_samples = len(Samples)
    iterations = number_samples//step_size
    
    for i in range(1,iterations + 1):
        number_subset = i * step_size
        subset = random.sample(Samples, k = number_subset)
        file_out = open(output_path + "/Batch-" + str(i) + ".txt", 'w')
        file_out.write("\n".join(subset) + "\n")
        file_out.close
        
    if iterations * step_size < number_samples:
        subset = random.sample(Samples, k = number_samples)
        file_out = open(output_path + "/Batch-" + str(iterations + 1) + ".txt", 'w')
        file_out.write("\n".join(subset) + "\n")
        file_out.close
    

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
