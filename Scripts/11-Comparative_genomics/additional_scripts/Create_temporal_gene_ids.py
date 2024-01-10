#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- CREATE TEMPORAL GENE IDS

This program creates temporal gene IDs and modifies FASTA and GTF files with this 
temporal gene IDs. Besides, it creates a database of temporal gene IDs and original 
gene IDS.

---------

Created on Thu Jan 19 12:30:00 2023

Last Modified:
    - Thu Jan 19 12:30:00 2023 --> Initial code. Two scripts.

@author: pvbermell - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse


# FUNCTIONS

def ModifyFasta (fasta_path, dicOLDNEW):
    
    dict_new_fasta = {}
    fasta_in = open(fasta_path, 'r')
    for line in fasta_in:
        line = line.strip()
        if line[0] == ">":
            old_id = line[1:]
            new_id = dicOLDNEW[old_id]
            dict_new_fasta[new_id] = ""
        else:
            value = line
            dict_new_fasta[new_id] += value
    fasta_in.close()
    
    return dict_new_fasta

def WriteFasta (dict_new_fasta, new_fasta_path):
    
    new_fasta_out = open(new_fasta_path, 'w')
    for key in list(dict_new_fasta.keys()):
        new_fasta_out.write(">" + key + "\n" + dict_new_fasta[key] + "\n")
    new_fasta_out.close()

def ModifyGtf (gtf_path, dicOLDNEW):
    
    list_new_gtf = []
    gtf_in = open(gtf_path, 'r')
    for line in gtf_in:
        line = line.strip().split("\t")
        transcript_id_old = line[8].split(";")[0].split(" ")[1].replace('"', "")
        transcript_id_new = dicOLDNEW[transcript_id_old]
        gene_id_new = transcript_id_new + "-gene"
        new_line = line[:8] + ['transcript_id "' + transcript_id_new + '"; gene_id "' + gene_id_new + '";']
        list_new_gtf.append(new_line)
    
    return list_new_gtf

def WriteGtf (list_new_gtf, new_gtf_path):
    
    new_gtf_out = open(new_gtf_path, 'w')
    for element in list_new_gtf:
        new_gtf_out.write("\t".join(element) + "\n")
    new_gtf_out.close()
    

# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
        prog='CreateTemporalGeneIDs', 
        description='''This program creates temporal gene IDs and modifies FASTA  \
            and GTF files with this temporal gene IDs. Besides, it creates a database \
            of temporal gene IDs and original gene IDS.''',  
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    
    parser.add_argument(
        "--spe", 
        type=str, 
        nargs=1,
        help="Specie."
        )
    parser.add_argument(
        "--ids", 
        type=str, 
        nargs=1,
        help="Absolute path of gene ids list (TXT file)."
        )
    parser.add_argument(
        "--gtf", 
        type=str, 
        nargs='?',
        help="Absolute path of genes GTF file."
        )
    parser.add_argument(
        "--fasta", 
        type=str, 
        nargs='?',
        help="Absolute path of FASTA file with the sequences coming from the \
            GTF/GFF3 file. You can use gffread or rsem-prepare-reference to\
            extract the sequences."
        )
    parser.add_argument(
        "--new-ids", 
        type=str, 
        nargs=1,
        help="Absolute path of new gene IDs list (TXT file)."
        )
    parser.add_argument(
        "--new-gtf", 
        type=str, 
        nargs='?',
        help="Absolute path of genes GTF file with the new IDs."
        )
    parser.add_argument(
        "--new-fasta", 
        type=str, 
        nargs='?',
        help="Absolute path of genes FASTA file with the new IDs."
        )
    parser.add_argument(
        "--database", 
        type=str, 
        nargs=1,
        help="Absolute path of new gene IDs database. key: temporal gene id, \
            value: previous gene id"
        )
    parser.add_argument(
        "--mode", 
        type=str, 
        nargs=1, 
        choices=['Blastn', 'OrthoFinder', 'OrthoMCL'],
        help="Absolute path of output table."
        )
    parser.add_argument(
        '--version', 
        action='version', 
        version='%(prog)s 3.0'
        )
    
    args = parser.parse_args()
    
    try:
        spe = args.spe[0]
        ids = args.ids[0]
        gtf = args.gtf
        fasta = args.fasta
        new_ids = args.new_ids[0]
        new_gtf = args.new_gtf
        new_fasta = args.new_fasta
        database = args.database[0]
        mode = args.mode[0]
    except:
        print("\nERROR: You have inserted a wrong argument or you are missing an argument.\n")
        parser.print_help()
        sys.exit()
    
    if gtf == None and fasta == None:
        print("\nERROR: You have to insert either gtf argument, gff3 argument or fasta argument. One of them are required.\n")
        parser.print_help()
        sys.exit()
    if (gtf == None and new_gtf != None) or (gtf != None and new_gtf == None):
        print("\nERROR: --gtf and --new-gtf always together.\n")
        parser.print_help()
        sys.exit()
    if (fasta == None and new_fasta != None) or (fasta != None and new_fasta == None):
        print("\nERROR: --fasta and --new-fasta always together.\n")
        parser.print_help()
        sys.exit()
    
    """
    spe = "cme"
    ids = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme/STEP-FINAL/Files/Genes/ORIGINAL_GENES_ids.txt"
    fasta = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme/STEP-FINAL/Files/Genes/ORIGINAL_GENES.fasta"
    gtf = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme/STEP-FINAL/Files/Genes/ORIGINAL_GENES.gtf"
    new_ids = "/home/pvbermell/Documentos/ORIGINAL_GENES_new_ids.txt"
    new_fasta = "/home/pvbermell/Documentos/ORIGINAL_GENES_new.fasta"
    new_gtf = "/home/pvbermell/Documentos/ORIGINAL_GENES_new.gtf"
    database = "/home/pvbermell/Documentos/cme.tsv"
    mode = "Blastn"
    """
    
    ### PIPELINE
    
    ## Create two dictionaries: 
    old_ids_list = [line for line in open(ids, 'r')]
    
    # dicOLDNEW['Old_gene_id'] = 'New_gene_id' and dicNEWOLD['New_gene_id'] = 'Old_gene_id'
    dicOLDNEW = {}
    dicNEWOLD = {}
    i = 1
    for old_id in old_ids_list:
        old_id = old_id.strip()
        if mode == "Blastn":
            new_id = spe + "_" + str(i)
        elif mode == "OrthoFinder" or mode == "OrthoMCL":
            new_id = spe + "_" + str(i) + "-" + spe
        dicOLDNEW[old_id] = new_id
        dicNEWOLD[new_id] = old_id
        i += 1
        
    # Save the database (dicNEWOLD)
    database_out = open(database, 'w')
    for key in list(dicNEWOLD.keys()):
        database_out.write(key + "\t" + dicNEWOLD[key] + "\n")
    database_out.close()
    
    # Save new_ids.
    new_ids_out = open(new_ids, 'w')
    for element in list(dicNEWOLD.keys()):
        new_ids_out.write(element + "\n")
    new_ids_out.close()
    
    
    ## Convert fasta to new_fasta.
    if fasta != None:
        if fasta.endswith(".fas") or fasta.endswith(".fa") or fasta.endswith(".fasta"):
            dict_new_fasta = ModifyFasta (fasta, dicOLDNEW)
            WriteFasta (dict_new_fasta, new_fasta)
        else:
            print("\nERROR (FASTA): Check the sequence file name.\n")
            sys.exit()
    
    ## Convert gtf to new_gtf.
    if gtf != None:
        if gtf.endswith("gtf"):
            list_new_gtf = ModifyGtf (gtf, dicOLDNEW)
            WriteGtf (list_new_gtf, new_gtf)
        else:
            print("\nERROR (GTF): Check the annotation file name.\n")
            sys.exit()
    
    
# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
