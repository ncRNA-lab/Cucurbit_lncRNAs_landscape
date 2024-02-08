#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- CONVERT GTF/GFF3 TO TABLE (TSV/CSV)

From an annotation file (GTF), coming from a genome-guided assembly approach 
(mode assembly) using Stringtie software followed by annotation with gffcompare
software, to a TSV or CSV table.

From an annotation file (GTF/GFF3) of genes to TSV or CSV table
(mode original). Original because it hasn't been assembled by our pipeline.

---------

Created on Wed Mar 23 12:30:00 2022

Last Modified:
    - Wed Mar 23 12:30:00 2022 --> Initial code. Two scripts.
    - Sun May 29 13:05:00 2022 --> Add arguments command line.
    - Sat Jun 11 17:25:00 2022 --> Merge the two scripts.

@author: pasviber - Pascual Villalba Bermell

"""

# MODULES

import sys
import argparse
import numpy as np
import pandas as pd


# FUNCTIONS

def GTFtoDICT (gtf, mode):
    """
    This function generates a dictionary from the information contained in 
    a GTF file.
    """
    
    Transcript = {}
    
    f_in = open(gtf, 'r')
    
    for line in f_in:
        line = line.strip().split("\t")
        
        ## Avoid commented lines.
        if line[0][:2] == "##" or line[0][0] == "#":
            continue
        
        elif line[2] == "transcript":
            
            chr_ = line[0]
            origin = line[1]
            pi = line[3]
            pf = line[4]
            strand = line[6]
            
            info = line[8].split(";")
            
            # REQUIRED: transcript_id and gene_id.
            if "transcript_id" in line[8] and "gene_id" in line[8]:
                trans_ID = [a for a in info if "transcript_id" in a][0].strip().split(" ")[1].replace('"', "")
                gene_ID = [a for a in info if "gene_id" in a and not "ref_" in a][0].strip().split(" ")[1].replace('"', "")
            else:
                print("\nERROR: Your GTF file requires gene_id and transcript_id attributes for the transcript feature. Check it.\n")
                sys.exit()
            # OPTIONAL: class_code.
            if "class_code" in line[8] and mode == "assembly":
                class_code = [a for a in info if "class_code" in a][0].strip().split(" ")[1].replace('"', "")
            else:
                class_code = np.nan
                
            exons = 0
            cdss = 0
            
            Transcript[trans_ID] = [chr_, origin, pi, pf, strand, trans_ID, gene_ID, class_code, exons, cdss]
            
        elif line[2] == "exon":
            
            info = line[8].split(";")
            
            # REQUIRED: transcript_id and gene_id.
            if "transcript_id" in line[8] and "gene_id" in line[8]:
                trans_ID = [a for a in info if "transcript_id" in a][0].strip().split(" ")[1].replace('"', "")
                gene_ID = [a for a in info if "gene_id" in a and not "ref_" in a][0].strip().split(" ")[1].replace('"', "")
            else:
                print("\nERROR: Your GTF file requires gene_id and transcript_id attributes for the exon feature. Check it.\n")
                sys.exit()
            
            Transcript[trans_ID][8] += 1
        
        elif line[2] == "CDS":
            
            info = line[8].split(";")
            
            # REQUIRED: transcript_id and gene_id.
            if "transcript_id" in line[8] and "gene_id" in line[8]:
                trans_ID = [a for a in info if "transcript_id" in a][0].strip().split(" ")[1].replace('"', "")
                gene_ID = [a for a in info if "gene_id" in a and not "ref_" in a][0].strip().split(" ")[1].replace('"', "")
            else:
                print("\nERROR: Your GTF file requires gene_id and transcript_id attributes for the cds feature. Check it.\n")
                sys.exit()
            
            Transcript[trans_ID][9] += 1
    
    for trans_ID in list(Transcript.keys()):
        if Transcript[trans_ID][8] > 0 and Transcript[trans_ID][9] == 0:
            Transcript[trans_ID] = Transcript[trans_ID][:8] + [str(Transcript[trans_ID][8])]
        elif Transcript[trans_ID][8] == 0 and Transcript[trans_ID][9] > 0:
            Transcript[trans_ID] = Transcript[trans_ID][:8] + [str(Transcript[trans_ID][9])]
        else:
            Transcript[trans_ID] = Transcript[trans_ID][:8] + [str(Transcript[trans_ID][8])]
    
    return Transcript


def GFF3toDICT (gff3, mode):
    """
    This function generates a dictionary from the information contained in 
    a GFF3 file.
    """
    
    Gene = {}
    Transcript = {}
    Conversion = {}
    
    f_in = open(gff3, 'r')
    
    for line in f_in:
        line = line.strip().split("\t")
        
        ## Avoid commented lines.
        if line[0][:2] == "##" or line[0][0] == "#":
            continue
                
        elif line[2] == "gene":
            chr_ = line[0]
            origin = line[1]
            pi = line[3]
            pf = line[4]
            strand = line[6]
            
            info = line[8].split(";")
            
            # REQUIRED: ID.
            if "ID" in line[8]:
                gene_ID = [a for a in info if "ID" in a][0].strip().split("=")[1]
            else:
                print("\nERROR: Your GFF3 file requires ID attribute for the gene feature. Check it.\n")
                sys.exit()
            # OPTIONAL: class_code
            if "class_code" in line[8] and mode == "assembly":
                class_code = [a for a in info if "class_code" in a][0].strip().split("=")[1]
            else:
                class_code = np.nan
            
            Gene[gene_ID] = [chr_, origin, pi, pf, strand, gene_ID, class_code]
            
        elif line[2] == "mRNA":
            
            info = line[8].split(";")
            
            # REQUIRED: ID.
            if "ID" in line[8]:
                trans_ID = [a for a in info if "ID" in a][0].strip().split("=")[1]
            else:
                print("\nERROR: Your GFF3 file requires ID attribute for the mRNA feature. Check it.\n")
                sys.exit()
            # REQUIRED: Parent.
            if "Parent" in line[8]:
                gene_ID = [a for a in info if "Parent" in a][0].strip().split("=")[1]
            else:
                print("\nERROR: Your GFF3 file requires Parent attribute for the mRNA feature. Check it.\n")
                sys.exit()
            
            exons = 0
            cdss = 0
            
            Transcript[trans_ID] = [gene_ID, trans_ID, exons, cdss]
            Conversion[gene_ID] = trans_ID
            
        elif line[2] == "exon":
            
            info = line[8].split(";")
            
            # REQUIRED: Parent.
            if "Parent" in line[8]:
                trans_ID = [a for a in info if "Parent" in a][0].strip().split("=")[1]
            else:
                print("\nERROR: Your GFF3 file requires Parent attribute for the exon feature. Check it.\n")
                sys.exit()
            
            Transcript[trans_ID][2] += 1
        
        elif line[2] == "CDS":
            
            info = line[8].split(";")
            
            # REQUIRED: Parent.
            if "Parent" in line[8]:
                trans_ID = [a for a in info if "Parent" in a][0].strip().split("=")[1]
            else:
                print("\nERROR: Your GFF3 file requires Parent attribute for the cds feature. Check it.\n")
                sys.exit()
            
            Transcript[trans_ID][3] += 1
    
    Transcript_def = {}
    for gene_ID in list(Gene.keys()):
        trans_ID = Conversion[gene_ID]
        if Transcript[trans_ID][2] > 0 and Transcript[trans_ID][3] == 0:
            Transcript_def[trans_ID] = Gene[gene_ID][:5] + [Transcript[trans_ID][1]] + Gene[gene_ID][5:] + [str(Transcript[trans_ID][2])]
        elif Transcript[trans_ID][2] == 0 and Transcript[trans_ID][3] > 0:
            Transcript_def[trans_ID] = Gene[gene_ID][:5] + [Transcript[trans_ID][1]] + Gene[gene_ID][5:] + [str(Transcript[trans_ID][3])]
        else:
            Transcript_def[trans_ID] = Gene[gene_ID][:5] + [Transcript[trans_ID][1]] + Gene[gene_ID][5:] + [str(Transcript[trans_ID][2])]
    
    return Transcript_def


def FASTAtoDICT (fasta):
    """
    This function generates a dictionary from the information contained in 
    a FASTA file.
    """
    
    Transcript = {}
    
    f_in = open(fasta, 'r')
    
    for line in f_in:
        line = line.strip()
        
        if line[0] == ">":
            header = line[1:]
        else:
            seq = line.upper()
            Length = len(seq)
            G = seq.count("G")
            C = seq.count("C")
            GC = round(((G+C)*100)/Length, 2)
            
            Transcript[header] = [str(Length), str(GC)]
    
    return Transcript


def CombineDICTSsameKEY (dict_annot, dict_fasta):
    """
    This function checks that the ids of the two dictionaries match in 
    nomenclature and none of them are missing, and then combines the two
    dictionaries (dict_annot, dict_fasta).
    """
    
    keys1 = list(dict_annot.keys())
    keys2 = list(dict_fasta.keys())
    
    check =  all(item in keys2 for item in keys1)
    
    if check == True:
        dict_combined = {}
        for k in keys1:
            dict_combined[k] = dict_annot[k] + dict_fasta[k]
    else:
        print("\nERROR: Your FASTA file doesn't correspond to the annotation file. Check it.\n")
        sys.exit()
    
    return dict_combined
    

def Write_table (tab, dict_, sep, flag):
    """
    This function writes the output table.
    """
    
    TAB = pd.DataFrame.from_dict(dict_, orient='index')
    
    if flag == "NON-FASTA":
        header = ["Chr", "Origin", "Start", "End", "Strand", "ID_transcript", "ID_Gene", "Class_code", "Exons"]
    elif flag == "FASTA":
        header = ["Chr", "Origin", "Start", "End", "Strand", "ID_transcript", "ID_Gene", "Class_code", "Exons", "Length", "GC"]
    
    TAB.columns = header
    
    TAB.to_csv(tab, sep = sep, header = True, index = None)
    

# MAIN PROGRAM
def main():
    """
    Main program.
    """
    
    parser = argparse.ArgumentParser(
        prog='GTF/GFF3toTable', 
        description='''This program creates a table with all the information \
            contained in the annotation file (GTF/GFF3).''',  
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    
    parser.add_argument(
        "--gtf", 
        type=str, 
        nargs='?',
        help="Absolute path of genes GTF file."
        )
    parser.add_argument(
        "--gff3", 
        type=str, 
        nargs='?',
        help="Absolute path of genes GFF3 file."
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
        "--tab", 
        type=str, 
        nargs=1,
        help="Absolute path of output table."
        )
    parser.add_argument(
        "--mode", 
        type=str, 
        nargs=1, 
        choices=['original', 'assembly'],
        help="Absolute path of output table."
        )
    parser.add_argument(
        "--sep", 
        type=str, 
        nargs='?', 
        choices=['\t', ','], 
        default='\t',
        help="Separator for the output table. (default: %(default)s)"
        )
    parser.add_argument(
        '--version', 
        action='version', 
        version='%(prog)s 3.0'
        )
    
    args = parser.parse_args()
    
    try:
        gtf = args.gtf
        gff3 = args.gff3
        fasta = args.fasta
        tab = args.tab[0]
        mode = args.mode[0]
        sep = args.sep
    except:
        print("\nERROR: You have inserted a wrong argument or you are missing an argument.\n")
        parser.print_help()
        sys.exit()
    
    if gff3 == None and gtf == None:
        print("\nERROR: You have to insert either gtf argument or gff3 argument. One of them are required.\n")
        parser.print_help()
        sys.exit()
    if gff3 != None and gtf != None:
        print("\nERROR: You can't insert both annotation files: gtf and gff3 argument. Use only one of them\n")
        parser.print_help()
        sys.exit()
    
    # gtf = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP1/Original_genes/ORIGINAL_GENES.gtf"
    # gff3 = None
    # fasta = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP1/Original_genes/ORIGINAL_GENES.fasta"
    # tab = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP1/Original_genes/ORIGINAL_GENES.tsv"
    # mode = "original"
    # sep = "\t"
    
    
    ### PIPELINE
    ## Get the information from annotation file.
    if gtf != None:
        if gtf.endswith("gtf"):
            dict_annot = GTFtoDICT (gtf, mode)
        else:
            print("\nERROR (GTF): Check the annotation file name.\n")
            sys.exit()
            
    if gff3 != None:
        if gff3.endswith("gff") or gff3.endswith("gff3"):
            dict_annot = GFF3toDICT (gff3, mode)
        else:
            print("\nERROR (GFF3): Check the annotation file name.\n")
            sys.exit()
                
    ## Get the information from fasta file.
    if fasta != None:
        if fasta.endswith(".fas") or fasta.endswith(".fa") or fasta.endswith(".fasta"):
            dict_fasta = FASTAtoDICT (fasta)
            flag = "FASTA"
            dict_combined = CombineDICTSsameKEY (dict_annot, dict_fasta)
        else:
            print("\nERROR (FASTA): Check the sequence file name.\n")
            sys.exit()
    else:
        flag = "NON-FASTA"
    
    ## Build the final table.
    if flag == "FASTA":
        Write_table (tab, dict_combined, sep, flag)
    elif flag == "NON-FASTA":
        Write_table (tab, dict_annot, sep, flag)
    

# CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call the main program.
    """
    main()
