################################################################################
#
# EXTRACT TRANSCRIPT IDS FROM A GTF ANNOTATION FILE
#
# Extract the transcript ids from the GTF annotation file and create a list.
#
# @author: pasviber - Pascual Villalba Bermell
#
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(rtracklayer))

## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("At least 2 arguments must be supplied.", call.=FALSE)
} else {
  path_gtf_sorted = args[1]
  path_ids_sorted = args[2]
}

# path_gtf_sorted = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/03-gene_and_lncRNAs_lists/High/intergenic/cme_def_sorted.gtf"
# path_ids_sorted = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/03-gene_and_lncRNAs_lists/High/intergenic/cme_ids_sorted.txt"


## 2. CREATE LISTS

z = import(path_gtf_sorted)
z = as.data.frame(z)
transcript_id = z[z$type == "transcript", c("transcript_id")]
write.table(transcript_id, path_ids_sorted, sep = "\n", row.names = F, col.names = F, quote = F)

