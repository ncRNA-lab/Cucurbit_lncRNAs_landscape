################################################################################
#
# Extract transcript id from the GTF annotation file and create a list.
#
################################################################################

rm(list = ls())

## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Extract_transcript_ids_from_GTF.R gtf_sorted.gtf ids_sorted.txt

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

# path_gtf_sorted = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/03-gene_and_lncRNAs_lists/High/intergenic/cme_def_sorted.gtf"
# path_ids_sorted = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/03-gene_and_lncRNAs_lists/High/intergenic/cme_ids_sorted.txt"


## 2. CREATE LISTS

### 3.1 SPECIE 1
z = import(path_gtf_sorted)
z = as.data.frame(z)
transcript_id = z[z$type == "transcript", c("transcript_id")]
write.table(transcript_id, path_ids_sorted, sep = "\n", row.names = F, col.names = F, quote = F)

