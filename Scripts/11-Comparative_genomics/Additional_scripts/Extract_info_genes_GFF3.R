################################################################################
#
# Extract gene name, strand and chromosome info from the GFF3 annotation file and
# create a list for each chromosome with the gene name and strand info.
#
################################################################################

rm(list = ls())

## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Extract_info_genes_GFF3.R spe1 spe2 path

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(rtracklayer))

## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied.", call.=FALSE)
} else {
  spe1 = args[1]
  spe2 = args[2]
  path = args[3]
}

# spe1 = "cma"
# spe2 = "cme"
# path = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sinteny/02-Adhore/PAIRWISE_SPECIES/cma-cme"

## 2. CREATE DIRECTORY

if (!dir.exists(paste0(path, "/Lists_", spe1))){
  dir.create(paste0(path, "/Lists_", spe1))
}
if (!dir.exists(paste0(path, "/Lists_", spe2))){
  dir.create(paste0(path, "/Lists_", spe2))
}

## 3. CREATE LISTS BY CHROMOSOME

### 3.1 SPECIE 1
z1 = import(paste0(path, "/", spe1, ".gff3"))
z1 = as.data.frame(z1)
dt1 = z1[z1$type == "mRNA", c("seqnames", "ID", "strand")]
rownames(dt1) = NULL
dt1$"Info" = paste0(dt1$ID, dt1$strand)
dt1$ID = NULL
dt1$strand = NULL

for (chr in unique(dt1$seqnames)) {
  dt1_chr = data.frame(Info = dt1[dt1$seqnames == chr, "Info"])
  write.table(dt1_chr, paste0(path, "/Lists_", spe1, "/", chr, ".lst"), row.names = F, col.names = F, quote = F)
}

### 3.2 SPECIE 2
z2 = import(paste0(path, "/", spe2, ".gff3"))
z2 = as.data.frame(z2)
dt2 = z2[z2$type == "mRNA", c("seqnames", "ID", "strand")]
rownames(dt2) = NULL
dt2$"Info" = paste0(dt2$ID, dt2$strand)
dt2$ID = NULL
dt2$strand = NULL

for (chr in unique(dt2$seqnames)) {
  dt2_chr = data.frame(Info = dt2[dt2$seqnames == chr, "Info"])
  write.table(dt2_chr, paste0(path, "/Lists_", spe2, "/", chr, ".lst"), row.names = F, col.names = F, quote = F)
}

