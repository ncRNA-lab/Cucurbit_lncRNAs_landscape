################################################################################
#
# Change headers proteome file in cme. Convert Name to ID.
#
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES
suppressMessages(library(rtracklayer))
suppressMessages(library(Biostrings))

## 1. VARIABLES
spe = "cme"
path_g = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Additional_info/GFF3_genes"
path_p = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Additional_info/Proteomes"

## 2. PIPELINE

### 2.1 Load annotation gff3 file.
z = import(paste0(path_g, "/", spe, ".gff3"))
z = as.data.frame(z)
dt = z[z$type=="mRNA", c("ID", "Name")]
rownames(dt) = NULL
dt$Name = gsub(".2.", ".", dt$Name, fixed = T)

### 2.2 Load the proteome fasta file.
s = readAAStringSet(paste0(path_p, "/security/", spe, ".fa"))
seq_name = names(s)
sequence = paste(s)
df = data.frame(seq_name, sequence)

### 2.3 Merge the info and change the header.
DT = merge(df, dt, by.x = "seq_name", by.y = "Name", all = T)
DT$"ID_mod" = paste0(">", DT$ID)
DT$seq_name = NULL
DT$ID = NULL
DT = DT[,c("ID_mod", "sequence")]
write.table(DT, paste0(path_p, "/", spe, ".fa"), col.names = F, row.names = F, sep = "\n", quote = F)
