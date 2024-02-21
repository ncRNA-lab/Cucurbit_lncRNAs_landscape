
######### MODULES

suppressMessages(library(dplyr))


######### VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments must be supplied.", call.=FALSE)
} else {
  WD_corr_S1 = args[1]
  WD_corr_S2 = args[2]
  WD_pred = args[3]
  spel = args[4]
  flag = args[5]
}



######### PIPELINE

## STEP 2: Create table with cis-targets.

# Load bedtools results: LncRNAs-Genes and Genes-Genes.
tab_cis = read.table(paste0(WD_corr_S1, "/LncRNA_Gene_cis_interactions.tsv"), sep = "\t", header = F, quote = "\"")
colnames(tab_cis) = c("Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", "Start.2", "End.2", "Strand.2", "Overlap")
tab_cis$Type = "LncRNAs_Genes"

# Modify and fill table.
tab_cis$Start.1 = tab_cis$Start.1 + 1
tab_cis$Start.2 = tab_cis$Start.2 + 1
tab_cis$Specie = spel
tab_cis$Overlap.Perc = round((tab_cis$Overlap * 100)/((tab_cis$End.2 - tab_cis$Start.2) + 1), 2)

# Sort columns.
tab_cis = tab_cis[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", "Start.2", 
                      "End.2", "Strand.2", "Overlap.Perc", "Overlap", "Type")]

# Save.
write.table(tab_cis, paste0(WD_corr_S2, "/TAB_CIS.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

