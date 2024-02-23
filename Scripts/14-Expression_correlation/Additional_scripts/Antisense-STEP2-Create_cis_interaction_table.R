################################################################################
#
# STEP2: CREATE CIS INTERACTION TABLE (CLOSEST)
#
# Modify the table previously generated in step 1 to create the final table with 
# the cis-interactions of all NAT-lncRNAs. The table only shows lncRNA-PCG 
# interactions. We only search for the closest gene.
#
# @author: pasviber - Pascual Villalba Bermell
# 
################################################################################


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

# WD_corr_S1 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/antisense/nr/STEP1"
# WD_corr_S2 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/antisense/nr/STEP2"
# WD_pred = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme"
# spel = "C. melo"
# flag = "nr"


######### PIPELINE

## CLOSEST: Build TAB_CIS_closest table.
cat(paste0("-Closest...\n"))

# Load bedtools results: LncRNAs-Genes.
tab_closest = read.table(paste0(WD_corr_S1, "/LncRNA_Gene_cis_interactions.tsv"), sep = "\t", header = F, quote = "\"")
colnames(tab_closest) = c("Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", "Start.2", "End.2", "Strand.2", "Overlap")
tab_closest$Type = "LncRNAs_Genes"

# Modify and fill table.
tab_closest$Start.1 = tab_closest$Start.1 + 1
tab_closest$Start.2 = tab_closest$Start.2 + 1
tab_closest$Specie = spel
tab_closest$Overlap.Perc = round((tab_closest$Overlap * 100)/((tab_closest$End.2 - tab_closest$Start.2) + 1), 2)

# Sort columns.
TAB_CIS_closest = tab_closest[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", "Start.2", 
                                  "End.2", "Strand.2", "Overlap.Perc", "Overlap", "Type")]

rm(list = c("tab_closest"))

write.table(TAB_CIS_closest, paste0(WD_corr_S2, "/TAB_CIS_closest.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

