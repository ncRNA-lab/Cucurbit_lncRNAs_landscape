################################################################################
#
# STEP2: CREATE CIS INTERACTION TABLE (CLOSEST AND RANGE)
#
# Modify the table previously generated in step 1 to create the final table with 
# the cis-interactions of all lincRNAs. In this case, the table shows lncRNA-PCG 
# and PCG-PCG interactions. In addition, two strategies are shown: (1) closest 
# gene and (2) distance ranges.
#
# NOTE: It may occur that two PCGs have mutually closest gene to each other. Then 
# the repeated pairs are eliminated leaving unique pairs.
#
# @author: pasviber - Pascual Villalba Bermell
# 
################################################################################


######### MODULES

suppressMessages(library(dplyr))


######### VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("At least 6 arguments must be supplied.", call.=FALSE)
} else {
  WD_corr_S1 = args[1]
  WD_corr_S2 = args[2]
  WD_pred = args[3]
  distances = as.numeric(unlist(strsplit(args[4], " ")))
  spel = args[5]
  flag = args[6]
}

# WD_corr_S1 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/intergenic/nr/STEP1"
# WD_corr_S2 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/intergenic/nr/STEP2"
# WD_pred = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme"
# distances = c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
# spel = "C. melo"
# flag = "nr"


######### PIPELINE

################################################################################

## Load lncRNAs database.
DB_LncRNAs = read.table(paste0(WD_pred, "/STEP-FINAL/Database/Database_LncRNAs_", toupper(flag), ".tsv"), sep = "\t", header = T, quote = "\"")
DB_LncRNAs = DB_LncRNAs[, c("ID_transcript", "Start", "End")]

## Load Genes database.
DB_Genes = read.table(paste0(WD_pred, "/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv"), sep = "\t", header = T, quote = "\"")
DB_Genes = DB_Genes[, c("ID_transcript", "Start", "End")]

################################################################################

## CLOSEST: Build tab_closest table.
cat(paste0("-Closest...\n"))

# Load bedtools results: LncRNAs-Genes and Genes-Genes.
tab_closest_LG = read.table(paste0(WD_corr_S1, "/LncRNA_Gene_closest.tsv"), sep = "\t", header = F, quote = "\"")
colnames(tab_closest_LG) = c("Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", "Start.2", "End.2", "Strand.2", "Distance")
tab_closest_LG$Type = "LncRNAs_Genes"
tab_closest_GG = read.table(paste0(WD_corr_S1, "/Gene_Gene_closest.tsv"), sep = "\t", header = F, quote = "\"")
colnames(tab_closest_GG) = c("Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", "Start.2", "End.2", "Strand.2", "Distance")
tab_closest_GG$Type = "Genes_Genes"
tab_closest = rbind(tab_closest_LG, tab_closest_GG)

# Sometimes if gene 2 is the closest to gene 1, gene 1 is also the closest to gene 2.
# So, it's necessary remove redundant pairs.
tab_closest = tab_closest %>%
  group_by(grp = paste(pmax(ID_transcript.1, ID_transcript.2), pmin(ID_transcript.1, ID_transcript.2), sep = "_")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-grp)
tab_closest = as.data.frame(tab_closest)

# Fill table.
tab_closest$Start.1 = tab_closest$Start.1 + 1
tab_closest$Start.2 = tab_closest$Start.2 + 1
tab_closest$"Specie" = spel
tab_closest = tab_closest[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", "Start.2", "End.2", "Strand.2", "Distance", "Type")]

TAB_CIS_closest = tab_closest

rm(list = c("tab_closest", "tab_closest_LG", "tab_closest_GG"))

write.table(TAB_CIS_closest, paste0(WD_corr_S2, "/TAB_CIS_closest.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

################################################################################

## RANGE: Build tab_range table.
cat(paste0("-Range...\n"))

TAB_CIS_range = data.frame()

for (dist in distances) {
  
  dist = as.integer(format(dist, scientific = F))
  
  cat(paste0("\t-", dist, "...\n"))
  
  # Load bedtools results: LncRNAs-Genes and Genes-Genes.
  tab_range_LG = read.table(paste0(WD_corr_S1, "/LncRNA_Gene_cis_interactions_range_", dist, ".tsv"), sep = "\t", header = F, quote = "\"")
  colnames(tab_range_LG) = c("Chr", "ID_transcript.1", "Start_range", "End_range", "Strand.1", "ID_transcript.2", "Start.2", "End.2", "Strand.2", "Overlap_range")
  tab_range_LG = merge(tab_range_LG, DB_LncRNAs, by.x = "ID_transcript.1", by.y = "ID_transcript", all = F)
  tab_range_LG$Type = "LncRNAs_Genes"
  tab_range_GG = read.table(paste0(WD_corr_S1, "/Gene_Gene_cis_interactions_range_", dist, ".tsv"), sep = "\t", header = F, quote = "\"")
  colnames(tab_range_GG) = c("Chr", "ID_transcript.1", "Start_range", "End_range", "Strand.1", "ID_transcript.2", "Start.2", "End.2", "Strand.2", "Overlap_range")
  tab_range_GG = merge(tab_range_GG, DB_Genes, by.x = "ID_transcript.1", by.y = "ID_transcript", all = F)
  tab_range_GG$Type = "Genes_Genes"
  tab_range = rbind(tab_range_LG, tab_range_GG)
  
  # Sometimes if gene 2 is the closest to gene 1, gene 1 is also the closest to gene 2.
  # So, it's necessary remove redundant pairs.
  tab_range = tab_range %>%
    group_by(grp = paste(pmax(ID_transcript.1, ID_transcript.2), pmin(ID_transcript.1, ID_transcript.2), sep = "_")) %>%
    slice(1) %>%
    ungroup() %>%
    select(-grp)
  tab_range = as.data.frame(tab_range)
  
  # Fill table.
  tab_range$Start_range = tab_range$Start_range + 1
  tab_range$Start.2 = tab_range$Start.2 + 1
  tab_range$Range = dist
  tab_range$Specie = spel
  tab_range$Overlap_range.Perc.2 = round((tab_range$Overlap_range * 100)/((tab_range$End.2 - tab_range$Start.2) + 1), 2)
  
  tab_range = tab_range[, c("Specie", "Chr", "ID_transcript.1", "Start", "End", "Strand.1", "ID_transcript.2", "Start.2", "End.2", 
                            "Strand.2", "Overlap_range.Perc.2", "Overlap_range", "Range", "Start_range", "End_range", "Type")]
  colnames(tab_range) = c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", "Start.2", "End.2", 
                          "Strand.2", "Overlap_range.Perc.2", "Overlap_range", "Range", "Start_range", "End_range", "Type")
  
  TAB_CIS_range = rbind(TAB_CIS_range, tab_range)
  
  rm(list = c("tab_range", "tab_range_LG", "tab_range_GG"))
}

rm(list = c("dist", "DB_LncRNAs", "DB_Genes"))

write.table(TAB_CIS_range, paste0(WD_corr_S2, "/TAB_CIS_range.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

