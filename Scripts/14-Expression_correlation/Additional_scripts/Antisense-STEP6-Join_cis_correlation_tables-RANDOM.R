################################################################################
#
# STEP6: JOIN CIS CORRELATION TABLES (CLOSEST-RANDOM)
#
# For NAT-lncRNAs, join the cis-correlation tables generated for each subexperiment 
# in step 5.
#
# NOTE: As they are random pairs, it may be the case that there are repeated pairs. 
# For example, gene1-gene2 and gene2-gene1. Therefore, we eliminate them. Repeated 
# pairs such as gen1-gen2 and gen1-gen2 are eliminated in step 5.
#
# @author: pasviber - Pascual Villalba Bermell
# 
################################################################################


######### MODULES

suppressMessages(library(dplyr))


######### VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("At least 2 arguments must be supplied.", call.=FALSE)
} else {
  WD_corr_S5 = args[1]
  WD_corr_S6 = args[2]
}

# WD_corr_S5 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/antisense/nr/STEP5"
# WD_corr_S6 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/antisense/nr/STEP6"


######### PIPELINE

cat(paste0("-Closest...\n"))

files = list.files(path = WD_corr_S5, pattern = "_closest-random.tsv")
TAB = data.frame()
for (file in files) {
  ID = unlist(strsplit(file, "-"))[2]
  cat(paste0("\t-", ID, "...\n"))
  tab = read.table(paste0(WD_corr_S5, "/", file), header = T, sep = "\t", quote = "\"")
  # Sometimes if gene 2 is the closest to gene 1, gene 1 is also the closest to gene 2.
  # So, it's necessary remove redundant pairs.
  tab = tab %>%
    group_by(grp = paste(pmax(ID_transcript.1, ID_transcript.2), pmin(ID_transcript.1, ID_transcript.2), sep = "_")) %>%
    slice(1) %>%
    ungroup() %>%
    select(-grp)
  TAB = rbind(TAB, tab)
}

write.table(TAB, paste0(WD_corr_S6, "/TAB_CIS-PEARSON_closest-random.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("files", "file", "ID", "tab", "TAB"))

