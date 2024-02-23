################################################################################
#
# STEP4: JOIN CIS CORRELATION TABLES (CLOSEST)
#
# For NAT-lncRNAs, join the cis-correlation tables generated for each subexperiment 
# in step 3.
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
  WD_corr_S3 = args[1]
  WD_corr_S4 = args[2]
}

# WD_corr_S3 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/antisense/nr/STEP3"
# WD_corr_S4 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/antisense/nr/STEP4"


######### PIPELINE

cat(paste0("-Closest...\n"))

files = list.files(path = WD_corr_S3, pattern = "_closest.tsv")
TAB = data.frame()
for (file in files) {
  ID = unlist(strsplit(file, "-"))[2]
  cat(paste0("\t-", ID, "...\n"))
  tab = read.table(paste0(WD_corr_S3, "/", file), header = T, sep = "\t", quote = "\"")
  TAB = rbind(TAB, tab)
}

write.table(TAB, paste0(WD_corr_S4, "/TAB_CIS-PEARSON_closest.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("files", "file", "ID", "tab", "TAB"))

