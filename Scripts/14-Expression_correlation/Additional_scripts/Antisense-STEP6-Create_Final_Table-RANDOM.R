
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


######### PIPELINE

## STEP 6: Create final table.

files = list.files(path = WD_corr_S5, pattern = "-random.tsv")
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

write.table(TAB, paste0(WD_corr_S6, "/TAB_CIS-PEARSON-random.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("files", "file", "ID", "tab", "TAB"))

