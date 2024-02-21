
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


######### PIPELINE

## STEP 4: Create final table.

files = list.files(path = WD_corr_S3, pattern = "-PEARSON.tsv")
TAB = data.frame()
for (file in files) {
  ID = unlist(strsplit(file, "-"))[2]
  cat(paste0(ID, "...\n"))
  tab = read.table(paste0(WD_corr_S3, "/", file), header = T, sep = "\t", quote = "\"")
  TAB = rbind(TAB, tab)
}

write.table(TAB, paste0(WD_corr_S4, "/TAB_CIS-PEARSON.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("files", "file", "ID", "tab", "TAB"))

