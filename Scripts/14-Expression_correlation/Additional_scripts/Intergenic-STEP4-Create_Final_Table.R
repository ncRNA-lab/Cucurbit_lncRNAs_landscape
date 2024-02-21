
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

################################################################################

## STEP 4 (Closest): Create final table.

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


################################################################################

## STEP 4 (Range): Create final table.

cat(paste0("-Range...\n"))

files = list.files(path = WD_corr_S3, pattern = "_range-")
TAB = data.frame()
for (file in files) {
  ID = unlist(strsplit(file, "-"))[2]
  dist = gsub(".tsv", "", unlist(strsplit(file, "-"))[4])
  cat(paste0("\t-", ID, "-", dist, "...\n"))
  tab = read.table(paste0(WD_corr_S3, "/", file), header = T, sep = "\t", quote = "\"")
  TAB = rbind(TAB, tab)
}

TAB$Range = as.integer(format(TAB$Range, scientific = F))

write.table(TAB, paste0(WD_corr_S4, "/TAB_CIS-PEARSON_range.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("files", "file", "ID", "dist", "tab", "TAB"))

