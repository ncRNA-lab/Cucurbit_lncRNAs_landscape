################################################################################
#
# CREATE THE BATCH TABLES COMING FROM EACH PERMUTATION
#
################################################################################


## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Batch_tables.R $WD2 $WD3 $n_batches $range $permutation $n_samples_total

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("dplyr"))


## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("At least 6 arguments must be supplied.", call.=FALSE)
} else {
  WD2 = args[1]
  WD3 = args[2]
  n_batches = as.integer(args[3])
  range = as.integer(args[4])
  permutation = as.integer(args[5])
  n_samples_total = as.integer(args[6])
}

# WD2 = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/13-saturation_plots/cme/Predictions/range_5-1"
# WD3 = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/13-saturation_plots/cme/Tables_CS"
# n_batches = 77
# range = 5
# permutation = 1
# n_samples_total = 383


## 2. TABLES

cat("TABLES...\n")

Batch_tab = data.frame()
for (i in 1:n_batches) {
  DB = read.table(paste0(WD2, "/04-Predict_lncRNAs/STEP-FINAL/Batch-", i, "/Database/Database_LncRNAs_NR.tsv"), sep = "\t", header = T, quote = "\"")
  DB$"Batch" = i
  if (i == n_batches) {
    DB$"Size" = n_samples_total
  } else {
    DB$"Size" = i * range
  }
  DB$"Permutation" = permutation
  DB = DB[,c("ID_transcript", "Class_code", "Significance_level", "Batch", "Size", "Permutation")]
  colnames(DB) = c("ID_transcript", "Class_code", "Confidence_level", "Batch", "Size", "Permutation")
  DB[DB == "u"] = "intergenic (u)"
  DB[DB == "x"] = "antisense (x)"
  DB[DB == "i"] = "intronic (i)"
  DB[DB == "o"] = "sense (o/e)"
  DB[DB == "e"] = "sense (o/e)"
  Batch_tab = rbind(Batch_tab, DB)
}

Batch_tab$Class_code = factor(Batch_tab$Class_code, levels = c("intergenic (u)", "antisense (x)", "intronic (i)", "sense (o/e)"))
Batch_tab$Confidence_level = factor(Batch_tab$Confidence_level, levels = c("Low", "Medium", "High"))

cat("\nBatch table 1...\n")
Batch_tab_1 = Batch_tab %>% group_by(Permutation, Batch, Size, Confidence_level, Class_code, .drop=FALSE) %>% summarise(Counts = n_distinct(ID_transcript))
Batch_tab_1 = as.data.frame(Batch_tab_1)
Batch_tab_1 = Batch_tab_1[order(Batch_tab_1$Batch),]
write.table(Batch_tab_1, paste0(WD3, "/Batch_tab_1-", range, "-", permutation, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

cat("\nBatch table 2...\n")
Batch_tab_2 = Batch_tab %>% group_by(Permutation, Batch, Size, Class_code, .drop=FALSE) %>% summarise(Counts = n_distinct(ID_transcript))
Batch_tab_2 = as.data.frame(Batch_tab_2)
Batch_tab_2 = Batch_tab_2[order(Batch_tab_2$Batch),]
write.table(Batch_tab_2, paste0(WD3, "/Batch_tab_2-", range, "-", permutation, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("DB", "i"))


