################################################################################
#
# JOIN THE BATCH TABLES COMING FROM EACH PERMUTATION
#
################################################################################


## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Batch_tables_mean_permutations.R $WD3 $range $permutations

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("dplyr"))


## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied.", call.=FALSE)
} else {
  WD3 = args[1]
  range = as.integer(args[2])
  permutations = as.integer(args[3])
}

# WD3 = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/13-saturation_plots/car/Tables_CS"
# range = 2
# permutations = 4


## 2. TABLES

cat("FINAL TABLES...\n")

cat("\nMean permutations: FINAL Batch table 1...\n")

Batch_tab_1_ALL = data.frame()
for (permutation in 1:permutations) {
  Batch_tab_1 = read.table(paste0(WD3, "/Batch_tab_1-", range, "-", permutation, ".tsv"), sep = "\t", header = T, quote = "\"")
  Batch_tab_1_ALL = rbind(Batch_tab_1_ALL, Batch_tab_1)
}

rm(list = c("Batch_tab_1", "permutation"))

Batch_tab_1_ALL$Class_code = factor(Batch_tab_1_ALL$Class_code, levels = c("intergenic (u)", "antisense (x)", "intronic (i)", "sense (o/e)"))
Batch_tab_1_ALL$Confidence_level = factor(Batch_tab_1_ALL$Confidence_level, levels = c("Low", "Medium", "High"))

Batch_tab_1_RED = Batch_tab_1_ALL %>% group_by(Batch, Size, Confidence_level, Class_code, .drop=FALSE) %>% summarise(Mean.Counts = mean(Counts))
Batch_tab_1_RED = as.data.frame(Batch_tab_1_RED)
Batch_tab_1_RED = Batch_tab_1_RED[order(Batch_tab_1_RED$Batch),]
write.table(Batch_tab_1_RED, paste0(WD3, "/FINAL-Batch_tab_1-", range, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


cat("\nMean permutations: FINAL Batch table 2...\n")

Batch_tab_2_ALL = data.frame()
for (permutation in 1:permutations) {
  Batch_tab_2 = read.table(paste0(WD3, "/Batch_tab_2-", range, "-", permutation, ".tsv"), sep = "\t", header = T, quote = "\"")
  Batch_tab_2_ALL = rbind(Batch_tab_2_ALL, Batch_tab_2)
}

rm(list = c("Batch_tab_2", "permutation"))

Batch_tab_2_ALL$Class_code = factor(Batch_tab_2_ALL$Class_code, levels = c("intergenic (u)", "antisense (x)", "intronic (i)", "sense (o/e)"))

Batch_tab_2_RED = Batch_tab_2_ALL %>% group_by(Batch, Size, Class_code, .drop=FALSE) %>% summarise(Mean.Counts = mean(Counts))
Batch_tab_2_RED = as.data.frame(Batch_tab_2_RED)
Batch_tab_2_RED = Batch_tab_2_RED[order(Batch_tab_2_RED$Batch),]
write.table(Batch_tab_2_RED, paste0(WD3, "/FINAL-Batch_tab_2-", range, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

