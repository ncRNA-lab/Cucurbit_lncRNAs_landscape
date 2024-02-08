################################################################################
#
# STEP 1: CREATE INITIAL TABLES
#
# This script creates two tables which show the number of bases in each category 
# (LncRNAs, PCGs and IR) that overlap with any type of repetitive element (Table 1) 
# or only transposons (Table 2). There will be sequences that will not overlap 
# with anything and others that may overlap with several different repetitive 
# elements.
#
# @author: pasviber - Pascual Villalba Bermell
#
################################################################################

rm(list = ls())


## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggExtra))
suppressMessages(options(bitmapType='cairo'))
options(stringsAsFactors = F)

## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("At least 6 arguments must be supplied.", call.=FALSE)
} else {
  spes = args[1]
  path_res = args[2]
  path_RepMask= args[3]
  path_pred = args[4]
  flag = args[5]
  confidence = args[6]
}

# spes = "cme"
# path_res = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-TEs_and_genomic_repeats/cme/02-Intersection"
# path_RepMask = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-TEs_and_genomic_repeats/cme/01-Repeat_calling/02-RepeatMasker"
# path_pred = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme"
# flag = "nr"
# confidence = "High"


## 2. PIPELINE

### 2.1 CREATE FINAL TABLE

### Directory.
if (!dir.exists(paste0(path_res, "/Final_tables"))){
  dir.create(paste0(path_res, "/Final_tables"))
}

### Genes, LncRNAs and IR intersected with Repeat regions.
res_genes = read.table(paste0(path_res, "/Intersection/ORIGINAL_GENES_intersect_Rep.tsv"), header = F, sep = "\t", quote = "\"")
colnames(res_genes) = c("Chr", "Start", "End", "ID_transcript", "Strand", "Start_rep", "End_rep", "Repeat_id", "Repeat_type_1", "Repeat_type_2", "Overlap")
res_lncRNAs = read.table(paste0(path_res, "/Intersection/POTENTIAL_LNCRNAS_intersect_Rep_", flag, ".tsv"), header = F, sep = "\t", quote = "\"")
colnames(res_lncRNAs) = c("Chr", "Start", "End", "ID_transcript", "Strand", "Start_rep", "End_rep", "Repeat_id", "Repeat_type_1", "Repeat_type_2", "Overlap")
res_IR = read.table(paste0(path_res, "/Intersection/Random_IR_intersect_Rep.tsv"), header = F, sep = "\t", quote = "\"")
colnames(res_IR) = c("Chr", "Start", "End", "ID_transcript", "Strand", "Start_rep", "End_rep", "Repeat_id", "Repeat_type_1", "Repeat_type_2", "Overlap")
res_IR$ID_transcript = paste0(res_IR$Chr, ":", res_IR$Start, "-", res_IR$End, "(", res_IR$Strand, ")")
RES = rbind(res_genes, res_lncRNAs, res_IR)

### Remove intersections with small RNAs.
RES = RES[!(RES$Repeat_type_2 %in% c("rRNA", "tRNA", "snRNA", "snoRNA")),]

### All the genes, lncRNAs and IR.
All_genes = read.table(paste0(path_res, "/Intersection/ORIGINAL_GENES.bed"), header = F, sep = "\t", quote = "\"")
All_genes = All_genes[,c(1:4,6)]
colnames(All_genes) = c("Chr", "Start", "End", "ID_transcript", "Strand")
All_genes$"Type" = "Protein Coding"

All_lncRNAs = read.table(paste0(path_res, "/Intersection/POTENTIAL_LNCRNAS_pred_", flag, ".bed"), header = F, sep = "\t", quote = "\"")
All_lncRNAs = All_lncRNAs[,c(1:4,6)]
colnames(All_lncRNAs) = c("Chr", "Start", "End", "ID_transcript", "Strand")
All_lncRNAs$"Type" = "LncRNA"
LncRNAs_db = read.table(paste0(path_pred, "/STEP-FINAL/Database/Database_LncRNAs_", flag, ".tsv"), header = T, sep = "\t", quote = "\"")
LncRNAs_db = LncRNAs_db[LncRNAs_db$Confidence == confidence,]
All_lncRNAs = All_lncRNAs[All_lncRNAs$ID_transcript %in% LncRNAs_db$ID_transcript,]

All_IR = read.table(paste0(path_res, "/Intersection/Random_IR.bed"), header = F, sep = "\t", quote = "\"")
All_IR = All_IR[,c(1:4,6)]
colnames(All_IR) = c("Chr", "Start", "End", "ID_transcript", "Strand")
All_IR$ID_transcript = paste0(All_IR$Chr, ":", All_IR$Start, "-", All_IR$End, "(", All_IR$Strand, ")")
All_IR$"Type" = "Intergenic Region"
ALL = rbind(All_genes, All_lncRNAs, All_IR)

### Join all the info ALL and RES. 
TAB = merge(ALL, RES[,c("ID_transcript", "Repeat_id", "Start_rep", "End_rep", "Repeat_type_1", "Repeat_type_2", "Overlap")], by = "ID_transcript", all.x = T, all.y = F)

### Convert all NA values to No_overlap except the overlap (by 0). They are the genes, lncRNAs and IR which haven't overlap on some repeat region.
TAB$Overlap = ifelse(is.na(TAB$Overlap), 0, TAB$Overlap)
TAB$Repeat_type_1 = ifelse(is.na(TAB$Repeat_type_1), "No repeat", TAB$Repeat_type_1)
TAB$Repeat_type_2 = ifelse(is.na(TAB$Repeat_type_2), "No repeat", TAB$Repeat_type_2)

rm(list = c("ALL", "All_genes", "All_IR", "All_lncRNAs", "RES", "res_genes", "res_IR", "res_lncRNAs"))

### Modify the Repeat_type_2 values to make easier the downstream analysis.
TAB1 = TAB
TAB1$"Repeat_type_2_mod" = ifelse(grepl("LINE", TAB1$Repeat_type_2, fixed = T), "LINE/SINE",
                                  ifelse(grepl("SINE", TAB1$Repeat_type_2, fixed = T), "LINE/SINE",
                                         ifelse(grepl("LTR", TAB1$Repeat_type_2, fixed = T), "LTR",
                                                ifelse(grepl("DNA/", TAB1$Repeat_type_2, fixed = T), "DNA transposon",
                                                       ifelse(grepl("RC/", TAB1$Repeat_type_2, fixed = T), "DNA transposon",
                                                              ifelse(TAB1$Repeat_type_2 == "Simple_repeat", "Simple repeat",
                                                                     ifelse(TAB1$Repeat_type_2 == "Low_complexity", "Low complexity repeat",
                                                                            ifelse(TAB1$Repeat_type_2 == "No repeat", "No repeat",
                                                                                   "Unknown"))))))))

TAB2 = TAB
TAB2$"Repeat_type_2_mod" = ifelse(grepl("LINE", TAB2$Repeat_type_2, fixed = T), "Transposon",
                                  ifelse(grepl("SINE", TAB2$Repeat_type_2, fixed = T), "Transposon",
                                         ifelse(grepl("LTR", TAB2$Repeat_type_2, fixed = T), "Transposon",
                                                ifelse(grepl("DNA/", TAB2$Repeat_type_2, fixed = T), "Transposon",
                                                       ifelse(grepl("RC/", TAB2$Repeat_type_2, fixed = T), "Transposon",
                                                              ifelse(grepl("Retroposon", TAB2$Repeat_type_2, fixed = T), "Transposon",
                                                                     ifelse(grepl("PLE/", TAB2$Repeat_type_2, fixed = T), "Transposon",
                                                                            "No transposon")))))))
TAB2 = TAB2[TAB2$Repeat_type_2_mod == "Transposon",]

### Add class code info to the table. u (intergenic), x (antisense), i (intronic), o/e (sense), pc (protein coding) and ir (intergenic region).
TAB1 = merge(TAB1, LncRNAs_db[,c("ID_transcript", "Class_code")], by = "ID_transcript", all.x = T, all.y = F)
TAB1$Class_code = ifelse(TAB1$Type == "Protein Coding", "pc", ifelse(TAB1$Type == "Intergenic Region", "ir", TAB1$Class_code))
TAB1[TAB1 == "o"] = "o/e"
TAB1[TAB1 == "e"] = "o/e"

TAB2 = merge(TAB2, LncRNAs_db[,c("ID_transcript", "Class_code")], by = "ID_transcript", all.x = T, all.y = F)
TAB2$Class_code = ifelse(TAB2$Type == "Protein Coding", "pc", ifelse(TAB2$Type == "Intergenic Region", "ir", TAB2$Class_code))
TAB2[TAB2 == "o"] = "o/e"
TAB2[TAB2 == "e"] = "o/e"

### Create factor classes.
TAB1$Class_code = factor(TAB1$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB1$Repeat_type_2_mod = factor(TAB1$Repeat_type_2_mod, levels = c("Simple repeat", "Low complexity repeat", "LINE/SINE", "LTR", "DNA transposon", "Unknown", "No repeat"))
TAB1 = TAB1[order(TAB1$Class_code),]

TAB2$Class_code = factor(TAB2$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB2$Repeat_type_2_mod = factor(TAB2$Repeat_type_2_mod, levels = c("Transposon"))
TAB2 = TAB2[order(TAB2$Class_code),]

### Add specie info.
TAB1$"Spe" = spes
TAB2$"Spe" = spes

### Order columns.
TAB1 = TAB1[,c("Spe", "Class_code", "Chr", "Strand", "ID_transcript", "Start", "End", "Repeat_id", "Start_rep", 
               "End_rep",  "Repeat_type_1", "Repeat_type_2", "Repeat_type_2_mod", "Overlap")]
TAB2 = TAB2[,c("Spe", "Class_code", "Chr", "Strand", "ID_transcript", "Start", "End", "Repeat_id", "Start_rep", 
               "End_rep",  "Repeat_type_1", "Repeat_type_2", "Repeat_type_2_mod", "Overlap")]

rownames(TAB1) = NULL
rownames(TAB2) = NULL

### Save TAB_FINAL.
write.table(TAB1, paste0(path_res, "/Final_tables/Final_tab-Repeat-", flag, "-", confidence, ".tsv"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(TAB2, paste0(path_res, "/Final_tables/Final_tab-Transposon-", flag, "-", confidence, ".tsv"), col.names = T, row.names = F, sep = "\t", quote = F)

rm(list = c("LncRNAs_db"))

