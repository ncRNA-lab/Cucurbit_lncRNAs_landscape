################################################################################
#
# Percentage of nucleotides by transcript coming from some kind of transposon or 
# repetitive element that can be found in the genome of each specie. We take always 
# into account the class code of each transcript.
#
# STEP 1
#
################################################################################

rm(list = ls())

## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Create_Tables_and_Figures_C-STEP1.R spes path_res path_RepMask path_pred flag confidence

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

# spes = "vvi"
# path_res = "/mnt/doctorado/3-lncRNAs/Vitis_Tom/Results/08-TEs_and_genomic_repeats/vvi/02-Comparison_PCGs_LncRNAs"
# path_RepMask = "/mnt/doctorado/3-lncRNAs/Vitis_Tom/Results/08-TEs_and_genomic_repeats/vvi/01-Repeat_calling/02-RepeatMasker"
# path_pred = "/mnt/doctorado/3-lncRNAs/Vitis_Tom/Results/05-LncRNAs_prediction/vvi"
# flag = "NR"
# confidence = "High"

# spes = "vvi"
# path_res = "/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results/08-TEs_and_genomic_repeats/vvi/02-Comparison_PCGs_LncRNAs"
# path_RepMask = "/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results/08-TEs_and_genomic_repeats/vvi/01-Repeat_calling/02-RepeatMasker"
# path_pred = "/storage/ncRNA/Projects/lncRNAs/Vitis_Tom/Results/05-LncRNAs_prediction/vvi"
# flag = "NR"
# confidence = "High"


## 2. PIPELINE

### 2.1 CREATE FINAL TABLE

### Directory.
if (!dir.exists(path_res)){
  dir.create(path_res)
}
if (!dir.exists(paste0(path_res, "/Final_tables"))){
  dir.create(paste0(path_res, "/Final_tables"))
}

### Genes, LncRNAs and IR intersected with Repeat regions.
res_genes = read.table(paste0(path_res, "/ORIGINAL_GENES_intersect_Rep_19.tsv"), header = F, sep = "\t", quote = "\"")
colnames(res_genes) = c("chr", "start", "end", "transcript_id", "strand", "start_rep", "end_rep", "Repeat_id", "Repeat_type_1", "Repeat_type_2", "overlap")
if (flag == "NR") {
  res_lncRNAs = read.table(paste0(path_res, "/POTENTIAL_LNCRNAS_intersect_Rep_NR_19.tsv"), header = F, sep = "\t", quote = "\"")
  colnames(res_lncRNAs) = c("chr", "start", "end", "transcript_id", "strand", "start_rep", "end_rep", "Repeat_id", "Repeat_type_1", "Repeat_type_2", "overlap")
}
if (flag == "R") {
  res_lncRNAs = read.table(paste0(path_res, "/POTENTIAL_LNCRNAS_intersect_Rep_R_19.tsv"), header = F, sep = "\t", quote = "\"")
  colnames(res_lncRNAs) = c("chr", "start", "end", "transcript_id", "strand", "start_rep", "end_rep", "Repeat_id", "Repeat_type_1", "Repeat_type_2", "overlap")
}
res_IR = read.table(paste0(path_res, "/Random_IR_intersect_Rep_15.tsv"), header = F, sep = "\t", quote = "\"")
colnames(res_IR) = c("chr", "start", "end", "transcript_id", "strand", "start_rep", "end_rep", "Repeat_id", "Repeat_type_1", "Repeat_type_2", "overlap")
res_IR$transcript_id = paste0(res_IR$chr, ":", res_IR$start, "-", res_IR$end, "(", res_IR$strand, ")")
RES = rbind(res_genes, res_lncRNAs, res_IR)

### Remove intersections with small RNAs.
RES = RES[!(RES$Repeat_type_2 %in% c("rRNA", "tRNA", "snRNA", "snoRNA")),]

### All the genes, lncRNAs and IR.
All_genes = read.table(paste0(path_res, "/ORIGINAL_GENES.bed"), header = F, sep = "\t", quote = "\"")
All_genes = All_genes[,c(1:4,6)]
colnames(All_genes) = c("chr", "start", "end", "transcript_id", "strand")
All_genes$"Type" = "Protein Coding"
if (flag == "NR") {
  All_lncRNAs = read.table(paste0(path_res, "/POTENTIAL_LNCRNAS_pred_NR.bed"), header = F, sep = "\t", quote = "\"")
  All_lncRNAs = All_lncRNAs[,c(1:4,6)]
  colnames(All_lncRNAs) = c("chr", "start", "end", "transcript_id", "strand")
  All_lncRNAs$"Type" = "LncRNA"
  LncRNAs_db = read.table(paste0(path_pred, "/STEP-FINAL/Database/Database_LncRNAs_NR.tsv"), header = T, sep = "\t", quote = "\"")
  LncRNAs_db = LncRNAs_db[LncRNAs_db$Confidence_level == confidence,]
  All_lncRNAs = All_lncRNAs[All_lncRNAs$transcript_id %in% LncRNAs_db$ID_transcript,]
}
if (flag == "R") {
  All_lncRNAs = read.table(paste0(path_res, "/POTENTIAL_LNCRNAS_pred_R.bed"), header = F, sep = "\t", quote = "\"")
  All_lncRNAs = All_lncRNAs[,c(1:4,6)]
  colnames(All_lncRNAs) = c("chr", "start", "end", "transcript_id", "strand")
  All_lncRNAs$"Type" = "LncRNA"
  LncRNAs_db = read.table(paste0(path_pred, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), header = T, sep = "\t", quote = "\"")
  LncRNAs_db = LncRNAs_db[LncRNAs_db$Confidence_level == confidence,]
  All_lncRNAs = All_lncRNAs[All_lncRNAs$transcript_id %in% LncRNAs_db$ID_transcript,]
}
All_IR = read.table(paste0(path_res, "/Random_IR.bed"), header = F, sep = "\t", quote = "\"")
All_IR = All_IR[,c(1:4,6)]
colnames(All_IR) = c("chr", "start", "end", "transcript_id", "strand")
All_IR$transcript_id = paste0(All_IR$chr, ":", All_IR$start, "-", All_IR$end, "(", All_IR$strand, ")")
All_IR$"Type" = "Intergenic Region"
ALL = rbind(All_genes, All_lncRNAs, All_IR)

### Join all the info ALL and RES. 
TAB = merge(ALL, RES[,c("transcript_id", "Repeat_id", "start_rep", "end_rep", "Repeat_type_1", "Repeat_type_2", "overlap")], by = "transcript_id", all.x = T, all.y = F)

### Convert all NA values to No_overlap except the overlap (by 0). They are the genes, lncRNAs and IR which haven't overlap on some repeat region.
TAB$overlap = ifelse(is.na(TAB$overlap), 0, TAB$overlap)
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
TAB1 = merge(TAB1, LncRNAs_db[,c("ID_transcript", "Class_code")], by.x = "transcript_id", by.y = "ID_transcript", all.x = T, all.y = F)
TAB1$Class_code = ifelse(TAB1$Type == "Protein Coding", "pc", ifelse(TAB1$Type == "Intergenic Region", "ir", TAB1$Class_code))
TAB1[TAB1 == "o"] = "o/e"
TAB1[TAB1 == "e"] = "o/e"

TAB2 = merge(TAB2, LncRNAs_db[,c("ID_transcript", "Class_code")], by.x = "transcript_id", by.y = "ID_transcript", all.x = T, all.y = F)
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
TAB1$"spe" = spes
TAB2$"spe" = spes

### Order columns
TAB1 = TAB1[,c("spe", "Class_code", "chr", "strand", "transcript_id", "start", "end", "Repeat_id", "start_rep", 
               "end_rep",  "Repeat_type_1", "Repeat_type_2", "Repeat_type_2_mod", "overlap")]
TAB2 = TAB2[,c("spe", "Class_code", "chr", "strand", "transcript_id", "start", "end", "Repeat_id", "start_rep", 
               "end_rep",  "Repeat_type_1", "Repeat_type_2", "Repeat_type_2_mod", "overlap")]

rownames(TAB1) = NULL
rownames(TAB2) = NULL

### Save TAB_FINAL.
write.table(TAB1, paste0(path_res, "/Final_tables/Final_tab-Repeat-", flag, "-", confidence, ".tsv"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(TAB2, paste0(path_res, "/Final_tables/Final_tab-Transposon-", flag, "-", confidence, ".tsv"), col.names = T, row.names = F, sep = "\t", quote = F)

rm(list = c("LncRNAs_db"))
