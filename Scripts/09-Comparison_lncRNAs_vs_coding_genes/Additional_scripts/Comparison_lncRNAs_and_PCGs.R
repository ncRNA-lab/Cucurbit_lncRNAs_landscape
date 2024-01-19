################################################################################
#
# FIGURES
#
# Create plots comparing genes and lncRNAs: GC content, Length, Exon number and 
# TPMs.
#
################################################################################

rm(list = ls())

## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Comparison_lncRNAs_and_mRNAs.R $WD1 $WD2 $WD3

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggExtra))
suppressMessages(options(bitmapType='cairo'))


## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("At least 4 arguments must be supplied.", call.=FALSE)
} else {
  pred = args[1]
  quant = args[2]
  comp = args[3]
  flag = args[4]
}

# pred = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/car"
# quant = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/06-quantification/car"
# comp = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/07-comparison_lncRNAs_vs_coding_genes/car"
# flag = "nr"

if (flag == "nr") {
  DB_path = paste0(pred, "/STEP-FINAL/Database/Database_LncRNAs_NR.tsv")
}
if (flag == "r") {
  DB_path = paste0(pred, "/STEP-FINAL/Database/Database_LncRNAs.tsv")
}

Genes_tab = paste0(pred, "/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv")
LncRNAs_tab = paste0(pred, "/STEP-FINAL/Files/LncRNAs/", flag, "/POTENTIAL_LNCRNAS_pred.tsv")
Quant_tab = paste0(quant, "/Salmon/ALL/", flag, "/04-Table/TPMs_summary.tsv")
IR_tab = paste0(comp, "/Random_IR.bed")
IR_fasta = paste0(comp, "/Random_IR.fasta")


## 2. LOAD INFO AND CREATES TABLES

# LncRNAs
DB = read.table(DB_path, sep = "\t", header = T, quote = "\"")
DB = DB[,c("ID_transcript", "Significance_level")]
TPMs = read.table(Quant_tab, sep = "\t", header = T, quote = "\"")
TAB_L_FINAL = read.table(LncRNAs_tab, sep = "\t", header = T, quote = "\"")

TAB_L_FINAL = merge(TAB_L_FINAL, DB, by = "ID_transcript", all = T)
TAB_L_FINAL = merge(TAB_L_FINAL, TPMs, by = "ID_transcript", all.x = T, all.y = F)
names(TAB_L_FINAL)[names(TAB_L_FINAL) == "Significance_level"] = "Type"
TAB_L_FINAL[TAB_L_FINAL == "High"] = "LncRNAs high"
TAB_L_FINAL[TAB_L_FINAL == "Medium"] = "LncRNAs medium"
TAB_L_FINAL[TAB_L_FINAL == "Low"] = "LncRNAs low"
TAB_L_FINAL = TAB_L_FINAL[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", 
                             "Class_code", "Exons", "Length", "GC", "Type", "TPMs.mean", 
                             "log2.TPMs.1.mean", "log2.TPMs.mean")]
TAB_L_FINAL[TAB_L_FINAL == "u"] = "intergenic (u)"
TAB_L_FINAL[TAB_L_FINAL == "x"] = "antisense (x)"
TAB_L_FINAL[TAB_L_FINAL == "i"] = "intronic (i)"
TAB_L_FINAL[TAB_L_FINAL == "o"] = "sense (o/e)"
TAB_L_FINAL[TAB_L_FINAL == "e"] = "sense (o/e)"
TAB_L_FINAL$Class_code = factor(TAB_L_FINAL$Class_code, levels = c("intergenic (u)", "antisense (x)", "intronic (i)", "sense (o/e)"))
TAB_L_FINAL$Type = factor(TAB_L_FINAL$Type, levels = c("LncRNAs low", "LncRNAs medium", "LncRNAs high"))

# Genes
TPMs = read.table(Quant_tab, sep = "\t", header = T, quote = "\"")
TAB_G_FINAL = read.table(Genes_tab, sep = "\t", header = T)

TAB_G_FINAL = merge(TAB_G_FINAL, TPMs, by = "ID_transcript", all.x = T, all.y = F)
TAB_G_FINAL$"Type" = "Genes"
TAB_G_FINAL$Class_code = "gene (=)"
TAB_G_FINAL = TAB_G_FINAL[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", 
                             "Class_code", "Exons", "Length", "GC", "Type", "TPMs.mean",
                             "log2.TPMs.1.mean", "log2.TPMs.mean")]

# IR
TAB_IR = read.table(IR_tab, sep = "\t", header = F)
TAB_IR$V5 = NULL
colnames(TAB_IR) = c("Chr", "Start", "End", "ID_transcript", "Strand")
TAB_IR$"Origin" = "bedtools"
TAB_IR$ID_transcript = paste0(TAB_IR$Chr, ":", TAB_IR$Start, "-", TAB_IR$End, "(", TAB_IR$Strand, ")")
FASTA_IR = read.table(IR_fasta, sep = "\t", header = F)
id = FASTA_IR$V1[grepl(">", FASTA_IR$V1, fixed = T)]
id = gsub(">", "", id)
seq = FASTA_IR$V1[!grepl(">", FASTA_IR$V1, fixed = T)]
tab = data.frame(ID_transcript = id, Seq = seq, Length = nchar(seq))
tab$"GC" = (nchar(gsub("[AT]", "", tab$Seq)) * 100)/tab$Length
tab$Seq = NULL
TAB_IR_FINAL = merge(TAB_IR, tab, by = "ID_transcript", all = T)
TAB_IR_FINAL$"Type" = "Intergenic Regions"
TAB_IR_FINAL = TAB_IR_FINAL[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", 
                               "Length", "GC", "Type")]

rm(list = c("TAB_IR", "FASTA_IR", "id", "seq", "tab"))


# Joined: LncRNAs and Genes
TAB_FINAL_JOIN_1 = rbind(TAB_G_FINAL, TAB_L_FINAL)
TAB_FINAL_JOIN_1$Type = factor(TAB_FINAL_JOIN_1$Type, levels = c("Genes", "LncRNAs low", "LncRNAs medium", "LncRNAs high"))
write.table(TAB_FINAL_JOIN_1, paste0(comp, "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

# Joined: LncRNAs, Genes and IR
TAB_FINAL_JOIN_2 = rbind(TAB_FINAL_JOIN_1[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Length", "GC", "Type")], TAB_IR_FINAL)
TAB_FINAL_JOIN_2$Type = factor(TAB_FINAL_JOIN_2$Type, levels = c("Genes", "LncRNAs low", "LncRNAs medium", "LncRNAs high", "Intergenic Regions"))
write.table(TAB_FINAL_JOIN_2, paste0(comp, "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("Genes_tab", "LncRNAs_tab", "DB", "IR_fasta", "IR_tab", "DB_path", "pred", "quant", "Quant_tab", "TPMs"))




## 3. GC CONTENT

gg = ggplot(TAB_FINAL_JOIN_2, aes(x = Type, y = GC, fill = Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  coord_cartesian(ylim = c(25, 55)) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5)) + 
  theme(legend.position = "none")

ggsave(paste0(comp, "/GC_content_by_type_", toupper(flag), ".png"), height = 7, width = 8, dpi = 600)

gg = ggplot(TAB_L_FINAL, aes(x = Class_code, y = GC, fill = Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#decfa9", "#dcb8b8", "#b5cf9b")) +
  facet_wrap(~Type, labeller = labeller(Type = unique(TAB_L_FINAL$Type)), scales="free") + 
  coord_cartesian(ylim = c(25, 55)) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "none")

ggsave(paste0(comp, "/GC_content_by_type_facet_and_class_code_", toupper(flag), ".png"), height = 7, width = 10, dpi = 600)

rm(list = c("gg"))




## 4. LENGTH

gg = ggplot(TAB_FINAL_JOIN_1, aes(x = Length, y = ..density.., fill = Type)) + 
  geom_density(alpha=0.4) + 
  scale_fill_manual(values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  coord_cartesian(xlim = c(0, 10000)) +
  xlab("") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5)) + 
  theme(legend.position = "top")

ggsave(paste0(comp, "/Length_by_type_", toupper(flag), ".png"), height = 7, width = 10, dpi = 600)


gg = ggplot(TAB_FINAL_JOIN_1, aes(x = Length, y = ..density.., fill = Type)) + 
  geom_density(alpha=0.4) + 
  scale_fill_manual(values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  coord_cartesian(xlim = c(0, 10000)) +
  facet_wrap(~Type, labeller = labeller(Type = unique(TAB_FINAL_JOIN_1$Type)), scales = "free") + 
  xlab("") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5)) + 
  theme(legend.position = "top")

ggsave(paste0(comp, "/Length_by_type_facet_", toupper(flag), ".png"), height = 7, width = 15, dpi = 600)

rm(list = c("gg"))




## 5. EXONS

gg = ggplot(TAB_FINAL_JOIN_1, aes(x = Exons, y = ..density.., fill = Type)) + 
  geom_density(alpha=0.4) + 
  scale_fill_manual(values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  facet_wrap(~Type, labeller = labeller(Type = unique(TAB_FINAL_JOIN_1$Type)), scales = "free") + 
  xlab("") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5)) + 
  theme(legend.position = "top")

ggsave(paste0(comp, "/Exons_by_type_facet_", toupper(flag), ".png"), height = 7, width = 15, dpi = 600)

rm(list = c("gg"))




## 6. TPMS

gg = ggplot(TAB_FINAL_JOIN_1, aes(x = Type, y = TPMs.mean, fill = Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b")) +
  coord_cartesian(ylim = c(0, 100)) +
  xlab("") +
  ylab("TPMs") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5)) + 
  theme(legend.position = "none")

ggsave(paste0(comp, "/TPMs_by_type_", toupper(flag), ".png"), height = 7, width = 8, dpi = 600)

rm(list = c("gg"))


