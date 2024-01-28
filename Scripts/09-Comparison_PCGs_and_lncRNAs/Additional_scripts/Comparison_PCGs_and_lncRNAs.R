################################################################################
#
# FIGURES
#
# Create plots comparing PCGs and lncRNAs: GC content, Length, Exon number and 
# TPM.
#
################################################################################

rm(list = ls())


## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggExtra))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(scales))
suppressMessages(library(ggpubr))

suppressMessages(options(bitmapType='cairo'))


## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop("At least 7 arguments must be supplied.", call.=FALSE)
} else {
  pred_path = args[1]
  quant_path = args[2]
  IR_path = args[3]
  repeat_path = args[4]
  comp_path = args[5]
  flag = args[6]
  spel = args[7]
}

# pred_path = "/mnt/doctorado/3-lncRNAs/Vitis_Tom/Results/05-LncRNAs_prediction/vvi"
# quant_path = "/mnt/doctorado/3-lncRNAs/Vitis_Tom/Results/06-Quantification/vvi"
# IR_path = "/mnt/doctorado/3-lncRNAs/Vitis_Tom/Results/07-Get_intergenic_regions/vvi"
# repeat_path = "/mnt/doctorado/3-lncRNAs/Vitis_Tom/Results/08-TEs_and_genomic_repeats/vvi"
# comp_path = "/mnt/doctorado/3-lncRNAs/Vitis_Tom/Results/09-Comparison_lncRNAs_vs_PCGs/vvi"
# flag = "nr"
# spel = "V. vinifera"


Genes_tab = paste0(pred_path, "/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv")
LncRNAs_tab = paste0(pred_path, "/STEP-FINAL/Files/LncRNAs/", flag, "/POTENTIAL_LNCRNAS_pred.tsv")
Quant_tab = paste0(quant_path, "/ALL/", flag, "/04-Table/TPMs_summary.tsv")
IR_tab = paste0(IR_path, "/Random_IR.bed")
IR_fasta = paste0(IR_path, "/Random_IR.fasta")
DB_tab = paste0(pred_path, "/STEP-FINAL/Database/Database_LncRNAs_", toupper(flag), ".tsv")
Types = c("LncRNAs low", "LncRNAs medium", "LncRNAs high")













## 2. CREATE GLOBAL TABLES (LONG AND WIDE)

cat(paste0("\n\n\nCREATE GLOBAL TABLES (LONG AND WIDE)..."))

# Load info
DB = read.table(DB_tab, sep = "\t", header = T, quote = "\"")
DB = DB[, c("ID_transcript", "Confidence")]
TPMs = read.table(Quant_tab, sep = "\t", header = T, quote = "\"")
TPMs = TPMs[, c("ID_transcript", "TPMs.mean")]
Repcon = data.frame()
for (le in Types) {
  le_mod = str_to_title(unlist(strsplit(le, " "))[2])
  Repcon_temp = read.table(paste0(repeat_path, "/02-Intersection/Final_tables/Final_tab-Repeat-", toupper(flag), "-", le_mod, "-COLLAPSED_REPEAT.tsv"), sep = "\t", header = T, quote = "\"")
  Repcon = rbind(Repcon, Repcon_temp)
}
Repcon = Repcon[!duplicated(Repcon),]
Repcon = Repcon[, c("ID_transcript", "Overlap_per")]

rm(list = c("Repcon_temp", "le_mod", "le"))

# LncRNAs
TAB_L_FINAL = read.table(LncRNAs_tab, sep = "\t", header = T, quote = "\"")

TAB_L_FINAL = merge(TAB_L_FINAL, DB, by = "ID_transcript", all = T)
TAB_L_FINAL = merge(TAB_L_FINAL, TPMs, by = "ID_transcript", all.x = T, all.y = F)
TAB_L_FINAL = merge(TAB_L_FINAL, Repcon, by = "ID_transcript", all.x = T, all.y = F)
names(TAB_L_FINAL)[names(TAB_L_FINAL) == "Confidence"] = "Type.2"
TAB_L_FINAL$"Type.1" = "LncRNAs"
TAB_L_FINAL$"Specie" = spel
TAB_L_FINAL = TAB_L_FINAL[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Specie", "Type.1", "Type.2", "Class_code", "Exons", "Length", "GC", "TPMs.mean", "Overlap_per")]

TAB_L_FINAL[TAB_L_FINAL == "High"] = "LncRNAs high"
TAB_L_FINAL[TAB_L_FINAL == "Medium"] = "LncRNAs medium"
TAB_L_FINAL[TAB_L_FINAL == "Low"] = "LncRNAs low"
TAB_L_FINAL[TAB_L_FINAL == "o"] = "o/e"
TAB_L_FINAL[TAB_L_FINAL == "e"] = "o/e"

# Genes
TAB_G_FINAL = read.table(Genes_tab, sep = "\t", header = T)

TAB_G_FINAL = merge(TAB_G_FINAL, TPMs, by = "ID_transcript", all.x = T, all.y = F)
TAB_G_FINAL = merge(TAB_G_FINAL, Repcon, by = "ID_transcript", all.x = T, all.y = F)
TAB_G_FINAL$"Type.1" = "Genes"
TAB_G_FINAL$"Type.2" = "Genes"
TAB_G_FINAL$Class_code = "pc"
TAB_G_FINAL$"Specie" = spel
TAB_G_FINAL = TAB_G_FINAL[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Specie", "Type.1", "Type.2", "Class_code", "Exons", "Length", "GC", "TPMs.mean", "Overlap_per")]

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
TAB_IR_FINAL = merge(TAB_IR_FINAL, Repcon, by = "ID_transcript", all.x = T, all.y = F)
TAB_IR_FINAL$"Type.1" = "Intergenic Regions"
TAB_IR_FINAL$"Type.2" = "Intergenic Regions"
TAB_IR_FINAL$"Class_code" = "ir"
TAB_IR_FINAL$"Exons" = NA
TAB_IR_FINAL$"Length" = NA
TAB_IR_FINAL$"TPMs.mean" = NA
TAB_IR_FINAL$"Specie" = spel
TAB_IR_FINAL = TAB_IR_FINAL[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Specie", "Type.1", "Type.2", "Class_code", "Exons", "Length", "GC", "TPMs.mean", "Overlap_per")]

rm(list = c("TAB_IR", "FASTA_IR", "id", "seq", "tab"))

# Joined: LncRNAs, Genes and Intergenic Regions
TAB_WIDE = rbind(TAB_G_FINAL, TAB_L_FINAL, TAB_IR_FINAL)
TAB_WIDE$Exons = as.numeric(TAB_WIDE$Exons)
TAB_WIDE$Length = as.numeric(TAB_WIDE$Length)
TAB_LONG = melt(setDT(TAB_WIDE), id.vars = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Specie", "Type.1", "Type.2", "Class_code"), variable.name = "Feature", value.name = "Value")
TAB_LONG$Feature = as.character(TAB_LONG$Feature)
TAB_LONG[TAB_LONG == "GC"] = "GC content"
TAB_LONG[TAB_LONG == "Exons"] = "Exon number"
TAB_LONG[TAB_LONG == "TPMs.mean"] = "Expression"
TAB_LONG[TAB_LONG == "Overlap_per"] = "Repeat content"

TAB_WIDE$Type.1 = factor(TAB_WIDE$Type.1, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_WIDE$Type.2 = factor(TAB_WIDE$Type.2, levels = c("Genes", Types, "Intergenic Regions"))
TAB_WIDE$Class_code = factor(TAB_WIDE$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG$Type.1 = factor(TAB_LONG$Type.1, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_LONG$Type.2 = factor(TAB_LONG$Type.2, levels = c("Genes", Types, "Intergenic Regions"))
TAB_LONG$Class_code = factor(TAB_LONG$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG$Feature = factor(TAB_LONG$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

TAB_WIDE = TAB_WIDE[order(TAB_WIDE$Type.2, TAB_WIDE$Class_code), ]
TAB_LONG = TAB_LONG[order(TAB_LONG$Feature, TAB_LONG$Type.2, TAB_LONG$Class_code), ]

write.table(TAB_LONG, paste0(comp_path, "/TAB_LONG-", toupper(flag), ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE, paste0(comp_path, "/TAB_WIDE-", toupper(flag), ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("Genes_tab", "LncRNAs_tab", "DB", "IR_fasta", "IR_tab", "pred_path", "quant_path", "Quant_tab", "TPMs", "TAB_G_FINAL", 
            "TAB_L_FINAL", "TAB_IR_FINAL", "Repcon"))










## 3. MEAN AND MEDIAN TABLES COMING FROM GLOBAL TABLES

cat(paste0("\n\n\nMEAN AND MEDIAN TABLES COMING FROM GLOBAL TABLES..."))

#### 3.1 (A) Without taking into account the class codes but taking into account Type.2 (Genes, LncRNAs low, LncRNAs medium, LncRNAs high, Intergenic Regions)

cat(paste0("\n\n\t-(A) Without taking into account the class codes but taking into account Type.2 (Genes, LncRNAs low, LncRNAs medium, LncRNAs high, Intergenic Regions)..."))

TAB_LONG_SUMMARY_A = TAB_LONG %>% 
  drop_na() %>%
  group_by(Specie, Type.2, Feature) %>% 
  summarise(MEAN = mean(Value),
            MEDIAN = median(Value))
TAB_WIDE_SUMMARY_A = pivot_wider(TAB_LONG_SUMMARY_A, 
                               id_cols = c("Specie", "Type.2"),
                               names_from = "Feature",
                               values_from = c("MEAN", "MEDIAN"),
                               names_glue = "{Feature}.{.value}")
colnames(TAB_WIDE_SUMMARY_A) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY_A))

TAB_LONG_SUMMARY_A$Type.2 = factor(TAB_LONG_SUMMARY_A$Type.2, levels = c("Genes", Types, "Intergenic Regions"))
TAB_LONG_SUMMARY_A$Feature = factor(TAB_LONG_SUMMARY_A$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

TAB_WIDE_SUMMARY_A$Type.2 = factor(TAB_WIDE_SUMMARY_A$Type.2, levels = c("Genes", Types, "Intergenic Regions"))

write.table(TAB_LONG_SUMMARY_A, paste0(comp_path, "/TAB_LONG_A-", toupper(flag), "-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_A, paste0(comp_path, "/TAB_WIDE_A-", toupper(flag), "-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


#### 3.2 (B) Taking into account the class codes and Type.2 (Genes, LncRNAs low, LncRNAs medium, LncRNAs high, Intergenic Regions)

cat(paste0("\n\n\t-(B) Taking into account the class codes and Type.2 (Genes, LncRNAs low, LncRNAs medium, LncRNAs high, Intergenic Regions)..."))

TAB_LONG_SUMMARY_B = TAB_LONG %>% 
  drop_na() %>%
  group_by(Specie, Type.2, Class_code, Feature) %>% 
  summarise(MEAN = mean(Value),
            MEDIAN = median(Value))
TAB_WIDE_SUMMARY_B = pivot_wider(TAB_LONG_SUMMARY_B, 
                               id_cols = c("Specie", "Type.2", "Class_code"),
                               names_from = "Feature",
                               values_from = c("MEAN", "MEDIAN"),
                               names_glue = "{Feature}.{.value}")
colnames(TAB_WIDE_SUMMARY_B) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY_B))

TAB_LONG_SUMMARY_B$Type.2 = factor(TAB_LONG_SUMMARY_B$Type.2, levels = c("Genes", Types, "Intergenic Regions"))
TAB_LONG_SUMMARY_B$Class_code = factor(TAB_LONG_SUMMARY_B$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG_SUMMARY_B$Feature = factor(TAB_LONG_SUMMARY_B$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

TAB_WIDE_SUMMARY_B$Type.2 = factor(TAB_WIDE_SUMMARY_B$Type.2, levels = c("Genes", Types, "Intergenic Regions"))
TAB_WIDE_SUMMARY_B$Class_code = factor(TAB_WIDE_SUMMARY_B$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))

write.table(TAB_LONG_SUMMARY_B, paste0(comp_path, "/TAB_LONG_B-", toupper(flag), "-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_B, paste0(comp_path, "/TAB_WIDE_B-", toupper(flag), "-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)










## 4. PLOT FIGURES BY CONFIDENCE LEVEL (ALL)

cat(paste0("\n\n\nPLOT FIGURES BY CONFIDENCE LEVEL (ALL)..."))

my_mean <- function(x) {
  log10(mean(10^x))
}

#### 4.1 (A) Without taking into account the class codes but taking into account Type.2 (Genes, LncRNAs low, LncRNAs medium, LncRNAs high, Intergenic Regions)

cat(paste0("\n\n\t-(A) Without taking into account the class codes but taking into account Type.2 (Genes, LncRNAs low, LncRNAs medium, LncRNAs high, Intergenic Regions)..."))

#### GC CONTENT

cat(paste0("\n\t\t-GC CONTENT..."))

GC_content = TAB_LONG[TAB_LONG$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Type.2, y = Value, fill = Type.2)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-LncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"), 
    values = c("#7cc1cf", c("#ce4fd3", "#eb92ef", "#edd0ee"), "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 7, 5.5, 5.5))

rm(list = c("GC_content"))

#### EXON NUMBER

cat(paste0("\n\t\t-EXON NUMBER..."))

Exon_number = TAB_LONG[TAB_LONG$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Type.2, y = Value, fill = Type.2)) + 
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-LncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"), 
    values = c("#7cc1cf", c("#ce4fd3", "#eb92ef", "#edd0ee"), "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 7, 5.5, 7))

rm(list = c("Exon_number"))

#### LENGTH

cat(paste0("\n\t\t-LENGTH..."))

Length = TAB_LONG[TAB_LONG$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Type.2, y = Value, fill = Type.2)) + 
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-LncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"), 
    values = c("#7cc1cf", c("#ce4fd3", "#eb92ef", "#edd0ee"), "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 7, 5.5, 7))

rm(list = c("Length"))

#### EXPRESSION

cat(paste0("\n\t\t-EXPRESSION..."))

Expression = TAB_LONG[TAB_LONG$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Type.2, y = Value, fill = Type.2)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-LncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"), 
    values = c("#7cc1cf", c("#ce4fd3", "#eb92ef", "#edd0ee"), "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 7, 5.5, 7))

rm(list = c("Expression"))

#### REPEAT CONTENT

cat(paste0("\n\t\t-REPEAT CONTENT..."))

Repeat_content = TAB_LONG[TAB_LONG$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Type.2, y = Value, fill = Type.2)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-LncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"), 
    values = c("#7cc1cf", c("#ce4fd3", "#eb92ef", "#edd0ee"), "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  xlab("") +
  ylab("Repeat content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 7))

rm(list = c("Repeat_content"))

#### FIGURE FEATURES

gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 5, common.legend = T, legend = "bottom")

ggsave(paste0(comp_path, "/FIGURE-A-", toupper(flag), ".png"), height = 6, width = 25, dpi = 600, bg = "white")
ggsave(paste0(comp_path, "/FIGURE-A-", toupper(flag), ".pdf"), height = 6, width = 25, dpi = 600, bg = "white")

rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))


#### 4.2 (B) Taking into account the class codes and Type.2 (Genes, LncRNAs low, LncRNAs medium, LncRNAs high, Intergenic Regions)

cat(paste0("\n\n\t-(B) Taking into account the class codes and Type.2 (Genes, LncRNAs low, LncRNAs medium, LncRNAs high, Intergenic Regions)..."))

G = TAB_LONG[TAB_LONG$Type.2 == "Genes",]
G$Type.1 = as.character(G$Type.1)
G$Type.1 = as.character(G$Type.2)
GL = G
GL[GL == "Genes"] = "LncRNAs low"
GM = G
GM[GM == "Genes"] = "LncRNAs medium"
GH = G
GH[GH == "Genes"] = "LncRNAs high"

IR = TAB_LONG[TAB_LONG$Type.2 == "Intergenic Regions",]
IR$Type.1 = as.character(IR$Type.1)
IR$Type.1 = as.character(IR$Type.2)
IRL = IR
IRL[IRL == "Intergenic Regions"] = "LncRNAs low"
IRM = IR
IRM[IRM == "Intergenic Regions"] = "LncRNAs medium"
IRH = IR
IRH[IRH == "Intergenic Regions"] = "LncRNAs high"

L = TAB_LONG[TAB_LONG$Type.2 == "LncRNAs low",]
L$Type.1 = as.character(L$Type.1)
L$Type.1 = as.character(L$Type.2)
L[L == "LncRNAs"] = "LncRNAs low"
M = TAB_LONG[TAB_LONG$Type.2 == "LncRNAs medium",]
M$Type.1 = as.character(M$Type.1)
M$Type.1 = as.character(M$Type.2)
M[M == "LncRNAs"] = "LncRNAs medium"
H = TAB_LONG[TAB_LONG$Type.2 == "LncRNAs high",]
H$Type.1 = as.character(H$Type.1)
H$Type.1 = as.character(H$Type.2)
H[H == "LncRNAs"] = "LncRNAs high"

FL = rbind(GL, L, IRL)
FM = rbind(GM, M, IRM)
FH = rbind(GH, H, IRH)

Final = rbind(FL, FM, FH)

Final$Type.2 = factor(Final$Type.2, levels = c("LncRNAs low", "LncRNAs medium", "LncRNAs high"))

rm(list = c("G", "GL", "GM", "GH", "IR", "IRL", "IRM", "IRH", "L", "M", "H", "FL", "FM", "FH"))

#### GC CONTENT

cat(paste0("\n\t\t-GC CONTENT..."))

GC_content = Final[Final$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
  ) +
  facet_grid(Type.2 ~ Feature) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 7, 5.5, 5.5))

rm(list = c("GC_content"))

#### EXON NUMBER

cat(paste0("\n\t\t-EXON NUMBER..."))

Exon_number = Final[Final$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
  ) +
  facet_grid(Type.2 ~ Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 7, 5.5, 7))

rm(list = c("Exon_number"))

#### LENGTH

cat(paste0("\n\t\t-LENGTH..."))

Length = Final[Final$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
  ) +
  facet_grid(Type.2 ~ Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 7, 5.5, 7))

rm(list = c("Length"))

#### EXPRESSION

cat(paste0("\n\t\t-EXPRESSION..."))

Expression = Final[Final$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Class_code, y = Value, fill = Class_code)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
  ) +
  facet_grid(Type.2 ~ Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 7, 5.5, 7))

rm(list = c("Expression"))

#### REPEAT CONTENT

cat(paste0("\n\t\t-REPEAT CONTENT..."))

Repeat_content = Final[Final$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
  ) +
  facet_grid(Type.2 ~ Feature) +
  xlab("") +
  ylab("Repeat content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text = element_text(size = 22, face = "bold"),
        strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 7))

rm(list = c("Repeat_content"))

#### FIGURE FEATURES

gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 5, common.legend = T, legend = "bottom")

ggsave(paste0(comp_path, "/FIGURE-B-", toupper(flag), ".png"), height = 18, width = 25, dpi = 600, bg = "white")
ggsave(paste0(comp_path, "/FIGURE-B-", toupper(flag), ".pdf"), height = 18, width = 25, dpi = 600, bg = "white")

rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5", "Final"))










## 5. PLOT FIGURES BY CONFIDENCE LEVEL (INDIVIDUAL)

cat(paste0("\n\n\nPLOT FIGURES BY CONFIDENCE LEVEL (INDIVIDUAL)..."))

TAB_STATISTICS_A = data.frame()
TAB_STATISTICS_B = data.frame()

for (i in 1:length(Types)) {
  
  type = Types[i]
  
  cat(paste0("\n\n\nType: ", type, "..."))
  
  ################################################################################
  ## 5.1 CREATE TABLES (LONG AND WIDE)
  
  cat(paste0("\n\n\t-CREATE TABLES (LONG AND WIDE)..."))
  
  # Select by confidence-level.
  tab_long = TAB_LONG[TAB_LONG$Type.2 %in% c("Genes", type, "Intergenic Regions"),]
  tab_wide = TAB_WIDE[TAB_WIDE$Type.2 %in% c("Genes", type, "Intergenic Regions"),]
  
  # Create factors
  tab_long$Type.1 = factor(tab_long$Type.1, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
  tab_long$Type.2 = factor(tab_long$Type.2, levels = c("Genes", type, "Intergenic Regions"))
  tab_long$Class_code = factor(tab_long$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  tab_long$Feature = factor(tab_long$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
  
  tab_wide$Type.1 = factor(tab_wide$Type.1, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
  tab_wide$Type.2 = factor(tab_wide$Type.2, levels = c("Genes", type, "Intergenic Regions"))
  tab_wide$Class_code = factor(tab_wide$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  
  # Save the tables
  write.table(tab_long, paste0(comp_path, "/TAB_LONG-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(tab_wide, paste0(comp_path, "/TAB_WIDE-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  ################################################################################
  ## 5.2 (A) CREATE SUMMARY TABLES (LONG AND WIDE)
  
  cat(paste0("\n\n\t-(A) CREATE SUMMARY TABLES (LONG AND WIDE)..."))
  
  # Select by confidence-level.
  tab_long_summary_A = TAB_LONG_SUMMARY_A[TAB_LONG_SUMMARY_A$Type.2 %in% c("Genes", type, "Intergenic Regions"),]
  tab_wide_summary_A = TAB_WIDE_SUMMARY_A[TAB_WIDE_SUMMARY_A$Type.2 %in% c("Genes", type, "Intergenic Regions"),]
  
  # Create factors
  tab_long_summary_A$Type.2 = factor(tab_long_summary_A$Type.2, levels = c("Genes", type, "Intergenic Regions"))
  tab_long_summary_A$Feature = factor(tab_long_summary_A$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
  
  tab_wide_summary_A$Type.2 = factor(tab_wide_summary_A$Type.2, levels = c("Genes", type, "Intergenic Regions"))
  
  write.table(tab_long_summary_A, paste0(comp_path, "/TAB_LONG_A-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), "-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(tab_wide_summary_A, paste0(comp_path, "/TAB_WIDE_A-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), "-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(list = c("tab_long_summary_A", "tab_wide_summary_A"))
  
  ################################################################################
  ## 5.3 (A) STATISTICAL ANALYSIS
  
  cat(paste0("\n\n\t-(A) STATISTICAL ANALYSIS..."))
  
  # Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
  # necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como 
  # el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
  # indipendientes. Algunos de estos papers son:
  # - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
  # - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)
  
  combinations = as.data.frame(t(combn(c("Genes", type, "Intergenic Regions"), 2)))
  rownames(combinations) = NULL
  colnames(combinations) = c("ty1", "ty2")
  
  tab_statistics_A = data.frame()
  for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
    for (j in 1:nrow(combinations)) {
      ty1 = combinations[j, "ty1"]
      ty2 = combinations[j, "ty2"]
      subset1 = tab_long[tab_long$Feature == feature & tab_long$Type.2 == ty1 & !is.na(tab_long$Value), "Value"]
      subset2 = tab_long[tab_long$Feature == feature & tab_long$Type.2 == ty2 & !is.na(tab_long$Value), "Value"]
      L1 = length(subset1$Value)
      L2 = length(subset2$Value)
      if (L1 > 0 & L2 > 0) {
        test = wilcox.test(subset1$Value, subset2$Value, paired = FALSE)
        Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
        row = data.frame(Feature = feature, Specie = spel, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
      } 
      else{
        row = data.frame(Feature = feature, Specie = spel, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
      }
      tab_statistics_A = rbind(tab_statistics_A, row)
    }
  }
  
  rownames(tab_statistics_A) = NULL
  
  write.table(tab_statistics_A, paste0(comp_path, "/TAB_STATISTICS_A-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  TAB_STATISTICS_A = rbind(TAB_STATISTICS_A, tab_statistics_A)
  
  rm(list = c("feature", "j", "ty1", "ty2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations", "tab_statistics_A"))
  
  ################################################################################
  ## 5.4 (A) FIGURES
  
  cat(paste0("\n\n\t-(A) FIGURES..."))
  
  my_mean <- function(x) {
    log10(mean(10^x))
  }
  
  colors = c("#ce4fd3", "#eb92ef", "#edd0ee")
  CLs = c("LC-LncRNAs", "MC-lncRNAs", "HC-lncRNAs")
  
  color = colors[i]
  CL = CLs[i]
  
  #### GC CONTENT
  
  cat(paste0("\n\t\t-GC CONTENT..."))
  
  GC_content = tab_long[tab_long$Feature == "GC content",]
  gg1 = ggplot(GC_content, aes(x = Type.2, y = Value, fill = Type.2)) + 
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", CL, "Intergenic regions"), 
      values = c("#7cc1cf", color, "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    xlab("") +
    ylab("GC content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 7, 5.5, 5.5))
  
  rm(list = c("GC_content"))
  
  #### EXON NUMBER
  
  cat(paste0("\n\t\t-EXON NUMBER..."))
  
  Exon_number = tab_long[tab_long$Feature == "Exon number",]
  gg2 = ggplot(Exon_number, aes(x = Type.2, y = Value, fill = Type.2)) + 
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", CL, "Intergenic regions"), 
      values = c("#7cc1cf", color, "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Exon number") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 7, 5.5, 7))
  
  rm(list = c("Exon_number"))
  
  #### LENGTH
  
  cat(paste0("\n\t\t-LENGTH..."))
  
  Length = tab_long[tab_long$Feature == "Length",]
  gg3 = ggplot(Length, aes(x = Type.2, y = Value, fill = Type.2)) + 
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", CL, "Intergenic regions"), 
      values = c("#7cc1cf", color, "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Length") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 7, 5.5, 7))
  
  rm(list = c("Length"))
  
  #### EXPRESSION
  
  cat(paste0("\n\t\t-EXPRESSION..."))
  
  Expression = tab_long[tab_long$Feature == "Expression",]
  gg4 = ggplot(Expression, aes(x = Type.2, y = Value, fill = Type.2)) + 
    geom_boxplot() + 
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", CL, "Intergenic regions"), 
      values = c("#7cc1cf", color, "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("TPM") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 7, 5.5, 7))
  
  rm(list = c("Expression"))
  
  #### REPEAT CONTENT
  
  cat(paste0("\n\t\t-REPEAT CONTENT..."))
  
  Repeat_content = tab_long[tab_long$Feature == "Repeat content",]
  gg5 = ggplot(Repeat_content, aes(x = Type.2, y = Value, fill = Type.2)) + 
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", CL, "Intergenic regions"), 
      values = c("#7cc1cf", color, "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    xlab("") +
    ylab("Repeat content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 7))
  
  rm(list = c("Repeat_content"))
  
  #### FIGURE FEATURES
  
  gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 5, common.legend = T, legend = "bottom")
  
  ggsave(paste0(comp_path, "/FIGURE-A-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), ".png"), height = 6, width = 25, dpi = 600, bg = "white")
  ggsave(paste0(comp_path, "/FIGURE-A-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), ".pdf"), height = 6, width = 25, dpi = 600, bg = "white")
  
  rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5", "color", "CL"))
  
  
  ################################################################################
  ################################################################################
  ################################################################################
  ## 5.5 (B) CREATE SUMMARY TABLES (LONG AND WIDE)
  
  cat(paste0("\n\n\t-(B) CREATE SUMMARY TABLES (LONG AND WIDE)..."))
  
  # Select by confidence-level.
  tab_long_summary_B = TAB_LONG_SUMMARY_B[TAB_LONG_SUMMARY_B$Type.2 %in% c("Genes", type, "Intergenic Regions"),]
  tab_wide_summary_B = TAB_WIDE_SUMMARY_B[TAB_WIDE_SUMMARY_B$Type.2 %in% c("Genes", type, "Intergenic Regions"),]
  
  # Create factors
  tab_long_summary_B$Type.2 = factor(tab_long_summary_B$Type.2, levels = c("Genes", type, "Intergenic Regions"))
  tab_long_summary_B$Class_code = factor(tab_long_summary_B$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  tab_long_summary_B$Feature = factor(tab_long_summary_B$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
  
  tab_wide_summary_B$Type.2 = factor(tab_wide_summary_B$Type.2, levels = c("Genes", type, "Intergenic Regions"))
  tab_wide_summary_B$Class_code = factor(tab_wide_summary_B$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  
  write.table(tab_long_summary_B, paste0(comp_path, "/TAB_LONG_B-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), "-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(tab_wide_summary_B, paste0(comp_path, "/TAB_WIDE_B-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), "-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(list = c("tab_long_summary_B", "tab_wide_summary_B"))
  
  ################################################################################
  ## 5.6 (B) STATISTICAL ANALYSIS
  
  cat(paste0("\n\n\t-(B) STATISTICAL ANALYSIS..."))
  
  # Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
  # necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como 
  # el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
  # indipendientes. Algunos de estos papers son:
  # - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
  # - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)
  
  combinations = as.data.frame(t(combn(c("pc", "u", "x", "i", "o/e", "ir"), 2)))
  rownames(combinations) = NULL
  colnames(combinations) = c("cl1", "cl2")
  
  tab_statistics_B = data.frame()
  for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
    for (j in 1:nrow(combinations)) {
      cl1 = combinations[j, "cl1"]
      cl2 = combinations[j, "cl2"]
      subset1 = tab_long[tab_long$Feature == feature & tab_long$Class_code == cl1 & !is.na(tab_long$Value), "Value"]
      subset2 = tab_long[tab_long$Feature == feature & tab_long$Class_code == cl2 & !is.na(tab_long$Value), "Value"]
      L1 = length(subset1$Value)
      L2 = length(subset2$Value)
      if (L1 > 0 & L2 > 0) {
        test = wilcox.test(subset1$Value, subset2$Value, paired = FALSE)
        Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
        row = data.frame(Type.2 = type, Specie = spel, Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
      } 
      else{
        row = data.frame(Type.2 = type, Specie = spel, Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
      }
      tab_statistics_B = rbind(tab_statistics_B, row)
    }
  }
  
  rownames(tab_statistics_B) = NULL
  
  write.table(tab_statistics_B, paste0(comp_path, "/TAB_STATISTICS_B-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  TAB_STATISTICS_B = rbind(TAB_STATISTICS_B, tab_statistics_B)
  
  rm(list = c("feature", "j", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations", "tab_statistics_B"))
    
  ################################################################################
  ## 5.7 (B) FIGURES
  
  cat(paste0("\n\n\t-(B) FIGURES..."))
  
  my_mean = function(x) {
    log10(mean(10^x))
  }
  
  #### GC CONTENT
  
  cat(paste0("\n\t\t-GC CONTENT..."))
  
  GC_content = tab_long[tab_long$Feature == "GC content",]
  gg1 = ggplot(GC_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    xlab("") +
    ylab("GC content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 7, 5.5, 5.5))
  
  rm(list = c("GC_content"))
  
  #### EXON NUMBER
  
  cat(paste0("\n\t\t-EXON NUMBER..."))
  
  Exon_number = tab_long[tab_long$Feature == "Exon number",]
  gg2 = ggplot(Exon_number, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Exon number") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 7, 5.5, 7))
  
  rm(list = c("Exon_number"))
  
  #### LENGTH
  
  cat(paste0("\n\t\t-LENGTH..."))
  
  Length = tab_long[tab_long$Feature == "Length",]
  gg3 = ggplot(Length, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Length") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 7, 5.5, 7))
  
  rm(list = c("Length"))
  
  #### EXPRESSION
  
  cat(paste0("\n\t\t-EXPRESSION..."))
  
  Expression = tab_long[tab_long$Feature == "Expression",]
  gg4 = ggplot(Expression, aes(x = Class_code, y = Value, fill = Class_code)) + 
    geom_boxplot() + 
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("TPM") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 7, 5.5, 7))
  
  rm(list = c("Expression"))
  
  #### REPEAT CONTENT
  
  cat(paste0("\n\t\t-REPEAT CONTENT..."))
  
  Repeat_content = tab_long[tab_long$Feature == "Repeat content",]
  gg5 = ggplot(Repeat_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    xlab("") +
    ylab("Repeat content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold"),
          strip.background.x = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 7))
  
  rm(list = c("Repeat_content"))
  
  #### FIGURE FEATURES
  
  gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 5, common.legend = T, legend = "bottom")
  
  ggsave(paste0(comp_path, "/FIGURE-B-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), ".png"), height = 6, width = 25, dpi = 600, bg = "white")
  ggsave(paste0(comp_path, "/FIGURE-B-", toupper(flag), "-", str_to_title(unlist(strsplit(type, " "))[2]), ".pdf"), height = 6, width = 25, dpi = 600, bg = "white")
  
  rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))
  
  rm(list = c("tab_wide", "tab_long"))
}

rm(list = c("CLs", "colors", "type", "i"))

write.table(TAB_STATISTICS_A, paste0(comp_path, "/TAB_STATISTICS_A-", toupper(flag), ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_STATISTICS_B, paste0(comp_path, "/TAB_STATISTICS_B-", toupper(flag), ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)



