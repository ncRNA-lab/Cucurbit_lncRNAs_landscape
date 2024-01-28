################################################################################
#
# FIGURES AND TABLES: COVERAGE GENOME
#
################################################################################

rm(list = ls())



## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("dplyr"))
suppressMessages(library("scales"))
suppressMessages(library("tibble"))



## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("At least 4 arguments must be supplied.", call.=FALSE)
} else {
  spes = args[1]
  spel = args[2]
  WD = args[3]
  AI = args[4]
  confidences = unlist(strsplit(args[5], " "))
}

# spes = "vvi"
# spel = "V. vinifera"
# WD = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/10-Genomic_distribution/vvi/nr"
# AI = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"




## 2. COVERAGE TRANSCRIPTOME.

# Load chromosomes file.
chrs = read.table(paste0(AI, "/Chromosomes/", spes, "_chrs.txt"), header = F, sep = "\t", quote = "\"")

# Load the tables.
tab = read.table(paste0(WD, "/TRANS-coverage.tsv"), sep = "\t", header = F, quote = "\"")
tab = tab[tab$V2 == 0, c(1,3:5)]
colnames(tab) = c("Chr", "nt_uncovered", "nt_total", "NoCov")
rownames(tab) = NULL

# Create table: all chromosomes.
tab_all = tab[tab$Chr == "genome",]
tab_all$"nt_covered" = tab_all$nt_total-tab_all$nt_uncovered
tab_all$"Cov" = (1-tab_all$NoCov)*100
tab_all$"Specie" = spes
tab_all$"Type" = "All chromosomes"
tab_all = tab_all[, c("Type", "Specie", "nt_covered", "nt_total", "Cov")]

# Create table: main chromosomes.
tab_filt = tab[tab$Chr %in% chrs$V1,]
tab_filt_red = tab_filt %>% summarise(nt_uncovered = sum(nt_uncovered), nt_total = sum(nt_total))
tab_filt_red$"nt_covered" = tab_filt_red$nt_total-tab_filt_red$nt_uncovered
tab_filt_red$"Cov" = (tab_filt_red$nt_covered*100)/tab_filt_red$nt_total
tab_filt_red$"Specie" = spes
tab_filt_red$"Type" = "Main chromosomes"
tab_filt_red = tab_filt_red[, c("Type", "Specie", "nt_covered", "nt_total", "Cov")]

# Join all and main tables.
TAB_TRANS = rbind(tab_all, tab_filt_red)
rownames(TAB_TRANS) = NULL

rm(list = c("tab", "tab_all", "tab_filt", "tab_filt_red", "chrs"))

# Save.
write.table(TAB_TRANS, paste0(WD, "/TABLE-TRANSCRIPTOME.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)



## 3. COVERAGE LNCRNAS AND GENES.

### 3.1 TABLE

TAB_GL = data.frame()

for (co in confidences) {
  
  # Load chromosomes file.
  chrs = read.table(paste0(AI, "/Chromosomes/", spes, "_chrs.txt"), header = F, sep = "\t", quote = "\"")
  
  # Load gene table.
  tab_G = read.table(paste0(WD, "/GENES-coverage.tsv"), sep = "\t", header = F, quote = "\"")
  tab_G = tab_G[tab_G$V2 == 0, c(1,3:5)]
  colnames(tab_G) = c("Chr", "nt_uncovered", "nt_total", "NoCov")
  rownames(tab_G) = NULL
  
  # Create table: all chromosomes.
  tab_G_all = tab_G[tab_G$Chr == "genome",]
  tab_G_all$"nt_covered" = tab_G_all$nt_total-tab_G_all$nt_uncovered
  tab_G_all$"Cov" = (1-tab_G_all$NoCov)*100
  tab_G_all$"Specie" = spes
  tab_G_all$"Type.1" = "All chromosomes"
  tab_G_all$"Type.2" = "PC genes"
  tab_G_all$"Confidence" = co
  tab_G_all = tab_G_all[, c("Type.1", "Type.2", "Confidence", "Specie", "nt_covered", "nt_total", "Cov")]
  
  # Create table: main chromosomes.
  tab_G_filt = tab_G[tab_G$Chr %in% chrs$V1,]
  tab_G_filt_red = tab_G_filt %>% summarise(nt_uncovered = sum(nt_uncovered), nt_total = sum(nt_total))
  tab_G_filt_red$"nt_covered" = tab_G_filt_red$nt_total-tab_G_filt_red$nt_uncovered
  tab_G_filt_red$"Cov" = (tab_G_filt_red$nt_covered*100)/tab_G_filt_red$nt_total
  tab_G_filt_red$"Specie" = spes
  tab_G_filt_red$"Type.1" = "Main chromosomes"
  tab_G_filt_red$"Type.2" = "PC genes"
  tab_G_filt_red$"Confidence" = co
  tab_G_filt_red = tab_G_filt_red[, c("Type.1", "Type.2", "Confidence", "Specie", "nt_covered", "nt_total", "Cov")]
  
  # Join all and main tables.
  tab_G = rbind(tab_G_all, tab_G_filt_red)
  
  
  # Load lncRNA table.
  tab_L = read.table(paste0(WD, "/LNCRNAS-", co, "-coverage.tsv"), sep = "\t", header = F, quote = "\"")
  tab_L = tab_L[tab_L$V2 == 0, c(1,3:5)]
  colnames(tab_L) = c("Chr", "nt_uncovered", "nt_total", "NoCov")
  rownames(tab_L) = NULL
  
  # Create table: all chromosomes.
  tab_L_all = tab_L[tab_L$Chr == "genome",]
  tab_L_all$"nt_covered" = tab_L_all$nt_total-tab_L_all$nt_uncovered
  tab_L_all$"Cov" = (1-tab_L_all$NoCov)*100
  tab_L_all$"Specie" = spes
  tab_L_all$"Type.1" = "All chromosomes"
  tab_L_all$"Type.2" = "LncRNAs"
  tab_L_all$"Confidence" = co
  tab_L_all = tab_L_all[, c("Type.1", "Type.2", "Confidence", "Specie", "nt_covered", "nt_total", "Cov")]
  
  # Create table: main chromosomes.
  tab_L_filt = tab_L[tab_L$Chr %in% chrs$V1,]
  tab_L_filt_red = tab_L_filt %>% summarise(nt_uncovered = sum(nt_uncovered), nt_total = sum(nt_total))
  tab_L_filt_red$"nt_covered" = tab_L_filt_red$nt_total-tab_L_filt_red$nt_uncovered
  tab_L_filt_red$"Cov" = (tab_L_filt_red$nt_covered*100)/tab_L_filt_red$nt_total
  tab_L_filt_red$"Specie" = spes
  tab_L_filt_red$"Type.1" = "Main chromosomes"
  tab_L_filt_red$"Type.2" = "LncRNAs"
  tab_L_filt_red$"Confidence" = co
  tab_L_filt_red = tab_L_filt_red[, c("Type.1", "Type.2", "Confidence", "Specie", "nt_covered", "nt_total", "Cov")]
  
  # Join all and main tables.
  tab_L = rbind(tab_L_all, tab_L_filt_red)
  
  # Join gene and lncRNA tables.
  tab_GL = rbind(tab_G[, c("Type.1", "Type.2", "Confidence", "Specie", "nt_covered", "nt_total", "Cov")], tab_L[, c("Type.1", "Type.2", "Confidence", "Specie", "nt_covered", "nt_total", "Cov")]) 
  
  # Join final table.
  TAB_GL = rbind(TAB_GL, tab_GL)
  
  rm(list = c("tab_L", "tab_G", "tab_GL", "tab_G_all", "tab_G_filt", "tab_G_filt_red", "tab_L_all", "tab_L_filt", "tab_L_filt_red", "chrs"))
}

rm(list = c("co"))

rownames(TAB_GL) = NULL

# Save.
write.table(TAB_GL, paste0(WD, "/TABLE-GL-BY_CONFIDENCE.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

### 3.2 FIGURES

## MAIN CHROMOSOMES

# Table
TAB_GL_mod = TAB_GL
TAB_GL_mod$Type_mod = ifelse(TAB_GL_mod$Type.2 == "LncRNAs" & TAB_GL_mod$Confidence == "Low", "LC-lncRNAs",
                             ifelse(TAB_GL_mod$Type.2 == "LncRNAs" & TAB_GL_mod$Confidence == "Medium", "MC-lncRNAs",
                                    ifelse(TAB_GL_mod$Type.2 == "LncRNAs" & TAB_GL_mod$Confidence == "High", "HC-lncRNAs",
                                           ifelse(TAB_GL_mod$Type.2 == "PC genes", "PC genes",
                                                  NA))))
TAB_GL_mod = TAB_GL_mod[TAB_GL_mod$Type.1 == "Main chromosomes", c("Type.1", "Type_mod", "Specie", "nt_covered", "nt_total", "Cov")]
TAB_GL_mod = TAB_GL_mod[!duplicated(TAB_GL_mod),]
rownames(TAB_GL_mod) = NULL

# Factors
TAB_GL_mod$Type_mod = factor(TAB_GL_mod$Type_mod, levels = c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))

# Figure
gg = ggplot(TAB_GL_mod, aes(x = Type_mod, y = Cov, fill = Type_mod)) +
  geom_bar(aes(y = Cov), stat="identity", color = "black", position=position_dodge()) +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee")) +
  scale_y_continuous(limits = c(0, 40)) +
  xlab("") +
  ylab("Genome coverage (%)") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.98, hjust=1, size = 16, angle = 45),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18)) + 
  geom_text(aes(label = round(Cov, 2)), nudge_y = 3, size = 6)

ggsave(paste0(WD, "/FIGURE-GL-MAIN.png"), height = 7, width = 9, dpi = 600)
ggsave(paste0(WD, "/FIGURE-GL-MAIN.pdf"), height = 7, width = 9, dpi = 600)

rm(list = c("gg", "TAB_GL_mod"))

