################################################################################
#
# FIGURES
#
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggbreak"))

options(bitmapType='cairo')
options(stringsAsFactors = F)

## 1. VARIABLES

flag = "nr"
WD1 = paste0("/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sequence_level/OrthoMCL/", flag, "/Prueba_genes/05-Families")
WD2 = paste0("/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sequence_level/OrthoMCL/", flag, "/Prueba_genes/06-Figures")
WD3 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-quantification"
species = c("car", "cla", "cma", "cme", "cmo", "cpe", "csa", "lsi", "mch")

if (!dir.exists(WD2)){
  dir.create(WD2)
}

# In OrthoMCL, it's important to check if the file exists. If there are some 
# blast table without hits OrthoMCL fails and doesn't generate the all_orthomcl.out 
# table. This table is necessary to get the families with the script get_families_from_OrthoMCL.py.
# So it's important to use file.exists() with Famtab and Gentab

# In OrthoFinder, this error doesn't appear because it always generates the Orthogroups.txt
# table and then the script get_families_from_OrthoFinder.py doesn't fail. However, if there are 
# few hits we can get some error with upsetR plots. So we use tryCatch to capture the error.

################################################################################
## 4. FAMILIES: CONSERVATION PERCENTAGE TABLE

# Aqui mostramos tanto una tabla como dos figuras globales con el porcentaje de familias 
# conservadas (Number of species by family > 1) y no conservadas (Number of species by 
# family = 1).

cat("\nFamilies: Building the conservation percentage table...\n")

TAB_1_FAM = data.frame()
TAB_2_FAM = data.frame()
TAB_3_FAM = data.frame()
TAB_ALL_FAM = data.frame()
if (file.exists(paste0(WD1, "/gen.tsv"))) {
  gen = read.table(paste0(WD1, "/gen.tsv"), sep = "\t", header = T, quote = "\"")
  gen$"Number_species_by_family" = rowSums(gen[,species])
  gen$"Type" = ifelse(gen$Number_species_by_family > 1, "Conserved", "Non-conserved")
  gen$"Conserved_level" = ifelse(gen$Number_species_by_family == 2 | gen$Number_species_by_family == 3, "Low-conserved",
                                 ifelse(gen$Number_species_by_family == 4 | gen$Number_species_by_family == 5 | gen$Number_species_by_family == 6, "Medium-conserved",
                                        ifelse(gen$Number_species_by_family == 7 | gen$Number_species_by_family == 8 | gen$Number_species_by_family == 9, "High-conserved", "Non-conserved")))
  gen$"Specie" = sapply(strsplit(gen$Member, "-"), "[[", 2)
  fam = gen[,c("Family", "Type", "Conserved_level", "Number_species_by_family", "Specie")]
  fam = fam[!duplicated(fam),]
  
  fam$Type = factor(fam$Type, levels = c("Non-conserved", "Conserved"))
  fam$Conserved_level = factor(fam$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
  fam$Number_species_by_family = factor(fam$Number_species_by_family, levels = 1:length(species))
  fam$Specie = factor(fam$Specie, levels = species)
  
  TAB_ALL_FAM = fam
  
  fam_red_1 = fam %>% 
    group_by(Specie, Number_species_by_family, .drop=FALSE) %>%
    summarise(
      Counts = n_distinct(Family)) %>%
    mutate(Total = sum(Counts),
           perc = round(Counts/sum(Counts) * 100, 2))
  fam_red_1 = as.data.frame(fam_red_1)
  colnames(fam_red_1) = c("Specie", "Number_species_by_family", "Counts_Families", "Total_Families", "Percentage_Families")
  TAB_1_FAM = fam_red_1
  
  fam_red_2 = fam %>% 
    group_by(Specie, Conserved_level, .drop=FALSE) %>%
    summarise(
      Counts = n_distinct(Family)) %>%
    mutate(Total = sum(Counts),
           perc = round(Counts/sum(Counts) * 100, 2))
  fam_red_2 = as.data.frame(fam_red_2)
  colnames(fam_red_2) = c("Specie", "Conserved_level", "Counts_Families", "Total_Families", "Percentage_Families")
  TAB_2_FAM = fam_red_2
  
  fam_red_3 = fam %>% 
    group_by(Specie, Type, .drop=FALSE) %>%
    summarise(
      Counts = n_distinct(Family)) %>%
    mutate(Total = sum(Counts),
           perc = round(Counts/sum(Counts) * 100, 2))
  fam_red_3 = as.data.frame(fam_red_3)
  colnames(fam_red_3) = c("Specie", "Type", "Counts_Families", "Total_Families", "Percentage_Families")
  TAB_3_FAM = fam_red_3
}

write.table(TAB_1_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_1.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_2_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_2.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_3_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_3.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_ALL_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_ALL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

gg1 = ggplot(TAB_1_FAM, aes(x = Specie, y = Percentage_Families, fill = Number_species_by_family)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_1.png"), height = 5, width = 10, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_1.pdf"), height = 5, width = 10, dpi = 600)

gg2 = ggplot(TAB_2_FAM, aes(x = Specie, y = Percentage_Families, fill = Conserved_level)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_2.png"), height = 5, width = 10, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_2.pdf"), height = 5, width = 10, dpi = 600)

gg3 = ggplot(TAB_3_FAM, aes(x = Specie, y = Percentage_Families, fill = Type)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_3.png"), height = 5, width = 10, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_3.pdf"), height = 5, width = 10, dpi = 600)


################################################################################
## 5. GENES: CONSERVATION PERCENTAGE TABLE

# Aqui mostramos tanto una tabla como dos figuras globales con el porcentaje de genes 
# conservados (Number of species by family > 1, de la familia donde se encuentra el
# gen) y no conservados (Number of species by family = 1, de la familia donde se 
# encuentra el gen).

cat("\nGenes: Building the conservation percentage table...\n")

TAB_1_GENES = data.frame()
TAB_2_GENES = data.frame()
TAB_3_GENES = data.frame()
TAB_ALL_GENES = data.frame()
if (file.exists(paste0(WD1, "/gen.tsv"))) {
  gen = read.table(paste0(WD1, "/gen.tsv"), sep = "\t", header = T, quote = "\"")
  gen$"Number_species_by_family" = rowSums(gen[,species])
  gen$"Type" = ifelse(gen$Number_species_by_family > 1, "Conserved", "Non-conserved")
  gen$"Conserved_level" = ifelse(gen$Number_species_by_family == 2 | gen$Number_species_by_family == 3, "Low-conserved",
                                 ifelse(gen$Number_species_by_family == 4 | gen$Number_species_by_family == 5 | gen$Number_species_by_family == 6, "Medium-conserved",
                                        ifelse(gen$Number_species_by_family == 7 | gen$Number_species_by_family == 8 | gen$Number_species_by_family == 9, "High-conserved", "Non-conserved")))
  gen$"Specie" = sapply(strsplit(gen$Member, "-"), "[[", 2)
  gen = gen[,c("Member", "Type", "Conserved_level", "Number_species_by_family", "Specie")]
  
  gen$Type = factor(gen$Type, levels = c("Non-conserved", "Conserved"))
  gen$Conserved_level = factor(gen$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
  gen$Number_species_by_family = factor(gen$Number_species_by_family, levels = 1:length(species))
  gen$Specie = factor(gen$Specie, levels = species)
  
  TAB_ALL_GENES = gen
  
  gen_red_1 = gen %>% 
    group_by(Specie, Number_species_by_family, .drop=FALSE) %>%
    summarise(
      Counts = n_distinct(Member)) %>%
    mutate(Total = sum(Counts),
           perc = round(Counts/sum(Counts) * 100, 2))
  gen_red_1 = as.data.frame(gen_red_1)
  colnames(gen_red_1) = c("Specie", "Number_species_by_family", "Counts_Genes", "Total_Genes", "Percentage_Genes")
  TAB_1_GENES = gen_red_1
  
  gen_red_2 = gen %>% 
    group_by(Specie, Conserved_level, .drop=FALSE) %>%
    summarise(
      Counts = n_distinct(Member)) %>%
    mutate(Total = sum(Counts),
           perc = round(Counts/sum(Counts) * 100, 2))
  gen_red_2 = as.data.frame(gen_red_2)
  colnames(gen_red_2) = c("Specie", "Conserved_level", "Counts_Genes", "Total_Genes", "Percentage_Genes")
  TAB_2_GENES = gen_red_2
  
  gen_red_3 = gen %>% 
    group_by(Specie, Type, .drop=FALSE) %>%
    summarise(
      Counts = n_distinct(Member)) %>%
    mutate(Total = sum(Counts),
           perc = round(Counts/sum(Counts) * 100, 2))
  gen_red_3 = as.data.frame(gen_red_3)
  colnames(gen_red_3) = c("Specie", "Type", "Counts_Genes", "Total_Genes", "Percentage_Genes")
  TAB_3_GENES = gen_red_3
}

write.table(TAB_1_GENES, paste0(WD2, "/TABLE_GENES_PERCENTAGE_1.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_2_GENES, paste0(WD2, "/TABLE_GENES_PERCENTAGE_2.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_3_GENES, paste0(WD2, "/TABLE_GENES_PERCENTAGE_3.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_ALL_GENES, paste0(WD2, "/TABLE_GENES_PERCENTAGE_ALL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

gg1 = ggplot(TAB_1_GENES, aes(x = Specie, y = Percentage_Genes, fill = Number_species_by_family)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Genes), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Genes_1.png"), height = 5, width = 10, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Genes_1.pdf"), height = 5, width = 10, dpi = 600)

gg2 = ggplot(TAB_2_GENES, aes(x = Specie, y = Percentage_Genes, fill = Conserved_level)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Genes), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Genes_2.png"), height = 5, width = 10, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Genes_2.pdf"), height = 5, width = 10, dpi = 600)

gg3 = ggplot(TAB_3_GENES, aes(x = Specie, y = Percentage_Genes, fill = Type)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Genes), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Genes_3.png"), height = 5, width = 10, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Genes_3.pdf"), height = 5, width = 10, dpi = 600)

