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

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("At least 4 arguments must be supplied.", call.=FALSE)
} else {
  WD1 = args[1]
  WD2 = args[2]
  species = unlist(strsplit(args[3], " "))
  classes = unlist(strsplit(args[4], " "))
  confidences = unlist(strsplit(args[5], " "))
}

# In OrthoFinder, this error doesn't appear because it always generates the Orthogroups.txt
# table and then the script get_families_from_OrthoFinder.py doesn't fail. However, if there are 
# few hits we can get some error with upsetR plots. So we use tryCatch to capture the error.



################################################################################
## 2. FAMILIES: CONSERVATION PERCENTAGE TABLE

# Aqui mostramos tanto una tabla como dos figuras globales con el porcentaje de familias 
# conservadas (Number of species by family > 1) y no conservadas (Number of species by 
# family = 1).

# Se representan dos figuras teniendo en cuenta:
#     -Fig 1: confidence level, class, specie y Number_species_by_family.
#     -Fig 2: confidence level, class, specie y type.

cat("\nFamilies: Building the conservation percentage table...\n")

TAB_1_FAM = data.frame()
TAB_2_FAM = data.frame()
TAB_3_FAM = data.frame()
TAB_ALL_FAM = data.frame()
for (confidence in confidences) {
  for (class in classes) {
    if (file.exists(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"))) {
      gen = read.table(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"), sep = "\t", header = T, quote = "\"")
      gen$"Number_species_by_family" = rowSums(gen[,species])
      gen$"Type" = ifelse(gen$Number_species_by_family > 1, "Conserved", "Non-conserved")
      gen$"Conserved_level" = ifelse(gen$Number_species_by_family == 2 | gen$Number_species_by_family == 3, "Low-conserved",
                                     ifelse(gen$Number_species_by_family == 4 | gen$Number_species_by_family == 5 | gen$Number_species_by_family == 6, "Medium-conserved",
                                            ifelse(gen$Number_species_by_family == 7 | gen$Number_species_by_family == 8 | gen$Number_species_by_family == 9, "High-conserved", "Non-conserved")))
      gen$"Specie" = sapply(strsplit(gen$Member, "-"), "[[", 2)
      fam = gen[,c("Family", "Type", "Conserved_level", "Number_species_by_family", "Specie")]
      fam = fam[!duplicated(fam),]
      fam$"Confidence" = confidence
      fam$"Class" = class
      
      fam$Type = factor(fam$Type, levels = c("Non-conserved", "Conserved"))
      fam$Conserved_level = factor(fam$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      fam$Number_species_by_family = factor(fam$Number_species_by_family, levels = 1:length(species))
      fam$Specie = factor(fam$Specie, levels = species)
      
      TAB_ALL_FAM = rbind(TAB_ALL_FAM, fam)
      
      fam_red_1 = suppressMessages(fam %>% 
        group_by(Confidence, Class, Specie, Number_species_by_family, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Family)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2)))
      fam_red_1 = as.data.frame(fam_red_1)
      colnames(fam_red_1) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_Families", "Total_Families", "Percentage_Families")
      TAB_1_FAM = rbind(TAB_1_FAM, fam_red_1)
      
      fam_red_2 = suppressMessages(fam %>% 
        group_by(Confidence, Class, Specie, Conserved_level, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Family)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2)))
      fam_red_2 = as.data.frame(fam_red_2)
      colnames(fam_red_2) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts_Families", "Total_Families", "Percentage_Families")
      TAB_2_FAM = rbind(TAB_2_FAM, fam_red_2)
      
      fam_red_3 = suppressMessages(fam %>% 
        group_by(Confidence, Class, Specie, Type, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Family)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2)))
      fam_red_3 = as.data.frame(fam_red_3)
      colnames(fam_red_3) = c("Confidence", "Class", "Specie", "Type", "Counts_Families", "Total_Families", "Percentage_Families")
      TAB_3_FAM = rbind(TAB_3_FAM, fam_red_3)
    }
  }
}

write.table(TAB_1_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_1.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_2_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_2.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_3_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_3.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_ALL_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_ALL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

TAB_1_FAM$Confidence = factor(TAB_1_FAM$Confidence, levels = confidences)
TAB_1_FAM$Class = factor(TAB_1_FAM$Class, levels = classes)

gg1 = ggplot(TAB_1_FAM, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_1.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_1.pdf"), height = 18, width = 16, dpi = 600)

TAB_2_FAM$Confidence = factor(TAB_2_FAM$Confidence, levels = confidences)
TAB_2_FAM$Class = factor(TAB_2_FAM$Class, levels = classes)

gg2 = ggplot(TAB_2_FAM, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_2.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_2.pdf"), height = 18, width = 16, dpi = 600)

TAB_3_FAM$Confidence = factor(TAB_3_FAM$Confidence, levels = confidences)
TAB_3_FAM$Class = factor(TAB_3_FAM$Class, levels = classes)

gg3 = ggplot(TAB_3_FAM, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_3.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_3.pdf"), height = 18, width = 16, dpi = 600)

################################################################################
## 3. LNCRNAS: CONSERVATION PERCENTAGE TABLE

# Aqui mostramos tanto una tabla como dos figuras globales con el porcentaje de lncRNAs 
# conservados (Number of species by family > 1, de la familia donde se encuentra el
# lncRNA) y no conservados (Number of species by family = 1, de la familia donde se 
# encuentra el lncRNA).

# Se representan dos figuras teniendo en cuenta:
#     -Fig 1: confidence level, class, specie y Number_species_by_family.
#     -Fig 2: confidence level, class, specie y type.

cat("\nLncRNAs: Building the conservation percentage table...\n")

TAB_1_LNCRNAS = data.frame()
TAB_2_LNCRNAS = data.frame()
TAB_3_LNCRNAS = data.frame()
TAB_ALL_LNCRNAS = data.frame()
for (confidence in confidences) {
  for (class in classes) {
    if (file.exists(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"))) {
      gen = read.table(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"), sep = "\t", header = T, quote = "\"")
      gen$"Number_species_by_family" = rowSums(gen[,species])
      gen$"Type" = ifelse(gen$Number_species_by_family > 1, "Conserved", "Non-conserved")
      gen$"Conserved_level" = ifelse(gen$Number_species_by_family == 2 | gen$Number_species_by_family == 3, "Low-conserved",
                                     ifelse(gen$Number_species_by_family == 4 | gen$Number_species_by_family == 5 | gen$Number_species_by_family == 6, "Medium-conserved",
                                            ifelse(gen$Number_species_by_family == 7 | gen$Number_species_by_family == 8 | gen$Number_species_by_family == 9, "High-conserved", "Non-conserved")))
      gen$"Specie" = sapply(strsplit(gen$Member, "-"), "[[", 2)
      gen = gen[,c("Member", "Type", "Conserved_level", "Number_species_by_family", "Specie")]
      gen$"Confidence" = confidence
      gen$"Class" = class
      
      gen$Type = factor(gen$Type, levels = c("Non-conserved", "Conserved"))
      gen$Conserved_level = factor(gen$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gen$Number_species_by_family = factor(gen$Number_species_by_family, levels = 1:length(species))
      gen$Specie = factor(gen$Specie, levels = species)
      
      TAB_ALL_LNCRNAS = rbind(TAB_ALL_LNCRNAS, gen)
      
      gen_red_1 = suppressMessages(gen %>% 
        group_by(Confidence, Class, Specie, Number_species_by_family, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Member)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2)))
      gen_red_1 = as.data.frame(gen_red_1)
      colnames(gen_red_1) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_LncRNAs", "Total_LncRNAs", "Percentage_LncRNAs")
      TAB_1_LNCRNAS = rbind(TAB_1_LNCRNAS, gen_red_1)
      
      gen_red_2 = suppressMessages(gen %>% 
        group_by(Confidence, Class, Specie, Conserved_level, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Member)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2)))
      gen_red_2 = as.data.frame(gen_red_2)
      colnames(gen_red_2) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts_LncRNAs", "Total_LncRNAs", "Percentage_LncRNAs")
      TAB_2_LNCRNAS = rbind(TAB_2_LNCRNAS, gen_red_2)
      
      gen_red_3 = suppressMessages(gen %>% 
        group_by(Confidence, Class, Specie, Type, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Member)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2)))
      gen_red_3 = as.data.frame(gen_red_3)
      colnames(gen_red_3) = c("Confidence", "Class", "Specie", "Type", "Counts_LncRNAs", "Total_LncRNAs", "Percentage_LncRNAs")
      TAB_3_LNCRNAS = rbind(TAB_3_LNCRNAS, gen_red_3)
    }
  }
}

write.table(TAB_1_LNCRNAS, paste0(WD2, "/TABLE_LNCRNAS_PERCENTAGE_1.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_2_LNCRNAS, paste0(WD2, "/TABLE_LNCRNAS_PERCENTAGE_2.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_3_LNCRNAS, paste0(WD2, "/TABLE_LNCRNAS_PERCENTAGE_3.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_ALL_LNCRNAS, paste0(WD2, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

TAB_1_LNCRNAS$Confidence = factor(TAB_1_LNCRNAS$Confidence, levels = confidences)
TAB_1_LNCRNAS$Class = factor(TAB_1_LNCRNAS$Class, levels = classes)

gg1 = ggplot(TAB_1_LNCRNAS, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_1.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_1.pdf"), height = 18, width = 16, dpi = 600)

TAB_2_LNCRNAS$Confidence = factor(TAB_2_LNCRNAS$Confidence, levels = confidences)
TAB_2_LNCRNAS$Class = factor(TAB_2_LNCRNAS$Class, levels = classes)

gg2 = ggplot(TAB_2_LNCRNAS, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_2.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_2.pdf"), height = 18, width = 16, dpi = 600)

TAB_3_LNCRNAS$Confidence = factor(TAB_3_LNCRNAS$Confidence, levels = confidences)
TAB_3_LNCRNAS$Class = factor(TAB_3_LNCRNAS$Class, levels = classes)

gg3 = ggplot(TAB_3_LNCRNAS, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_3.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_3.pdf"), height = 18, width = 16, dpi = 600)

