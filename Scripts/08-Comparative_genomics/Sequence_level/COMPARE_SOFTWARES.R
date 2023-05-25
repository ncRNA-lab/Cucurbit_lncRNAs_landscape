################################################################################
#
# FIGURES: COMPARISON
#
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggvenn"))

options(bitmapType='cairo')
options(stringsAsFactors = F)

## 1. VARIABLES

WD1 = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sequence_level"
species = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
classes = c("ALL", "intergenic", "antisense", "intronic", "sense")
confidences = c("Low", "Medium", "High")


## 2. DIRECTORIES

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats"))
}

################################################################################
## 3. NON-REDUNDANT AND REDUNDANT COMPARISON

cat("\n\n\nNON-REDUNDANT AND REDUNDANT COMPARISON")

### 3.1 Blastn (Families)

cat("\n\nBlastn (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_Blastn_NR_def = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_Blastn_R_def = read.table(paste0(WD1, "/Blastn/r/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_Blastn_NR_def = TAB_FAM_Blastn_NR_def[TAB_FAM_Blastn_NR_def$Class == class,]
    TAB_FAM_Blastn_R_def = TAB_FAM_Blastn_R_def[TAB_FAM_Blastn_R_def$Class == class,]
    
    # Add label.
    TAB_FAM_Blastn_NR_def$"Label" = "Blastn non-redundant"
    TAB_FAM_Blastn_R_def$"Label" = "Blastn redundant"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_Blastn_NR_def, TAB_FAM_Blastn_R_def)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("Blastn non-redundant", "Blastn redundant"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 3.2 Blastn (LncRNAs)

cat("\n\nBlastn (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_Blastn_NR_def = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_Blastn_R_def = read.table(paste0(WD1, "/Blastn/r/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_Blastn_NR_def = TAB_LNCRNAS_Blastn_NR_def[TAB_LNCRNAS_Blastn_NR_def$Class == class,]
    TAB_LNCRNAS_Blastn_R_def = TAB_LNCRNAS_Blastn_R_def[TAB_LNCRNAS_Blastn_R_def$Class == class,]
    
    # Add label.
    TAB_LNCRNAS_Blastn_NR_def$"Label" = "Blastn non-redundant"
    TAB_LNCRNAS_Blastn_R_def$"Label" = "Blastn redundant"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_Blastn_NR_def, TAB_LNCRNAS_Blastn_R_def)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("Blastn non-redundant", "Blastn redundant"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 3.3 OrthoFinder (Families)

cat("\n\nOrthoFinder (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_OrthoFinder_NR_def = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoFinder_R_def = read.table(paste0(WD1, "/OrthoFinder/r/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_OrthoFinder_NR_def = TAB_FAM_OrthoFinder_NR_def[TAB_FAM_OrthoFinder_NR_def$Class == class,]
    TAB_FAM_OrthoFinder_R_def = TAB_FAM_OrthoFinder_R_def[TAB_FAM_OrthoFinder_R_def$Class == class,]
    
    # Add label.
    TAB_FAM_OrthoFinder_NR_def$"Label" = "OrthoFinder non-redundant"
    TAB_FAM_OrthoFinder_R_def$"Label" = "OrthoFinder redundant"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_OrthoFinder_NR_def, TAB_FAM_OrthoFinder_R_def)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoFinder non-redundant", "OrthoFinder redundant"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 3.4 OrthoFinder (LncRNAs)

cat("\n\nOrthoFinder (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_OrthoFinder_NR_def = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoFinder_R_def = read.table(paste0(WD1, "/OrthoFinder/r/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_OrthoFinder_NR_def = TAB_LNCRNAS_OrthoFinder_NR_def[TAB_LNCRNAS_OrthoFinder_NR_def$Class == class,]
    TAB_LNCRNAS_OrthoFinder_R_def = TAB_LNCRNAS_OrthoFinder_R_def[TAB_LNCRNAS_OrthoFinder_R_def$Class == class,]
    
    # Add label.
    TAB_LNCRNAS_OrthoFinder_NR_def$"Label" = "OrthoFinder non-redundant"
    TAB_LNCRNAS_OrthoFinder_R_def$"Label" = "OrthoFinder redundant"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_OrthoFinder_NR_def, TAB_LNCRNAS_OrthoFinder_R_def)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoFinder non-redundant", "OrthoFinder redundant"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 3.5 OrthoMCL (Families)

cat("\n\nOrthoMCL (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_OrthoMCL_NR_def = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoMCL_R_def = read.table(paste0(WD1, "/OrthoMCL/r/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_OrthoMCL_NR_def = TAB_FAM_OrthoMCL_NR_def[TAB_FAM_OrthoMCL_NR_def$Class == class,]
    TAB_FAM_OrthoMCL_R_def = TAB_FAM_OrthoMCL_R_def[TAB_FAM_OrthoMCL_R_def$Class == class,]
    
    # Add label.
    TAB_FAM_OrthoMCL_NR_def$"Label" = "OrthoMCL non-redundant"
    TAB_FAM_OrthoMCL_R_def$"Label" = "OrthoMCL redundant"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_OrthoMCL_NR_def, TAB_FAM_OrthoMCL_R_def)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoMCL non-redundant", "OrthoMCL redundant"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 3.6 OrthoMCL (LncRNAs)

cat("\n\nOrthoMCL (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_OrthoMCL_NR_def = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoMCL_R_def = read.table(paste0(WD1, "/OrthoMCL/r/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_OrthoMCL_NR_def = TAB_LNCRNAS_OrthoMCL_NR_def[TAB_LNCRNAS_OrthoMCL_NR_def$Class == class,]
    TAB_LNCRNAS_OrthoMCL_R_def = TAB_LNCRNAS_OrthoMCL_R_def[TAB_LNCRNAS_OrthoMCL_R_def$Class == class,]
    
    # Add label.
    TAB_LNCRNAS_OrthoMCL_NR_def$"Label" = "OrthoMCL non-redundant"
    TAB_LNCRNAS_OrthoMCL_R_def$"Label" = "OrthoMCL redundant"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_OrthoMCL_NR_def, TAB_LNCRNAS_OrthoMCL_R_def)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoMCL non-redundant", "OrthoMCL redundant"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/NR_and_R_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}


################################################################################
## 4. PARAMETERS COMPARISON

cat("\n\n\nPARAMETERS COMPARISON")

### 4.1 Blastn (Families)

cat("\n\nBlastn (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_Blastn_NR_def = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_Blastn_NR_P1 = read.table(paste0(WD1, "/Blastn/nr/Prueba_e-value/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_Blastn_NR_P2 = read.table(paste0(WD1, "/Blastn/nr/Prueba_filters_identity_percentage_and_alignment_length/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_Blastn_NR_P3 = read.table(paste0(WD1, "/Blastn/nr/Prueba_max_5/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_Blastn_NR_def = TAB_FAM_Blastn_NR_def[TAB_FAM_Blastn_NR_def$Class == class,]
    TAB_FAM_Blastn_NR_P1 = TAB_FAM_Blastn_NR_P1[TAB_FAM_Blastn_NR_P1$Class == class,]
    TAB_FAM_Blastn_NR_P2 = TAB_FAM_Blastn_NR_P2[TAB_FAM_Blastn_NR_P2$Class == class,]
    TAB_FAM_Blastn_NR_P3 = TAB_FAM_Blastn_NR_P3[TAB_FAM_Blastn_NR_P3$Class == class,]
    
    # Add label.
    TAB_FAM_Blastn_NR_def$"Label" = "Blastn (Default)"
    TAB_FAM_Blastn_NR_P1$"Label" = "Blastn (Test 1)"
    TAB_FAM_Blastn_NR_P2$"Label" = "Blastn (Test 2)"
    TAB_FAM_Blastn_NR_P3$"Label" = "Blastn (Test 3)"
    
    # Add legend.
    TAB_FAM_Blastn_NR_def$"Info" = "-evalue 1e-3 -max_target_seqs 1 -max_hsps 1"
    TAB_FAM_Blastn_NR_P1$"Info" = "-evalue 1e-5"
    TAB_FAM_Blastn_NR_P2$"Info" = "Identity percentage >= 20 & Alignment length >= 50"
    TAB_FAM_Blastn_NR_P3$"Info" = "-max_target_seqs 5 -max_hsps 5"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_Blastn_NR_def, TAB_FAM_Blastn_NR_P1, TAB_FAM_Blastn_NR_P2, TAB_FAM_Blastn_NR_P3)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("Blastn (Default)", "Blastn (Test 1)", "Blastn (Test 2)", "Blastn (Test 3)"))
    TAB_FINAL$Info = factor(TAB_FINAL$Info, levels = c("-evalue 1e-3 -max_target_seqs 1 -max_hsps 1", "-evalue 1e-5", "Identity percentage >= 20 & Alignment length >= 50", "-max_target_seqs 5 -max_hsps 5"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
  }
}

### 4.2 Blastn (LncRNAs)

cat("\n\nBlastn (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_Blastn_NR_def = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_Blastn_NR_P1 = read.table(paste0(WD1, "/Blastn/nr/Prueba_e-value/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_Blastn_NR_P2 = read.table(paste0(WD1, "/Blastn/nr/Prueba_filters_identity_percentage_and_alignment_length/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_Blastn_NR_P3 = read.table(paste0(WD1, "/Blastn/nr/Prueba_max_5/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_Blastn_NR_def = TAB_LNCRNAS_Blastn_NR_def[TAB_LNCRNAS_Blastn_NR_def$Class == class,]
    TAB_LNCRNAS_Blastn_NR_P1 = TAB_LNCRNAS_Blastn_NR_P1[TAB_LNCRNAS_Blastn_NR_P1$Class == class,]
    TAB_LNCRNAS_Blastn_NR_P2 = TAB_LNCRNAS_Blastn_NR_P2[TAB_LNCRNAS_Blastn_NR_P2$Class == class,]
    TAB_LNCRNAS_Blastn_NR_P3 = TAB_LNCRNAS_Blastn_NR_P3[TAB_LNCRNAS_Blastn_NR_P3$Class == class,]
    
    # Add label.
    TAB_LNCRNAS_Blastn_NR_def$"Label" = "Blastn (Default)"
    TAB_LNCRNAS_Blastn_NR_P1$"Label" = "Blastn (Test 1)"
    TAB_LNCRNAS_Blastn_NR_P2$"Label" = "Blastn (Test 2)"
    TAB_LNCRNAS_Blastn_NR_P3$"Label" = "Blastn (Test 3)"
    
    # Add legend.
    TAB_LNCRNAS_Blastn_NR_def$"Info" = "-evalue 1e-3 -max_target_seqs 1 -max_hsps 1"
    TAB_LNCRNAS_Blastn_NR_P1$"Info" = "-evalue 1e-5"
    TAB_LNCRNAS_Blastn_NR_P2$"Info" = "Identity percentage >= 20 & Alignment length >= 50"
    TAB_LNCRNAS_Blastn_NR_P3$"Info" = "-max_target_seqs 5 -max_hsps 5"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_Blastn_NR_def, TAB_LNCRNAS_Blastn_NR_P1, TAB_LNCRNAS_Blastn_NR_P2, TAB_LNCRNAS_Blastn_NR_P3)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("Blastn (Default)", "Blastn (Test 1)", "Blastn (Test 2)", "Blastn (Test 3)"))
    TAB_FINAL$Info = factor(TAB_FINAL$Info, levels = c("-evalue 1e-3 -max_target_seqs 1 -max_hsps 1", "-evalue 1e-5", "Identity percentage >= 20 & Alignment length >= 50", "-max_target_seqs 5 -max_hsps 5"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
  }
}


### 4.3 OrthoFinder (Families)

cat("\n\nOrthoFinder (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_OrthoFinder_NR_def = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoFinder_NR_P1 = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_e-value/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoFinder_NR_P2 = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_filters_identity_percentage_and_alignment_length/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoFinder_NR_P3 = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_max_5/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_OrthoFinder_NR_def = TAB_FAM_OrthoFinder_NR_def[TAB_FAM_OrthoFinder_NR_def$Class == class,]
    TAB_FAM_OrthoFinder_NR_P1 = TAB_FAM_OrthoFinder_NR_P1[TAB_FAM_OrthoFinder_NR_P1$Class == class,]
    TAB_FAM_OrthoFinder_NR_P2 = TAB_FAM_OrthoFinder_NR_P2[TAB_FAM_OrthoFinder_NR_P2$Class == class,]
    TAB_FAM_OrthoFinder_NR_P3 = TAB_FAM_OrthoFinder_NR_P3[TAB_FAM_OrthoFinder_NR_P3$Class == class,]
    
    # Add label.
    TAB_FAM_OrthoFinder_NR_def$"Label" = "OrthoFinder (Default)"
    TAB_FAM_OrthoFinder_NR_P1$"Label" = "OrthoFinder (Test 1)"
    TAB_FAM_OrthoFinder_NR_P2$"Label" = "OrthoFinder (Test 2)"
    TAB_FAM_OrthoFinder_NR_P3$"Label" = "OrthoFinder (Test 3)"
    
    # Add legend.
    TAB_FAM_OrthoFinder_NR_def$"Info" = "-evalue 1e-3 -max_target_seqs 1 -max_hsps 1"
    TAB_FAM_OrthoFinder_NR_P1$"Info" = "-evalue 1e-5"
    TAB_FAM_OrthoFinder_NR_P2$"Info" = "Identity percentage >= 20 & Alignment length >= 50"
    TAB_FAM_OrthoFinder_NR_P3$"Info" = "-max_target_seqs 5 -max_hsps 5"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_OrthoFinder_NR_def, TAB_FAM_OrthoFinder_NR_P1, TAB_FAM_OrthoFinder_NR_P2, TAB_FAM_OrthoFinder_NR_P3)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoFinder (Default)", "OrthoFinder (Test 1)", "OrthoFinder (Test 2)", "OrthoFinder (Test 3)"))
    TAB_FINAL$Info = factor(TAB_FINAL$Info, levels = c("-evalue 1e-3 -max_target_seqs 1 -max_hsps 1", "-evalue 1e-5", "Identity percentage >= 20 & Alignment length >= 50", "-max_target_seqs 5 -max_hsps 5"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
  }
}

### 4.4 OrthoFinder (LncRNAs)

cat("\n\nOrthoFinder (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_OrthoFinder_NR_def = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoFinder_NR_P1 = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_e-value/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoFinder_NR_P2 = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_filters_identity_percentage_and_alignment_length/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoFinder_NR_P3 = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_max_5/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_OrthoFinder_NR_def = TAB_LNCRNAS_OrthoFinder_NR_def[TAB_LNCRNAS_OrthoFinder_NR_def$Class == class,]
    TAB_LNCRNAS_OrthoFinder_NR_P1 = TAB_LNCRNAS_OrthoFinder_NR_P1[TAB_LNCRNAS_OrthoFinder_NR_P1$Class == class,]
    TAB_LNCRNAS_OrthoFinder_NR_P2 = TAB_LNCRNAS_OrthoFinder_NR_P2[TAB_LNCRNAS_OrthoFinder_NR_P2$Class == class,]
    TAB_LNCRNAS_OrthoFinder_NR_P3 = TAB_LNCRNAS_OrthoFinder_NR_P3[TAB_LNCRNAS_OrthoFinder_NR_P3$Class == class,]
    
    # Add label.
    TAB_LNCRNAS_OrthoFinder_NR_def$"Label" = "OrthoFinder (Default)"
    TAB_LNCRNAS_OrthoFinder_NR_P1$"Label" = "OrthoFinder (Test 1)"
    TAB_LNCRNAS_OrthoFinder_NR_P2$"Label" = "OrthoFinder (Test 2)"
    TAB_LNCRNAS_OrthoFinder_NR_P3$"Label" = "OrthoFinder (Test 3)"
    
    # Add legend.
    TAB_LNCRNAS_OrthoFinder_NR_def$"Info" = "-evalue 1e-3 -max_target_seqs 1 -max_hsps 1"
    TAB_LNCRNAS_OrthoFinder_NR_P1$"Info" = "-evalue 1e-5"
    TAB_LNCRNAS_OrthoFinder_NR_P2$"Info" = "Identity percentage >= 20 & Alignment length >= 50"
    TAB_LNCRNAS_OrthoFinder_NR_P3$"Info" = "-max_target_seqs 5 -max_hsps 5"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_OrthoFinder_NR_def, TAB_LNCRNAS_OrthoFinder_NR_P1, TAB_LNCRNAS_OrthoFinder_NR_P2, TAB_LNCRNAS_OrthoFinder_NR_P3)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoFinder (Default)", "OrthoFinder (Test 1)", "OrthoFinder (Test 2)", "OrthoFinder (Test 3)"))
    TAB_FINAL$Info = factor(TAB_FINAL$Info, levels = c("-evalue 1e-3 -max_target_seqs 1 -max_hsps 1", "-evalue 1e-5", "Identity percentage >= 20 & Alignment length >= 50", "-max_target_seqs 5 -max_hsps 5"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
  }
}

### 4.5 OrthoMCL (Families)

cat("\n\nOrthoMCL (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_OrthoMCL_NR_def = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoMCL_NR_P1 = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_e-value/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoMCL_NR_P2 = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_filters_identity_percentage_and_alignment_length/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoMCL_NR_P3 = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_max_5/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_OrthoMCL_NR_def = TAB_FAM_OrthoMCL_NR_def[TAB_FAM_OrthoMCL_NR_def$Class == class,]
    TAB_FAM_OrthoMCL_NR_P1 = TAB_FAM_OrthoMCL_NR_P1[TAB_FAM_OrthoMCL_NR_P1$Class == class,]
    TAB_FAM_OrthoMCL_NR_P2 = TAB_FAM_OrthoMCL_NR_P2[TAB_FAM_OrthoMCL_NR_P2$Class == class,]
    TAB_FAM_OrthoMCL_NR_P3 = TAB_FAM_OrthoMCL_NR_P3[TAB_FAM_OrthoMCL_NR_P3$Class == class,]
    
    # Add label.
    TAB_FAM_OrthoMCL_NR_def$"Label" = "OrthoMCL (Default)"
    TAB_FAM_OrthoMCL_NR_P1$"Label" = "OrthoMCL (Test 1)"
    TAB_FAM_OrthoMCL_NR_P2$"Label" = "OrthoMCL (Test 2)"
    TAB_FAM_OrthoMCL_NR_P3$"Label" = "OrthoMCL (Test 3)"
    
    # Add legend.
    TAB_FAM_OrthoMCL_NR_def$"Info" = "-evalue 1e-3 -max_target_seqs 1 -max_hsps 1"
    TAB_FAM_OrthoMCL_NR_P1$"Info" = "-evalue 1e-5"
    TAB_FAM_OrthoMCL_NR_P2$"Info" = "Identity percentage >= 20 & Alignment length >= 50"
    TAB_FAM_OrthoMCL_NR_P3$"Info" = "-max_target_seqs 5 -max_hsps 5"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_OrthoMCL_NR_def, TAB_FAM_OrthoMCL_NR_P1, TAB_FAM_OrthoMCL_NR_P2, TAB_FAM_OrthoMCL_NR_P3)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoMCL (Default)", "OrthoMCL (Test 1)", "OrthoMCL (Test 2)", "OrthoMCL (Test 3)"))
    TAB_FINAL$Info = factor(TAB_FINAL$Info, levels = c("-evalue 1e-3 -max_target_seqs 1 -max_hsps 1", "-evalue 1e-5", "Identity percentage >= 20 & Alignment length >= 50", "-max_target_seqs 5 -max_hsps 5"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
  }
}

### 4.6 OrthoMCL (LncRNAs)

cat("\n\nOrthoMCL (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_OrthoMCL_NR_def = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoMCL_NR_P1 = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_e-value/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoMCL_NR_P2 = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_filters_identity_percentage_and_alignment_length/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoMCL_NR_P3 = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_max_5/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_OrthoMCL_NR_def = TAB_LNCRNAS_OrthoMCL_NR_def[TAB_LNCRNAS_OrthoMCL_NR_def$Class == class,]
    TAB_LNCRNAS_OrthoMCL_NR_P1 = TAB_LNCRNAS_OrthoMCL_NR_P1[TAB_LNCRNAS_OrthoMCL_NR_P1$Class == class,]
    TAB_LNCRNAS_OrthoMCL_NR_P2 = TAB_LNCRNAS_OrthoMCL_NR_P2[TAB_LNCRNAS_OrthoMCL_NR_P2$Class == class,]
    TAB_LNCRNAS_OrthoMCL_NR_P3 = TAB_LNCRNAS_OrthoMCL_NR_P3[TAB_LNCRNAS_OrthoMCL_NR_P3$Class == class,]
    
    # Add label.
    TAB_LNCRNAS_OrthoMCL_NR_def$"Label" = "OrthoMCL (Default)"
    TAB_LNCRNAS_OrthoMCL_NR_P1$"Label" = "OrthoMCL (Test 1)"
    TAB_LNCRNAS_OrthoMCL_NR_P2$"Label" = "OrthoMCL (Test 2)"
    TAB_LNCRNAS_OrthoMCL_NR_P3$"Label" = "OrthoMCL (Test 3)"
    
    # Add legend.
    TAB_LNCRNAS_OrthoMCL_NR_def$"Info" = "-evalue 1e-3 -max_target_seqs 1 -max_hsps 1"
    TAB_LNCRNAS_OrthoMCL_NR_P1$"Info" = "-evalue 1e-5"
    TAB_LNCRNAS_OrthoMCL_NR_P2$"Info" = "Identity percentage >= 20 & Alignment length >= 50"
    TAB_LNCRNAS_OrthoMCL_NR_P3$"Info" = "-max_target_seqs 5 -max_hsps 5"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_OrthoMCL_NR_def, TAB_LNCRNAS_OrthoMCL_NR_P1, TAB_LNCRNAS_OrthoMCL_NR_P2, TAB_LNCRNAS_OrthoMCL_NR_P3)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoMCL (Default)", "OrthoMCL (Test 1)", "OrthoMCL (Test 2)", "OrthoMCL (Test 3)"))
    TAB_FINAL$Info = factor(TAB_FINAL$Info, levels = c("-evalue 1e-3 -max_target_seqs 1 -max_hsps 1", "-evalue 1e-5", "Identity percentage >= 20 & Alignment length >= 50", "-max_target_seqs 5 -max_hsps 5"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 18, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Parameters_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 18, width = 16, dpi = 600)
    }
  }
}


################################################################################
## 5. PIPELINES COMPARISON

cat("\n\n\nPIPELINES COMPARISON")

### 5.1 Families

cat("\n\nFamilies")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_Blastn_NR_def = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoFinder_NR_def = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoMCL_NR_def = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_Blastn_NR_def = TAB_FAM_Blastn_NR_def[TAB_FAM_Blastn_NR_def$Class == class,]
    TAB_FAM_OrthoFinder_NR_def = TAB_FAM_OrthoFinder_NR_def[TAB_FAM_OrthoFinder_NR_def$Class == class,]
    TAB_FAM_OrthoMCL_NR_def = TAB_FAM_OrthoMCL_NR_def[TAB_FAM_OrthoMCL_NR_def$Class == class,]
    
    # Add label.
    TAB_FAM_Blastn_NR_def$"Label" = "Blastn"
    TAB_FAM_OrthoFinder_NR_def$"Label" = "OrthoFinder"
    TAB_FAM_OrthoMCL_NR_def$"Label" = "OrthoMCL"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_Blastn_NR_def, TAB_FAM_OrthoFinder_NR_def, TAB_FAM_OrthoMCL_NR_def)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("Blastn", "OrthoFinder", "OrthoMCL"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 12, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 12, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 12, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 12, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 12, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 12, width = 16, dpi = 600)
    }
  }
}


### 5.2 LncRNAs

cat("\n\nLncRNAs")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_Blastn_NR_def = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoFinder_NR_def = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_LNCRNAS_OrthoMCL_NR_def = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_Blastn_NR_def = TAB_LNCRNAS_Blastn_NR_def[TAB_LNCRNAS_Blastn_NR_def$Class == class,]
    TAB_LNCRNAS_OrthoFinder_NR_def = TAB_LNCRNAS_OrthoFinder_NR_def[TAB_LNCRNAS_OrthoFinder_NR_def$Class == class,]
    TAB_LNCRNAS_OrthoMCL_NR_def = TAB_LNCRNAS_OrthoMCL_NR_def[TAB_LNCRNAS_OrthoMCL_NR_def$Class == class,]
    
    # Add label.
    TAB_LNCRNAS_Blastn_NR_def$"Label" = "Blastn"
    TAB_LNCRNAS_OrthoFinder_NR_def$"Label" = "OrthoFinder"
    TAB_LNCRNAS_OrthoMCL_NR_def$"Label" = "OrthoMCL"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_Blastn_NR_def, TAB_LNCRNAS_OrthoFinder_NR_def, TAB_LNCRNAS_OrthoMCL_NR_def)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("Blastn", "OrthoFinder", "OrthoMCL"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 12, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 12, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 12, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 12, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 12, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Pipelines_comparison/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 12, width = 16, dpi = 600)
    }
  }
}


################################################################################
## 6. GENES AND LNCRNAS COMPARISON

cat("\n\n\nGENES AND LNCRNAS COMPARISON")

### 6.1 Blastn (Families)

cat("\n\nBlastn (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_Blastn_NR_def = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_Blastn_genes = read.table(paste0(WD1, "/Blastn/nr/Prueba_genes/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_Blastn_NR_def = TAB_FAM_Blastn_NR_def[TAB_FAM_Blastn_NR_def$Class == class,]
    
    # Modify genes table.
    TAB_FAM_Blastn_genes_Low = TAB_FAM_Blastn_genes
    TAB_FAM_Blastn_genes_Low$"Confidence" = "Low"
    TAB_FAM_Blastn_genes_Medium = TAB_FAM_Blastn_genes
    TAB_FAM_Blastn_genes_Medium$"Confidence" = "Medium"
    TAB_FAM_Blastn_genes_High = TAB_FAM_Blastn_genes
    TAB_FAM_Blastn_genes_High$"Confidence" = "High"
    TAB_FAM_Blastn_genes = rbind(TAB_FAM_Blastn_genes_Low, TAB_FAM_Blastn_genes_Medium, TAB_FAM_Blastn_genes_High)
    TAB_FAM_Blastn_genes$"Class" = class
    if (i == 1) {
      TAB_FAM_Blastn_genes = TAB_FAM_Blastn_genes[,c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    if (i == 2) {
      TAB_FAM_Blastn_genes = TAB_FAM_Blastn_genes[,c("Confidence", "Class", "Specie", "Conserved_level", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    if (i == 3) {
      TAB_FAM_Blastn_genes = TAB_FAM_Blastn_genes[,c("Confidence", "Class", "Specie", "Type", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    
    # Add label.
    TAB_FAM_Blastn_NR_def$"Label" = "Blastn LncRNAs"
    TAB_FAM_Blastn_genes$"Label" = "Blastn Genes"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_Blastn_NR_def, TAB_FAM_Blastn_genes)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("Blastn LncRNAs", "Blastn Genes"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 6.2 Blastn (LncRNAs)

cat("\n\nBlastn (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_Blastn_NR_def = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_GENES_Blastn = read.table(paste0(WD1, "/Blastn/nr/Prueba_genes/06-Figures/TABLE_GENES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_Blastn_NR_def = TAB_LNCRNAS_Blastn_NR_def[TAB_LNCRNAS_Blastn_NR_def$Class == class,]
    
    # Modify genes table.
    TAB_GENES_Blastn_Low = TAB_GENES_Blastn
    TAB_GENES_Blastn_Low$"Confidence" = "Low"
    TAB_GENES_Blastn_Medium = TAB_GENES_Blastn
    TAB_GENES_Blastn_Medium$"Confidence" = "Medium"
    TAB_GENES_Blastn_High = TAB_GENES_Blastn
    TAB_GENES_Blastn_High$"Confidence" = "High"
    TAB_GENES_Blastn = rbind(TAB_GENES_Blastn_Low, TAB_GENES_Blastn_Medium, TAB_GENES_Blastn_High)
    TAB_GENES_Blastn$"Class" = class
    if (i == 1) {
      TAB_GENES_Blastn = TAB_GENES_Blastn[,c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    if (i == 2) {
      TAB_GENES_Blastn = TAB_GENES_Blastn[,c("Confidence", "Class", "Specie", "Conserved_level", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    if (i == 3) {
      TAB_GENES_Blastn = TAB_GENES_Blastn[,c("Confidence", "Class", "Specie", "Type", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    
    # Rename columns.
    if (i == 1) {
      colnames(TAB_LNCRNAS_Blastn_NR_def) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_Blastn) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts", "Total", "Percentage")
    }
    if (i == 2) {
      colnames(TAB_LNCRNAS_Blastn_NR_def) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_Blastn) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts", "Total", "Percentage")
    }
    if (i == 3) {
      colnames(TAB_LNCRNAS_Blastn_NR_def) = c("Confidence", "Class", "Specie", "Type", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_Blastn) = c("Confidence", "Class", "Specie", "Type", "Counts", "Total", "Percentage")
    }
    
    # Add label.
    TAB_LNCRNAS_Blastn_NR_def$"Label" = "Blastn LncRNAs"
    TAB_GENES_Blastn$"Label" = "Blastn Genes"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_Blastn_NR_def, TAB_GENES_Blastn)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("Blastn LncRNAs", "Blastn Genes"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/Blastn/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 6.3 OrthoFinder (Families)

cat("\n\nOrthoFinder (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_OrthoFinder_NR_def = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoFinder_genes = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_genes/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_OrthoFinder_NR_def = TAB_FAM_OrthoFinder_NR_def[TAB_FAM_OrthoFinder_NR_def$Class == class,]
    
    # Modify genes table.
    TAB_FAM_OrthoFinder_genes_Low = TAB_FAM_OrthoFinder_genes
    TAB_FAM_OrthoFinder_genes_Low$"Confidence" = "Low"
    TAB_FAM_OrthoFinder_genes_Medium = TAB_FAM_OrthoFinder_genes
    TAB_FAM_OrthoFinder_genes_Medium$"Confidence" = "Medium"
    TAB_FAM_OrthoFinder_genes_High = TAB_FAM_OrthoFinder_genes
    TAB_FAM_OrthoFinder_genes_High$"Confidence" = "High"
    TAB_FAM_OrthoFinder_genes = rbind(TAB_FAM_OrthoFinder_genes_Low, TAB_FAM_OrthoFinder_genes_Medium, TAB_FAM_OrthoFinder_genes_High)
    TAB_FAM_OrthoFinder_genes$"Class" = class
    if (i == 1) {
      TAB_FAM_OrthoFinder_genes = TAB_FAM_OrthoFinder_genes[,c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    if (i == 2) {
      TAB_FAM_OrthoFinder_genes = TAB_FAM_OrthoFinder_genes[,c("Confidence", "Class", "Specie", "Conserved_level", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    if (i == 3) {
      TAB_FAM_OrthoFinder_genes = TAB_FAM_OrthoFinder_genes[,c("Confidence", "Class", "Specie", "Type", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    
    # Add label.
    TAB_FAM_OrthoFinder_NR_def$"Label" = "OrthoFinder LncRNAs"
    TAB_FAM_OrthoFinder_genes$"Label" = "OrthoFinder Genes"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_OrthoFinder_NR_def, TAB_FAM_OrthoFinder_genes)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoFinder LncRNAs", "OrthoFinder Genes"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 6.4 OrthoFinder (LncRNAs)

cat("\n\nOrthoFinder (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_OrthoFinder_NR_def = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_GENES_OrthoFinder = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_genes/06-Figures/TABLE_GENES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_OrthoFinder_NR_def = TAB_LNCRNAS_OrthoFinder_NR_def[TAB_LNCRNAS_OrthoFinder_NR_def$Class == class,]
    
    # Modify genes table.
    TAB_GENES_OrthoFinder_Low = TAB_GENES_OrthoFinder
    TAB_GENES_OrthoFinder_Low$"Confidence" = "Low"
    TAB_GENES_OrthoFinder_Medium = TAB_GENES_OrthoFinder
    TAB_GENES_OrthoFinder_Medium$"Confidence" = "Medium"
    TAB_GENES_OrthoFinder_High = TAB_GENES_OrthoFinder
    TAB_GENES_OrthoFinder_High$"Confidence" = "High"
    TAB_GENES_OrthoFinder = rbind(TAB_GENES_OrthoFinder_Low, TAB_GENES_OrthoFinder_Medium, TAB_GENES_OrthoFinder_High)
    TAB_GENES_OrthoFinder$"Class" = class
    if (i == 1) {
      TAB_GENES_OrthoFinder = TAB_GENES_OrthoFinder[,c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    if (i == 2) {
      TAB_GENES_OrthoFinder = TAB_GENES_OrthoFinder[,c("Confidence", "Class", "Specie", "Conserved_level", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    if (i == 3) {
      TAB_GENES_OrthoFinder = TAB_GENES_OrthoFinder[,c("Confidence", "Class", "Specie", "Type", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    
    # Rename columns.
    if (i == 1) {
      colnames(TAB_LNCRNAS_OrthoFinder_NR_def) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_OrthoFinder) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts", "Total", "Percentage")
    }
    if (i == 2) {
      colnames(TAB_LNCRNAS_OrthoFinder_NR_def) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_OrthoFinder) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts", "Total", "Percentage")
    }
    if (i == 3) {
      colnames(TAB_LNCRNAS_OrthoFinder_NR_def) = c("Confidence", "Class", "Specie", "Type", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_OrthoFinder) = c("Confidence", "Class", "Specie", "Type", "Counts", "Total", "Percentage")
    }
    
    # Add label.
    TAB_LNCRNAS_OrthoFinder_NR_def$"Label" = "OrthoFinder LncRNAs"
    TAB_GENES_OrthoFinder$"Label" = "OrthoFinder Genes"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_OrthoFinder_NR_def, TAB_GENES_OrthoFinder)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoFinder LncRNAs", "OrthoFinder Genes"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoFinder/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 6.5 OrthoMCL (Families)

cat("\n\nOrthoMCL (Families)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/Families"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_FAM_OrthoMCL_NR_def = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_FAM_OrthoMCL_genes = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_genes/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_FAM_OrthoMCL_NR_def = TAB_FAM_OrthoMCL_NR_def[TAB_FAM_OrthoMCL_NR_def$Class == class,]
    
    # Modify genes table.
    TAB_FAM_OrthoMCL_genes_Low = TAB_FAM_OrthoMCL_genes
    TAB_FAM_OrthoMCL_genes_Low$"Confidence" = "Low"
    TAB_FAM_OrthoMCL_genes_Medium = TAB_FAM_OrthoMCL_genes
    TAB_FAM_OrthoMCL_genes_Medium$"Confidence" = "Medium"
    TAB_FAM_OrthoMCL_genes_High = TAB_FAM_OrthoMCL_genes
    TAB_FAM_OrthoMCL_genes_High$"Confidence" = "High"
    TAB_FAM_OrthoMCL_genes = rbind(TAB_FAM_OrthoMCL_genes_Low, TAB_FAM_OrthoMCL_genes_Medium, TAB_FAM_OrthoMCL_genes_High)
    TAB_FAM_OrthoMCL_genes$"Class" = class
    if (i == 1) {
      TAB_FAM_OrthoMCL_genes = TAB_FAM_OrthoMCL_genes[,c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    if (i == 2) {
      TAB_FAM_OrthoMCL_genes = TAB_FAM_OrthoMCL_genes[,c("Confidence", "Class", "Specie", "Conserved_level", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    if (i == 3) {
      TAB_FAM_OrthoMCL_genes = TAB_FAM_OrthoMCL_genes[,c("Confidence", "Class", "Specie", "Type", "Counts_Families", "Total_Families", "Percentage_Families")]
    }
    
    # Add label.
    TAB_FAM_OrthoMCL_NR_def$"Label" = "OrthoMCL LncRNAs"
    TAB_FAM_OrthoMCL_genes$"Label" = "OrthoMCL Genes"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_FAM_OrthoMCL_NR_def, TAB_FAM_OrthoMCL_genes)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoMCL LncRNAs", "OrthoMCL Genes"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/Families/", class, "-stacked_barplot_Families_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}

### 6.6 OrthoMCL (LncRNAs)

cat("\n\nOrthoMCL (LncRNAs)")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL"))
}
if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/LncRNAs"))
}

for (class in classes) {
  for (i in c(1,2,3)) {
    
    cat(paste0("\nClass: ", class, " Table: ", i))
    
    # Load the tables.
    TAB_LNCRNAS_OrthoMCL_NR_def = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    TAB_GENES_OrthoMCL = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_genes/06-Figures/TABLE_GENES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Select by class code.
    TAB_LNCRNAS_OrthoMCL_NR_def = TAB_LNCRNAS_OrthoMCL_NR_def[TAB_LNCRNAS_OrthoMCL_NR_def$Class == class,]
    
    # Modify genes table.
    TAB_GENES_OrthoMCL_Low = TAB_GENES_OrthoMCL
    TAB_GENES_OrthoMCL_Low$"Confidence" = "Low"
    TAB_GENES_OrthoMCL_Medium = TAB_GENES_OrthoMCL
    TAB_GENES_OrthoMCL_Medium$"Confidence" = "Medium"
    TAB_GENES_OrthoMCL_High = TAB_GENES_OrthoMCL
    TAB_GENES_OrthoMCL_High$"Confidence" = "High"
    TAB_GENES_OrthoMCL = rbind(TAB_GENES_OrthoMCL_Low, TAB_GENES_OrthoMCL_Medium, TAB_GENES_OrthoMCL_High)
    TAB_GENES_OrthoMCL$"Class" = class
    if (i == 1) {
      TAB_GENES_OrthoMCL = TAB_GENES_OrthoMCL[,c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    if (i == 2) {
      TAB_GENES_OrthoMCL = TAB_GENES_OrthoMCL[,c("Confidence", "Class", "Specie", "Conserved_level", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    if (i == 3) {
      TAB_GENES_OrthoMCL = TAB_GENES_OrthoMCL[,c("Confidence", "Class", "Specie", "Type", "Counts_Genes", "Total_Genes", "Percentage_Genes")]
    }
    
    # Rename columns.
    if (i == 1) {
      colnames(TAB_LNCRNAS_OrthoMCL_NR_def) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_OrthoMCL) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts", "Total", "Percentage")
    }
    if (i == 2) {
      colnames(TAB_LNCRNAS_OrthoMCL_NR_def) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_OrthoMCL) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts", "Total", "Percentage")
    }
    if (i == 3) {
      colnames(TAB_LNCRNAS_OrthoMCL_NR_def) = c("Confidence", "Class", "Specie", "Type", "Counts", "Total", "Percentage")
      colnames(TAB_GENES_OrthoMCL) = c("Confidence", "Class", "Specie", "Type", "Counts", "Total", "Percentage")
    }
    
    # Add label.
    TAB_LNCRNAS_OrthoMCL_NR_def$"Label" = "OrthoMCL LncRNAs"
    TAB_GENES_OrthoMCL$"Label" = "OrthoMCL Genes"
    
    # Join the tables.
    TAB_FINAL = rbind(TAB_LNCRNAS_OrthoMCL_NR_def, TAB_GENES_OrthoMCL)
    
    # Factors
    TAB_FINAL$Confidence = factor(TAB_FINAL$Confidence, levels = confidences)
    TAB_FINAL$Specie = factor(TAB_FINAL$Specie, levels = species)
    TAB_FINAL$Label = factor(TAB_FINAL$Label, levels = c("OrthoMCL LncRNAs", "OrthoMCL Genes"))
    
    # Figures
    if (i == 1) {
      TAB_FINAL$Number_species_by_family = factor(TAB_FINAL$Number_species_by_family, levels = 1:length(species))
      gg1 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Number_species_by_family)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 2) {
      TAB_FINAL$Conserved_level = factor(TAB_FINAL$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gg2 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Conserved_level)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
    
    if (i == 3) {
      TAB_FINAL$Type = factor(TAB_FINAL$Type, levels = c("Non-conserved", "Conserved"))
      gg3 = ggplot(TAB_FINAL, aes(x = Confidence, y = Percentage, fill = Type)) +
        geom_bar(colour = "black", position="fill", stat="identity") +
        facet_grid(Label ~ Specie) +
        scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
        scale_y_continuous(expand = c(0.02, 0.05)) +
        xlab("") +
        ylab("Percentage (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(nrow = 1)) +
        geom_text(aes(y = 1, label = Total), vjust = -0.3)
      
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".png"), height = 8, width = 16, dpi = 600)
      ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Genes_and_LncRNAs_comparison/OrthoMCL/LncRNAs/", class, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 8, width = 16, dpi = 600)
    }
  }
}


################################################################################
## 7. ACCURACY BETWEEN PIPELINES COMPARISON

cat("\n\n\nACCURACY BETWEEN PIPELINES COMPARISON")

### 7.1 LncRNAs: -max_target_seqs 1 and -max_hsps 1

cat("\n\nLncRNAs: -max_target_seqs 1 and -max_hsps 1")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison/LncRNAs"))
}

for (class in classes) {
    
  cat(paste0("\nClass: ", class))
  
  # Load the tables.
  TAB_LNCRNAS_BLASTN = read.table(paste0(WD1, "/Blastn/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
  TAB_LNCRNAS_ORTHOFINDER = read.table(paste0(WD1, "/OrthoFinder/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
  TAB_LNCRNAS_ORTHOMCL = read.table(paste0(WD1, "/OrthoMCL/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
  
  # Select by class code.
  TAB_LNCRNAS_Blastn_red = TAB_LNCRNAS_BLASTN[TAB_LNCRNAS_BLASTN$Class == class,]
  TAB_LNCRNAS_OrthoFinder_red = TAB_LNCRNAS_ORTHOFINDER[TAB_LNCRNAS_ORTHOFINDER$Class == class,]
  TAB_LNCRNAS_OrthoMCL_red = TAB_LNCRNAS_ORTHOMCL[TAB_LNCRNAS_ORTHOMCL$Class == class,]
  
  # Figures
  ids_B = TAB_LNCRNAS_Blastn_red[TAB_LNCRNAS_Blastn_red$Conserved_level != "Non-conserved", "Member"]
  ids_OF = TAB_LNCRNAS_OrthoFinder_red[TAB_LNCRNAS_OrthoFinder_red$Conserved_level != "Non-conserved", "Member"]
  ids_OM = TAB_LNCRNAS_OrthoMCL_red[TAB_LNCRNAS_OrthoMCL_red$Conserved_level != "Non-conserved", "Member"]
    
  x = list(Blastn = ids_B, OrthoFinder = ids_OF, OrthoMCL = ids_OM)
    
  png(filename = paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison/LncRNAs/1-", class, ".png"), width = 3000, height = 3000, res = 500)
  print(ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#abd6bf"), stroke_size = 0.5, set_name_size = 4))
  invisible(dev.off())
  
  # Figures by Conserved_level
  for (level in c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved")) {
    
    ids_B = TAB_LNCRNAS_Blastn_red[TAB_LNCRNAS_Blastn_red$Conserved_level == level, "Member"]
    ids_OF = TAB_LNCRNAS_OrthoFinder_red[TAB_LNCRNAS_OrthoFinder_red$Conserved_level == level, "Member"]
    ids_OM = TAB_LNCRNAS_OrthoMCL_red[TAB_LNCRNAS_OrthoMCL_red$Conserved_level == level, "Member"]
    
    x = list(Blastn = ids_B, OrthoFinder = ids_OF, OrthoMCL = ids_OM)
    
    png(filename = paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison/LncRNAs/1-", class, "-", level, ".png"), width = 3000, height = 3000, res = 500)
    print(ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#abd6bf"), stroke_size = 0.5, set_name_size = 4))
    invisible(dev.off())
  }
}

### 7.2 LncRNAs: -max_target_seqs 5 and -max_hsps 5

cat("\n\nLncRNAs: -max_target_seqs 5 and -max_hsps 5")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison/LncRNAs"))
}

for (class in classes) {
  
  cat(paste0("\nClass: ", class))
  
  # Load the tables.
  TAB_LNCRNAS_BLASTN = read.table(paste0(WD1, "/Blastn/nr/Prueba_e-value/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
  TAB_LNCRNAS_ORTHOFINDER = read.table(paste0(WD1, "/OrthoFinder/nr/Prueba_e-value/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
  TAB_LNCRNAS_ORTHOMCL = read.table(paste0(WD1, "/OrthoMCL/nr/Prueba_e-value/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
  
  # Select by class code.
  TAB_LNCRNAS_Blastn_red = TAB_LNCRNAS_BLASTN[TAB_LNCRNAS_BLASTN$Class == class,]
  TAB_LNCRNAS_OrthoFinder_red = TAB_LNCRNAS_ORTHOFINDER[TAB_LNCRNAS_ORTHOFINDER$Class == class,]
  TAB_LNCRNAS_OrthoMCL_red = TAB_LNCRNAS_ORTHOMCL[TAB_LNCRNAS_ORTHOMCL$Class == class,]
  
  # Figures by Conserved_level
  for (level in c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved")) {
    
    ids_B = TAB_LNCRNAS_Blastn_red[TAB_LNCRNAS_Blastn_red$Conserved_level == level, "Member"]
    ids_OF = TAB_LNCRNAS_OrthoFinder_red[TAB_LNCRNAS_OrthoFinder_red$Conserved_level == level, "Member"]
    ids_OM = TAB_LNCRNAS_OrthoMCL_red[TAB_LNCRNAS_OrthoMCL_red$Conserved_level == level, "Member"]
    
    x = list(Blastn = ids_B, OrthoFinder = ids_OF, OrthoMCL = ids_OM)
    
    png(filename = paste0(WD1, "/COMPARISON_SOFTWARES/Accuracy_comparison/LncRNAs/5-", class, "-", level, ".png"), width = 3000, height = 3000, res = 500)
    print(ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#abd6bf"), stroke_size = 0.5, set_name_size = 4))
    invisible(dev.off())
  }
}

cat("\n")


################################################################################
## 8. OTHER FORMAT FIGURES

cat("\n\n\nOTHER FORMAT FIGURES")

### 8.1 Families

cat("\n\nFamilies")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/Families"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/Families"))
}

for (pipeline in c("Blastn", "OrthoFinder", "OrthoMCL")) {
  for (confidence in confidences) {
    for (i in c(1,2,3)) {
      
      cat(paste0("\nPipeline: ", pipeline, " Confidence: ", confidence, " Table: ", i))
      
      # Load the tables.
      TAB_FAM_NR_def = read.table(paste0(WD1, "/", pipeline, "/nr/Definitive/06-Figures/TABLE_FAMILIES_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
      
      # Select by confidence level.
      TAB_FAM_NR_def = TAB_FAM_NR_def[TAB_FAM_NR_def$Confidence == confidence,]
      
      # Factors
      TAB_FAM_NR_def$Class = factor(TAB_FAM_NR_def$Class, levels = classes)
      TAB_FAM_NR_def$Specie = factor(TAB_FAM_NR_def$Specie, levels = species)
      
      # Figures
      if (i == 1) {
        TAB_FAM_NR_def$Number_species_by_family = factor(TAB_FAM_NR_def$Number_species_by_family, levels = 1:length(species))
        gg1 = ggplot(TAB_FAM_NR_def, aes(x = Specie, y = Percentage_Families, fill = Number_species_by_family)) +
          geom_bar(colour = "black", position="fill", stat="identity") +
          facet_grid(Class ~ .) +
          scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
          scale_y_continuous(expand = c(0.02, 0.08)) +
          xlab("") +
          ylab("Percentage (%)") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
          theme(legend.position = "top") +
          guides(fill = guide_legend(nrow = 2, title = "Conservation level")) +
          geom_text(aes(y = 1, label = Total_Families), vjust = -0.3, size = 3)
        
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/Families/", pipeline, "-", confidence, "-stacked_barplot_Families_", i, ".png"), height = 12, width = 6, dpi = 700)
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/Families/", pipeline, "-", confidence, "-stacked_barplot_Families_", i, ".pdf"), height = 12, width = 6, dpi = 700)
      }
      
      if (i == 2) {
        TAB_FAM_NR_def$Conserved_level = factor(TAB_FAM_NR_def$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
        gg2 = ggplot(TAB_FAM_NR_def, aes(x = Specie, y = Percentage_Families, fill = Conserved_level)) +
          geom_bar(colour = "black", position="fill", stat="identity") +
          facet_grid(Class ~ .) +
          scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
          scale_y_continuous(expand = c(0.02, 0.08)) +
          xlab("") +
          ylab("Percentage (%)") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
          theme(legend.position = "top") +
          guides(fill = guide_legend(nrow = 2, title = "Conservation level")) +
          geom_text(aes(y = 1, label = Total_Families), vjust = -0.3, size = 3)
        
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/Families/", pipeline, "-", confidence, "-stacked_barplot_Families_", i, ".png"), height = 12, width = 6, dpi = 700)
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/Families/", pipeline, "-", confidence, "-stacked_barplot_Families_", i, ".pdf"), height = 12, width = 6, dpi = 700)
      }
      
      if (i == 3) {
        TAB_FAM_NR_def$Type = factor(TAB_FAM_NR_def$Type, levels = c("Non-conserved", "Conserved"))
        gg3 = ggplot(TAB_FAM_NR_def, aes(x = Specie, y = Percentage_Families, fill = Type)) +
          geom_bar(colour = "black", position="fill", stat="identity") +
          facet_grid(Class ~ .) +
          scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
          scale_y_continuous(expand = c(0.02, 0.08)) +
          xlab("") +
          ylab("Percentage (%)") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
          theme(legend.position = "top") +
          guides(fill = guide_legend(nrow = 2, title = "Conservation level")) +
          geom_text(aes(y = 1, label = Total_Families), vjust = -0.3, size = 3)
        
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/Families/", pipeline, "-", confidence, "-stacked_barplot_Families_", i, ".png"), height = 12, width = 6, dpi = 700)
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/Families/", pipeline, "-", confidence, "-stacked_barplot_Families_", i, ".pdf"), height = 12, width = 6, dpi = 700)
      }
    }
  }
}


### 8.2 LncRNAs

cat("\n\nLncRNAs")

if (!dir.exists(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/LncRNAs"))){
  dir.create(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/LncRNAs"))
}

for (pipeline in c("Blastn", "OrthoFinder", "OrthoMCL")) {
  for (confidence in confidences) {
    for (i in c(1,2,3)) {
      
      cat(paste0("\nPipeline: ", pipeline, " Confidence: ", confidence, " Table: ", i))
      
      # Load the tables.
      TAB_LNCRNAS_NR_def = read.table(paste0(WD1, "/", pipeline, "/nr/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_", i, ".tsv"), sep = "\t", header = T, quote = "\"")
      
      # Select by confidence level.
      TAB_LNCRNAS_NR_def = TAB_LNCRNAS_NR_def[TAB_LNCRNAS_NR_def$Confidence == confidence,]
      
      # Factors
      TAB_LNCRNAS_NR_def$Class = factor(TAB_LNCRNAS_NR_def$Class, levels = classes)
      TAB_LNCRNAS_NR_def$Specie = factor(TAB_LNCRNAS_NR_def$Specie, levels = species)
      
      # Figures
      if (i == 1) {
        TAB_LNCRNAS_NR_def$Number_species_by_family = factor(TAB_LNCRNAS_NR_def$Number_species_by_family, levels = 1:length(species))
        gg1 = ggplot(TAB_LNCRNAS_NR_def, aes(x = Specie, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
          geom_bar(colour = "black", position="fill", stat="identity") +
          facet_grid(Class ~ .) +
          scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
          scale_y_continuous(expand = c(0.02, 0.08)) +
          xlab("") +
          ylab("Percentage (%)") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
          theme(legend.position = "top") +
          guides(fill = guide_legend(nrow = 2, title = "Conservation level")) +
          geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3, size = 3)
        
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/LncRNAs/", pipeline, "-", confidence, "-stacked_barplot_LncRNAs_", i, ".png"), height = 12, width = 5, dpi = 700)
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/LncRNAs/", pipeline, "-", confidence, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 12, width = 5, dpi = 700)
      }
      
      if (i == 2) {
        TAB_LNCRNAS_NR_def$Conserved_level = factor(TAB_LNCRNAS_NR_def$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
        gg2 = ggplot(TAB_LNCRNAS_NR_def, aes(x = Specie, y = Percentage_LncRNAs, fill = Conserved_level)) +
          geom_bar(colour = "black", position="fill", stat="identity") +
          facet_grid(Class ~ .) +
          scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
          scale_y_continuous(expand = c(0.02, 0.08)) +
          xlab("") +
          ylab("Percentage (%)") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
          theme(legend.position = "top") +
          guides(fill = guide_legend(nrow = 2, title = "Conservation level")) +
          geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3, size = 3)
        
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/LncRNAs/", pipeline, "-", confidence, "-stacked_barplot_LncRNAs_", i, ".png"), height = 12, width = 5, dpi = 700)
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/LncRNAs/", pipeline, "-", confidence, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 12, width = 5, dpi = 700)
      }
      
      if (i == 3) {
        TAB_LNCRNAS_NR_def$Type = factor(TAB_LNCRNAS_NR_def$Type, levels = c("Non-conserved", "Conserved"))
        gg3 = ggplot(TAB_LNCRNAS_NR_def, aes(x = Specie, y = Percentage_LncRNAs, fill = Type)) +
          geom_bar(colour = "black", position="fill", stat="identity") +
          facet_grid(Class ~ .) +
          scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
          scale_y_continuous(expand = c(0.02, 0.08)) +
          xlab("") +
          ylab("Percentage (%)") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
          theme(legend.position = "top") +
          guides(fill = guide_legend(nrow = 2, title = "Conservation level")) +
          geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3, size = 3)
        
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/LncRNAs/", pipeline, "-", confidence, "-stacked_barplot_LncRNAs_", i, ".png"), height = 12, width = 5, dpi = 700)
        ggsave(paste0(WD1, "/COMPARISON_SOFTWARES/Other_formats/LncRNAs/", pipeline, "-", confidence, "-stacked_barplot_LncRNAs_", i, ".pdf"), height = 12, width = 5, dpi = 700)
      }
    }
  }
}

cat("\n")

