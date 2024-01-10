################################################################################
#
# FIGURES NON-REDUNDANT
#
# Join plots comparing genes and lncRNAs: GC content, Length, Exon number, TPMs 
# and Repeat content.
#
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggExtra))
suppressMessages(library(data.table))
suppressMessages(library(grid))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(options(bitmapType='cairo'))


## 1. VARIABLES

WD = "storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/07-comparison_lncRNAs_vs_coding_genes"
repeat_content_path = "storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-TEs_and_genomic_repeats/02-Comparison_Genes_LncRNAs/Figures_and_Tables/C"
species_short_name = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
species_long_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")
types = c("LncRNAs low", "LncRNAs medium", "LncRNAs high")

if (!dir.exists(paste0(WD, "/ALL"))){
  dir.create(paste0(WD, "/ALL"))
}
if (!dir.exists(paste0(WD, "/ALL/DEFINITIVE"))){
  dir.create(paste0(WD, "/ALL/DEFINITIVE"))
}










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## 2. GLOBAL TABLES

cat(paste0("\n\n\nGLOBAL TABLES..."))

TAB_LONG_FINAL = data.frame()
TAB_WIDE_FINAL = data.frame()

for (type in types) {
  cat(paste0("\n\n\nType: ", type, "..."))
  
  ## 2.1 CREATE LONG AND WIDE TABLE.
  
  ### 2.1.1 TAB LONG
  
  cat(paste0("\n\n\t-LONG TABLE..."))
  TAB_LONG = data.frame()
  
  #### GC CONTENT
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-GC CONTENT: Spe: ", species_short_name[i], "..."))
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "GC")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "GC")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "GC")]
    
    tab = rbind(L, G, IR)
    tab$"Feature" = "GC content"
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  #### EXON NUMBER
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-EXON NUMBER: Spe: ", species_short_name[i], "..."))
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Exons" = NA
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons")]
    
    tab = rbind(L, G, IR)
    tab$"Feature" = "Exon number"
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  #### LENGTH
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-LENGTH: Spe: ", species_short_name[i], "..."))
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Length")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Length")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Length" = NA
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Length")]
    
    tab = rbind(L, G, IR)
    tab$"Feature" = "Length"
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  #### EXPRESSION
  
  # All the species have some NA values due to the detection of duplicated transcripts in the quantification with salmon. These duplicated 
  # transcripts are removed and then they are not quantified. This happens in the same way with LncRNAs and genes. This kind of duplication 
  # isn't intra-locus, so they weren't removed in the redundancy filter. It exists other region in the genome where exists a transcript equal 
  # to this transcript. It happens few times. They will not be plotted because they are kept as NA values.
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-EXPRESSION: Spe: ", species_short_name[i], "..."))
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "TPMs.mean")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "TPMs.mean")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"TPMs.mean" = NA
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "TPMs.mean")]
    
    tab = rbind(L, G, IR)
    tab$"Feature" = "Expression"
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  #### REPEAT CONTENT
  
  # The table coming from repeat content analysis only contains transcripts with more than 0% of repeat content. So, when we merge
  # tables we find NA values which will be converted to 0 because they are transcripts with 0% of repeat content.
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-REPEAT CONTENT: Spe: ", species_short_name[i], "..."))
    
    tab_rep = read.table(paste0(repeat_content_path, "/Final_tab-Repeat-NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-COLLAPSED_REPEAT.tsv"), sep = "\t", header = T, quote = "\"")
    tab_rep = tab_rep[tab_rep$spe == species_short_name[i], c("transcript_id", "overlap_per")]
    colnames(tab_rep) = c("ID_transcript", "overlap_per")
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie")]
    
    tab = rbind(L, G, IR)
    tab = merge(tab, tab_rep, by = "ID_transcript", all = T)
    tab$"Feature" = "Repeat content"
    tab = tab[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code",  "Type", "Specie", "overlap_per", "Feature")]
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    tab[is.na(tab)] = 0
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  ### 2.1.2 TAB WIDE
  
  cat(paste0("\n\n\t-WIDE TABLE..."))
  TAB_WIDE = data.frame()
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-Spe: ", species_short_name[i], "..."))
    
    tab_rep = read.table(paste0(repeat_content_path, "/Final_tab-Repeat-NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-COLLAPSED_REPEAT.tsv"), sep = "\t", header = T, quote = "\"")
    tab_rep = tab_rep[tab_rep$spe == species_short_name[i], c("transcript_id", "overlap_per")]
    colnames(tab_rep) = c("ID_transcript", "overlap_per")
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons", "Length", "GC", "TPMs.mean")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons", "Length", "GC", "TPMs.mean")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Length" = NA
    IR$"Exons" = NA
    IR$"TPMs.mean" = NA
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons", "Length", "GC", "TPMs.mean")]
    
    tab = rbind(L, G, IR)
    tab = merge(tab, tab_rep, by = "ID_transcript", all = T)
    colnames(tab) = c("ID_transcript", "Chr", "Start", "End", "Strand", "Origin", "Class_code", "Type.1", "Specie", "Exons", "Length", "GC", "TPMs.mean", "overlap_per")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Exons", "Length", "GC", "TPMs.mean", "overlap_per")]
    
    TAB_WIDE = rbind(TAB_WIDE, tab)
  }
  
  
  # Modify the class codes.
  TAB_LONG[TAB_LONG == "intergenic (u)"] = "u"
  TAB_LONG[TAB_LONG == "antisense (x)"] = "x"
  TAB_LONG[TAB_LONG == "intronic (i)"] = "i"
  TAB_LONG[TAB_LONG == "sense (o/e)"] = "o/e"
  TAB_LONG[TAB_LONG == "gene (=)"] = "pc"
  
  TAB_WIDE[TAB_WIDE == "intergenic (u)"] = "u"
  TAB_WIDE[TAB_WIDE == "antisense (x)"] = "x"
  TAB_WIDE[TAB_WIDE == "intronic (i)"] = "i"
  TAB_WIDE[TAB_WIDE == "sense (o/e)"] = "o/e"
  TAB_WIDE[TAB_WIDE == "gene (=)"] = "pc"
  
  # Save.
  write.table(TAB_LONG, paste0(WD, "/ALL/DEFINITIVE/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-TAB_LONG.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(TAB_WIDE, paste0(WD, "/ALL/DEFINITIVE/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-TAB_WIDE.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(list = c("tab_1", "tab_2", "tab", "i", "L", "G", "IR", "tab_rep"))
  
  ## 2.2 FINAL TABLES
  
  cat(paste0("\n\n\t-JOIN TABLES TO FINAL TABLES..."))
  TAB_LONG_FINAL = rbind(TAB_LONG_FINAL, TAB_LONG)
  TAB_WIDE_FINAL = rbind(TAB_WIDE_FINAL, TAB_WIDE)
  
  rm(list = c("TAB_LONG", "TAB_WIDE"))
}

rm(list = c("type"))

# Remove duplicated rows. For example, genes and intergenic regions in each iteration (confidence level) are the same, so they are repeated.
# In summary tables, we don't remove duplicated rows (genes and intergenic regions) to make easier the comparison between genes, lncRNAs and 
# intergenic regions.
TAB_LONG_FINAL = TAB_LONG_FINAL[!duplicated(TAB_LONG_FINAL),]
TAB_WIDE_FINAL = TAB_WIDE_FINAL[!duplicated(TAB_WIDE_FINAL),]

write.table(TAB_LONG_FINAL, paste0(WD, "/ALL/DEFINITIVE/NR-TAB_LONG-FINAL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_FINAL, paste0(WD, "/ALL/DEFINITIVE/NR-TAB_WIDE-FINAL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("TAB_LONG_FINAL", "TAB_WIDE_FINAL"))










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## 3.A FIGURES, SUMMARY TABLES AND STATISTICS "A" (SPECIES AND CLASS_CODES)

cat(paste0("\n\n\nFIGURES, SUMMARY TABLES AND STATISTICS \"A\" (SPECIES AND CLASS_CODES)..."))

if (!dir.exists(paste0(WD, "/ALL/DEFINITIVE/A"))){
  dir.create(paste0(WD, "/ALL/DEFINITIVE/A"))
}

TAB_LONG_SUMMARY_FINAL = data.frame()
TAB_WIDE_SUMMARY_FINAL = data.frame()
TAB_STATISTICS_FINAL = data.frame()

for (type in types) {
  cat(paste0("\n\n\nType: ", type, "..."))
  
  ## 3.A.1 LOAD TABLES
  
  cat(paste0("\n\n\t-LOAD LONG AND WIDE TABLES..."))
  
  # Load tables.
  TAB_LONG = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-TAB_LONG.tsv"), sep = "\t", header = T, quote = "\"")
  TAB_WIDE = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-TAB_WIDE.tsv"), sep = "\t", header = T, quote = "\"")
  
  # Create factors
  TAB_LONG$Type.1 = factor(TAB_LONG$Type.1, levels = c("Genes", type, "Intergenic Regions"))
  TAB_LONG$Type.2 = factor(TAB_LONG$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
  TAB_LONG$Class_code = factor(TAB_LONG$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  TAB_LONG$Specie = factor(TAB_LONG$Specie, levels = species_long_name)
  TAB_LONG$Feature = factor(TAB_LONG$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
  
  TAB_WIDE$Type.1 = factor(TAB_WIDE$Type.1, levels = c("Genes", type, "Intergenic Regions"))
  TAB_WIDE$Type.2 = factor(TAB_WIDE$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
  TAB_WIDE$Class_code = factor(TAB_WIDE$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  TAB_WIDE$Specie = factor(TAB_WIDE$Specie, levels = species_long_name)
  
  ## 3.A.2 FIGURES
  
  cat(paste0("\n\n\t-FIGURES..."))
  
  my_mean = function(x) {
    log10(mean(10^x))
  }
  
  #### GC CONTENT
  
  cat(paste0("\n\t\t-GC CONTENT..."))
  
  # Grid
  GC_content = TAB_LONG[TAB_LONG$Feature == "GC content",]
  gg1 = ggplot(GC_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")) +
    facet_grid(Feature~Specie) +
    xlab("") +
    ylab("GC content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold.italic"),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 5.5, 0, 5.5))
  
  gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 2, widths = c(0.0046, 0.9954))
  
  rm(list = c("GC_content"))
  
  #### EXON NUMBER
  
  cat(paste0("\n\t\t-EXON NUMBER..."))
  
  # Grid
  Exon_number = TAB_LONG[TAB_LONG$Feature == "Exon number",]
  gg2 = ggplot(Exon_number, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")) +
    facet_grid(Feature~Specie) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Exon number") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))
  
  rm(list = c("Exon_number"))
  
  #### LENGTH
  
  cat(paste0("\n\t\t-LENGTH..."))
  
  # Grid
  Length = TAB_LONG[TAB_LONG$Feature == "Length",]
  gg3 = ggplot(Length, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")) +
    facet_grid(Feature~Specie) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Length") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))
  
  rm(list = c("Length"))
  
  #### EXPRESSION
  
  cat(paste0("\n\t\t-EXPRESSION..."))
  
  # Grid
  Expression = TAB_LONG[TAB_LONG$Feature == "Expression",]
  gg4 = ggplot(Expression, aes(x = Class_code, y = Value, fill = Class_code)) + 
    geom_boxplot() + 
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")) +
    facet_grid(Feature~Specie) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("TPM") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))
  
  gg4 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg4, ncol = 2, widths = c(0.002, 0.998))
  
  rm(list = c("Expression"))
  
  #### REPEAT CONTENT
  
  cat(paste0("\n\t\t-REPEAT CONTENT..."))
  
  # Grid
  Repeat_content = TAB_LONG[TAB_LONG$Feature == "Repeat content",]
  gg5 = ggplot(Repeat_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
    ) +
    facet_grid(Feature~Specie) +
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
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(0, 5.5, 5.5, 5.5))
  
  gg5 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg5, ncol = 2, widths = c(0.005, 0.995))
  
  rm(list = c("Repeat_content"))
  
  #### FIGURE FEATURES
  
  # Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
  # el siguiente error:
  # 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
  # 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
  # Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
  # -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
  # a logaritmo serian -Inf.
  # 5: Transformation introduced infinite values in continuous y-axis.
  
  cat(paste0("\n\t\t-ALL FEATURES..."))
  
  gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 1, heights = c(0.21, 0.19, 0.19, 0.19, 0.22))
  
  ggsave(paste0(WD, "/ALL/DEFINITIVE/A/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-A.png"), height = 25, width = 25, dpi = 600)
  ggsave(paste0(WD, "/ALL/DEFINITIVE/A/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-A.pdf"), height = 25, width = 25, dpi = 600)
  
  rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))
  
  ## 3.A.3 MEAN AND MEDIAN TABLES
  
  cat(paste0("\n\n\t-MEAN AND MEDIAN TABLE..."))
  
  TAB_LONG_SUMMARY = suppressMessages(TAB_LONG %>% 
                                        drop_na() %>%
                                        group_by(Specie, Type.1, Class_code, Feature) %>% 
                                        summarise(MEAN = mean(Value),
                                                  MEDIAN = median(Value)))
  TAB_WIDE_SUMMARY = pivot_wider(TAB_LONG_SUMMARY, 
                                 id_cols = c("Specie", "Type.1", "Class_code"),
                                 names_from = "Feature",
                                 values_from = c("MEAN", "MEDIAN"),
                                 names_glue = "{Feature}.{.value}")
  colnames(TAB_WIDE_SUMMARY) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY))
  
  TAB_LONG_SUMMARY = as.data.frame(TAB_LONG_SUMMARY)
  TAB_WIDE_SUMMARY = as.data.frame(TAB_WIDE_SUMMARY)
  
  write.table(TAB_LONG_SUMMARY, paste0(WD, "/ALL/DEFINITIVE/A/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-SUMMARY-TAB_LONG-A.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(TAB_WIDE_SUMMARY, paste0(WD, "/ALL/DEFINITIVE/A/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-SUMMARY-TAB_WIDE-A.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  ## 3.A.4 STATISTICAL ANALYSIS
  
  # Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
  # necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como 
  # el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
  # indipendientes. Algunos de estos papers son:
  # - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
  # - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)
  
  cat(paste0("\n\n\t-STATISTICAL ANALYSIS..."))
  
  combinations = as.data.frame(t(combn(c("pc", "u", "x", "i", "o/e", "ir"), 2)))
  rownames(combinations) = NULL
  colnames(combinations) = c("cl1", "cl2")
  
  TAB_STATISTICS = data.frame()
  for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
    for (spe in species_long_name) {
      for (i in 1:nrow(combinations)) {
        cl1 = combinations[i, "cl1"]
        cl2 = combinations[i, "cl2"]
        subset1 = TAB_LONG[TAB_LONG$Feature == feature & TAB_LONG$Specie == spe & TAB_LONG$Class_code == cl1 & !is.na(TAB_LONG$Value), "Value"]
        subset2 = TAB_LONG[TAB_LONG$Feature == feature & TAB_LONG$Specie == spe & TAB_LONG$Class_code == cl2 & !is.na(TAB_LONG$Value), "Value"]
        L1 = length(subset1)
        L2 = length(subset2)
        if (L1 > 0 & L2 > 0) {
          test = wilcox.test(subset1, subset2, paired = FALSE)
          Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
          row = data.frame(Type.1 = type, Specie = spe, Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
        } 
        else{
          row = data.frame(Type.1 = type, Specie = spe, Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
        }
        TAB_STATISTICS = rbind(TAB_STATISTICS, row)
      }
    }
  }
  
  rownames(TAB_STATISTICS) = NULL
  
  write.table(TAB_STATISTICS, paste0(WD, "/ALL/DEFINITIVE/A/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-TAB_STATISTICS-A.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(list = c("feature", "spe", "i", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))
  
  ## 3.A.5 FINAL TABLES
  
  cat(paste0("\n\n\t-JOIN TABLES TO FINAL TABLES..."))
  
  TAB_LONG_SUMMARY_FINAL = rbind(TAB_LONG_SUMMARY_FINAL, TAB_LONG_SUMMARY)
  TAB_WIDE_SUMMARY_FINAL = rbind(TAB_WIDE_SUMMARY_FINAL, TAB_WIDE_SUMMARY)
  TAB_STATISTICS_FINAL = rbind(TAB_STATISTICS_FINAL, TAB_STATISTICS)
  
  rm(list = c("TAB_LONG", "TAB_WIDE", "TAB_LONG_SUMMARY", "TAB_WIDE_SUMMARY", "TAB_STATISTICS"))
}

rm(list = c("type"))

# Remove duplicated rows. For example, genes and intergenic regions in each iteration (confidence level) are the same, so they are repeated.
TAB_LONG_SUMMARY_FINAL = TAB_LONG_SUMMARY_FINAL[!duplicated(TAB_LONG_SUMMARY_FINAL),]
TAB_WIDE_SUMMARY_FINAL = TAB_WIDE_SUMMARY_FINAL[!duplicated(TAB_WIDE_SUMMARY_FINAL),]

# Sort tables.
TAB_LONG_SUMMARY_FINAL = TAB_LONG_SUMMARY_FINAL[order(TAB_LONG_SUMMARY_FINAL$Type.1, TAB_LONG_SUMMARY_FINAL$Specie),]
TAB_WIDE_SUMMARY_FINAL = TAB_WIDE_SUMMARY_FINAL[order(TAB_WIDE_SUMMARY_FINAL$Type.1, TAB_WIDE_SUMMARY_FINAL$Specie),]

# Save.
write.table(TAB_LONG_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/A/NR-SUMMARY-TAB_LONG_FINAL-A.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/A/NR-SUMMARY-TAB_WIDE_FINAL-A.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_STATISTICS_FINAL, paste0(WD, "/ALL/DEFINITIVE/A/NR-TAB_STATISTICS_FINAL-A.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("TAB_LONG_SUMMARY_FINAL", "TAB_WIDE_SUMMARY_FINAL", "TAB_STATISTICS_FINAL", "my_mean"))










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## 3.B FIGURES, SUMMARY TABLES AND STATISTICS "B" (CLASS_CODES)

cat(paste0("\n\n\nFIGURES, SUMMARY TABLES AND STATISTICS \"B\" (CLASS_CODES)..."))

if (!dir.exists(paste0(WD, "/ALL/DEFINITIVE/B"))){
  dir.create(paste0(WD, "/ALL/DEFINITIVE/B"))
}

TAB_LONG_SUMMARY_FINAL = data.frame()
TAB_WIDE_SUMMARY_FINAL = data.frame()
TAB_STATISTICS_FINAL = data.frame()

for (type in types) {
  cat(paste0("\n\n\nType: ", type, "..."))
  
  ## 3.B.1 LOAD TABLES
  
  cat(paste0("\n\n\t-LOAD LONG AND WIDE TABLES..."))
  
  # Load tables.
  TAB_LONG = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-TAB_LONG.tsv"), sep = "\t", header = T, quote = "\"")
  TAB_WIDE = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-TAB_WIDE.tsv"), sep = "\t", header = T, quote = "\"")
  
  # Create factors
  TAB_LONG$Type.1 = factor(TAB_LONG$Type.1, levels = c("Genes", type, "Intergenic Regions"))
  TAB_LONG$Type.2 = factor(TAB_LONG$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
  TAB_LONG$Class_code = factor(TAB_LONG$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  TAB_LONG$Specie = factor(TAB_LONG$Specie, levels = species_long_name)
  TAB_LONG$Feature = factor(TAB_LONG$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
  
  TAB_WIDE$Type.1 = factor(TAB_WIDE$Type.1, levels = c("Genes", type, "Intergenic Regions"))
  TAB_WIDE$Type.2 = factor(TAB_WIDE$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
  TAB_WIDE$Class_code = factor(TAB_WIDE$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  TAB_WIDE$Specie = factor(TAB_WIDE$Specie, levels = species_long_name)
  
  ## 3.B.2 FIGURES
  
  cat(paste0("\n\n\t-FIGURES..."))
  
  my_mean = function(x) {
    log10(mean(10^x))
  }
  
  #### GC CONTENT
  
  cat(paste0("\n\t\t-GC CONTENT..."))
  
  # Grid
  GC_content = TAB_LONG[TAB_LONG$Feature == "GC content",]
  gg1 = ggplot(GC_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    xlab("") +
    ylab("GC content (%)") +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1)) +
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
  
  # Grid
  Exon_number = TAB_LONG[TAB_LONG$Feature == "Exon number",]
  gg2 = ggplot(Exon_number, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Exon number") +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1)) +
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
  
  # Grid
  Length = TAB_LONG[TAB_LONG$Feature == "Length",]
  gg3 = ggplot(Length, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Length") +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1)) +
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
  
  # Grid
  Expression = TAB_LONG[TAB_LONG$Feature == "Expression",]
  gg4 = ggplot(Expression, aes(x = Class_code, y = Value, fill = Class_code)) + 
    geom_boxplot() + 
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
    ) +
    facet_grid(.~Feature) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("TPM") +
    theme_bw() + 
    guides(fill = guide_legend(nrow = 1)) +
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
  
  # Grid
  Repeat_content = TAB_LONG[TAB_LONG$Feature == "Repeat content",]
  gg5 = ggplot(Repeat_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
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
  
  # Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
  # el siguiente error:
  # 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
  # 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
  # Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
  # -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
  # a logaritmo serian -Inf.
  # 5: Transformation introduced infinite values in continuous y-axis.
  
  cat(paste0("\n\t\t-ALL FEATURES..."))
  
  gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 5, common.legend = T, legend = "bottom")
  
  ggsave(paste0(WD, "/ALL/DEFINITIVE/B/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-B.png"), height = 6, width = 25, dpi = 600, bg = "white")
  ggsave(paste0(WD, "/ALL/DEFINITIVE/B/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-B.pdf"), height = 6, width = 25, dpi = 600, bg = "white")
  
  rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))
  
  ## 3.B.3 MEAN AND MEDIAN TABLES

  cat(paste0("\n\n\t-MEAN AND MEDIAN TABLE..."))

  TAB_LONG_SUMMARY = suppressMessages(TAB_LONG %>%
                                        drop_na() %>%
                                        group_by(Type.1, Class_code, Feature) %>%
                                        summarise(MEAN = mean(Value),
                                                  MEDIAN = median(Value)))
  TAB_WIDE_SUMMARY = pivot_wider(TAB_LONG_SUMMARY,
                                 id_cols = c("Type.1", "Class_code"),
                                 names_from = "Feature",
                                 values_from = c("MEAN", "MEDIAN"),
                                 names_glue = "{Feature}.{.value}")
  colnames(TAB_WIDE_SUMMARY) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY))
  
  TAB_LONG_SUMMARY = as.data.frame(TAB_LONG_SUMMARY)
  TAB_WIDE_SUMMARY = as.data.frame(TAB_WIDE_SUMMARY)

  write.table(TAB_LONG_SUMMARY, paste0(WD, "/ALL/DEFINITIVE/B/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-SUMMARY-TAB_LONG-B.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(TAB_WIDE_SUMMARY, paste0(WD, "/ALL/DEFINITIVE/B/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-SUMMARY-TAB_WIDE-B.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  ## 3.B.4 STATISTICAL ANALYSIS

  # Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
  # necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como
  # el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
  # indipendientes. Algunos de estos papers son:
  # - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
  # - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)

  cat(paste0("\n\n\t-STATISTICAL ANALYSIS..."))

  combinations = as.data.frame(t(combn(c("pc", "u", "x", "i", "o/e", "ir"), 2)))
  rownames(combinations) = NULL
  colnames(combinations) = c("cl1", "cl2")

  TAB_STATISTICS = data.frame()
  for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
    for (i in 1:nrow(combinations)) {
      cl1 = combinations[i, "cl1"]
      cl2 = combinations[i, "cl2"]
      subset1 = TAB_LONG[TAB_LONG$Feature == feature & TAB_LONG$Class_code == cl1 & !is.na(TAB_LONG$Value), "Value"]
      subset2 = TAB_LONG[TAB_LONG$Feature == feature & TAB_LONG$Class_code == cl2 & !is.na(TAB_LONG$Value), "Value"]
      L1 = length(subset1)
      L2 = length(subset2)
      if (L1 > 0 & L2 > 0) {
        test = wilcox.test(subset1, subset2, paired = FALSE)
        Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
        row = data.frame(Type.1 = type, Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
      }
      else{
        row = data.frame(Type.1 = type, Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
      }
      TAB_STATISTICS = rbind(TAB_STATISTICS, row)
    }
  }

  rownames(TAB_STATISTICS) = NULL

  write.table(TAB_STATISTICS, paste0(WD, "/ALL/DEFINITIVE/B/NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-TAB_STATISTICS-B.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

  rm(list = c("feature", "i", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

  ## 3.B.5 FINAL TABLES

  cat(paste0("\n\n\t-JOIN TABLES TO FINAL TABLES..."))
  
  TAB_LONG_SUMMARY_FINAL = rbind(TAB_LONG_SUMMARY_FINAL, TAB_LONG_SUMMARY)
  TAB_WIDE_SUMMARY_FINAL = rbind(TAB_WIDE_SUMMARY_FINAL, TAB_WIDE_SUMMARY)
  TAB_STATISTICS_FINAL = rbind(TAB_STATISTICS_FINAL, TAB_STATISTICS)

  rm(list = c("TAB_LONG", "TAB_WIDE", "TAB_LONG_SUMMARY", "TAB_WIDE_SUMMARY", "TAB_STATISTICS"))
}

rm(list = c("type"))

# Remove duplicated rows. For example, genes and intergenic regions in each iteration (confidence level) are the same, so they are repeated.
TAB_LONG_SUMMARY_FINAL = TAB_LONG_SUMMARY_FINAL[!duplicated(TAB_LONG_SUMMARY_FINAL),]
TAB_WIDE_SUMMARY_FINAL = TAB_WIDE_SUMMARY_FINAL[!duplicated(TAB_WIDE_SUMMARY_FINAL),]

# Sort tables.
TAB_LONG_SUMMARY_FINAL = TAB_LONG_SUMMARY_FINAL[order(TAB_LONG_SUMMARY_FINAL$Type.1),]
TAB_WIDE_SUMMARY_FINAL = TAB_WIDE_SUMMARY_FINAL[order(TAB_WIDE_SUMMARY_FINAL$Type.1),]

# Save.
write.table(TAB_LONG_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/B/NR-SUMMARY-TAB_LONG_FINAL-B.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/B/NR-SUMMARY-TAB_WIDE_FINAL-B.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_STATISTICS_FINAL, paste0(WD, "/ALL/DEFINITIVE/B/NR-TAB_STATISTICS_FINAL-B.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("TAB_LONG_SUMMARY_FINAL", "TAB_WIDE_SUMMARY_FINAL", "TAB_STATISTICS_FINAL", "my_mean"))










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## 3.C FIGURES, SUMMARY TABLES AND STATISTICS "C" (CLASS_CODES ALL)

cat(paste0("\n\n\nFIGURES, SUMMARY TABLES AND STATISTICS \"C\" (CLASS_CODES ALL)..."))

if (!dir.exists(paste0(WD, "/ALL/DEFINITIVE/C"))){
  dir.create(paste0(WD, "/ALL/DEFINITIVE/C"))
}

## 3.C.1 LOAD TABLES

cat(paste0("\n\n\t-LOAD LONG AND WIDE TABLES..."))

# Load tables.
TAB_LONG_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_LONG-FINAL.tsv"), sep = "\t", header = T, quote = "\"")
TAB_WIDE_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_WIDE-FINAL.tsv"), sep = "\t", header = T, quote = "\"")

# Create factors
TAB_LONG_FINAL$Type.1 = factor(TAB_LONG_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_LONG_FINAL$Type.2 = factor(TAB_LONG_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_LONG_FINAL$Class_code = factor(TAB_LONG_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG_FINAL$Specie = factor(TAB_LONG_FINAL$Specie, levels = species_long_name)
TAB_LONG_FINAL$Feature = factor(TAB_LONG_FINAL$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

TAB_WIDE_FINAL$Type.1 = factor(TAB_WIDE_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_WIDE_FINAL$Type.2 = factor(TAB_WIDE_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_WIDE_FINAL$Class_code = factor(TAB_WIDE_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_WIDE_FINAL$Specie = factor(TAB_WIDE_FINAL$Specie, levels = species_long_name)

## 3.C.2 FIGURES

cat(paste0("\n\n\t-FIGURES..."))

my_mean = function(x) {
  log10(mean(10^x))
}

#### GC CONTENT

cat(paste0("\n\t\t-GC CONTENT..."))

# Grid
GC_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Exon_number = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Length = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Expression = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Class_code, y = Value, fill = Class_code)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Repeat_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")
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

# Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
# el siguiente error:
# 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
# 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
# Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
# -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
# a logaritmo serian -Inf.
# 5: Transformation introduced infinite values in continuous y-axis.

cat(paste0("\n\t\t-ALL FEATURES..."))

gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 5, common.legend = T, legend = "bottom")

ggsave(paste0(WD, "/ALL/DEFINITIVE/C/NR-C.png"), height = 6, width = 25, dpi = 600, bg = "white")
ggsave(paste0(WD, "/ALL/DEFINITIVE/C/NR-C.pdf"), height = 6, width = 25, dpi = 600, bg = "white")

rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))

## 3.C.3 MEAN AND MEDIAN TABLES

cat(paste0("\n\n\t-MEAN AND MEDIAN TABLE..."))

TAB_LONG_SUMMARY_FINAL = suppressMessages(TAB_LONG_FINAL %>%
                                            drop_na() %>%
                                            group_by(Class_code, Feature) %>%
                                            summarise(MEAN = mean(Value),
                                                      MEDIAN = median(Value)))
TAB_WIDE_SUMMARY_FINAL = pivot_wider(TAB_LONG_SUMMARY_FINAL,
                                     id_cols = c("Class_code"),
                                     names_from = "Feature",
                                     values_from = c("MEAN", "MEDIAN"),
                                     names_glue = "{Feature}.{.value}")
colnames(TAB_WIDE_SUMMARY_FINAL) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY_FINAL))

TAB_LONG_SUMMARY_FINAL = as.data.frame(TAB_LONG_SUMMARY_FINAL)
TAB_WIDE_SUMMARY_FINAL = as.data.frame(TAB_WIDE_SUMMARY_FINAL)

write.table(TAB_LONG_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/C/NR-SUMMARY-TAB_LONG_FINAL-C.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/C/NR-SUMMARY-TAB_WIDE_FINAL-C.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

## 3.C.4 STATISTICAL ANALYSIS

# Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
# necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como
# el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
# indipendientes. Algunos de estos papers son:
# - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
# - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)

cat(paste0("\n\n\t-STATISTICAL ANALYSIS..."))

combinations = as.data.frame(t(combn(c("pc", "u", "x", "i", "o/e", "ir"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("cl1", "cl2")

TAB_STATISTICS_FINAL = data.frame()
for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
  for (i in 1:nrow(combinations)) {
    cl1 = combinations[i, "cl1"]
    cl2 = combinations[i, "cl2"]
    subset1 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Class_code == cl1 & !is.na(TAB_LONG_FINAL$Value), "Value"]
    subset2 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Class_code == cl2 & !is.na(TAB_LONG_FINAL$Value), "Value"]
    L1 = length(subset1)
    L2 = length(subset2)
    if (L1 > 0 & L2 > 0) {
      test = wilcox.test(subset1, subset2, paired = FALSE)
      Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
      row = data.frame(Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
    }
    else{
      row = data.frame(Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
    }
    TAB_STATISTICS_FINAL = rbind(TAB_STATISTICS_FINAL, row)
  }
}

rownames(TAB_STATISTICS_FINAL) = NULL

write.table(TAB_STATISTICS_FINAL, paste0(WD, "/ALL/DEFINITIVE/C/NR-TAB_STATISTICS_FINAL-C.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("feature", "i", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

rm(list = c("TAB_LONG_FINAL", "TAB_WIDE_FINAL", "TAB_LONG_SUMMARY_FINAL", "TAB_WIDE_SUMMARY_FINAL", "TAB_STATISTICS_FINAL", "my_mean"))










###################################################################################################################################################################################################
###################################################################################################################################################################################################










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## 3.D FIGURES, SUMMARY TABLES AND STATISTICS "D" (SPECIES AND CONFIDENCE_LEVEL)

cat(paste0("\n\n\nFIGURES, SUMMARY TABLES AND STATISTICS \"D\" (SPECIES AND CONFIDENCE_LEVEL)..."))

if (!dir.exists(paste0(WD, "/ALL/DEFINITIVE/D"))){
  dir.create(paste0(WD, "/ALL/DEFINITIVE/D"))
}

## 3.D.1 LOAD TABLES

cat(paste0("\n\n\t-LOAD LONG AND WIDE TABLES..."))

# Load tables.
TAB_LONG_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_LONG-FINAL.tsv"), sep = "\t", header = T, quote = "\"")
TAB_WIDE_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_WIDE-FINAL.tsv"), sep = "\t", header = T, quote = "\"")

# Create factors
TAB_LONG_FINAL$Type.1 = factor(TAB_LONG_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_LONG_FINAL$Type.2 = factor(TAB_LONG_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_LONG_FINAL$Class_code = factor(TAB_LONG_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG_FINAL$Specie = factor(TAB_LONG_FINAL$Specie, levels = species_long_name)
TAB_LONG_FINAL$Feature = factor(TAB_LONG_FINAL$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

TAB_WIDE_FINAL$Type.1 = factor(TAB_WIDE_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_WIDE_FINAL$Type.2 = factor(TAB_WIDE_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_WIDE_FINAL$Class_code = factor(TAB_WIDE_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_WIDE_FINAL$Specie = factor(TAB_WIDE_FINAL$Specie, levels = species_long_name)

## 3.D.2 FIGURES

cat(paste0("\n\n\t-FIGURES..."))

my_mean = function(x) {
  log10(mean(10^x))
}

#### GC CONTENT

cat(paste0("\n\t\t-GC CONTENT..."))

# Grid
GC_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Type.1, y = Value, fill = Type.1)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")) +
  facet_grid(Feature~Specie) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold.italic"),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 2, widths = c(0.0046, 0.9954))

rm(list = c("GC_content"))

#### EXON NUMBER

cat(paste0("\n\t\t-EXON NUMBER..."))

# Grid
Exon_number = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Type.1, y = Value, fill = Type.1)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Exon_number"))

#### LENGTH

cat(paste0("\n\t\t-LENGTH..."))

# Grid
Length = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Type.1, y = Value, fill = Type.1)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Length"))

#### EXPRESSION

cat(paste0("\n\t\t-EXPRESSION..."))

# Grid
Expression = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

gg4 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg4, ncol = 2, widths = c(0.002, 0.998))

rm(list = c("Expression"))

#### REPEAT CONTENT

cat(paste0("\n\t\t-REPEAT CONTENT..."))

# Grid
Repeat_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Type.1, y = Value, fill = Type.1)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")
  ) +
  facet_grid(Feature~Specie) +
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
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

gg5 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg5, ncol = 2, widths = c(0.005, 0.995))

rm(list = c("Repeat_content"))

#### FIGURE FEATURES

# Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
# el siguiente error:
# 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
# 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
# Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
# -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
# a logaritmo serian -Inf.
# 5: Transformation introduced infinite values in continuous y-axis.

cat(paste0("\n\t\t-ALL FEATURES..."))

gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 1, heights = c(0.21, 0.19, 0.19, 0.19, 0.22))

ggsave(paste0(WD, "/ALL/DEFINITIVE/D/NR-D.png"), height = 25, width = 25, dpi = 600)
ggsave(paste0(WD, "/ALL/DEFINITIVE/D/NR-D.pdf"), height = 25, width = 25, dpi = 600)

rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))

## 3.D.3 MEAN AND MEDIAN TABLES

cat(paste0("\n\n\t-MEAN AND MEDIAN TABLE..."))

TAB_LONG_SUMMARY_FINAL = suppressMessages(TAB_LONG_FINAL %>% 
                                      drop_na() %>%
                                      group_by(Specie, Type.1, Feature) %>% 
                                      summarise(MEAN = mean(Value),
                                                MEDIAN = median(Value)))
TAB_WIDE_SUMMARY_FINAL = pivot_wider(TAB_LONG_SUMMARY_FINAL, 
                               id_cols = c("Specie", "Type.1"),
                               names_from = "Feature",
                               values_from = c("MEAN", "MEDIAN"),
                               names_glue = "{Feature}.{.value}")
colnames(TAB_WIDE_SUMMARY_FINAL) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY_FINAL))

TAB_LONG_SUMMARY_FINAL = as.data.frame(TAB_LONG_SUMMARY_FINAL)
TAB_WIDE_SUMMARY_FINAL = as.data.frame(TAB_WIDE_SUMMARY_FINAL)

write.table(TAB_LONG_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/D/NR-SUMMARY-TAB_LONG_FINAL-D.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/D/NR-SUMMARY-TAB_WIDE_FINAL-D.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

## 3.D.4 STATISTICAL ANALYSIS

# Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
# necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como 
# el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
# indipendientes. Algunos de estos papers son:
# - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
# - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)

cat(paste0("\n\n\t-STATISTICAL ANALYSIS..."))

combinations = as.data.frame(t(combn(c("Genes", types, "Intergenic Regions"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("ty1", "ty2")

TAB_STATISTICS_FINAL = data.frame()
for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
  for (spe in species_long_name) {
    for (i in 1:nrow(combinations)) {
      ty1 = combinations[i, "ty1"]
      ty2 = combinations[i, "ty2"]
      subset1 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Specie == spe & TAB_LONG_FINAL$Type.1 == ty1 & !is.na(TAB_LONG_FINAL$Value), "Value"]
      subset2 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Specie == spe & TAB_LONG_FINAL$Type.1 == ty2 & !is.na(TAB_LONG_FINAL$Value), "Value"]
      L1 = length(subset1)
      L2 = length(subset2)
      if (L1 > 0 & L2 > 0) {
        test = wilcox.test(subset1, subset2, paired = FALSE)
        Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
        row = data.frame(Specie = spe, Feature = feature, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
      } 
      else{
        row = data.frame(Specie = spe, Feature = feature, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
      }
      TAB_STATISTICS_FINAL = rbind(TAB_STATISTICS_FINAL, row)
    }
  }
}

rownames(TAB_STATISTICS_FINAL) = NULL

write.table(TAB_STATISTICS_FINAL, paste0(WD, "/ALL/DEFINITIVE/D/NR-TAB_STATISTICS_FINAL-D.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("feature", "spe", "i", "ty1", "ty2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

rm(list = c("TAB_LONG_FINAL", "TAB_WIDE_FINAL", "TAB_LONG_SUMMARY_FINAL", "TAB_WIDE_SUMMARY_FINAL", "TAB_STATISTICS_FINAL", "my_mean"))












################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## 3.E FIGURES, SUMMARY TABLES AND STATISTICS "E" (CONFIDENCE_LEVEL)

cat(paste0("\n\n\nFIGURES, SUMMARY TABLES AND STATISTICS \"E\" (CONFIDENCE_LEVEL)..."))

if (!dir.exists(paste0(WD, "/ALL/DEFINITIVE/E"))){
  dir.create(paste0(WD, "/ALL/DEFINITIVE/E"))
}

## 3.E.1 LOAD TABLES

cat(paste0("\n\n\t-LOAD LONG AND WIDE TABLES..."))

# Load tables.
TAB_LONG_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_LONG-FINAL.tsv"), sep = "\t", header = T, quote = "\"")
TAB_WIDE_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_WIDE-FINAL.tsv"), sep = "\t", header = T, quote = "\"")

# Create factors
TAB_LONG_FINAL$Type.1 = factor(TAB_LONG_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_LONG_FINAL$Type.2 = factor(TAB_LONG_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_LONG_FINAL$Class_code = factor(TAB_LONG_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG_FINAL$Specie = factor(TAB_LONG_FINAL$Specie, levels = species_long_name)
TAB_LONG_FINAL$Feature = factor(TAB_LONG_FINAL$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

TAB_WIDE_FINAL$Type.1 = factor(TAB_WIDE_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_WIDE_FINAL$Type.2 = factor(TAB_WIDE_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_WIDE_FINAL$Class_code = factor(TAB_WIDE_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_WIDE_FINAL$Specie = factor(TAB_WIDE_FINAL$Specie, levels = species_long_name)

## 3.E.2 FIGURES

cat(paste0("\n\n\t-FIGURES..."))

my_mean = function(x) {
  log10(mean(10^x))
}

#### GC CONTENT

cat(paste0("\n\t\t-GC CONTENT..."))

# Grid
GC_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Type.1, y = Value, fill = Type.1)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")) +
  facet_grid(.~Feature) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Exon_number = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Type.1, y = Value, fill = Type.1)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Length = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Type.1, y = Value, fill = Type.1)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Expression = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Repeat_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Type.1, y = Value, fill = Type.1)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee", "#5d65b4")
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

# Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
# el siguiente error:
# 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
# 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
# Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
# -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
# a logaritmo serian -Inf.
# 5: Transformation introduced infinite values in continuous y-axis.

cat(paste0("\n\t\t-ALL FEATURES..."))

gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 5, common.legend = T, legend = "bottom")

ggsave(paste0(WD, "/ALL/DEFINITIVE/E/NR-E.png"), height = 6, width = 25, dpi = 600)
ggsave(paste0(WD, "/ALL/DEFINITIVE/E/NR-E.pdf"), height = 6, width = 25, dpi = 600)

rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))

## 3.E.3 MEAN AND MEDIAN TABLES

cat(paste0("\n\n\t-MEAN AND MEDIAN TABLE..."))

TAB_LONG_SUMMARY_FINAL = suppressMessages(TAB_LONG_FINAL %>% 
                                            drop_na() %>%
                                            group_by(Type.1, Feature) %>% 
                                            summarise(MEAN = mean(Value),
                                                      MEDIAN = median(Value)))
TAB_WIDE_SUMMARY_FINAL = pivot_wider(TAB_LONG_SUMMARY_FINAL, 
                                     id_cols = c("Type.1"),
                                     names_from = "Feature",
                                     values_from = c("MEAN", "MEDIAN"),
                                     names_glue = "{Feature}.{.value}")
colnames(TAB_WIDE_SUMMARY_FINAL) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY_FINAL))

TAB_LONG_SUMMARY_FINAL = as.data.frame(TAB_LONG_SUMMARY_FINAL)
TAB_WIDE_SUMMARY_FINAL = as.data.frame(TAB_WIDE_SUMMARY_FINAL)

write.table(TAB_LONG_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/E/NR-SUMMARY-TAB_LONG_FINAL-E.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/E/NR-SUMMARY-TAB_WIDE_FINAL-E.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

## 3.E.4 STATISTICAL ANALYSIS

# Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
# necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como 
# el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
# indipendientes. Algunos de estos papers son:
# - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
# - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)

cat(paste0("\n\n\t-STATISTICAL ANALYSIS..."))

combinations = as.data.frame(t(combn(c("Genes", types, "Intergenic Regions"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("ty1", "ty2")

TAB_STATISTICS_FINAL = data.frame()
for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
  for (i in 1:nrow(combinations)) {
    ty1 = combinations[i, "ty1"]
    ty2 = combinations[i, "ty2"]
    subset1 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Type.1 == ty1 & !is.na(TAB_LONG_FINAL$Value), "Value"]
    subset2 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Type.1 == ty2 & !is.na(TAB_LONG_FINAL$Value), "Value"]
    L1 = length(subset1)
    L2 = length(subset2)
    if (L1 > 0 & L2 > 0) {
      test = wilcox.test(subset1, subset2, paired = FALSE)
      Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
      row = data.frame(Feature = feature, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
    } 
    else{
      row = data.frame(Feature = feature, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
    }
    TAB_STATISTICS_FINAL = rbind(TAB_STATISTICS_FINAL, row)
  }
}

rownames(TAB_STATISTICS_FINAL) = NULL

write.table(TAB_STATISTICS_FINAL, paste0(WD, "/ALL/DEFINITIVE/E/NR-TAB_STATISTICS_FINAL-E.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("feature", "i", "ty1", "ty2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

rm(list = c("TAB_LONG_FINAL", "TAB_WIDE_FINAL", "TAB_LONG_SUMMARY_FINAL", "TAB_WIDE_SUMMARY_FINAL", "TAB_STATISTICS_FINAL", "my_mean"))










###################################################################################################################################################################################################
###################################################################################################################################################################################################










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## 3.F FIGURES, SUMMARY TABLES AND STATISTICS "F" (SPECIES AND -GENES, LNCRNAS AND INTERGINC REGIONS-)

cat(paste0("\n\n\nFIGURES, SUMMARY TABLES AND STATISTICS \"F\" (SPECIES AND -GENES, LNCRNAS AND INTERGINC REGIONS-)..."))

if (!dir.exists(paste0(WD, "/ALL/DEFINITIVE/F"))){
  dir.create(paste0(WD, "/ALL/DEFINITIVE/F"))
}

## 3.F.1 LOAD TABLES

cat(paste0("\n\n\t-LOAD LONG AND WIDE TABLES..."))

# Load tables.
TAB_LONG_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_LONG-FINAL.tsv"), sep = "\t", header = T, quote = "\"")
TAB_WIDE_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_WIDE-FINAL.tsv"), sep = "\t", header = T, quote = "\"")

# Create factors
TAB_LONG_FINAL$Type.1 = factor(TAB_LONG_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_LONG_FINAL$Type.2 = factor(TAB_LONG_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_LONG_FINAL$Class_code = factor(TAB_LONG_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG_FINAL$Specie = factor(TAB_LONG_FINAL$Specie, levels = species_long_name)
TAB_LONG_FINAL$Feature = factor(TAB_LONG_FINAL$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

TAB_WIDE_FINAL$Type.1 = factor(TAB_WIDE_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_WIDE_FINAL$Type.2 = factor(TAB_WIDE_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_WIDE_FINAL$Class_code = factor(TAB_WIDE_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_WIDE_FINAL$Specie = factor(TAB_WIDE_FINAL$Specie, levels = species_long_name)

## 3.F.2 FIGURES

cat(paste0("\n\n\t-FIGURES..."))

my_mean = function(x) {
  log10(mean(10^x))
}

#### GC CONTENT

cat(paste0("\n\t\t-GC CONTENT..."))

# Grid
GC_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Type.2, y = Value, fill = Type.2)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#0089b2", "#5d65b4")) +
  facet_grid(Feature~Specie) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold.italic"),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 2, widths = c(0.0046, 0.9954))

rm(list = c("GC_content"))

#### EXON NUMBER

cat(paste0("\n\t\t-EXON NUMBER..."))

# Grid
Exon_number = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Type.2, y = Value, fill = Type.2)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#0089b2", "#5d65b4")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Exon_number"))

#### LENGTH

cat(paste0("\n\t\t-LENGTH..."))

# Grid
Length = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Type.2, y = Value, fill = Type.2)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#0089b2", "#5d65b4")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Length"))

#### EXPRESSION

cat(paste0("\n\t\t-EXPRESSION..."))

# Grid
Expression = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Type.2, y = Value, fill = Type.2)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#7cc1cf", "#0089b2", "#5d65b4")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

gg4 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg4, ncol = 2, widths = c(0.002, 0.998))

rm(list = c("Expression"))

#### REPEAT CONTENT

cat(paste0("\n\t\t-REPEAT CONTENT..."))

# Grid
Repeat_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Type.2, y = Value, fill = Type.2)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#0089b2", "#5d65b4")
  ) +
  facet_grid(Feature~Specie) +
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
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

gg5 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg5, ncol = 2, widths = c(0.005, 0.995))

rm(list = c("Repeat_content"))

#### FIGURE FEATURES

# Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
# el siguiente error:
# 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
# 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
# Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
# -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
# a logaritmo serian -Inf.
# 5: Transformation introduced infinite values in continuous y-axis.

cat(paste0("\n\t\t-ALL FEATURES..."))

gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 1, heights = c(0.21, 0.19, 0.19, 0.19, 0.22))

ggsave(paste0(WD, "/ALL/DEFINITIVE/F/NR-F.png"), height = 25, width = 25, dpi = 600)
ggsave(paste0(WD, "/ALL/DEFINITIVE/F/NR-F.pdf"), height = 25, width = 25, dpi = 600)

rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))

## 3.F.3 MEAN AND MEDIAN TABLES

cat(paste0("\n\n\t-MEAN AND MEDIAN TABLE..."))

TAB_LONG_SUMMARY_FINAL = suppressMessages(TAB_LONG_FINAL %>% 
                                            drop_na() %>%
                                            group_by(Specie, Type.2, Feature) %>% 
                                            summarise(MEAN = mean(Value),
                                                      MEDIAN = median(Value)))
TAB_WIDE_SUMMARY_FINAL = pivot_wider(TAB_LONG_SUMMARY_FINAL, 
                                     id_cols = c("Specie", "Type.2"),
                                     names_from = "Feature",
                                     values_from = c("MEAN", "MEDIAN"),
                                     names_glue = "{Feature}.{.value}")
colnames(TAB_WIDE_SUMMARY_FINAL) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY_FINAL))

TAB_LONG_SUMMARY_FINAL = as.data.frame(TAB_LONG_SUMMARY_FINAL)
TAB_WIDE_SUMMARY_FINAL = as.data.frame(TAB_WIDE_SUMMARY_FINAL)

write.table(TAB_LONG_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/F/NR-SUMMARY-TAB_LONG_FINAL-F.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/F/NR-SUMMARY-TAB_WIDE_FINAL-F.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

## 3.F.4 STATISTICAL ANALYSIS

# Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
# necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como 
# el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
# indipendientes. Algunos de estos papers son:
# - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
# - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)

cat(paste0("\n\n\t-STATISTICAL ANALYSIS..."))

combinations = as.data.frame(t(combn(c("Genes", "LncRNAs", "Intergenic Regions"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("ty1", "ty2")

TAB_STATISTICS_FINAL = data.frame()
for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
  for (spe in species_long_name) {
    for (i in 1:nrow(combinations)) {
      ty1 = combinations[i, "ty1"]
      ty2 = combinations[i, "ty2"]
      subset1 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Specie == spe & TAB_LONG_FINAL$Type.2 == ty1 & !is.na(TAB_LONG_FINAL$Value), "Value"]
      subset2 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Specie == spe & TAB_LONG_FINAL$Type.2 == ty2 & !is.na(TAB_LONG_FINAL$Value), "Value"]
      L1 = length(subset1)
      L2 = length(subset2)
      if (L1 > 0 & L2 > 0) {
        test = wilcox.test(subset1, subset2, paired = FALSE)
        Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
        row = data.frame(Specie = spe, Feature = feature, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
      } 
      else{
        row = data.frame(Specie = spe, Feature = feature, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
      }
      TAB_STATISTICS_FINAL = rbind(TAB_STATISTICS_FINAL, row)
    }
  }
}

rownames(TAB_STATISTICS_FINAL) = NULL

write.table(TAB_STATISTICS_FINAL, paste0(WD, "/ALL/DEFINITIVE/F/NR-TAB_STATISTICS_FINAL-F.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("feature", "spe", "i", "ty1", "ty2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

rm(list = c("TAB_LONG_FINAL", "TAB_WIDE_FINAL", "TAB_LONG_SUMMARY_FINAL", "TAB_WIDE_SUMMARY_FINAL", "TAB_STATISTICS_FINAL", "my_mean"))










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## 3.G FIGURES, SUMMARY TABLES AND STATISTICS "G" (-GENES, LNCRNAS AND INTERGINC REGIONS-)

cat(paste0("\n\n\nFIGURES, SUMMARY TABLES AND STATISTICS \"G\" (-GENES, LNCRNAS AND INTERGINC REGIONS-)..."))

if (!dir.exists(paste0(WD, "/ALL/DEFINITIVE/G"))){
  dir.create(paste0(WD, "/ALL/DEFINITIVE/G"))
}

## 3.G.1 LOAD TABLES

cat(paste0("\n\n\t-LOAD LONG AND WIDE TABLES..."))

# Load tables.
TAB_LONG_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_LONG-FINAL.tsv"), sep = "\t", header = T, quote = "\"")
TAB_WIDE_FINAL = read.table(paste0(WD, "/ALL/DEFINITIVE/NR-TAB_WIDE-FINAL.tsv"), sep = "\t", header = T, quote = "\"")

# Create factors
TAB_LONG_FINAL$Type.1 = factor(TAB_LONG_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_LONG_FINAL$Type.2 = factor(TAB_LONG_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_LONG_FINAL$Class_code = factor(TAB_LONG_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG_FINAL$Specie = factor(TAB_LONG_FINAL$Specie, levels = species_long_name)
TAB_LONG_FINAL$Feature = factor(TAB_LONG_FINAL$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

TAB_WIDE_FINAL$Type.1 = factor(TAB_WIDE_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_WIDE_FINAL$Type.2 = factor(TAB_WIDE_FINAL$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
TAB_WIDE_FINAL$Class_code = factor(TAB_WIDE_FINAL$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_WIDE_FINAL$Specie = factor(TAB_WIDE_FINAL$Specie, levels = species_long_name)

## 3.G.2 FIGURES

cat(paste0("\n\n\t-FIGURES..."))

my_mean = function(x) {
  log10(mean(10^x))
}

#### GC CONTENT

cat(paste0("\n\t\t-GC CONTENT..."))

# Grid
GC_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Type.2, y = Value, fill = Type.2)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#0089b2", "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Exon_number = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Type.2, y = Value, fill = Type.2)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#0089b2", "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Length = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Type.2, y = Value, fill = Type.2)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#0089b2", "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Expression = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Type.2, y = Value, fill = Type.2)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#0089b2", "#5d65b4")
  ) +
  facet_grid(.~Feature) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  guides(fill = guide_legend(nrow = 1)) +
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

# Grid
Repeat_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Type.2, y = Value, fill = Type.2)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LncRNAs", "Intergenic regions"),
    values = c("#7cc1cf", "#0089b2", "#5d65b4")
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

# Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
# el siguiente error:
# 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
# 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
# Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
# -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
# a logaritmo serian -Inf.
# 5: Transformation introduced infinite values in continuous y-axis.

cat(paste0("\n\t\t-ALL FEATURES..."))

gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 5, common.legend = T, legend = "bottom")

ggsave(paste0(WD, "/ALL/DEFINITIVE/G/NR-G.png"), height = 6, width = 25, dpi = 600)
ggsave(paste0(WD, "/ALL/DEFINITIVE/G/NR-G.pdf"), height = 6, width = 25, dpi = 600)

rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))

## 3.G.3 MEAN AND MEDIAN TABLES

cat(paste0("\n\n\t-MEAN AND MEDIAN TABLE..."))

TAB_LONG_SUMMARY_FINAL = suppressMessages(TAB_LONG_FINAL %>% 
                                            drop_na() %>%
                                            group_by(Type.2, Feature) %>% 
                                            summarise(MEAN = mean(Value),
                                                      MEDIAN = median(Value)))
TAB_WIDE_SUMMARY_FINAL = pivot_wider(TAB_LONG_SUMMARY_FINAL, 
                                     id_cols = c("Type.2"),
                                     names_from = "Feature",
                                     values_from = c("MEAN", "MEDIAN"),
                                     names_glue = "{Feature}.{.value}")
colnames(TAB_WIDE_SUMMARY_FINAL) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY_FINAL))

TAB_LONG_SUMMARY_FINAL = as.data.frame(TAB_LONG_SUMMARY_FINAL)
TAB_WIDE_SUMMARY_FINAL = as.data.frame(TAB_WIDE_SUMMARY_FINAL)

write.table(TAB_LONG_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/G/NR-SUMMARY-TAB_LONG_FINAL-G.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_FINAL, paste0(WD, "/ALL/DEFINITIVE/G/NR-SUMMARY-TAB_WIDE_FINAL-G.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

## 3.G.4 STATISTICAL ANALYSIS

# Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
# necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como 
# el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
# indipendientes. Algunos de estos papers son:
# - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
# - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)

cat(paste0("\n\n\t-STATISTICAL ANALYSIS..."))

combinations = as.data.frame(t(combn(c("Genes", "LncRNAs", "Intergenic Regions"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("ty1", "ty2")

TAB_STATISTICS_FINAL = data.frame()
for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
  for (i in 1:nrow(combinations)) {
    ty1 = combinations[i, "ty1"]
    ty2 = combinations[i, "ty2"]
    subset1 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Type.2 == ty1 & !is.na(TAB_LONG_FINAL$Value), "Value"]
    subset2 = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == feature & TAB_LONG_FINAL$Type.2 == ty2 & !is.na(TAB_LONG_FINAL$Value), "Value"]
    L1 = length(subset1)
    L2 = length(subset2)
    if (L1 > 0 & L2 > 0) {
      test = wilcox.test(subset1, subset2, paired = FALSE)
      Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
      row = data.frame(Feature = feature, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
    } 
    else{
      row = data.frame(Feature = feature, TY1 = ty1, TY2 = ty2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
    }
    TAB_STATISTICS_FINAL = rbind(TAB_STATISTICS_FINAL, row)
  }
}

rownames(TAB_STATISTICS_FINAL) = NULL

write.table(TAB_STATISTICS_FINAL, paste0(WD, "/ALL/DEFINITIVE/G/NR-TAB_STATISTICS_FINAL-G.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("feature", "i", "ty1", "ty2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

rm(list = c("TAB_LONG_FINAL", "TAB_WIDE_FINAL", "TAB_LONG_SUMMARY_FINAL", "TAB_WIDE_SUMMARY_FINAL", "TAB_STATISTICS_FINAL", "my_mean"))

