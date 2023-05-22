################################################################################
#
# INDIVIDUAL: TISSUE SPECIFICITY STUDY: APPROACH 1 - STEP 5
#
# Analyse transcripts with TAU > 0.8.
# 
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggExtra))
suppressMessages(library(data.table))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(circlize))
suppressMessages(options(bitmapType='cairo'))

## 1. PATHS

# Own computer
path_tissue_specificity = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
path_comp = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/07-comparison_lncRNAs_vs_coding_genes"
flag = "nr"
species_short_name = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
species_long_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")

# # Garnatxa
# path_tissue_specificity = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
# path_comp = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/07-comparison_lncRNAs_vs_coding_genes"
# flag = "nr"
# species_short_name = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
# species_long_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")


## 2. MEAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEAN TAU (TABLE): \n"))

TAB_mean = data.frame()
for (i in 1:length(species_short_name)) {
  
  spe = species_short_name[i]
  spe_l = species_long_name[i]
  
  cat(paste0("------ Specie: ", spe, "\n"))
  
  # Directories.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5"))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5"))
  }
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/", spe))
  }
  
  # Get files with TAU metric values and log-transformed expression values.
  files = list.files(path = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe), pattern = ".tsv")
  
  if (length(files) > 0) {
    SRA.Studies = unique(unlist(lapply(strsplit(unique(unlist(lapply(strsplit(files, "_"), `[[`, 1))), "-"), `[[`, 2)))
    
    for (SRA.Study in SRA.Studies) {
      cat(paste0("--------- SRA.Study: ", SRA.Study, "\n"))
      
      # Load table with the tissue specificity results and the log-tranformed expression values. 
      tab_mean_TAU_G = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/GENES-", SRA.Study, "_mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU_L = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/LNCRNAS-", SRA.Study, "_mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU = rbind(tab_mean_TAU_G, tab_mean_TAU_L)
      tab_mean_TAU = tab_mean_TAU[, c("ID_transcript", "Confidence", "Class_code", "TAU")]
      
      # Add column to specify TAU > 0.8
      tab_mean_TAU$"Type" = ifelse(tab_mean_TAU$TAU >= 0.8, "TS", "Non-TS")
      
      # Load LncRNA and gene features info.
      if (flag == "nr") {
        tab_info = read.table(paste0(path_comp, "/ALL/TAB_WIDE_FINAL_ALL_SPECIES-NR.tsv"), sep = "\t", header = T, quote = "\"")
        tab_info = tab_info[tab_info$Specie == spe_l & tab_info$Type.1 != "Intergenic Regions", c("ID_transcript", "Type.1", "Class_code", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")]
        colnames(tab_info) = c("ID_transcript", "Confidence", "Class_code", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")
      } 
      if (flag == "r") {
        tab_info = read.table(paste0(path_comp, "/ALL/TAB_WIDE_FINAL_ALL_SPECIES-R.tsv"), sep = "\t", header = T, quote = "\"")
        tab_info = tab_info[tab_info$Specie == spe_l & tab_info$Type.1 != "Intergenic Regions", c("ID_transcript", "Type.1", "Class_code", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")]
        colnames(tab_info) = c("ID_transcript", "Confidence", "Class_code", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")
      }
      tab_info[tab_info == 'Genes'] = 'Protein-coding'
      tab_info[tab_info == 'LncRNAs low'] = 'Low-confidence lncRNA'
      tab_info[tab_info == 'LncRNAs medium'] = 'Medium-confidence lncRNA'
      tab_info[tab_info == 'LncRNAs high'] = 'High-confidence lncRNA'
      
      # Join TAU info with features info. NA values refer to lncRNAs and genes that had no more than 1 TPM in any of the tissues and those that were duplicated during quantification.
      tab_mean_merged = merge(tab_mean_TAU, tab_info, by = c("ID_transcript", "Confidence", "Class_code"), all = T)
      tab_mean_merged$Type = ifelse(is.na(tab_mean_merged$TAU), "Non-analysed", tab_mean_merged$Type)
      
      # Convert to factors.
      tab_mean_merged$Confidence = factor(tab_mean_merged$Confidence, levels = c("Protein-coding", "Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
      tab_mean_merged$Class_code = factor(tab_mean_merged$Class_code, levels = c("pc", "u", "x", "i", "o/e"))
      tab_mean_merged$Type = factor(tab_mean_merged$Type, levels = c("Non-analysed", "Non-TS", "TS"))
      
      # Save table.
      write.table(tab_mean_merged, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/", spe, "/", SRA.Study, "_mean-FEATURES.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      
      # Add info.
      tab_mean_merged$"Spe" = spe
      tab_mean_merged$"SRA.Study" = SRA.Study
      tab_mean_merged = tab_mean_merged[,c("Spe", "SRA.Study", "ID_transcript", "Confidence", "Class_code", "Type", "TAU", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")]
      TAB_mean = rbind(TAB_mean, tab_mean_merged)
      
      rm(list = c("tab_mean_TAU", "tab_info", "tab_mean_merged"))
    }
  }
}

write.table(TAB_mean, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/mean-FEATURES.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

rm(list = c("files", "spe", "spe_l", "SRA.Study", "SRA.Studies", "i"))

# Convert wide to long tab.
cols = names(which(sapply(TAB_mean, is.integer)))
TAB_mean[,cols]=lapply(cols,function(x) {tryCatch({as.numeric(TAB_mean[[x]])},warning = function(w) {TAB_mean[[x]]})} )
TAB_mean_long = melt(setDT(TAB_mean), id.vars = c("Spe", "SRA.Study", "ID_transcript", "Confidence", "Class_code", "TAU", "Type"), variable.name = "Feature", value.name = "Value")
TAB_mean_long = TAB_mean_long[TAB_mean_long$Feature != "Exons" & TAB_mean_long$Feature != "Length", ]
TAB_mean_long$Feature = ifelse(TAB_mean_long$Feature == "GC", "GC content", 
                               ifelse(TAB_mean_long$Feature == "Log.Exons.1", "Exon number",
                                      ifelse(TAB_mean_long$Feature == "Log.Length", "Length",
                                             ifelse(TAB_mean_long$Feature == "log2.TPMs.1.mean", "Expression",
                                                    ifelse(TAB_mean_long$Feature == "overlap_per", "Repeat content", TAB_mean_long$Feature)))))
TAB_mean_long$Feature = factor(TAB_mean_long$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
TAB_mean_long$Confidence = ifelse(TAB_mean_long$Confidence == "Protein-coding", "pc", 
                                  ifelse(TAB_mean_long$Confidence == "Low-confidence lncRNA", "LC-lncRNA", 
                                         ifelse(TAB_mean_long$Confidence == "Medium-confidence lncRNA", "MC-lncRNA", 
                                                ifelse(TAB_mean_long$Confidence == "High-confidence lncRNA", "HC-lncRNA", TAB_mean_long$Confidence))))
TAB_mean_long$Confidence = factor(TAB_mean_long$Confidence, levels = c("pc", "LC-lncRNA", "MC-lncRNA", "HC-lncRNA"))
TAB_mean_long$Spe = factor(TAB_mean_long$Spe, levels = species_short_name)

# Figure (all class code).
cat(paste0("\n\n--- MEAN TAU (FIGURE ALL CLASS CODES): \n"))

species_label = species_long_name
names(species_label) = species_short_name

label_y = paste0(
  "         ", "Repeat content (%)", "         ",
  "            ", "Log2(TPMs + 1)", "           ",
  "             ", "Log2(Length)", "           ",
  "     ", "Log2(Number of exons + 1)", "      ",
  "         ", "GC content (%)", "                "
)

gg = ggplot(TAB_mean_long, aes(x = Confidence, y = Value, fill = Confidence, alpha = Type)) + 
  geom_boxplot(aes(fill = Confidence, alpha = Type), colour = "black", outlier.size = 0.2, position = position_dodge2(width = 0.9, preserve = "single")) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2.5, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
  scale_fill_manual(
    labels = c("Protein-coding (pc)", "Low-Confidence lncRNA (LC-lncRNA)", "Medium-Confidence lncRNA (MC-lncRNA)", "High-Confidence lncRNA (HC-lncRNA)"), 
    values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b")
  ) +
  scale_alpha_manual(
    labels = c("Non-Analysed", "Non-Tissue Specific", "Tissue Specific"),
    values = c(0.2,0.6,1)
  ) +
  facet_grid(Feature~Spe, scales = "free_y", labeller = labeller(Spe = species_label)) +
  xlab("") +
  ylab(label_y) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 23, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 2),
         alpha = guide_legend(nrow = 2)) +
  theme(axis.text.x = element_text(hjust =1, size = 19, angle = 65),
        axis.text.y = element_text(size = 19),
        axis.title.y.left = element_text(size = 22, hjust = 0),
        strip.text.x = element_text(size = 19, face = "bold.italic"),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1))

ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/mean-FEATURES.png"), height = 25, width = 25, dpi = 600)
ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/mean-FEATURES.pdf"), height = 25, width = 25, dpi = 600)

rm(list = c("gg", "label_y"))

# Figure (by class code).
cat(paste0("\n\n--- MEAN TAU (FIGURE BY CLASS CODES): \n"))

species_label = species_long_name
names(species_label) = species_short_name

class_codes = as.character(unique(TAB_mean_long$Class_code))[2:length(as.character(unique(TAB_mean_long$Class_code)))]
for (class in class_codes) {
  cat(paste0("------ Class code: ", class, "\n"))
  
  label_y = paste0(
    "         ", "Repeat content (%)", "         ",
    "            ", "Log2(TPMs + 1)", "           ",
    "             ", "Log2(Length)", "           ",
    "     ", "Log2(Number of exons + 1)", "      ",
    "         ", "GC content (%)", "                "
  )
  
  subset = TAB_mean_long[TAB_mean_long$Class_code == "pc" | TAB_mean_long$Class_code == class, ]
  
  gg = ggplot(subset, aes(x = Confidence, y = Value, fill = Confidence, alpha = Type)) + 
    geom_boxplot(aes(fill = Confidence, alpha = Type), colour = "black", outlier.size = 0.2, position = position_dodge2(width = 0.9, preserve = "single")) + 
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2.5, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
    scale_fill_manual(
      labels = c("Protein-coding (pc)", "Low-Confidence lncRNA (LC-lncRNA)", "Medium-Confidence lncRNA (MC-lncRNA)", "High-Confidence lncRNA (HC-lncRNA)"), 
      values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b")
    ) +
    scale_alpha_manual(
      labels = c("Non-Analysed", "Non-Tissue Specific", "Tissue Specific"),
      values = c(0.2,0.6,1)
    ) +
    facet_grid(Feature~Spe, scales = "free_y", labeller = labeller(Spe = species_label)) +
    xlab("") +
    ylab(label_y) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 23, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 2),
           alpha = guide_legend(nrow = 2)) +
    theme(axis.text.x = element_text(hjust =1, size = 19, angle = 65),
          axis.text.y = element_text(size = 19),
          axis.title.y.left = element_text(size = 22, hjust = 0),
          strip.text.x = element_text(size = 19, face = "bold.italic"),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/mean-FEATURES-", gsub("/", "-", class), ".png"), height = 25, width = 25, dpi = 600)
  ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/mean-FEATURES-", gsub("/", "-", class), ".pdf"), height = 25, width = 25, dpi = 600)
  
  rm(list = c("gg", "label_y", "subset"))
}

rm(list = c("species_label"))


## 3. MEDIAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEDIAN TAU: \n"))

TAB_median = data.frame()
for (i in 1:length(species_short_name)) {
  
  spe = species_short_name[i]
  spe_l = species_long_name[i]
  
  cat(paste0("------ Specie: ", spe, "\n"))
  
  # Directories.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5"))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5"))
  }
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/", spe))
  }
  
  # Get files with TAU metric values and log-transformed expression values.
  files = list.files(path = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe), pattern = ".tsv")
  
  if (length(files) > 0) {
    SRA.Studies = unique(unlist(lapply(strsplit(unique(unlist(lapply(strsplit(files, "_"), `[[`, 1))), "-"), `[[`, 2)))
    
    for (SRA.Study in SRA.Studies) {
      cat(paste0("--------- SRA.Study: ", SRA.Study, "\n"))
      
      # Load table with the tissue specificity results and the log-tranformed expression values. 
      tab_median_TAU_G = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/GENES-", SRA.Study, "_median-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_median_TAU_L = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/LNCRNAS-", SRA.Study, "_median-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_median_TAU = rbind(tab_median_TAU_G, tab_median_TAU_L)
      tab_median_TAU = tab_median_TAU[, c("ID_transcript", "Confidence", "Class_code", "TAU")]
      
      # Add column to specify TAU > 0.8
      tab_median_TAU$"Type" = ifelse(tab_median_TAU$TAU >= 0.8, "TS", "Non-TS")
      
      # Load LncRNA and gene features info.
      if (flag == "nr") {
        tab_info = read.table(paste0(path_comp, "/ALL/TAB_WIDE_FINAL_ALL_SPECIES-NR.tsv"), sep = "\t", header = T, quote = "\"")
        tab_info = tab_info[tab_info$Specie == spe_l & tab_info$Type.1 != "Intergenic Regions", c("ID_transcript", "Type.1", "Class_code", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")]
        colnames(tab_info) = c("ID_transcript", "Confidence", "Class_code", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")
      } 
      if (flag == "r") {
        tab_info = read.table(paste0(path_comp, "/ALL/TAB_WIDE_FINAL_ALL_SPECIES-R.tsv"), sep = "\t", header = T, quote = "\"")
        tab_info = tab_info[tab_info$Specie == spe_l & tab_info$Type.1 != "Intergenic Regions", c("ID_transcript", "Type.1", "Class_code", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")]
        colnames(tab_info) = c("ID_transcript", "Confidence", "Class_code", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")
      }
      tab_info[tab_info == 'Genes'] = 'Protein-coding'
      tab_info[tab_info == 'LncRNAs low'] = 'Low-confidence lncRNA'
      tab_info[tab_info == 'LncRNAs medium'] = 'Medium-confidence lncRNA'
      tab_info[tab_info == 'LncRNAs high'] = 'High-confidence lncRNA'
      
      # Join TAU info with features info. NA values refer to lncRNAs and genes that had no more than 1 TPM in any of the tissues and those that were duplicated during quantification.
      tab_median_merged = merge(tab_median_TAU, tab_info, by = c("ID_transcript", "Confidence", "Class_code"), all = T)
      tab_median_merged$Type = ifelse(is.na(tab_median_merged$TAU), "Non-analysed", tab_median_merged$Type)
      
      # Convert to factors.
      tab_median_merged$Confidence = factor(tab_median_merged$Confidence, levels = c("Protein-coding", "Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
      tab_median_merged$Class_code = factor(tab_median_merged$Class_code, levels = c("pc", "u", "x", "i", "o/e"))
      tab_median_merged$Type = factor(tab_median_merged$Type, levels = c("Non-analysed", "Non-TS", "TS"))
      
      # Save table.
      write.table(tab_median_merged, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/", spe, "/", SRA.Study, "_median-FEATURES.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      
      # Add info.
      tab_median_merged$"Spe" = spe
      tab_median_merged$"SRA.Study" = SRA.Study
      tab_median_merged = tab_median_merged[,c("Spe", "SRA.Study", "ID_transcript", "Confidence", "Class_code", "Type", "TAU", "Exons", "Log.Exons.1", "Length", "Log.Length", "GC", "log2.TPMs.1.mean", "overlap_per")]
      TAB_median = rbind(TAB_median, tab_median_merged)
      
      rm(list = c("tab_median_TAU", "tab_info", "tab_median_merged"))
      
    }
  }
}

write.table(TAB_median, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/median-FEATURES.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

rm(list = c("files", "spe", "spe_l", "SRA.Study", "SRA.Studies", "i"))

# Convert wide to long tab.
cols = names(which(sapply(TAB_median, is.integer)))
TAB_median[,cols]=lapply(cols,function(x) {tryCatch({as.numeric(TAB_median[[x]])},warning = function(w) {TAB_median[[x]]})} )
TAB_median_long = melt(setDT(TAB_median), id.vars = c("Spe", "SRA.Study", "ID_transcript", "Confidence", "Class_code", "TAU", "Type"), variable.name = "Feature", value.name = "Value")
TAB_median_long = TAB_median_long[TAB_median_long$Feature != "Exons" & TAB_median_long$Feature != "Length", ]
TAB_median_long$Feature = ifelse(TAB_median_long$Feature == "GC", "GC content", 
                               ifelse(TAB_median_long$Feature == "Log.Exons.1", "Exon number",
                                      ifelse(TAB_median_long$Feature == "Log.Length", "Length",
                                             ifelse(TAB_median_long$Feature == "log2.TPMs.1.mean", "Expression",
                                                    ifelse(TAB_median_long$Feature == "overlap_per", "Repeat content", TAB_median_long$Feature)))))
TAB_median_long$Feature = factor(TAB_median_long$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
TAB_median_long$Confidence = ifelse(TAB_median_long$Confidence == "Protein-coding", "pc", 
                                  ifelse(TAB_median_long$Confidence == "Low-confidence lncRNA", "LC-lncRNA", 
                                         ifelse(TAB_median_long$Confidence == "Medium-confidence lncRNA", "MC-lncRNA", 
                                                ifelse(TAB_median_long$Confidence == "High-confidence lncRNA", "HC-lncRNA", TAB_median_long$Confidence))))
TAB_median_long$Confidence = factor(TAB_median_long$Confidence, levels = c("pc", "LC-lncRNA", "MC-lncRNA", "HC-lncRNA"))
TAB_median_long$Spe = factor(TAB_median_long$Spe, levels = species_short_name)

# Figure (all class code).
cat(paste0("\n\n--- MEDIAN TAU (FIGURE ALL CLASS CODES): \n"))

species_label = species_long_name
names(species_label) = species_short_name

label_y = paste0(
  "         ", "Repeat content (%)", "         ",
  "            ", "Log2(TPMs + 1)", "           ",
  "             ", "Log2(Length)", "           ",
  "     ", "Log2(Number of exons + 1)", "      ",
  "         ", "GC content (%)", "                "
)

gg = ggplot(TAB_median_long, aes(x = Confidence, y = Value, fill = Confidence, alpha = Type)) + 
  geom_boxplot(aes(fill = Confidence, alpha = Type), colour = "black", outlier.size = 0.2, position = position_dodge2(width = 0.9, preserve = "single")) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2.5, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
  scale_fill_manual(
    labels = c("Protein-coding (pc)", "Low-Confidence lncRNA (LC-lncRNA)", "Medium-Confidence lncRNA (MC-lncRNA)", "High-Confidence lncRNA (HC-lncRNA)"), 
    values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b")
  ) +
  scale_alpha_manual(
    labels = c("Non-Analysed", "Non-Tissue Specific", "Tissue Specific"),
    values = c(0.2,0.6,1)
  ) +
  facet_grid(Feature~Spe, scales = "free_y", labeller = labeller(Spe = species_label)) +
  xlab("") +
  ylab(label_y) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 23, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 2),
         alpha = guide_legend(nrow = 2)) +
  theme(axis.text.x = element_text(hjust =1, size = 19, angle = 65),
        axis.text.y = element_text(size = 19),
        axis.title.y.left = element_text(size = 22, hjust = 0),
        strip.text.x = element_text(size = 19, face = "bold.italic"),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1))

ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/median-FEATURES.png"), height = 25, width = 25, dpi = 600)
ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/median-FEATURES.pdf"), height = 25, width = 25, dpi = 600)

rm(list = c("gg", "label_y"))

# Figure (by class code).
cat(paste0("\n\n--- MEDIAN TAU (FIGURE BY CLASS CODES): \n"))

species_label = species_long_name
names(species_label) = species_short_name

class_codes = as.character(unique(TAB_median_long$Class_code))[2:length(as.character(unique(TAB_median_long$Class_code)))]
for (class in class_codes) {
  cat(paste0("------ Class code: ", class, "\n"))
  
  label_y = paste0(
    "         ", "Repeat content (%)", "         ",
    "            ", "Log2(TPMs + 1)", "           ",
    "             ", "Log2(Length)", "           ",
    "     ", "Log2(Number of exons + 1)", "      ",
    "         ", "GC content (%)", "                "
  )
  
  subset = TAB_median_long[TAB_median_long$Class_code == "pc" | TAB_median_long$Class_code == class, ]
  
  gg = ggplot(subset, aes(x = Confidence, y = Value, fill = Confidence, alpha = Type)) + 
    geom_boxplot(aes(fill = Confidence, alpha = Type), colour = "black", outlier.size = 0.2, position = position_dodge2(width = 0.9, preserve = "single")) + 
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2.5, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
    scale_fill_manual(
      labels = c("Protein-coding (pc)", "Low-Confidence lncRNA (LC-lncRNA)", "Medium-Confidence lncRNA (MC-lncRNA)", "High-Confidence lncRNA (HC-lncRNA)"), 
      values = c("#a2ded9", "#decfa9", "#dcb8b8", "#b5cf9b")
    ) +
    scale_alpha_manual(
      labels = c("Non-Analysed", "Non-Tissue Specific", "Tissue Specific"),
      values = c(0.2,0.6,1)
    ) +
    facet_grid(Feature~Spe, scales = "free_y", labeller = labeller(Spe = species_label)) +
    xlab("") +
    ylab(label_y) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 23, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 2),
           alpha = guide_legend(nrow = 2)) +
    theme(axis.text.x = element_text(hjust =1, size = 19, angle = 65),
          axis.text.y = element_text(size = 19),
          axis.title.y.left = element_text(size = 22, hjust = 0),
          strip.text.x = element_text(size = 19, face = "bold.italic"),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/median-FEATURES-", gsub("/", "-", class), ".png"), height = 25, width = 25, dpi = 600)
  ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP5/median-FEATURES-", gsub("/", "-", class), ".pdf"), height = 25, width = 25, dpi = 600)
  
  rm(list = c("gg", "label_y", "subset"))
}

rm(list = c("species_label"))

