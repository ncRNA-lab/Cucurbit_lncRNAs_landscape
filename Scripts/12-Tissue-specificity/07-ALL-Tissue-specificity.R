################################################################################
#
# ALL: TISSUE SPECIFICITY STUDY: APPROACH 1 - STEP 7
#
# Motif-level conservation, Positional-level conservation and Tissue-specificity
# analysis.
# 
################################################################################


#https://www.life-science-alliance.org/content/6/6/e202302002
#https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07313-3


rm(list = ls())

se = function(x) sqrt(var(x)/length(x)) 


## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggExtra))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(PupillometryR))
suppressMessages(options(bitmapType='cairo'))


## 1. PATHS

# Own computer
path_tissue_specificity = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
path_comp = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/07-comparison_lncRNAs_vs_coding_genes"
path_cons = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics"
flag = "nr"
species_short_name = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
species_long_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")
comparisons = c("ALL", "intergenic", "antisense", "intronic", "sense")
confidences = c("Low", "Medium", "High")

# # Garnatxa
# path_tissue_specificity = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
# path_comp = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/07-comparison_lncRNAs_vs_coding_genes"
# path_cons = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-comparative_genomics"
# flag = "nr"
# species_short_name = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
# species_long_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")



################################################################################
## 2. LOAD CONSERVATION INFO.

cat(paste0("\n\n--- LOAD CONSERVATION INFO... \n"))

# Load lncRNAs conservation info at sequence level (BLastn).
cat(paste0("------ Sequence Level (Blastn)\n"))
SeqlevB = read.table(paste0(path_cons, "/Sequence_level/Blastn/", flag, "/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), header = T, sep = "\t", quote = "\"")
SeqlevB = SeqlevB[, c("Member", "Type", "Conserved_level", "Number_species_by_family", "Class")]
colnames(SeqlevB) = c("ID_transcript", "SeqLev.Blastn.Conservation_1", "SeqLev.Blastn.Conservation_2", "SeqLev.Blastn.Conservation_3", "Comparison")
SeqB_extra_info = data.frame()
for (con in confidences) {
  for (com in comparisons) {
    tab = tryCatch({read.table(paste0(path_cons, "/Sequence_level/Blastn/", flag, "/Definitive/05-Families/", con, "/", com, "/gen.tsv"), header = T, sep = "\t", quote = "\"")}, 
                   error = function(e) {NULL}, 
                   warning = function(e) {message(paste("gen.tsv file doesn't exist for clustering method OrthoFinder, confidence-level (", co, ") and class (", com, "). So, we will not use this info."))}
    )
    if (!is.null(tab)) {
      tab = tab[, c("Member", "Family")]
      colnames(tab) = c("ID_transcript", "SeqLev.Blastn.Family")
      tab$"Comparison" = com
      SeqB_extra_info = rbind(SeqB_extra_info, tab)
    }
  }
}
SeqlevB = merge(SeqlevB, SeqB_extra_info, by = c("ID_transcript", "Comparison"), all = T)

# Load lncRNAs conservation info at sequence level (OrthoFinder).
cat(paste0("------ Sequence Level (OrthoFinder)\n"))
SeqlevOF = read.table(paste0(path_cons, "/Sequence_level/OrthoFinder/", flag, "/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), header = T, sep = "\t", quote = "\"")
SeqlevOF = SeqlevOF[, c("Member", "Type", "Conserved_level", "Number_species_by_family", "Class")]
colnames(SeqlevOF) = c("ID_transcript", "SeqLev.OF.Conservation_1", "SeqLev.OF.Conservation_2", "SeqLev.OF.Conservation_3", "Comparison")
SeqOF_extra_info = data.frame()
for (con in confidences) {
  for (com in comparisons) {
    tab = tryCatch({read.table(paste0(path_cons, "/Sequence_level/OrthoFinder/", flag, "/Definitive/05-Families/", con, "/", com, "/gen.tsv"), header = T, sep = "\t", quote = "\"")}, 
                   error = function(e) {NULL}, 
                   warning = function(e) {message(paste("gen.tsv file doesn't exist for clustering method OrthoFinder, confidence-level (", co, ") and class (", com, "). So, we will not use this info."))}
    )
    if (!is.null(tab)) {
      tab = tab[, c("Member", "Family")]
      colnames(tab) = c("ID_transcript", "SeqLev.OF.Family")
      tab$"Comparison" = com
      SeqOF_extra_info = rbind(SeqOF_extra_info, tab)
    }
  }
}
SeqlevOF = merge(SeqlevOF, SeqOF_extra_info, by = c("ID_transcript", "Comparison"), all = T)

# Load lncRNAs conservation info at sequence level (OrthoMCL).
cat(paste0("------ Sequence Level (OrthoMCL)\n"))
SeqlevOM = read.table(paste0(path_cons, "/Sequence_level/OrthoMCL/", flag, "/Definitive/06-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), header = T, sep = "\t", quote = "\"")
SeqlevOM = SeqlevOM[, c("Member", "Type", "Conserved_level", "Number_species_by_family", "Class")]
colnames(SeqlevOM) = c("ID_transcript", "SeqLev.OM.Conservation_1", "SeqLev.OM.Conservation_2", "SeqLev.OM.Conservation_3", "Comparison")
SeqOM_extra_info = data.frame()
for (con in confidences) {
  for (com in comparisons) {
    tab = tryCatch({read.table(paste0(path_cons, "/Sequence_level/OrthoMCL/", flag, "/Definitive/05-Families/", con, "/", com, "/gen.tsv"), header = T, sep = "\t", quote = "\"")}, 
             error = function(e) {NULL}, 
             warning = function(e) {message(paste("gen.tsv file doesn't exist for clustering method OrthoMCL, confidence-level (", con, ") and class (", com, "). So, we will not use this info."))}
    )
    if (!is.null(tab)) {
      tab = tab[, c("Member", "Family")]
      colnames(tab) = c("ID_transcript", "SeqLev.OM.Family")
      tab$"Comparison" = com
      SeqOM_extra_info = rbind(SeqOM_extra_info, tab)
    }
  }
}
SeqlevOM = merge(SeqlevOM, SeqOM_extra_info, by = c("ID_transcript", "Comparison"), all = T)

# Load lncRNAs conservation info at positional level.
cat(paste0("------ Positional Level\n"))
Poslev = read.table(paste0(path_cons, "/Positional_level/Approach_2/", flag, "/05-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), header = T, sep = "\t", quote = "\"")
Poslev = Poslev[Poslev$Strictness == "ORIGINAL" & Poslev$NonMatch == "no", c("Member", "Type", "Conserved_level", "Number_species_by_family", "Class")]
colnames(Poslev) = c("ID_transcript", "PosLev.Conservation_1", "PosLev.Conservation_2", "PosLev.Conservation_3", "Comparison")
Pos_extra_info = data.frame()
for (con in confidences) {
  for (com in comparisons) {
    tab = read.table(paste0(path_cons, "/Positional_level/Approach_2/", flag, "/04-Families/", con, "/", com, "/gen_ORIGINAL_no.tsv"), header = T, sep = "\t", quote = "\"")
    tab = tab[, c("Member", "Family")]
    colnames(tab) = c("ID_transcript", "PosLev.Family")
    tab$"Comparison" = com
    Pos_extra_info = rbind(Pos_extra_info, tab)
  }
}
Poslev = merge(Poslev, Pos_extra_info, by = c("ID_transcript", "Comparison"), all = T)

# Load lncRNAs conservation info at motif level.
cat(paste0("------ Motif Level\n"))
Motlev = read.table(paste0(path_cons, "/Motif_level/", flag, "/Positional_conserved/06-Figures_and_tables/Tables/ORIGINAL/no/GLOBAL_TABLE_MEME-GOMO-POSLEV.tsv"), header = T, sep = "\t", quote = "\"")
Motlev = Motlev[Motlev$Type == "REAL", c("LncRNA", "Comparison", "Mode", "Width", "Meme_Motif.Identifier", "Meme_Motif.E.value", "Meme_Motif.LncRNA", "Meme_Motif.P.value", "Meme_Motif.Width", "Gomo_GO.Term.Identifier", "Gomo_Q.value")]
colnames(Motlev) = c("ID_transcript", "Comparison", "MotLev.Mode", "MotLev.Width", "MotLev.Meme_Motif.Identifier", "MotLev.Meme_Motif.E.value", "MotLev.Meme_Motif.LncRNA", "MotLev.Meme_Motif.P.value", "MotLev.Meme_Motif.Width", "MotLev.Gomo_GO.Term.Identifier", "MotLev.Gomo_Q.value")

rm(list = c("con", "com", "SeqB_extra_info", "SeqOF_extra_info", "SeqOM_extra_info", "Pos_extra_info", "tab"))



################################################################################
## 3. MEAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEAN TAU (TABLE): \n"))

TAB_mean = data.frame()

for (i in 1:length(species_short_name)) {
  
  spe = species_short_name[i]
  spe_l = species_long_name[i]
  
  cat(paste0("------ Specie: ", spe, "\n"))
  
  # Directories.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7"))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7"))
  }
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe))
  }
  
  # Get files with TAU metric values and log-transformed expression values.
  files = list.files(path = paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP2/", spe), pattern = ".tsv")
  
  if (length(files) > 0) {
    SRA.Studies = unique(unlist(lapply(strsplit(files, "_"), `[[`, 1)))
    
    for (SRA.Study in SRA.Studies) {
      cat(paste0("--------- SRA.Study: ", SRA.Study, "\n"))
      
      # Load table with the tissue specificity results and the log-tranformed expression values. 
      tab_mean_TAU = read.table(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP2/", spe, "/", SRA.Study, "_mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU = tab_mean_TAU[tab_mean_TAU$Class_code != "pc", c("ID_transcript", "Confidence", "Class_code", "TAU")]
      tab_mean_TAU$ID_transcript = paste0(tab_mean_TAU$ID_transcript, "-", spe)
      
      # Add column to specify TAU > 0.8
      tab_mean_TAU$"Type" = ifelse(tab_mean_TAU$TAU >= 0.8, "TS", "Non-TS")
      
      # Load LncRNA features info.
      if (flag == "nr") {
        tab_info = read.table(paste0(path_comp, "/ALL/TAB_WIDE_FINAL_ALL_SPECIES-NR.tsv"), sep = "\t", header = T, quote = "\"")
        tab_info = tab_info[tab_info$Specie == spe_l & tab_info$Type.1 != "Intergenic Regions" & tab_info$Type.1 != "Genes", c("ID_transcript", "Type.1", "Class_code", "Exons", "Length", "GC", "TPMs.mean", "overlap_per")]
        tab_info$ID_transcript = paste0(tab_info$ID_transcript, "-", spe)
        colnames(tab_info) = c("ID_transcript", "Confidence", "Class_code", "Exons", "Length", "GC", "TPMs.mean", "overlap_per")
      } 
      if (flag == "r") {
        tab_info = read.table(paste0(path_comp, "/ALL/TAB_WIDE_FINAL_ALL_SPECIES-R.tsv"), sep = "\t", header = T, quote = "\"")
        tab_info = tab_info[tab_info$Specie == spe_l & tab_info$Type.1 != "Intergenic Regions" & tab_info$Type.1 != "Genes", c("ID_transcript", "Type.1", "Class_code", "Exons", "Length", "GC", "TPMs.mean", "overlap_per")]
        tab_info$ID_transcript = paste0(tab_info$ID_transcript, "-", spe)
        colnames(tab_info) = c("ID_transcript", "Confidence", "Class_code", "Exons", "Length", "GC", "TPMs.mean", "overlap_per")
      }
      
      tab_info[tab_info == 'LncRNAs low'] = 'Low-confidence lncRNA'
      tab_info[tab_info == 'LncRNAs medium'] = 'Medium-confidence lncRNA'
      tab_info[tab_info == 'LncRNAs high'] = 'High-confidence lncRNA'
      
      # Join TAU info with features info. NA values refer to lncRNAs and genes that had no more than 1 TPM in any of the tissues and those that were duplicated during quantification.
      tab_mean_merged = merge(tab_mean_TAU, tab_info, by = c("ID_transcript", "Confidence", "Class_code"), all = T)
      tab_mean_merged$Type = ifelse(is.na(tab_mean_merged$TAU), "Non-analysed", tab_mean_merged$Type)
      
      # Add conservation info.
      tab_mean_merged = merge(tab_mean_merged, SeqlevB, by = "ID_transcript", all.x = T, all.y = F)
      tab_mean_merged = merge(tab_mean_merged, SeqlevOF, by = c("ID_transcript", "Comparison"), all.x = T, all.y = F)
      tab_mean_merged = merge(tab_mean_merged, SeqlevOM, by = c("ID_transcript", "Comparison"), all.x = T, all.y = F)
      tab_mean_merged = merge(tab_mean_merged, Poslev, by = c("ID_transcript", "Comparison"), all.x = T, all.y = F)
      tab_mean_merged = merge(tab_mean_merged, Motlev, by = c("ID_transcript", "Comparison"), all.x = T, all.y = F)
      
      # Convert to factors.
      tab_mean_merged$Confidence = factor(tab_mean_merged$Confidence, levels = c("Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
      tab_mean_merged$Class_code = factor(tab_mean_merged$Class_code, levels = c("u", "x", "i", "o/e"))
      tab_mean_merged$Type = factor(tab_mean_merged$Type, levels = c("Non-analysed", "Non-TS", "TS"))
      tab_mean_merged$SeqLev.Blastn.Conservation_1 = factor(tab_mean_merged$SeqLev.Blastn.Conservation_1, levels = c("Non-conserved", "Conserved"))
      tab_mean_merged$SeqLev.Blastn.Conservation_2 = factor(tab_mean_merged$SeqLev.Blastn.Conservation_2, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      tab_mean_merged$SeqLev.Blastn.Conservation_3 = factor(tab_mean_merged$SeqLev.Blastn.Conservation_3, levels = 1:9)
      tab_mean_merged$SeqLev.OF.Conservation_1 = factor(tab_mean_merged$SeqLev.OF.Conservation_1, levels = c("Non-conserved", "Conserved"))
      tab_mean_merged$SeqLev.OF.Conservation_2 = factor(tab_mean_merged$SeqLev.OF.Conservation_2, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      tab_mean_merged$SeqLev.OF.Conservation_3 = factor(tab_mean_merged$SeqLev.OF.Conservation_3, levels = 1:9)
      tab_mean_merged$SeqLev.OM.Conservation_1 = factor(tab_mean_merged$SeqLev.OM.Conservation_1, levels = c("Non-conserved", "Conserved"))
      tab_mean_merged$SeqLev.OM.Conservation_2 = factor(tab_mean_merged$SeqLev.OM.Conservation_2, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      tab_mean_merged$SeqLev.OM.Conservation_3 = factor(tab_mean_merged$SeqLev.OM.Conservation_3, levels = 1:9)
      tab_mean_merged$PosLev.Conservation_1 = factor(tab_mean_merged$PosLev.Conservation_1, levels = c("Non-conserved", "Conserved"))
      tab_mean_merged$PosLev.Conservation_2 = factor(tab_mean_merged$PosLev.Conservation_2, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      tab_mean_merged$PosLev.Conservation_3 = factor(tab_mean_merged$PosLev.Conservation_3, levels = 1:9)
      
      # Save table.
      write.table(tab_mean_merged, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/", SRA.Study, "_mean.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      
      # Add info.
      tab_mean_merged = cbind(SRA.Study = SRA.Study, tab_mean_merged)
      tab_mean_merged = cbind(Spe = spe, tab_mean_merged)
      TAB_mean = rbind(TAB_mean, tab_mean_merged)
      
      rm(list = c("tab_mean_TAU", "tab_info", "tab_mean_merged"))
    }
  }
}

write.table(TAB_mean, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/mean.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

rm(list = c("files", "spe", "spe_l", "SRA.Study", "SRA.Studies", "i"))











#################################
cat(paste0("\n\n--- MEAN TAU (CORRELATION - FUSION PROJECTS): \n"))

# Directories.
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7"))
}
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures"))
}

# Convert to factors.
TAB_mean$Confidence = factor(TAB_mean$Confidence, levels = c("Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
TAB_mean$Type = factor(TAB_mean$Type, levels = c("Non-analysed", "Non-TS", "TS"))


############################################
## COMPARISON == ALL
cat(paste0("------ COMPARISON == ALL\n"))

if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL"))
}

for (con in confidences) {
  cat(paste0("--------- ", con, '-confidence lncRNA\n'))
  
  ############################################
  # OrthoFinder A
  tab_temp_OF = TAB_mean[TAB_mean$Type != "Non-analysed" & TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison == "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "SeqLev.OF.Conservation_3", "TAU")]
  tab_temp_OF = tab_temp_OF[!duplicated(tab_temp_OF),]
  tab_temp_OF$Class_code = as.character(tab_temp_OF$Class_code)
  tab_temp_OF[tab_temp_OF == "u"] = "lincRNAs"
  tab_temp_OF[tab_temp_OF == "x"] = "NAT-lncRNAs"
  tab_temp_OF[tab_temp_OF == "i"] = "int-lncRNAs"
  tab_temp_OF[tab_temp_OF == "o/e"] = "SOT-lncRNAs"
  tab_temp_OF$SeqLev.OF.Conservation_3 = factor(tab_temp_OF$SeqLev.OF.Conservation_3, levels = 9:1)
  tab_temp_OF$Class_code = factor(tab_temp_OF$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  
  gg1 = ggplot(tab_temp_OF, aes(x = SeqLev.OF.Conservation_3, y = TAU)) +
    geom_flat_violin(position = position_nudge(x = .18, y = 0), adjust = 1, trim = TRUE, alpha = 1, fill = "#f2c36a") +
    geom_point(aes(x = as.numeric(SeqLev.OF.Conservation_3)-.22, y = TAU),position = position_jitter(width = .05), size = .25, shape = 20) +
    geom_boxplot(aes(x = SeqLev.OF.Conservation_3, y = TAU), outlier.shape = NA, alpha = 1, width = .25, colour = "black") +
    stat_summary(aes(x = SeqLev.OF.Conservation_3), fun=mean, geom="point", shape=18, size=2, color="red") +
    scale_y_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
    theme_bw() +
    xlab("Conservation level (Sequence OF)") + ylab("TAU") +
    theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    coord_flip() +
    theme(plot.margin = margin(5.5, 0, 5.5, 5.5))
  
  gg1_facet = gg1 + facet_grid(. ~ Class_code)
  
  gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 1, heights = c(0.09, 0.91))
  
  gg1_facet = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1_facet, ncol = 1, heights = c(0.058, 0.942))
  
  # OrthoFinder B
  tab_temp_OF_all = TAB_mean[TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison == "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "SeqLev.OF.Conservation_3", "Type")]
  tab_temp_OF_all = tab_temp_OF_all[!duplicated(tab_temp_OF_all),]
  tab_temp_OF_all$Class_code = as.character(tab_temp_OF_all$Class_code)
  tab_temp_OF_all[tab_temp_OF_all == "u"] = "lincRNAs"
  tab_temp_OF_all[tab_temp_OF_all == "x"] = "NAT-lncRNAs"
  tab_temp_OF_all[tab_temp_OF_all == "i"] = "int-lncRNAs"
  tab_temp_OF_all[tab_temp_OF_all == "o/e"] = "SOT-lncRNAs"
  tab_temp_OF_all$"Count" = 1
  
  tab_temp_OF_all_red = tab_temp_OF_all %>%
    group_by(SRA.Study, SeqLev.OF.Conservation_3, Class_code, Type) %>%
    summarise(Counts = sum(Count)) %>%
    group_by(SRA.Study, SeqLev.OF.Conservation_3, Class_code) %>%
    mutate(Percent = (100*Counts)/sum(Counts)) %>%
    group_by(SeqLev.OF.Conservation_3, Class_code, Type) %>%
    summarise(Mean.Percent = mean(Percent), SE.Percent = se(Percent))
  tab_temp_OF_all_red = as.data.frame(tab_temp_OF_all_red)
  tab_temp_OF_all_red$SeqLev.OF.Conservation_3 = factor(tab_temp_OF_all_red$SeqLev.OF.Conservation_3, levels = 9:1)
  tab_temp_OF_all_red$Class_code = factor(tab_temp_OF_all_red$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  
  gg2 = ggplot(tab_temp_OF_all_red, aes(x = SeqLev.OF.Conservation_3, y = Mean.Percent, ymin=Mean.Percent-SE.Percent, ymax=Mean.Percent+SE.Percent, group = Type, fill = Type)) +
    geom_bar(aes(group = Type, fill = Type), position=position_dodge(preserve="single"), stat="identity", color = "black", width = 0.7) +
    geom_errorbar(position=position_dodge(width=0.7, preserve="single"), colour="black", alpha = 1, width=0.4) +
    scale_fill_manual(values = c("#9ED4F5", "#6E93AA", "#2A4758")) +
    xlab("") +
    ylab("LncRNAs (%)") +
    facet_grid(. ~ Class_code) +
    theme_bw() +
    theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 10, margin = margin(r = 0.5, unit = 'cm')),
          legend.key.size = unit(0.8, "lines"),
          legend.text.align = 0) +
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    coord_flip() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
  
  gg2 = ggarrange(gg2, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), ncol = 1, heights = c(0.997, 0.003))
  
  # Joinn figures.
  gg_final_OF = ggarrange(gg1, gg2, ncol = 2, widths = c(0.60, 0.40))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-SEQ_OF-", con, ".png"), height = 8, width = 16, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-SEQ_OF-", con, ".pdf"), height = 8, width = 16, dpi = 600, bg = "white")
  
  gg_final_OF_facet = ggarrange(gg1_facet, gg2, ncol = 2, widths = c(0.60, 0.40))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-SEQ_OF-", con, "_facet.png"), height = 8, width = 20, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-SEQ_OF-", con, "_facet.pdf"), height = 8, width = 20, dpi = 600, bg = "white")
  
  rm(list = c("tab_temp_OF", "tab_temp_OF_all", "tab_temp_OF_all_red", "gg1", "gg1_facet", "gg2"))
  
  
  ############################################
  # Positional A
  tab_temp_Pos = TAB_mean[TAB_mean$Type != "Non-analysed" & TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison == "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "PosLev.Conservation_3", "TAU")]
  tab_temp_Pos = tab_temp_Pos[!duplicated(tab_temp_Pos),]
  tab_temp_Pos$Class_code = as.character(tab_temp_Pos$Class_code)
  tab_temp_Pos[tab_temp_Pos == "u"] = "lincRNAs"
  tab_temp_Pos[tab_temp_Pos == "x"] = "NAT-lncRNAs"
  tab_temp_Pos[tab_temp_Pos == "i"] = "int-lncRNAs"
  tab_temp_Pos[tab_temp_Pos == "o/e"] = "SOT-lncRNAs"
  tab_temp_Pos$PosLev.Conservation_3 = factor(tab_temp_Pos$PosLev.Conservation_3, levels = 9:1)
  tab_temp_Pos$Class_code = factor(tab_temp_Pos$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  
  # tab_temp_Pos = TAB_mean[TAB_mean$Type != "Non-analysed" & TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison == "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "PosLev.Conservation_3", "TAU")]
  # tab_temp_Pos = tab_temp_Pos[!duplicated(tab_temp_Pos),]
  # tab_temp_Pos$Class_code = as.character(tab_temp_Pos$Class_code)
  # 
  # tab_temp_Mot = TAB_mean[TAB_mean$Type != "Non-analysed" & TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison == "ALL" & TAB_mean$PosLev.Conservation_3 != 1 & TAB_mean$MotLev.Width == "6-50", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "MotLev.Meme_Motif.Identifier")]
  # tab_temp_Mot$"MotLev.Meme_Motif.Presence_Ausence" = ifelse(!is.na(tab_temp_Mot$MotLev.Meme_Motif.Identifier), "Presence", "Ausence")
  # tab_temp_Mot = tab_temp_Mot[, c("Spe", "SRA.Study", "ID_transcript", "Class_code", "MotLev.Meme_Motif.Presence_Ausence")]
  # tab_temp_Mot = tab_temp_Mot[!duplicated(tab_temp_Mot),]
  # 
  # tab_temp_Pos = merge(tab_temp_Pos, tab_temp_Mot, by = c("Spe", "SRA.Study", "ID_transcript", "Class_code"), all = T)
  # tab_temp_Pos$MotLev.Meme_Motif.Presence_Ausence = ifelse(is.na(tab_temp_Pos$MotLev.Meme_Motif.Presence_Ausence), "Ausence", tab_temp_Pos$MotLev.Meme_Motif.Presence_Ausence)
  # tab_temp_Pos[tab_temp_Pos == "u"] = "lincRNAs"
  # tab_temp_Pos[tab_temp_Pos == "x"] = "NAT-lncRNAs"
  # tab_temp_Pos[tab_temp_Pos == "i"] = "int-lncRNAs"
  # tab_temp_Pos[tab_temp_Pos == "o/e"] = "SOT-lncRNAs"
  # tab_temp_Pos$PosLev.Conservation_3 = factor(tab_temp_Pos$PosLev.Conservation_3, levels = 9:1)
  # tab_temp_Pos$Class_code = factor(tab_temp_Pos$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  # tab_temp_Pos$MotLev.Meme_Motif.Presence_Ausence = factor(tab_temp_Pos$MotLev.Meme_Motif.Presence_Ausence, levels = c("Presence", "Ausence"))
  
  gg1 = ggplot(tab_temp_Pos, aes(x = PosLev.Conservation_3, y = TAU)) +
    geom_flat_violin(position = position_nudge(x = .18, y = 0), adjust = 1, trim = TRUE, alpha = 1, fill = "#f2c36a") +
    geom_point(aes(x = as.numeric(PosLev.Conservation_3)-.22, y = TAU), position = position_jitter(width = .05), size = .25, shape = 20) +
    geom_boxplot(aes(x = PosLev.Conservation_3, y = TAU), outlier.shape = NA, alpha = 1, width = .25, colour = "black") +
    stat_summary(aes(x = PosLev.Conservation_3), fun=mean, geom="point", shape=18, size=2, color="red") +
    scale_y_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
    theme_bw() +
    xlab("Conservation level (Positonal)") + ylab("TAU") +
    theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    coord_flip() +
    theme(plot.margin = margin(5.5, 0, 5.5, 5.5))
  
  gg1_facet = gg1 + facet_grid(. ~ Class_code)
  
  gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 1, heights = c(0.09, 0.91))
  
  gg1_facet = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1_facet, ncol = 1, heights = c(0.058, 0.942))
  
  # Positional B
  tab_temp_Pos_all = TAB_mean[TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison == "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "PosLev.Conservation_3", "Type")]
  tab_temp_Pos_all = tab_temp_Pos_all[!duplicated(tab_temp_Pos_all),]
  tab_temp_Pos_all$Class_code = as.character(tab_temp_Pos_all$Class_code)
  tab_temp_Pos_all[tab_temp_Pos_all == "u"] = "lincRNAs"
  tab_temp_Pos_all[tab_temp_Pos_all == "x"] = "NAT-lncRNAs"
  tab_temp_Pos_all[tab_temp_Pos_all == "i"] = "int-lncRNAs"
  tab_temp_Pos_all[tab_temp_Pos_all == "o/e"] = "SOT-lncRNAs"
  tab_temp_Pos_all$"Count" = 1
  
  tab_temp_Pos_all_red = tab_temp_Pos_all %>%
    group_by(SRA.Study, PosLev.Conservation_3, Class_code, Type) %>%
    summarise(Counts = sum(Count)) %>%
    group_by(SRA.Study, PosLev.Conservation_3, Class_code) %>%
    mutate(Percent = (100*Counts)/sum(Counts)) %>%
    group_by(PosLev.Conservation_3, Class_code, Type) %>%
    summarise(Mean.Percent = mean(Percent), SE.Percent = se(Percent))
  tab_temp_Pos_all_red = as.data.frame(tab_temp_Pos_all_red)
  tab_temp_Pos_all_red$PosLev.Conservation_3 = factor(tab_temp_Pos_all_red$PosLev.Conservation_3, levels = 9:1)
  tab_temp_Pos_all_red$Class_code = factor(tab_temp_Pos_all_red$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  
  gg2 = ggplot(tab_temp_Pos_all_red, aes(x = PosLev.Conservation_3, y = Mean.Percent, ymin=Mean.Percent-SE.Percent, ymax=Mean.Percent+SE.Percent, group = Type, fill = Type)) +
    geom_bar(aes(group = Type, fill = Type), position=position_dodge(preserve="single"), stat="identity", color = "black", width = 0.7) +
    geom_errorbar(position=position_dodge(width=0.7, preserve="single"), colour="black", alpha = 1, width=0.4) +
    scale_fill_manual(values = c("#9ED4F5", "#6E93AA", "#2A4758")) +
    xlab("") +
    ylab("LncRNAs (%)") +
    facet_grid(. ~ Class_code) +
    theme_bw() +
    theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 10, margin = margin(r = 0.5, unit = 'cm')),
          legend.key.size = unit(0.8, "lines"),
          legend.text.align = 0) +
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    coord_flip() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
  
  gg2 = ggarrange(gg2, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), ncol = 1, heights = c(0.997, 0.003))
  
  # Joinn figures.
  gg_final_Pos = ggarrange(gg1, gg2, ncol = 2, widths = c(0.60, 0.40))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-POS-", con, ".png"), height = 8, width = 16, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-POS-", con, ".pdf"), height = 8, width = 16, dpi = 600, bg = "white")
  
  gg_final_Pos_facet = ggarrange(gg1_facet, gg2, ncol = 2, widths = c(0.60, 0.40))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-POS-", con, "_facet.png"), height = 8, width = 20, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-POS-", con, "_facet.pdf"), height = 8, width = 20, dpi = 600, bg = "white")
  
  rm(list = c("tab_temp_Pos", "tab_temp_Pos_all", "tab_temp_Pos_all_red", "gg1", "gg1_facet", "gg2"))
  
  
  ############################################
  # Join Seq and Pos figures
  gg_final = ggarrange(gg_final_OF, gg_final_Pos, nrow = 2)
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-", con, ".png"), height = 16, width = 16, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-", con, ".pdf"), height = 16, width = 16, dpi = 600, bg = "white")
  
  gg_final_facet = ggarrange(gg_final_OF_facet, gg_final_Pos_facet, nrow = 2)
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-", con, "_facet.png"), height = 16, width = 20, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/ALL/mean-", con, "_facet.pdf"), height = 16, width = 20, dpi = 600, bg = "white")
  
  rm(list = c("gg_final_OF", "gg_final_Pos", "gg_final_OF_facet", "gg_final_Pos_facet", "gg_final", "gg_final_facet"))
}

rm(list = c("con"))


############################################
## COMPARISON != ALL
cat(paste0("------ COMPARISON != ALL\n"))

if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES"))
}

for (con in confidences) {
  cat(paste0("--------- ", con, '-confidence lncRNA\n'))
  
  ############################################
  # OrthoFinder A
  tab_temp_OF = TAB_mean[TAB_mean$Type != "Non-analysed" & TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison != "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "SeqLev.OF.Conservation_3", "TAU")]
  tab_temp_OF = tab_temp_OF[!duplicated(tab_temp_OF),]
  tab_temp_OF$Class_code = as.character(tab_temp_OF$Class_code)
  tab_temp_OF[tab_temp_OF == "u"] = "lincRNAs"
  tab_temp_OF[tab_temp_OF == "x"] = "NAT-lncRNAs"
  tab_temp_OF[tab_temp_OF == "i"] = "int-lncRNAs"
  tab_temp_OF[tab_temp_OF == "o/e"] = "SOT-lncRNAs"
  tab_temp_OF$SeqLev.OF.Conservation_3 = factor(tab_temp_OF$SeqLev.OF.Conservation_3, levels = 9:1)
  tab_temp_OF$Class_code = factor(tab_temp_OF$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  
  gg1 = ggplot(tab_temp_OF, aes(x = SeqLev.OF.Conservation_3, y = TAU)) +
    geom_flat_violin(position = position_nudge(x = .18, y = 0), adjust = 1, trim = TRUE, alpha = 1, fill = "#f2c36a") +
    geom_point(aes(x = as.numeric(SeqLev.OF.Conservation_3)-.22, y = TAU),position = position_jitter(width = .05), size = .25, shape = 20) +
    geom_boxplot(aes(x = SeqLev.OF.Conservation_3, y = TAU), outlier.shape = NA, alpha = 1, width = .25, colour = "black") +
    stat_summary(aes(x = SeqLev.OF.Conservation_3), fun=mean, geom="point", shape=18, size=2, color="red") +
    scale_y_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
    theme_bw() +
    xlab("Conservation level (Sequence OF)") + ylab("TAU") +
    theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    coord_flip() +
    theme(plot.margin = margin(5.5, 0, 5.5, 5.5))
  
  gg1_facet = gg1 + facet_grid(. ~ Class_code)
  
  gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 1, heights = c(0.09, 0.91))
  
  gg1_facet = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1_facet, ncol = 1, heights = c(0.058, 0.942))
  
  # OrthoFinder B
  tab_temp_OF_all = TAB_mean[TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison != "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "SeqLev.OF.Conservation_3", "Type")]
  tab_temp_OF_all = tab_temp_OF_all[!duplicated(tab_temp_OF_all),]
  tab_temp_OF_all$Class_code = as.character(tab_temp_OF_all$Class_code)
  tab_temp_OF_all[tab_temp_OF_all == "u"] = "lincRNAs"
  tab_temp_OF_all[tab_temp_OF_all == "x"] = "NAT-lncRNAs"
  tab_temp_OF_all[tab_temp_OF_all == "i"] = "int-lncRNAs"
  tab_temp_OF_all[tab_temp_OF_all == "o/e"] = "SOT-lncRNAs"
  tab_temp_OF_all$"Count" = 1
  
  tab_temp_OF_all_red = tab_temp_OF_all %>%
    group_by(SRA.Study, SeqLev.OF.Conservation_3, Class_code, Type) %>%
    summarise(Counts = sum(Count)) %>%
    group_by(SRA.Study, SeqLev.OF.Conservation_3, Class_code) %>%
    mutate(Percent = (100*Counts)/sum(Counts)) %>%
    group_by(SeqLev.OF.Conservation_3, Class_code, Type) %>%
    summarise(Mean.Percent = mean(Percent), SE.Percent = se(Percent))
  tab_temp_OF_all_red = as.data.frame(tab_temp_OF_all_red)
  tab_temp_OF_all_red$SeqLev.OF.Conservation_3 = factor(tab_temp_OF_all_red$SeqLev.OF.Conservation_3, levels = 9:1)
  tab_temp_OF_all_red$Class_code = factor(tab_temp_OF_all_red$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  
  gg2 = ggplot(tab_temp_OF_all_red, aes(x = SeqLev.OF.Conservation_3, y = Mean.Percent, ymin=Mean.Percent-SE.Percent, ymax=Mean.Percent+SE.Percent, group = Type, fill = Type)) +
    geom_bar(aes(group = Type, fill = Type), position=position_dodge(preserve="single"), stat="identity", color = "black", width = 0.7) +
    geom_errorbar(position=position_dodge(width=0.7, preserve="single"), colour="black", alpha = 1, width=0.4) +
    scale_fill_manual(values = c("#9ED4F5", "#6E93AA", "#2A4758")) +
    xlab("") +
    ylab("LncRNAs (%)") +
    facet_grid(. ~ Class_code) +
    theme_bw() +
    theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 10, margin = margin(r = 0.5, unit = 'cm')),
          legend.key.size = unit(0.8, "lines"),
          legend.text.align = 0) +
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    coord_flip() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
  
  gg2 = ggarrange(gg2, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), ncol = 1, heights = c(0.997, 0.003))
  
  # Joinn figures.
  gg_final_OF = ggarrange(gg1, gg2, ncol = 2, widths = c(0.60, 0.40))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean-SEQ_OF-", con, ".png"), height = 8, width = 16, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean-SEQ_OF-", con, ".pdf"), height = 8, width = 16, dpi = 600, bg = "white")
  
  gg_final_OF_facet = ggarrange(gg1_facet, gg2, ncol = 2, widths = c(0.60, 0.40))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean-SEQ_OF-", con, "_facet.png"), height = 8, width = 20, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean-SEQ_OF-", con, "_facet.pdf"), height = 8, width = 20, dpi = 600, bg = "white")
  
  rm(list = c("tab_temp_OF", "tab_temp_OF_all", "tab_temp_OF_all_red", "gg1", "gg1_facet", "gg2"))
  
  
  ############################################
  # Positional A
  tab_temp_Pos = TAB_mean[TAB_mean$Type != "Non-analysed" & TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison != "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "PosLev.Conservation_3", "TAU")]
  tab_temp_Pos = tab_temp_Pos[!duplicated(tab_temp_Pos),]
  tab_temp_Pos$Class_code = as.character(tab_temp_Pos$Class_code)
  tab_temp_Pos[tab_temp_Pos == "u"] = "lincRNAs"
  tab_temp_Pos[tab_temp_Pos == "x"] = "NAT-lncRNAs"
  tab_temp_Pos[tab_temp_Pos == "i"] = "int-lncRNAs"
  tab_temp_Pos[tab_temp_Pos == "o/e"] = "SOT-lncRNAs"
  tab_temp_Pos$PosLev.Conservation_3 = factor(tab_temp_Pos$PosLev.Conservation_3, levels = 9:1)
  tab_temp_Pos$Class_code = factor(tab_temp_Pos$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  
  gg1 = ggplot(tab_temp_Pos, aes(x = PosLev.Conservation_3, y = TAU)) +
    geom_flat_violin(position = position_nudge(x = .18, y = 0), adjust = 1, trim = TRUE, alpha = 1, fill = "#f2c36a") +
    geom_point(aes(x = as.numeric(PosLev.Conservation_3)-.22, y = TAU),position = position_jitter(width = .05), size = .25, shape = 20) +
    geom_boxplot(aes(x = PosLev.Conservation_3, y = TAU), outlier.shape = NA, alpha = 1, width = .25, colour = "black") +
    stat_summary(aes(x = PosLev.Conservation_3), fun=mean, geom="point", shape=18, size=2, color="red") +
    scale_y_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
    theme_bw() +
    xlab("Conservation level (Positonal)") + ylab("TAU") +
    theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    coord_flip() +
    theme(plot.margin = margin(5.5, 0, 5.5, 5.5))
  
  gg1_facet = gg1 + facet_grid(. ~ Class_code)
  
  gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 1, heights = c(0.09, 0.91))
  
  gg1_facet = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1_facet, ncol = 1, heights = c(0.058, 0.942))
  
  # Positional B
  tab_temp_Pos_all = TAB_mean[TAB_mean$Confidence == paste0(con, '-confidence lncRNA') & TAB_mean$Comparison != "ALL", c("Spe", "SRA.Study", "ID_transcript", "Class_code", "PosLev.Conservation_3", "Type")]
  tab_temp_Pos_all = tab_temp_Pos_all[!duplicated(tab_temp_Pos_all),]
  tab_temp_Pos_all$Class_code = as.character(tab_temp_Pos_all$Class_code)
  tab_temp_Pos_all[tab_temp_Pos_all == "u"] = "lincRNAs"
  tab_temp_Pos_all[tab_temp_Pos_all == "x"] = "NAT-lncRNAs"
  tab_temp_Pos_all[tab_temp_Pos_all == "i"] = "int-lncRNAs"
  tab_temp_Pos_all[tab_temp_Pos_all == "o/e"] = "SOT-lncRNAs"
  tab_temp_Pos_all$"Count" = 1
  
  tab_temp_Pos_all_red = tab_temp_Pos_all %>%
    group_by(SRA.Study, PosLev.Conservation_3, Class_code, Type) %>%
    summarise(Counts = sum(Count)) %>%
    group_by(SRA.Study, PosLev.Conservation_3, Class_code) %>%
    mutate(Percent = (100*Counts)/sum(Counts)) %>%
    group_by(PosLev.Conservation_3, Class_code, Type) %>%
    summarise(Mean.Percent = mean(Percent), SE.Percent = se(Percent))
  tab_temp_Pos_all_red = as.data.frame(tab_temp_Pos_all_red)
  tab_temp_Pos_all_red$PosLev.Conservation_3 = factor(tab_temp_Pos_all_red$PosLev.Conservation_3, levels = 9:1)
  tab_temp_Pos_all_red$Class_code = factor(tab_temp_Pos_all_red$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  
  gg2 = ggplot(tab_temp_Pos_all_red, aes(x = PosLev.Conservation_3, y = Mean.Percent, ymin=Mean.Percent-SE.Percent, ymax=Mean.Percent+SE.Percent, group = Type, fill = Type)) +
    geom_bar(aes(group = Type, fill = Type), position=position_dodge(preserve="single"), stat="identity", color = "black", width = 0.7) +
    geom_errorbar(position=position_dodge(width=0.7, preserve="single"), colour="black", alpha = 1, width=0.4) +
    scale_fill_manual(values = c("#9ED4F5", "#6E93AA", "#2A4758")) +
    xlab("") +
    ylab("LncRNAs (%)") +
    facet_grid(. ~ Class_code) +
    theme_bw() +
    theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 10, margin = margin(r = 0.5, unit = 'cm')),
          legend.key.size = unit(0.8, "lines"),
          legend.text.align = 0) +
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    coord_flip() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
  
  gg2 = ggarrange(gg2, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), ncol = 1, heights = c(0.997, 0.003))
  
  # Joinn figures.
  gg_final_Pos = ggarrange(gg1, gg2, ncol = 2, widths = c(0.60, 0.40))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean-POS-", con, ".png"), height = 8, width = 16, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean-POS-", con, ".pdf"), height = 8, width = 16, dpi = 600, bg = "white")
  
  gg_final_Pos_facet = ggarrange(gg1_facet, gg2, ncol = 2, widths = c(0.60, 0.40))
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean-POS-", con, "_facet.png"), height = 8, width = 20, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean-POS-", con, "_facet.pdf"), height = 8, width = 20, dpi = 600, bg = "white")
  
  rm(list = c("tab_temp_Pos", "tab_temp_Pos_all", "tab_temp_Pos_all_red", "gg1", "gg1_facet", "gg2"))
  
  
  ############################################
  # Join Seq and Pos figures
  gg_final = ggarrange(gg_final_OF, gg_final_Pos, nrow = 2)
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean", con, ".png"), height = 16, width = 16, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean", con, ".pdf"), height = 16, width = 16, dpi = 600, bg = "white")
  
  gg_final_facet = ggarrange(gg_final_OF_facet, gg_final_Pos_facet, nrow = 2)
  
  ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean", con, "_facet.png"), height = 16, width = 20, dpi = 600, bg = "white")
  #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/Figures/CLASS_CODES/mean", con, "_facet.pdf"), height = 16, width = 20, dpi = 600, bg = "white")
  
  rm(list = c("gg_final_OF", "gg_final_Pos", "gg_final_OF_facet", "gg_final_Pos_facet", "gg_final", "gg_final_facet"))
}

rm(list = c("con"))











#################################
cat(paste0("\n\n--- MEAN TAU (CORRELATION BY PROJECT): \n"))

for (i in 1:length(species_short_name)) {
  
  spe = species_short_name[i]
  spe_l = species_long_name[i]
  
  cat(paste0("------ Specie: ", spe, "\n"))
  
  # Directories.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7"))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7"))
  }
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe))
  }
  
  # Get files with TAU metric values and log-transformed expression values.
  files = list.files(path = paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP2/", spe), pattern = ".tsv")
  
  if (length(files) > 0) {
    SRA.Studies = unique(unlist(lapply(strsplit(files, "_"), `[[`, 1)))
    
    for (SRA.Study in SRA.Studies) {
      cat(paste0("--------- SRA.Study: ", SRA.Study, "\n"))
      
      # Table
      tab_mean = TAB_mean[TAB_mean$SRA.Study == SRA.Study,]
      
      # Convert to factors.
      tab_mean$Confidence = factor(tab_mean$Confidence, levels = c("Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
      tab_mean$Type = factor(tab_mean$Type, levels = c("Non-analysed", "Non-TS", "TS"))
      
      
      ############################################
      ## COMPARISON == ALL
      cat(paste0("------------ COMPARISON == ALL\n"))
      
      if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL"))){
        dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL"))
      }
      
      for (con in confidences) {
        cat(paste0("--------------- ", con, '-confidence lncRNA\n'))
        
        ############################################
        # OrthoFinder A
        tab_temp_OF = tab_mean[tab_mean$Type != "Non-analysed" & tab_mean$Confidence == paste0(con, '-confidence lncRNA') & tab_mean$Comparison == "ALL", c("ID_transcript", "Class_code", "SeqLev.OF.Conservation_3", "TAU")]
        tab_temp_OF = tab_temp_OF[!duplicated(tab_temp_OF),]
        tab_temp_OF$Class_code = as.character(tab_temp_OF$Class_code)
        tab_temp_OF[tab_temp_OF == "u"] = "lincRNAs"
        tab_temp_OF[tab_temp_OF == "x"] = "NAT-lncRNAs"
        tab_temp_OF[tab_temp_OF == "i"] = "int-lncRNAs"
        tab_temp_OF[tab_temp_OF == "o/e"] = "SOT-lncRNAs"
        tab_temp_OF$SeqLev.OF.Conservation_3 = factor(tab_temp_OF$SeqLev.OF.Conservation_3, levels = 9:1)
        tab_temp_OF$Class_code = factor(tab_temp_OF$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
        
        gg1 = ggplot(tab_temp_OF, aes(x = SeqLev.OF.Conservation_3, y = TAU)) +
          geom_flat_violin(position = position_nudge(x = .18, y = 0), adjust = 1, trim = TRUE, alpha = 1, fill = "#f2c36a") +
          geom_point(aes(x = as.numeric(SeqLev.OF.Conservation_3)-.22, y = TAU),position = position_jitter(width = .05), size = .25, shape = 20) +
          geom_boxplot(aes(x = SeqLev.OF.Conservation_3, y = TAU), outlier.shape = NA, alpha = 1, width = .25, colour = "black") +
          stat_summary(aes(x = SeqLev.OF.Conservation_3), fun=mean, geom="point", shape=18, size=2, color="red") +
          scale_y_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
          theme_bw() +
          xlab("Conservation level (Sequence OF)") + ylab("TAU") +
          theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank()) +
          theme(axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14)) +
          coord_flip() +
          theme(plot.margin = margin(5.5, 0, 5.5, 5.5))
        
        gg1_facet = gg1 + facet_grid(. ~ Class_code)
        
        gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 1, heights = c(0.09, 0.91))
        
        gg1_facet = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1_facet, ncol = 1, heights = c(0.058, 0.942))
        
        # OrthoFinder B
        tab_temp_OF_all = tab_mean[tab_mean$Confidence == paste0(con, '-confidence lncRNA') & tab_mean$Comparison == "ALL", c("ID_transcript", "Class_code", "SeqLev.OF.Conservation_3", "Type")]
        tab_temp_OF_all = tab_temp_OF_all[!duplicated(tab_temp_OF_all),]
        tab_temp_OF_all$Class_code = as.character(tab_temp_OF_all$Class_code)
        tab_temp_OF_all[tab_temp_OF_all == "u"] = "lincRNAs"
        tab_temp_OF_all[tab_temp_OF_all == "x"] = "NAT-lncRNAs"
        tab_temp_OF_all[tab_temp_OF_all == "i"] = "int-lncRNAs"
        tab_temp_OF_all[tab_temp_OF_all == "o/e"] = "SOT-lncRNAs"
        tab_temp_OF_all$"Count" = 1
        
        tab_temp_OF_all_red = tab_temp_OF_all %>%
          group_by(SeqLev.OF.Conservation_3, Class_code, Type) %>%
          summarise(Counts = sum(Count)) %>%
          group_by(SeqLev.OF.Conservation_3, Class_code) %>%
          mutate(Percent = (100*Counts)/sum(Counts))
        tab_temp_OF_all_red = as.data.frame(tab_temp_OF_all_red)
        tab_temp_OF_all_red$SeqLev.OF.Conservation_3 = factor(tab_temp_OF_all_red$SeqLev.OF.Conservation_3, levels = 9:1)
        tab_temp_OF_all_red$Class_code = factor(tab_temp_OF_all_red$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
        
        gg2 = ggplot(tab_temp_OF_all_red, aes(x = SeqLev.OF.Conservation_3, y = Percent, group = Type, fill = Type)) +
          geom_bar(aes(group = Type, fill = Type), position=position_dodge(preserve="single"), stat="identity", color = "black", width = 0.7) +
          scale_fill_manual(values = c("#9ED4F5", "#6E93AA", "#2A4758")) +
          xlab("") +
          ylab("LncRNAs (%)") +
          facet_grid(. ~ Class_code) +
          theme_bw() +
          theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank()) +
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) +
          theme(legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size = 10, margin = margin(r = 0.5, unit = 'cm')),
                legend.key.size = unit(0.8, "lines"),
                legend.text.align = 0) +
          theme(axis.title.x = element_text(size = 14),
                axis.text.x = element_text(size = 12)) +
          coord_flip() +
          theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
        
        gg2 = ggarrange(gg2, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), ncol = 1, heights = c(0.997, 0.003))
        
        # Joinn figures.
        gg_final_OF = ggarrange(gg1, gg2, ncol = 2, widths = c(0.60, 0.40))
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_meanSEQ_OF-", con, ".png"), height = 8, width = 16, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_meanSEQ_OF-", con, ".pdf"), height = 8, width = 16, dpi = 600, bg = "white")
        
        gg_final_OF_facet = ggarrange(gg1_facet, gg2, ncol = 2, widths = c(0.60, 0.40))
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-SEQ_OF-", con, "_facet.png"), height = 8, width = 20, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-SEQ_OF-", con, "_facet.pdf"), height = 8, width = 20, dpi = 600, bg = "white")
        
        rm(list = c("tab_temp_OF", "tab_temp_OF_all", "tab_temp_OF_all_red", "gg1", "gg1_facet", "gg2"))
        
        
        ############################################
        # Positional A
        tab_temp_Pos = tab_mean[tab_mean$Type != "Non-analysed" & tab_mean$Confidence == paste0(con, '-confidence lncRNA') & tab_mean$Comparison == "ALL", c("ID_transcript", "Class_code", "PosLev.Conservation_3", "TAU")]
        tab_temp_Pos = tab_temp_Pos[!duplicated(tab_temp_Pos),]
        tab_temp_Pos$Class_code = as.character(tab_temp_Pos$Class_code)
        tab_temp_Pos[tab_temp_Pos == "u"] = "lincRNAs"
        tab_temp_Pos[tab_temp_Pos == "x"] = "NAT-lncRNAs"
        tab_temp_Pos[tab_temp_Pos == "i"] = "int-lncRNAs"
        tab_temp_Pos[tab_temp_Pos == "o/e"] = "SOT-lncRNAs"
        tab_temp_Pos$PosLev.Conservation_3 = factor(tab_temp_Pos$PosLev.Conservation_3, levels = 9:1)
        tab_temp_Pos$Class_code = factor(tab_temp_Pos$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
        
        gg1 = ggplot(tab_temp_Pos, aes(x = PosLev.Conservation_3, y = TAU)) +
          geom_flat_violin(position = position_nudge(x = .18, y = 0), adjust = 1, trim = TRUE, alpha = 1, fill = "#f2c36a") +
          geom_point(aes(x = as.numeric(PosLev.Conservation_3)-.22, y = TAU),position = position_jitter(width = .05), size = .25, shape = 20) +
          geom_boxplot(aes(x = PosLev.Conservation_3, y = TAU), outlier.shape = NA, alpha = 1, width = .25, colour = "black") +
          stat_summary(aes(x = PosLev.Conservation_3), fun=mean, geom="point", shape=18, size=2, color="red") +
          scale_y_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
          theme_bw() +
          xlab("Conservation level (Positonal)") + ylab("TAU") +
          theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank()) +
          theme(axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14)) +
          coord_flip() +
          theme(plot.margin = margin(5.5, 0, 5.5, 5.5))
        
        gg1_facet = gg1 + facet_grid(. ~ Class_code)
        
        gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 1, heights = c(0.09, 0.91))
        
        gg1_facet = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1_facet, ncol = 1, heights = c(0.058, 0.942))
        
        # Positional B
        tab_temp_Pos_all = tab_mean[tab_mean$Confidence == paste0(con, '-confidence lncRNA') & tab_mean$Comparison == "ALL", c("ID_transcript", "Class_code", "PosLev.Conservation_3", "Type")]
        tab_temp_Pos_all = tab_temp_Pos_all[!duplicated(tab_temp_Pos_all),]
        tab_temp_Pos_all$Class_code = as.character(tab_temp_Pos_all$Class_code)
        tab_temp_Pos_all[tab_temp_Pos_all == "u"] = "lincRNAs"
        tab_temp_Pos_all[tab_temp_Pos_all == "x"] = "NAT-lncRNAs"
        tab_temp_Pos_all[tab_temp_Pos_all == "i"] = "int-lncRNAs"
        tab_temp_Pos_all[tab_temp_Pos_all == "o/e"] = "SOT-lncRNAs"
        tab_temp_Pos_all$"Count" = 1
        
        tab_temp_Pos_all_red = tab_temp_Pos_all %>%
          group_by(PosLev.Conservation_3, Class_code, Type) %>%
          summarise(Counts = sum(Count)) %>%
          group_by(PosLev.Conservation_3, Class_code) %>%
          mutate(Percent = (100*Counts)/sum(Counts))
        tab_temp_Pos_all_red = as.data.frame(tab_temp_Pos_all_red)
        tab_temp_Pos_all_red$PosLev.Conservation_3 = factor(tab_temp_Pos_all_red$PosLev.Conservation_3, levels = 9:1)
        tab_temp_Pos_all_red$Class_code = factor(tab_temp_Pos_all_red$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
        
        gg2 = ggplot(tab_temp_Pos_all_red, aes(x = PosLev.Conservation_3, y = Percent, group = Type, fill = Type)) +
          geom_bar(aes(group = Type, fill = Type), position=position_dodge(preserve="single"), stat="identity", color = "black", width = 0.7) +
          scale_fill_manual(values = c("#9ED4F5", "#6E93AA", "#2A4758")) +
          xlab("") +
          ylab("LncRNAs (%)") +
          facet_grid(. ~ Class_code) +
          theme_bw() +
          theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank()) +
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) +
          theme(legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size = 10, margin = margin(r = 0.5, unit = 'cm')),
                legend.key.size = unit(0.8, "lines"),
                legend.text.align = 0) +
          theme(axis.title.x = element_text(size = 14),
                axis.text.x = element_text(size = 12)) +
          coord_flip() +
          theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
        
        gg2 = ggarrange(gg2, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), ncol = 1, heights = c(0.997, 0.003))
        
        # Joinn figures.
        gg_final_Pos = ggarrange(gg1, gg2, ncol = 2, widths = c(0.60, 0.40))
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-POS-", con, ".png"), height = 8, width = 16, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-POS-", con, ".pdf"), height = 8, width = 16, dpi = 600, bg = "white")
        
        gg_final_Pos_facet = ggarrange(gg1_facet, gg2, ncol = 2, widths = c(0.60, 0.40))
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-POS-", con, "_facet.png"), height = 8, width = 20, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-POS-", con, "_facet.pdf"), height = 8, width = 20, dpi = 600, bg = "white")
        
        rm(list = c("tab_temp_Pos", "tab_temp_Pos_all", "tab_temp_Pos_all_red", "gg1", "gg1_facet", "gg2"))
        
        
        ############################################
        # Join Seq and Pos figures
        gg_final = ggarrange(gg_final_OF, gg_final_Pos, nrow = 2)
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-", con, ".png"), height = 16, width = 16, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-", con, ".pdf"), height = 16, width = 16, dpi = 600, bg = "white")
        
        gg_final_facet = ggarrange(gg_final_OF_facet, gg_final_Pos_facet, nrow = 2)
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-", con, "_facet.png"), height = 16, width = 20, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/ALL/", SRA.Study, "_mean-", con, "_facet.pdf"), height = 16, width = 20, dpi = 600, bg = "white")
        
        rm(list = c("gg_final_OF", "gg_final_Pos", "gg_final_OF_facet", "gg_final_Pos_facet", "gg_final", "gg_final_facet"))
      }
      
      rm(list = c("con"))
      
      
      ############################################
      ## COMPARISON != ALL
      cat(paste0("------------ COMPARISON != ALL\n"))
      
      if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES"))){
        dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES"))
      }
      
      for (con in confidences) {
        cat(paste0("--------------- ", con, '-confidence lncRNA\n'))
        
        ############################################
        # OrthoFinder A
        tab_temp_OF = tab_mean[tab_mean$Type != "Non-analysed" & tab_mean$Confidence == paste0(con, '-confidence lncRNA') & tab_mean$Comparison != "ALL", c("ID_transcript", "Class_code", "SeqLev.OF.Conservation_3", "TAU")]
        tab_temp_OF = tab_temp_OF[!duplicated(tab_temp_OF),]
        tab_temp_OF$Class_code = as.character(tab_temp_OF$Class_code)
        tab_temp_OF[tab_temp_OF == "u"] = "lincRNAs"
        tab_temp_OF[tab_temp_OF == "x"] = "NAT-lncRNAs"
        tab_temp_OF[tab_temp_OF == "i"] = "int-lncRNAs"
        tab_temp_OF[tab_temp_OF == "o/e"] = "SOT-lncRNAs"
        tab_temp_OF$SeqLev.OF.Conservation_3 = factor(tab_temp_OF$SeqLev.OF.Conservation_3, levels = 9:1)
        tab_temp_OF$Class_code = factor(tab_temp_OF$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
        
        gg1 = ggplot(tab_temp_OF, aes(x = SeqLev.OF.Conservation_3, y = TAU)) +
          geom_flat_violin(position = position_nudge(x = .18, y = 0), adjust = 1, trim = TRUE, alpha = 1, fill = "#f2c36a") +
          geom_point(aes(x = as.numeric(SeqLev.OF.Conservation_3)-.22, y = TAU),position = position_jitter(width = .05), size = .25, shape = 20) +
          geom_boxplot(aes(x = SeqLev.OF.Conservation_3, y = TAU), outlier.shape = NA, alpha = 1, width = .25, colour = "black") +
          stat_summary(aes(x = SeqLev.OF.Conservation_3), fun=mean, geom="point", shape=18, size=2, color="red") +
          scale_y_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
          theme_bw() +
          xlab("Conservation level (Sequence OF)") + ylab("TAU") +
          theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank()) +
          theme(axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14)) +
          coord_flip() +
          theme(plot.margin = margin(5.5, 0, 5.5, 5.5))
        
        gg1_facet = gg1 + facet_grid(. ~ Class_code)
        
        gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 1, heights = c(0.09, 0.91))
        
        gg1_facet = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1_facet, ncol = 1, heights = c(0.058, 0.942))
        
        # OrthoFinder B
        tab_temp_OF_all = tab_mean[tab_mean$Confidence == paste0(con, '-confidence lncRNA') & tab_mean$Comparison != "ALL", c("ID_transcript", "Class_code", "SeqLev.OF.Conservation_3", "Type")]
        tab_temp_OF_all = tab_temp_OF_all[!duplicated(tab_temp_OF_all),]
        tab_temp_OF_all$Class_code = as.character(tab_temp_OF_all$Class_code)
        tab_temp_OF_all[tab_temp_OF_all == "u"] = "lincRNAs"
        tab_temp_OF_all[tab_temp_OF_all == "x"] = "NAT-lncRNAs"
        tab_temp_OF_all[tab_temp_OF_all == "i"] = "int-lncRNAs"
        tab_temp_OF_all[tab_temp_OF_all == "o/e"] = "SOT-lncRNAs"
        tab_temp_OF_all$"Count" = 1
        
        tab_temp_OF_all_red = tab_temp_OF_all %>%
          group_by(SeqLev.OF.Conservation_3, Class_code, Type) %>%
          summarise(Counts = sum(Count)) %>%
          group_by(SeqLev.OF.Conservation_3, Class_code) %>%
          mutate(Percent = (100*Counts)/sum(Counts))
        tab_temp_OF_all_red = as.data.frame(tab_temp_OF_all_red)
        tab_temp_OF_all_red$SeqLev.OF.Conservation_3 = factor(tab_temp_OF_all_red$SeqLev.OF.Conservation_3, levels = 9:1)
        tab_temp_OF_all_red$Class_code = factor(tab_temp_OF_all_red$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
        
        gg2 = ggplot(tab_temp_OF_all_red, aes(x = SeqLev.OF.Conservation_3, y = Percent, group = Type, fill = Type)) +
          geom_bar(aes(group = Type, fill = Type), position=position_dodge(preserve="single"), stat="identity", color = "black", width = 0.7) +
          scale_fill_manual(values = c("#9ED4F5", "#6E93AA", "#2A4758")) +
          xlab("") +
          ylab("LncRNAs (%)") +
          facet_grid(. ~ Class_code) +
          theme_bw() +
          theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank()) +
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) +
          theme(legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size = 10, margin = margin(r = 0.5, unit = 'cm')),
                legend.key.size = unit(0.8, "lines"),
                legend.text.align = 0) +
          theme(axis.title.x = element_text(size = 14),
                axis.text.x = element_text(size = 12)) +
          coord_flip() +
          theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
        
        gg2 = ggarrange(gg2, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), ncol = 1, heights = c(0.997, 0.003))
        
        # Joinn figures.
        gg_final_OF = ggarrange(gg1, gg2, ncol = 2, widths = c(0.60, 0.40))
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-SEQ_OF-", con, ".png"), height = 8, width = 16, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-SEQ_OF-", con, ".pdf"), height = 8, width = 16, dpi = 600, bg = "white")
        
        gg_final_OF_facet = ggarrange(gg1_facet, gg2, ncol = 2, widths = c(0.60, 0.40))
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-SEQ_OF-", con, "_facet.png"), height = 8, width = 20, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-SEQ_OF-", con, "_facet.pdf"), height = 8, width = 20, dpi = 600, bg = "white")
        
        rm(list = c("tab_temp_OF", "tab_temp_OF_all", "tab_temp_OF_all_red", "gg1", "gg1_facet", "gg2"))
        
        
        ############################################
        # Positional A
        tab_temp_Pos = tab_mean[tab_mean$Type != "Non-analysed" & tab_mean$Confidence == paste0(con, '-confidence lncRNA') & tab_mean$Comparison != "ALL", c("ID_transcript", "Class_code", "PosLev.Conservation_3", "TAU")]
        tab_temp_Pos = tab_temp_Pos[!duplicated(tab_temp_Pos),]
        tab_temp_Pos$Class_code = as.character(tab_temp_Pos$Class_code)
        tab_temp_Pos[tab_temp_Pos == "u"] = "lincRNAs"
        tab_temp_Pos[tab_temp_Pos == "x"] = "NAT-lncRNAs"
        tab_temp_Pos[tab_temp_Pos == "i"] = "int-lncRNAs"
        tab_temp_Pos[tab_temp_Pos == "o/e"] = "SOT-lncRNAs"
        tab_temp_Pos$PosLev.Conservation_3 = factor(tab_temp_Pos$PosLev.Conservation_3, levels = 9:1)
        tab_temp_Pos$Class_code = factor(tab_temp_Pos$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
        
        gg1 = ggplot(tab_temp_Pos, aes(x = PosLev.Conservation_3, y = TAU)) +
          geom_flat_violin(position = position_nudge(x = .18, y = 0), adjust = 1, trim = TRUE, alpha = 1, fill = "#f2c36a") +
          geom_point(aes(x = as.numeric(PosLev.Conservation_3)-.22, y = TAU),position = position_jitter(width = .05), size = .25, shape = 20) +
          geom_boxplot(aes(x = PosLev.Conservation_3, y = TAU), outlier.shape = NA, alpha = 1, width = .25, colour = "black") +
          stat_summary(aes(x = PosLev.Conservation_3), fun=mean, geom="point", shape=18, size=2, color="red") +
          scale_y_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
          theme_bw() +
          xlab("Conservation level (Positonal)") + ylab("TAU") +
          theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank()) +
          theme(axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14)) +
          coord_flip() +
          theme(plot.margin = margin(5.5, 0, 5.5, 5.5))
        
        gg1_facet = gg1 + facet_grid(. ~ Class_code)
        
        gg1 = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1, ncol = 1, heights = c(0.09, 0.91))
        
        gg1_facet = ggarrange(ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg1_facet, ncol = 1, heights = c(0.058, 0.942))
        
        # Positional B
        tab_temp_Pos_all = tab_mean[tab_mean$Confidence == paste0(con, '-confidence lncRNA') & tab_mean$Comparison != "ALL", c("ID_transcript", "Class_code", "PosLev.Conservation_3", "Type")]
        tab_temp_Pos_all = tab_temp_Pos_all[!duplicated(tab_temp_Pos_all),]
        tab_temp_Pos_all$Class_code = as.character(tab_temp_Pos_all$Class_code)
        tab_temp_Pos_all[tab_temp_Pos_all == "u"] = "lincRNAs"
        tab_temp_Pos_all[tab_temp_Pos_all == "x"] = "NAT-lncRNAs"
        tab_temp_Pos_all[tab_temp_Pos_all == "i"] = "int-lncRNAs"
        tab_temp_Pos_all[tab_temp_Pos_all == "o/e"] = "SOT-lncRNAs"
        tab_temp_Pos_all$"Count" = 1
        
        tab_temp_Pos_all_red = tab_temp_Pos_all %>%
          group_by(PosLev.Conservation_3, Class_code, Type) %>%
          summarise(Counts = sum(Count)) %>%
          group_by(PosLev.Conservation_3, Class_code) %>%
          mutate(Percent = (100*Counts)/sum(Counts))
        tab_temp_Pos_all_red = as.data.frame(tab_temp_Pos_all_red)
        tab_temp_Pos_all_red$PosLev.Conservation_3 = factor(tab_temp_Pos_all_red$PosLev.Conservation_3, levels = 9:1)
        tab_temp_Pos_all_red$Class_code = factor(tab_temp_Pos_all_red$Class_code, levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
        
        gg2 = ggplot(tab_temp_Pos_all_red, aes(x = PosLev.Conservation_3, y = Percent, group = Type, fill = Type)) +
          geom_bar(aes(group = Type, fill = Type), position=position_dodge(preserve="single"), stat="identity", color = "black", width = 0.7) +
          scale_fill_manual(values = c("#9ED4F5", "#6E93AA", "#2A4758")) +
          xlab("") +
          ylab("LncRNAs (%)") +
          facet_grid(. ~ Class_code) +
          theme_bw() +
          theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank()) +
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) +
          theme(legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size = 10, margin = margin(r = 0.5, unit = 'cm')),
                legend.key.size = unit(0.8, "lines"),
                legend.text.align = 0) +
          theme(axis.title.x = element_text(size = 14),
                axis.text.x = element_text(size = 12)) +
          coord_flip() +
          theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
        
        gg2 = ggarrange(gg2, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), ncol = 1, heights = c(0.997, 0.003))
        
        # Joinn figures.
        gg_final_Pos = ggarrange(gg1, gg2, ncol = 2, widths = c(0.60, 0.40))
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-POS-", con, ".png"), height = 8, width = 16, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-POS-", con, ".pdf"), height = 8, width = 16, dpi = 600, bg = "white")
        
        gg_final_Pos_facet = ggarrange(gg1_facet, gg2, ncol = 2, widths = c(0.60, 0.40))
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-POS-", con, "_facet.png"), height = 8, width = 20, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-POS-", con, "_facet.pdf"), height = 8, width = 20, dpi = 600, bg = "white")
        
        rm(list = c("tab_temp_Pos", "tab_temp_Pos_all", "tab_temp_Pos_all_red", "gg1", "gg1_facet", "gg2"))
        
        
        ############################################
        # Join Seq and Pos figures
        gg_final = ggarrange(gg_final_OF, gg_final_Pos, nrow = 2)
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-", con, ".png"), height = 16, width = 16, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-", con, ".pdf"), height = 16, width = 16, dpi = 600, bg = "white")
        
        gg_final_facet = ggarrange(gg_final_OF_facet, gg_final_Pos_facet, nrow = 2)
        
        ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-", con, "_facet.png"), height = 16, width = 20, dpi = 600, bg = "white")
        #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-", con, "_facet.pdf"), height = 16, width = 20, dpi = 600, bg = "white")
        
        rm(list = c("gg_final_OF", "gg_final_Pos", "gg_final_OF_facet", "gg_final_Pos_facet", "gg_final", "gg_final_facet"))
      }
      
      rm(list = c("con", "tab_mean"))
    }
  }
}

rm(list = c("files", "spe", "spe_l", "SRA.Study", "SRA.Studies", "i"))


# # Motifs
# tab_mean_red_Mot = tab_mean_merged[tab_mean_merged$Comparison != "ALL" & tab_mean_merged$PosLev.Conservation_3 != 1, c("ID_transcript", "Comparison", "Confidence", "Class_code", "Type", "MotLev.Width", "MotLev.Meme_Motif.Identifier")]
# tab_mean_red_Mot = tab_mean_red_Mot[tab_mean_red_Mot$MotLev.Width == "6-50", c("ID_transcript", "Comparison", "Confidence", "Class_code", "Type", "MotLev.Meme_Motif.Identifier")]
# tab_mean_red_Mot$"MotLev.Meme_Motif.Presence_Ausence" = ifelse(!is.na(tab_mean_red_Mot$MotLev.Meme_Motif.Identifier), "Presence", "Ausence")
# tab_mean_red_Mot = tab_mean_red_Mot[, c("ID_transcript", "Comparison", "Confidence", "Class_code", "Type", "MotLev.Meme_Motif.Presence_Ausence")]
# tab_mean_red_Mot = tab_mean_red_Mot[!duplicated(tab_mean_red_Mot),]
# tab_mean_red_Mot$MotLev.Meme_Motif.Presence_Ausence = factor(tab_mean_red_Mot$MotLev.Meme_Motif.Presence_Ausence, levels = c("Presence", "Ausence"))
# 
# tab_mean_collapsed_Mot = tab_mean_red_Mot %>%
#   group_by(Comparison, Confidence, Class_code, MotLev.Meme_Motif.Presence_Ausence, Type) %>%
#   summarise(Counts = n_distinct(ID_transcript)) %>%
#   group_by( Comparison, Confidence, Class_code, MotLev.Meme_Motif.Presence_Ausence) %>%
#   mutate(Percent = (100*Counts)/sum(Counts))
# tab_mean_collapsed_Mot = as.data.frame(tab_mean_collapsed_Mot)
# 
# gg = ggplot(tab_mean_collapsed_Mot, aes(x = MotLev.Meme_Motif.Presence_Ausence, y = Percent, group = Type, alpha = Type)) +
#   geom_bar(aes(group = Type, alpha = Type), position=position_dodge(preserve = "single"), stat="identity", fill = "blue", color = "black") +
#   scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
#   facet_grid(Confidence ~ Class_code) +
#   xlab("") +
#   ylab("LncRNAs (%)") +
#   theme_bw() +
#   theme(legend.position = "top") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
#   theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, angle = 60)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.3)
# 
# ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-MOT.png"), height = 10, width = 15, dpi = 600)
# #ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP7/", spe, "/CLASS_CODES/", SRA.Study, "_mean-MOT.pdf"), height = 10, width = 15, dpi = 600)







## 4. MEDIAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEDIAN TAU: \n"))

















