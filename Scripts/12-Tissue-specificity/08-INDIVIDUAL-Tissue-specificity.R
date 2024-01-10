################################################################################
#
# ALL: TISSUE SPECIFICITY STUDY: APPROACH 1 - STEP 8
# 
################################################################################


rm(list = ls())


## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))


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
## 2. MEAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEAN TAU (TABLE): \n"))

TAB_mean = data.frame()

for (i in 1:length(species_short_name)) {
  
  spe = species_short_name[i]
  spe_l = species_long_name[i]
  
  cat(paste0("------ Specie: ", spe, "\n"))
  
  # Directories.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8"))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8"))
  }
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/", spe))
  }
  
  # Get files with TAU metric values and log-transformed expression values.
  files = list.files(path = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe), pattern = ".tsv")
  
  if (length(files) > 0) {
    SRA.Studies = unique(unlist(lapply(strsplit(unique(unlist(lapply(strsplit(files, "_"), `[[`, 1))), "-"), `[[`, 2)))
    
    for (SRA.Study in SRA.Studies) {
      cat(paste0("--------- SRA.Study: ", SRA.Study, "\n"))
      
      # Load table with expression values. 
      tab_mean_TAU_1_G = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP1/", spe, "/GENES-", SRA.Study, "_mean.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU_1_L = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP1/", spe, "/LNCRNAS-", SRA.Study, "_mean.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU_1 = rbind(tab_mean_TAU_1_G, tab_mean_TAU_1_L)
      tab_mean_TAU_1_long = melt(setDT(tab_mean_TAU_1), id.vars = c("ID_transcript", "Confidence", "Class_code"), variable.name = "Tissue", value.name = "TPMs.mean")
      tab_mean_TAU_1_long_filt = tab_mean_TAU_1_long %>% group_by(ID_transcript) %>% filter(TPMs.mean == max(TPMs.mean))
      tab_mean_TAU_1_long_filt$ID_transcript = paste0(tab_mean_TAU_1_long_filt$ID_transcript, "-", spe)
      
      # Load table with the tissue specificity results and the log-tranformed expression values.
      tab_mean_TAU_2_G = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/GENES-", SRA.Study, "_mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU_2_L = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/LNCRNAS-", SRA.Study, "_mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU_2 = rbind(tab_mean_TAU_2_G, tab_mean_TAU_2_L)
      tab_mean_TAU_2 = tab_mean_TAU_2[, c("ID_transcript", "Confidence", "Class_code", "TAU")]
      tab_mean_TAU_2$ID_transcript = paste0(tab_mean_TAU_2$ID_transcript, "-", spe)
      
      # Load table STEP 7 with conservation data.
      tab_mean_TAU_3 = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP7/", spe, "/", SRA.Study, "_mean.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU_3 = tab_mean_TAU_3[, c("ID_transcript", "Confidence", "Class_code", "SeqLev.OF.Family", "SeqLev.OF.Conservation_1", "SeqLev.OF.Conservation_2", "SeqLev.OF.Conservation_3", 
                                          "PosLev.Family", "PosLev.Conservation_1", "PosLev.Conservation_2", "PosLev.Conservation_3", "MotLev.Mode", "MotLev.Width", "MotLev.Meme_Motif.Identifier",
                                          "MotLev.Meme_Motif.E.value","MotLev.Meme_Motif.LncRNA", "MotLev.Meme_Motif.P.value", "MotLev.Meme_Motif.Width", "MotLev.Gomo_GO.Term.Identifier",
                                          "MotLev.Gomo_Q.value")]
      
      # Join tables
      tab_mean_TAU = merge(tab_mean_TAU_1_long_filt, tab_mean_TAU_2, by = c("ID_transcript", "Confidence", "Class_code"), all = T)
      tab_mean_TAU = merge(tab_mean_TAU, tab_mean_TAU_3, by = c("ID_transcript", "Confidence", "Class_code"), all = T)
      
      rm(list = c("tab_mean_TAU_1", "tab_mean_TAU_1_G", "tab_mean_TAU_1_L", "tab_mean_TAU_1_long", "tab_mean_TAU_1_long_filt", "tab_mean_TAU_2", "tab_mean_TAU_2_G", "tab_mean_TAU_2_L", "tab_mean_TAU_3"))
      
      # Add column to specify TAU > 0.8. NA values refer to lncRNAs and genes that had no more than 1 TPM in any of the tissues and those that were duplicated during quantification.
      tab_mean_TAU$"Type" = ifelse(is.na(tab_mean_TAU$TAU), "Non-analysed", ifelse(tab_mean_TAU$TAU >= 0.8, "TS", "Non-TS"))
      
      # Change values
      tab_mean_TAU[tab_mean_TAU == "u"] = "lincRNAs"
      tab_mean_TAU[tab_mean_TAU == "x"] = "NAT-lncRNAs"
      tab_mean_TAU[tab_mean_TAU == "i"] = "int-lncRNAs"
      tab_mean_TAU[tab_mean_TAU == "o/e"] = "SOT-lncRNAs"
      tab_mean_TAU[tab_mean_TAU == "pc"] = "PC genes"
      
      # Save table.
      write.table(tab_mean_TAU, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/", spe, "/", SRA.Study, "_mean-ALL.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      
      # Add info.
      tab_mean_TAU = cbind(SRA.Study = SRA.Study, tab_mean_TAU)
      tab_mean_TAU = cbind(Spe = spe, tab_mean_TAU)
      TAB_mean = rbind(TAB_mean, tab_mean_TAU)
      
      rm(list = c("tab_mean_TAU"))
    }
  }
}

write.table(TAB_mean, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/mean-ALL.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

rm(list = c("files", "spe", "spe_l", "SRA.Study", "SRA.Studies", "i"))













#################################

cat(paste0("\n\n--- MEAN TAU (CORRELATION - FUSION PROJECTS): \n"))

# Directories.
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8"))
}
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/Figures"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/Figures"))
}

for (con in confidences) {
  
  cat(paste0("------ ", con, '-confidence lncRNA\n'))
  
  # Filter
  TAB_mean_temp = tryCatch({TAB_mean[TAB_mean$Type != "Non-analysed" & TAB_mean$Confidence %in% c(paste0(con, '-confidence lncRNA'), "Protein-coding"), c("Spe", "SRA.Study", "ID_transcript", "Class_code", "Tissue", "TPMs.mean", "TAU", "Type")]},
                     error = function(e) {NULL},
                     warning = function(e) {message(paste0("There is no rows for ", con, " confidence-level."))}
  )
  
  if (!is.null(TAB_mean_temp)) {
    
    # Remove duplicated rows.
    TAB_mean_temp = TAB_mean_temp[!duplicated(TAB_mean_temp),]
    
    # Convert to factors.
    TAB_mean_temp$Type = factor(TAB_mean_temp$Type, levels = c("Non-analysed", "Non-TS", "TS"))
    TAB_mean_temp$Class_code = factor(TAB_mean_temp$Class_code, levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
    TAB_mean_temp$Tissue = factor(TAB_mean_temp$Tissue, levels = sort(unique(as.character(TAB_mean_temp$Tissue))))
    
    # Filter by number of transcript ids >= 50.
    TAB_mean_temp_summary = TAB_mean_temp %>% group_by(Class_code, Tissue) %>% summarize(Counts = n()) %>% filter(Counts >= 50)
    TAB_mean_temp = TAB_mean_temp[paste0(TAB_mean_temp$Class_code, "|", TAB_mean_temp$Tissue) %in% paste0(TAB_mean_temp_summary$Class_code, "|", TAB_mean_temp_summary$Tissue), ]
    
    # Figure
    gg = ggplot(TAB_mean_temp, aes(x = TAU, y = TPMs.mean)) +
      geom_point(position = position_jitter(width = .05), size = 0.1, shape = 20) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 4)) +
      facet_grid(Tissue ~ Class_code) +
      scale_x_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000)) +
      theme_bw() +
      xlab("TAU") + ylab("Max TPM") +
      theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank()) +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14))
      
    ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/Figures/mean-", con, ".png"), height = 14, width = 16, dpi = 600, bg = "white")
    #ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/Figures/mean-", con, ".pdf"), height = 14, width = 16, dpi = 600, bg = "white")
  }
}

rm(list = c("con", "TAB_mean_temp", "TAB_mean_temp_summary", "gg"))









#################################

cat(paste0("\n\n--- MEAN TAU (CORRELATION - BY SPECIE): \n"))

for (i in 1:length(species_short_name)) {
  
  spe = species_short_name[i]
  
  cat(paste0("------ Specie: ", spe, "\n"))

  # Directories.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8"))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8"))
  }
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/", spe))
  }
  
  for (con in confidences) {
    
    cat(paste0("--------- ", con, '-confidence lncRNA\n'))
    
    # Filter
    TAB_mean_temp = tryCatch({TAB_mean[TAB_mean$Spe == spe & TAB_mean$Type != "Non-analysed" & TAB_mean$Confidence %in% c(paste0(con, '-confidence lncRNA'), "Protein-coding"), c("Spe", "SRA.Study", "ID_transcript", "Class_code", "Tissue", "TPMs.mean", "TAU", "Type")]},
                             error = function(e) {NULL},
                             warning = function(e) {message(paste0("There is no rows for specie ", spe, " and ", con, " confidence-level."))}
    )
    
    if (!is.null(TAB_mean_temp)) {
    
      TAB_mean_temp = TAB_mean_temp[!duplicated(TAB_mean_temp),]
      
      # Convert to factors.
      TAB_mean_temp$Type = factor(TAB_mean_temp$Type, levels = c("Non-analysed", "Non-TS", "TS"))
      TAB_mean_temp$Class_code = factor(TAB_mean_temp$Class_code, levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
      TAB_mean_temp$Tissue = factor(TAB_mean_temp$Tissue, levels = sort(unique(as.character(TAB_mean_temp$Tissue))))
      
      # Filter by number of transcript ids >= 50.
      TAB_mean_temp_summary = TAB_mean_temp %>% group_by(Class_code, Tissue) %>% summarize(Counts = n()) %>% filter(Counts >= 50)
      TAB_mean_temp = TAB_mean_temp[paste0(TAB_mean_temp$Class_code, "|", TAB_mean_temp$Tissue) %in% paste0(TAB_mean_temp_summary$Class_code, "|", TAB_mean_temp_summary$Tissue), ]
      
      # Figure
      gg = ggplot(TAB_mean_temp, aes(x = TAU, y = TPMs.mean)) +
        geom_point(position = position_jitter(width = .05), size = 0.1, shape = 20) +
        stat_smooth(method = "lm", formula = y ~ poly(x, 4)) +
        facet_grid(Tissue ~ Class_code) +
        scale_x_continuous(labels = c(0,0.25,0.50,0.75,1), breaks = c(0,0.25,0.50,0.75,1), expand = c(0.015, 0.015)) +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000)) +
        theme_bw() +
        xlab("TAU") + ylab("Max TPM") +
        theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype = "dashed"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.background = element_blank()) +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))
      
      ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/", spe, "/mean-", con, ".png"), height = 14, width = 16, dpi = 600, bg = "white")
      #ggsave(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP8/", spe, "/mean-", con, ".pdf"), height = 14, width = 16, dpi = 600, bg = "white")
    }
  }
}

rm(list = c("con", "TAB_mean_temp", "TAB_mean_temp_summary", "gg", "i", "spe"))



















# # Convert to factors.
# tab_mean_TAU$Confidence = factor(tab_mean_TAU$Confidence, levels = c("Protein-coding", "Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
# tab_mean_TAU$Class_code = factor(tab_mean_TAU$Class_code, levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
# tab_mean_TAU$Type = factor(tab_mean_TAU$Type, levels = c("Non-analysed", "Non-TS", "TS"))
# tab_mean_TAU$SeqLev.OF.Conservation_1 = factor(tab_mean_TAU$SeqLev.OF.Conservation_1, levels = c("Non-conserved", "Conserved"))
# tab_mean_TAU$SeqLev.OF.Conservation_2 = factor(tab_mean_TAU$SeqLev.OF.Conservation_2, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
# tab_mean_TAU$SeqLev.OF.Conservation_3 = factor(tab_mean_TAU$SeqLev.OF.Conservation_3, levels = 1:9)
# tab_mean_TAU$PosLev.Conservation_1 = factor(tab_mean_TAU$PosLev.Conservation_1, levels = c("Non-conserved", "Conserved"))
# tab_mean_TAU$PosLev.Conservation_2 = factor(tab_mean_TAU$PosLev.Conservation_2, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
# tab_mean_TAU$PosLev.Conservation_3 = factor(tab_mean_TAU$PosLev.Conservation_3, levels = 1:9)
# 













## 3. MEDIAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEDIAN TAU: \n"))

















