################################################################################
#
# ALL: TISSUE SPECIFICITY STUDY: STEP 3
#
# Draw Violin plots to compare tissue-specificity in PCGs and lncRNAs.
#
# @author: pasviber - Pascual Villalba Bermell
# 
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggExtra))
suppressMessages(library(data.table))
suppressMessages(options(bitmapType='cairo'))

## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied.", call.=FALSE)
} else {
  path_tissue_specificity = args[1]
  flag = args[2]
  specie = args[3]
}

# path_tissue_specificity = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/12-Tissue-specificity/cme"
# flag = "nr"
# specie = "cme"

## 2. MEAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEAN TAU: \n"))

TAB_mean_TAU_mod = data.frame()

# Get files with TAU metric values and log-transformed expression values.
files = list.files(path = paste0(path_tissue_specificity, "/ALL/", flag, "/STEP2/Studies"), pattern = ".tsv")

if (length(files) > 0) {
  SRA.Studies = unique(unlist(lapply(strsplit(files, "_"), `[[`, 1)))
  
  for (SRA.Study in SRA.Studies) {
    cat(paste0("--------- SRA.Study: ", SRA.Study, "\n"))
    
    # Load table with the tissue specificity results and the log-tranformed expression values. 
    tab_mean_TAU = read.table(paste0(path_tissue_specificity, "/ALL/", flag, "/STEP2/Studies/", SRA.Study, "_mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")
    
    # Convert to factors.
    tab_mean_TAU$Confidence = factor(tab_mean_TAU$Confidence, levels = c("Protein-coding", "Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
    tab_mean_TAU$Class_code = factor(tab_mean_TAU$Class_code, levels = c("pc", "u", "x", "i", "o/e"))
    
    # Violin plot
    gg = ggplot(tab_mean_TAU, aes(x = Confidence, y = TAU, fill = Confidence)) + 
      geom_violin(color = "grey", position = position_dodge(1)) +
      geom_jitter(width = 0.4, color = "#8c8f8f", size = 0.1, alpha = 0.3) +
      stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
      geom_boxplot(width = 0.1, color = "black", alpha = 0.2, outlier.shape = NA, position = position_dodge(1)) +
      scale_y_continuous(limits = c(-0.001,1.001), breaks = seq(0,1,0.1), expand = c(0.01, 0.01)) + 
      theme_classic() + 
      scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee")) +
      scale_x_discrete(labels=c("Protein-coding" = "Protein-coding", 
                                "Low-confidence lncRNA" = "Low-confidence\nlncRNA", 
                                "Medium-confidence lncRNA" = "Medium-confidence\nlncRNA",
                                "High-confidence lncRNA" = "High-confidence\nlncRNA")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            panel.border = element_blank(),
            legend.position="none",
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16)) +
      xlab("") + ylab("tau")
    
    ggsave(paste0(path_tissue_specificity, "/ALL/", flag, "/STEP3/Studies/", SRA.Study, "_mean-TAU-Violin.png"), height = 10, width = 10, dpi = 600)
    ggsave(paste0(path_tissue_specificity, "/ALL/", flag, "/STEP3/Studies/", SRA.Study, "_mean-TAU-Violin.pdf"), height = 10, width = 10, dpi = 600)
    
    rm(list = c("gg"))
    
    # Violin plot by class code.
    tab_mean_TAU_mod = tab_mean_TAU[tab_mean_TAU$Class_code != "pc", c("ID_transcript", "Confidence", "Class_code", "TAU")]
    for (c in c("u", "x", "i", "o/e")) {
      pc = tab_mean_TAU[tab_mean_TAU$Class_code == "pc", c("ID_transcript", "Confidence", "Class_code", "TAU")]
      pc[pc == "pc"] = c
      tab_mean_TAU_mod = rbind(tab_mean_TAU_mod, pc)
    }
    tab_mean_TAU_mod$Class_code = ifelse(tab_mean_TAU_mod$Class_code == "u", "Intergenic lncRNA (u)", 
                                         ifelse(tab_mean_TAU_mod$Class_code == "x", "Antisense lncRNA (x)",
                                                ifelse(tab_mean_TAU_mod$Class_code == "i", "Intronic lncRNA (i)",
                                                       ifelse(tab_mean_TAU_mod$Class_code == "o/e", "Sense lncRNA (o/e)", tab_mean_TAU_mod$Class_code))))
    tab_mean_TAU_mod$Class_code = factor(tab_mean_TAU_mod$Class_code, levels = c("Intergenic lncRNA (u)", "Antisense lncRNA (x)", "Intronic lncRNA (i)", "Sense lncRNA (o/e)"))
    tab_mean_TAU_mod$Confidence = factor(tab_mean_TAU_mod$Confidence, levels = c("Protein-coding", "Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
    
    tab_mean_TAU_mod$"SRA.Study" = SRA.Study
    tab_mean_TAU_mod$"Specie" = specie
    tab_mean_TAU_mod = tab_mean_TAU_mod[,c("Specie", "SRA.Study", "ID_transcript", "Confidence", "Class_code", "TAU")]
    TAB_mean_TAU_mod = rbind(TAB_mean_TAU_mod, tab_mean_TAU_mod)
    
    gg = ggplot(tab_mean_TAU_mod, aes(x = Confidence, y = TAU, fill = Confidence)) + 
      geom_violin(color = "grey", position = position_dodge(1)) +
      geom_jitter(width = 0.4, color = "#8c8f8f", size = 0.1, alpha = 0.3) +
      stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
      geom_boxplot(width = 0.1, color = "black", alpha = 0.2, outlier.shape = NA, position = position_dodge(1)) +
      scale_y_continuous(limits = c(-0.001,1.001), breaks = seq(0,1,0.1), expand = c(0.01, 0.01)) + 
      facet_grid(. ~ Class_code) +
      theme_classic() + 
      scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee")) +
      scale_x_discrete(labels=c("Protein-coding" = "Protein-coding", 
                                "Low-confidence lncRNA" = "Low-confidence\nlncRNA", 
                                "Medium-confidence lncRNA" = "Medium-confidence\nlncRNA",
                                "High-confidence lncRNA" = "High-confidence\nlncRNA")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            panel.border = element_blank(),
            legend.position="none",
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 16, angle = 45, vjust = 0.5, hjust = 0.5),
            axis.text.y = element_text(size = 16),
            strip.text = element_text(size = 16)) +
      xlab("") + ylab("tau")
    
    ggsave(paste0(path_tissue_specificity, "/ALL/", flag, "/STEP3/Studies/", SRA.Study, "_mean-TAU-Violin_by_class_code.png"), height = 10, width = 20, dpi = 600)
    ggsave(paste0(path_tissue_specificity, "/ALL/", flag, "/STEP3/Studies/", SRA.Study, "_mean-TAU-Violin_by_class_code.pdf"), height = 10, width = 20, dpi = 600)
    
    rm(list = c("gg", "c", "pc", "tab_mean_TAU_mod"))
  }
}

rm(list = c("files", "specie", "SRA.Study", "SRA.Studies"))

write.table(TAB_mean_TAU_mod, paste0(path_tissue_specificity, "/ALL/", flag, "/STEP3/mean-TAU.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
