
######### MODULES
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))


######### VARIABLES

WD_pred = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs"
WD_corr = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/17-Correlation_DEFINITIVE"
species_short = c("cme", "csa", "cla")
species_long = c("C. melo", "C. sativus", "C. lanatus")


######### PIPELINE

if (!dir.exists(paste0(WD_corr, "/ALL"))){
  dir.create(paste0(WD_corr, "/ALL"))
}

## STEP 7 (Closest)

cat(paste0("-Closest...\n"))


################################################################################
## Load the lncRNA and PCG database.
cat(paste0("\t-Load the lncRNA database...\n"))

DB_lnc = data.frame()

for (i in 1:length(species_short)) {
  
  spes = species_short[i]
  spel = species_long[i]
  
  cat(paste0("\t\t-Specie: ", spes, "...\n"))
  
  # Load lncRNA database.
  db_lnc = read.table(paste0(WD_pred, "/", spes, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), sep = "\t", header = T, quote = "\"")
  db_lnc = db_lnc[, c("ID_transcript", "Significance_level")]
  colnames(db_lnc) = c("ID_transcript", "Confidence")
  db_lnc$Specie = spel
  db_lnc = db_lnc[, c("Specie", "ID_transcript", "Confidence")]
  
  DB_lnc = rbind(DB_lnc, db_lnc)
  
  rm(list = c("spes", "spel", "db_lnc"))
}

rm(list = c("i"))

################################################################################
## Load correlation tables with all the information about:
##    - Intergenic (TAB_u_corr) --> LncRNAs-PCGs and PCGs-PCGs.
##    - Antisense (TAB_x_corr) --> LncRNAs-PCGs.
cat(paste0("\t-Load correlation tables about intergenic and antisense lncRNAs...\n"))

TAB_u_corr = data.frame()
TAB_x_corr = data.frame()

for (i in 1:length(species_short)) {
  
  spes = species_short[i]
  spel = species_long[i]
  
  cat(paste0("\t\t-Specie: ", spes, "...\n"))
  
  # Load lncRNAs database.
  db_lnc = read.table(paste0(WD_pred, "/", spes, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), sep = "\t", header = T, quote = "\"")
  db_lnc = db_lnc[, c("ID_transcript", "Significance_level")]
  db_lnc_1 = db_lnc
  colnames(db_lnc_1) = c("ID_transcript.1", "Confidence.1")
  db_lnc_2 = db_lnc
  colnames(db_lnc_2) = c("ID_transcript.2", "Confidence.2")
  
  # Load the table closest from intergenic folder.
  tab_u = read.table(paste0(WD_corr, "/", spes, "/intergenic/nr/STEP4/TAB_CIS-PEARSON_closest.tsv"), header = T, sep = "\t", quote = "\"")
  tab_u = merge(tab_u, db_lnc_1, by = "ID_transcript.1", all.x = T, all.y = F)
  tab_u = merge(tab_u, db_lnc_2, by = "ID_transcript.2", all.x = T, all.y = F)
  tab_u$Confidence.1 = ifelse(is.na(tab_u$Confidence.1), "PC gene", tab_u$Confidence.1)
  tab_u$Confidence.2 = ifelse(is.na(tab_u$Confidence.2), "PC gene", tab_u$Confidence.2)
  tab_u = tab_u[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "Confidence.1", "ID_transcript.2", "Start.2", 
                    "End.2", "Strand.2", "Confidence.2", "Distance", "Type", "Experiment", "Cor.Deseq.Vst")]
  TAB_u_corr = rbind(TAB_u_corr, tab_u)
  
  rm(list = c("tab_u"))
  
  # Load the table closest from antisense folder.
  tab_x = read.table(paste0(WD_corr, "/", spes, "/antisense/nr/STEP4/TAB_CIS-PEARSON.tsv"), header = T, sep = "\t", quote = "\"")
  tab_x = merge(tab_x, db_lnc_1, by = "ID_transcript.1", all.x = T, all.y = F)
  tab_x = merge(tab_x, db_lnc_2, by = "ID_transcript.2", all.x = T, all.y = F)
  tab_x$Confidence.1 = ifelse(is.na(tab_x$Confidence.1), "PC gene", tab_x$Confidence.1)
  tab_x$Confidence.2 = ifelse(is.na(tab_x$Confidence.2), "PC gene", tab_x$Confidence.2)
  tab_x = tab_x[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "Confidence.1", "ID_transcript.2", "Start.2", 
                    "End.2", "Strand.2", "Confidence.2", "Overlap.Perc", "Overlap", "Type", "Experiment", "Cor.Deseq.Vst")]
  TAB_x_corr = rbind(TAB_x_corr, tab_x)
  
  rm(list = c("tab_x"))
  
  rm(list = c("spes", "spel", "db_lnc", "db_lnc_1", "db_lnc_2"))
}

rm(list = c("i"))


################################################################################
## Load random pairs.
cat(paste0("\t-Load random pairs...\n"))

TAB_u_corr_random = data.frame()
TAB_x_corr_random = data.frame()

for (i in 1:length(species_short)) {
  
  spes = species_short[i]
  spel = species_long[i]
  
  cat(paste0("\t\t-Specie: ", spes, "...\n"))
  
  # Load the table from intergenic folder.
  tab_u_random = read.table(paste0(WD_corr, "/", spes, "/intergenic/nr/STEP6/TAB_CIS-PEARSON_closest-random.tsv"), header = T, sep = "\t", quote = "\"")
  tab_u_random = merge(tab_u_random, DB_lnc, by.x = c("Specie", "ID_transcript.1"), by.y = c("Specie", "ID_transcript"), all.x = T, all.y = F)
  names(tab_u_random)[names(tab_u_random) == "Confidence"] = "Confidence.1"
  tab_u_random$Confidence.1 = ifelse(is.na(tab_u_random$Confidence.1), "PC gene", tab_u_random$Confidence.1)
  tab_u_random = merge(tab_u_random, DB_lnc, by.x = c("Specie", "ID_transcript.2"), by.y = c("Specie", "ID_transcript"), all.x = T, all.y = F)
  names(tab_u_random)[names(tab_u_random) == "Confidence"] = "Confidence.2"
  tab_u_random$Confidence.2 = ifelse(is.na(tab_u_random$Confidence.2), "PC gene", tab_u_random$Confidence.2)
  tab_u_random = tab_u_random[, c("Specie", "ID_transcript.1", "Confidence.1", "ID_transcript.2", "Confidence.2", "Type", "Experiment", "Cor.Deseq.Vst")]
  TAB_u_corr_random = rbind(TAB_u_corr_random, tab_u_random)
  
  rm(list = c("tab_u_random"))
  
  # Load the table from antisense folder.
  tab_x_random = read.table(paste0(WD_corr, "/", spes, "/antisense/nr/STEP6/TAB_CIS-PEARSON-random.tsv"), header = T, sep = "\t", quote = "\"")
  tab_x_random = merge(tab_x_random, DB_lnc, by.x = c("Specie", "ID_transcript.1"), by.y = c("Specie", "ID_transcript"), all.x = T, all.y = F)
  names(tab_x_random)[names(tab_x_random) == "Confidence"] = "Confidence.1"
  tab_x_random$Confidence.1 = ifelse(is.na(tab_x_random$Confidence.1), "PC gene", tab_x_random$Confidence.1)
  tab_x_random = merge(tab_x_random, DB_lnc, by.x = c("Specie", "ID_transcript.2"), by.y = c("Specie", "ID_transcript"), all.x = T, all.y = F)
  names(tab_x_random)[names(tab_x_random) == "Confidence"] = "Confidence.2"
  tab_x_random$Confidence.2 = ifelse(is.na(tab_x_random$Confidence.2), "PC gene", tab_x_random$Confidence.2)
  tab_x_random = tab_x_random[, c("Specie", "ID_transcript.1", "Confidence.1", "ID_transcript.2", "Confidence.2", "Type", "Experiment", "Cor.Deseq.Vst")]
  TAB_x_corr_random = rbind(TAB_x_corr_random, tab_x_random)
  
  rm(list = c("tab_x_random", "spes", "spel"))
}

rm(list = c("i"))

################################################################################
## Reduce the number of columns.
cat(paste0("\t-Reduce the number of columns...\n"))

TAB_u_corr_red = TAB_u_corr[TAB_u_corr$Confidence.1 %in% c("High", "PC gene"), c("Specie", "ID_transcript.1", "Confidence.1", "ID_transcript.2", "Confidence.2", "Type", "Experiment", "Cor.Deseq.Vst")]
TAB_u_corr_random_red = TAB_u_corr_random[(TAB_u_corr_random$Confidence.1 %in% c("High", "PC gene")) & (TAB_u_corr_random$Confidence.2 %in% c("High", "PC gene")), ]
TAB_u_corr_red_all = rbind(TAB_u_corr_red, TAB_u_corr_random_red)
TAB_u_corr_red_all$Specie = factor(TAB_u_corr_red_all$Specie, levels = species_long)
TAB_u_corr_red_all$Type = factor(TAB_u_corr_red_all$Type, levels = c("LncRNAs_Genes", "Genes_Genes", "Random pairs"))

TAB_x_corr_red = TAB_x_corr[TAB_x_corr$Confidence.1 %in% c("High"), c("Specie", "ID_transcript.1", "Confidence.1", "ID_transcript.2", "Confidence.2", "Type", "Experiment", "Cor.Deseq.Vst")]
TAB_x_corr_random_red = TAB_x_corr_random[(TAB_x_corr_random$Confidence.1 %in% c("High", "PC gene")) & (TAB_x_corr_random$Confidence.2 %in% c("High", "PC gene")), ]
TAB_x_corr_red_all = rbind(TAB_x_corr_red, TAB_x_corr_random_red)
TAB_x_corr_red_all$Specie = factor(TAB_x_corr_red_all$Specie, levels = species_long)
TAB_x_corr_red_all$Type = factor(TAB_x_corr_red_all$Type, levels = c("LncRNAs_Genes", "Random pairs"))

################################################################################
## Mean.
TAB_u_corr_red_all_sum = TAB_u_corr_red_all %>% group_by(Specie, Type) %>% summarize(mean = mean(Cor.Deseq.Vst))
TAB_x_corr_red_all_sum = TAB_x_corr_red_all %>% group_by(Specie, Type) %>% summarize(mean = mean(Cor.Deseq.Vst))

################################################################################
## Global correlation DENSITY.
cat(paste0("\t-Figure (DENSITY)...\n"))

ggU = ggplot(TAB_u_corr_red_all, aes(x = Cor.Deseq.Vst, color = Type)) +
  geom_density() +
  geom_vline(data = TAB_u_corr_red_all_sum, aes(xintercept = mean, color = Type), linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("#e5dd6c", "#7cc1cf", "#C1BFBE")) + 
  facet_grid(Specie~., scales = "free_y") +
  xlab("correlation coefficient (lincRNAs)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top",
        legend.title = element_blank())

ggX = ggplot(TAB_x_corr_red_all, aes(x = Cor.Deseq.Vst, color = Type)) +
  geom_density() +
  geom_vline(data = TAB_x_corr_red_all_sum, aes(xintercept = mean, color = Type), linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("#e1ad60", "#C1BFBE")) + 
  facet_grid(Specie~., scales = "free_y") +
  xlab("correlation coefficient (NAT-lncRNAs)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top",
        legend.title = element_blank())

gg_final = ggarrange(ggU,
                     ggX,
                     ncol = 2,
                     common.legend = F,
                     legend = "top")

ggsave(paste0(WD_corr, "/ALL/Global_correlation-DENSITY.png"), height = 8, width = 12, dpi = 600)
ggsave(paste0(WD_corr, "/ALL/Global_correlation-DENSITY.pdf"), height = 8, width = 12, dpi = 600)

rm(list = c("ggU", "ggX", "gg_final"))

################################################################################
## Global correlation BOXPLOT.
cat(paste0("\t-Figure (BOXPLOT)...\n"))

ggU = ggplot(TAB_u_corr_red_all, aes(x = Type, y = Cor.Deseq.Vst, fill = Type)) +
  geom_boxplot() +
  #stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#e5dd6c", "#7cc1cf", "#C1BFBE")) + 
  facet_grid(Specie~., scales = "free_y") +
  ylab("correlation coefficient (lincRNAs)") + xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top",
        legend.title = element_blank())

ggX = ggplot(TAB_x_corr_red_all, aes(x = Type, y = Cor.Deseq.Vst, fill = Type)) +
  geom_boxplot() +
  #stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#e1ad60", "#C1BFBE")) + 
  facet_grid(Specie~., scales = "free_y") +
  ylab("correlation coefficient (NAT-lncRNAs)") + xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top",
        legend.title = element_blank())

gg_final = ggarrange(ggU,
                     ggX,
                     ncol = 2,
                     common.legend = F,
                     legend = "top")

ggsave(paste0(WD_corr, "/ALL/Global_correlation-BOXPLOT.png"), height = 8, width = 12, dpi = 600)
ggsave(paste0(WD_corr, "/ALL/Global_correlation-BOXPLOT.pdf"), height = 8, width = 12, dpi = 600)

rm(list = c("ggU", "ggX", "gg_final"))

################################################################################
## Global correlation VIOLIN.
cat(paste0("\t-Figure (VIOLIN)...\n"))

ggU = ggplot(TAB_u_corr_red_all, aes(x = Type, y = Cor.Deseq.Vst, fill = Type)) +
  geom_violin() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#e5dd6c", "#7cc1cf", "#C1BFBE")) + 
  facet_grid(Specie~., scales = "free_y") +
  ylab("correlation coefficient (lincRNAs)") + xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top",
        legend.title = element_blank())

ggX = ggplot(TAB_x_corr_red_all, aes(x = Type, y = Cor.Deseq.Vst, fill = Type)) +
  geom_violin() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#e1ad60", "#C1BFBE")) + 
  facet_grid(Specie~., scales = "free_y") +
  ylab("correlation coefficient (NAT-lncRNAs)") + xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top",
        legend.title = element_blank())

gg_final = ggarrange(ggU,
                     ggX,
                     ncol = 2,
                     common.legend = F,
                     legend = "top")

ggsave(paste0(WD_corr, "/ALL/Global_correlation-VIOLIN.png"), height = 8, width = 12, dpi = 600)
ggsave(paste0(WD_corr, "/ALL/Global_correlation-VIOLIN.pdf"), height = 8, width = 12, dpi = 600)

rm(list = c("ggU", "ggX", "gg_final"))

################################################################################
## Intergenic: Statistical analysis.
cat(paste0("\t-Statistics (Intergenic)...\n"))

combinations = as.data.frame(t(combn(c("LncRNAs_Genes", "Genes_Genes", "Random pairs"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("cl1", "cl2")

TAB_u_STATISTICS = data.frame()

for (spe in c("C. lanatus", "C. melo", "C. sativus")) {
  for (i in 1:nrow(combinations)) {
    cl1 = combinations[i, "cl1"]
    cl2 = combinations[i, "cl2"]
    subset1 = TAB_u_corr_red_all[TAB_u_corr_red_all$Specie == spe & TAB_u_corr_red_all$Type == cl1 & !is.na(TAB_u_corr_red_all$Cor.Deseq.Vst), "Cor.Deseq.Vst"]
    subset2 = TAB_u_corr_red_all[TAB_u_corr_red_all$Specie == spe & TAB_u_corr_red_all$Type == cl2 & !is.na(TAB_u_corr_red_all$Cor.Deseq.Vst), "Cor.Deseq.Vst"]
    L1 = length(subset1)
    L2 = length(subset2)
    if (L1 > 0 & L2 > 0) {
      test = wilcox.test(subset1, subset2, paired = FALSE)
      Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
      row = data.frame(Specie = spe, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
    } 
    else{
      row = data.frame(Specie = spe, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
    }
    TAB_u_STATISTICS = rbind(TAB_u_STATISTICS, row)
  }
}

rownames(TAB_u_STATISTICS) = NULL

write.table(TAB_u_STATISTICS, paste0(WD_corr, "/ALL/Global_correlation-STATISTICS-U.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("spe", "i", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

## Antisense: Statistical analysis.
cat(paste0("\t-Statistics (Antisense)...\n"))

combinations = as.data.frame(t(combn(c("LncRNAs_Genes", "Random pairs"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("cl1", "cl2")

TAB_x_STATISTICS = data.frame()

for (spe in c("C. lanatus", "C. melo", "C. sativus")) {
  for (i in 1:nrow(combinations)) {
    cl1 = combinations[i, "cl1"]
    cl2 = combinations[i, "cl2"]
    subset1 = TAB_x_corr_red_all[TAB_x_corr_red_all$Specie == spe & TAB_x_corr_red_all$Type == cl1 & !is.na(TAB_x_corr_red_all$Cor.Deseq.Vst), "Cor.Deseq.Vst"]
    subset2 = TAB_x_corr_red_all[TAB_x_corr_red_all$Specie == spe & TAB_x_corr_red_all$Type == cl2 & !is.na(TAB_x_corr_red_all$Cor.Deseq.Vst), "Cor.Deseq.Vst"]
    L1 = length(subset1)
    L2 = length(subset2)
    if (L1 > 0 & L2 > 0) {
      test = wilcox.test(subset1, subset2, paired = FALSE)
      Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
      row = data.frame(Specie = spe, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
    } 
    else{
      row = data.frame(Specie = spe, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
    }
    TAB_x_STATISTICS = rbind(TAB_x_STATISTICS, row)
  }
}

rownames(TAB_x_STATISTICS) = NULL

write.table(TAB_x_STATISTICS, paste0(WD_corr, "/ALL/Global_correlation-STATISTICS-X.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("spe", "i", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

rm(list = c("TAB_u_corr_random", "TAB_x_corr_random", "TAB_u_corr_random_red", "TAB_x_corr_random_red", "TAB_u_corr_red", "TAB_x_corr_red", 
            "TAB_u_corr_red_all_sum", "TAB_x_corr_red_all_sum", "TAB_u_corr_red_all", "TAB_x_corr_red_all", "TAB_u_STATISTICS", "TAB_x_STATISTICS"))










################################################################################
################################################################################
################################################################################











# ################################################################################
# ## Load interaction tables with all the information about:
# ##    - Intergenic (TAB_u_int) --> LncRNAs-PCGs.
# ##    - Antisense (TAB_x_int) --> LncRNAs-PCGs.
# cat(paste0("\t-Load interaction tables about intergenic and antisense lncRNAs...\n"))
# 
# TAB_u_int = data.frame()
# TAB_x_int = data.frame()
# 
# for (i in 1:length(species_short)) {
#   
#   spes = species_short[i]
#   spel = species_long[i]
#   
#   cat(paste0("\t\t-Specie: ", spes, "...\n"))
#   
#   # Load lncRNAs database.
#   db_lnc = read.table(paste0(WD_pred, "/", spes, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), sep = "\t", header = T, quote = "\"")
#   db_lnc = db_lnc[, c("ID_transcript", "Significance_level")]
#   db_lnc_1 = db_lnc
#   colnames(db_lnc_1) = c("ID_transcript.1", "Confidence.1")
#   db_lnc_2 = db_lnc
#   colnames(db_lnc_2) = c("ID_transcript.2", "Confidence.2")
#   
#   # Load the table closest from intergenic folder.
#   tab_u = read.table(paste0(WD_corr, "/", spes, "/intergenic/nr/STEP2/TAB_CIS_closest.tsv"), header = T, sep = "\t", quote = "\"")
#   tab_u = merge(tab_u, db_lnc_1, by = "ID_transcript.1", all.x = T, all.y = F)
#   tab_u = merge(tab_u, db_lnc_2, by = "ID_transcript.2", all.x = T, all.y = F)
#   tab_u$Confidence.1 = ifelse(is.na(tab_u$Confidence.1), "PC gene", tab_u$Confidence.1)
#   tab_u$Confidence.2 = ifelse(is.na(tab_u$Confidence.2), "PC gene", tab_u$Confidence.2)
#   tab_u = tab_u[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "Confidence.1", "ID_transcript.2", "Start.2", 
#                     "End.2", "Strand.2", "Confidence.2", "Distance", "Type")]
#   TAB_u_int = rbind(TAB_u_int, tab_u)
#   
#   rm(list = c("tab_u"))
#   
#   # Load the table closest from antisense folder.
#   tab_x = read.table(paste0(WD_corr, "/", spes, "/antisense/nr/STEP2/TAB_CIS.tsv"), header = T, sep = "\t", quote = "\"")
#   tab_x = merge(tab_x, db_lnc_1, by = "ID_transcript.1", all.x = T, all.y = F)
#   tab_x = merge(tab_x, db_lnc_2, by = "ID_transcript.2", all.x = T, all.y = F)
#   tab_x$Confidence.1 = ifelse(is.na(tab_x$Confidence.1), "PC gene", tab_x$Confidence.1)
#   tab_x$Confidence.2 = ifelse(is.na(tab_x$Confidence.2), "PC gene", tab_x$Confidence.2)
#   tab_x = tab_x[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "Confidence.1", "ID_transcript.2", "Start.2", 
#                     "End.2", "Strand.2", "Confidence.2", "Overlap.Perc", "Overlap", "Type")]
#   TAB_x_int = rbind(TAB_x_int, tab_x)
#   
#   rm(list = c("tab_x"))
#   
#   rm(list = c("spes", "spel", "db_lnc", "db_lnc_1", "db_lnc_2"))
# }
# 
# rm(list = c("i"))
# 
# ################################################################################
# ## Correlation tables with LncRNAs-PCGs only (Intergenic and Antisense lncRNAs)
# cat(paste0("\t-Filter correlation tables and keep only lncRNAs-Genes...\n"))
# 
# TAB_u_corr_lnc = TAB_u_corr[TAB_u_corr$Type == "LncRNAs_Genes",]
# TAB_u_int_lnc = TAB_u_int[TAB_u_int$Type == "LncRNAs_Genes",]
# 
# 
# 
# TAB_u_corr_lnc$"Class.1" = "lincRNA"
# TAB_u_corr_lnc = TAB_u_corr_lnc[, c("Specie", "ID_transcript.1", "Confidence.1", "Class.1", "ID_transcript.2", "Confidence.2", "Distance", "Experiment", "Cor.Deseq.Vst")]
# TAB_x_corr_lnc = TAB_x_corr
# TAB_x_corr_lnc$"Class.1" = "NAT-lncRNA" 
# TAB_x_corr_lnc = TAB_x_corr_lnc[, c("Specie", "ID_transcript.1", "Confidence.1", "Class.1", "ID_transcript.2", "Confidence.2", "Overlap.Perc", "Overlap", "Experiment", "Cor.Deseq.Vst")]










################################################################################
########################## ANTISENSE CLASSIFICATION ############################
################################################################################











################################################################################
## Load interaction tables with all the information about:
##    - Antisense (TAB_x_int) --> LncRNAs-PCGs.
cat(paste0("\t-Load interaction tables about intergenic and antisense lncRNAs...\n"))

TAB_x_int = data.frame()

for (i in 1:length(species_short)) {
  
  spes = species_short[i]
  spel = species_long[i]
  
  cat(paste0("\t\t-Specie: ", spes, "...\n"))
  
  # Load lncRNAs database.
  db_lnc = read.table(paste0(WD_pred, "/", spes, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), sep = "\t", header = T, quote = "\"")
  db_lnc = db_lnc[, c("ID_transcript", "Significance_level")]
  db_lnc_1 = db_lnc
  colnames(db_lnc_1) = c("ID_transcript.1", "Confidence.1")
  db_lnc_2 = db_lnc
  colnames(db_lnc_2) = c("ID_transcript.2", "Confidence.2")
  
  # Load the table closest from antisense folder.
  tab_x = read.table(paste0(WD_corr, "/", spes, "/antisense/nr/STEP2/TAB_CIS.tsv"), header = T, sep = "\t", quote = "\"")
  tab_x = merge(tab_x, db_lnc_1, by = "ID_transcript.1", all.x = T, all.y = F)
  tab_x = merge(tab_x, db_lnc_2, by = "ID_transcript.2", all.x = T, all.y = F)
  tab_x$Confidence.1 = ifelse(is.na(tab_x$Confidence.1), "PC gene", tab_x$Confidence.1)
  tab_x$Confidence.2 = ifelse(is.na(tab_x$Confidence.2), "PC gene", tab_x$Confidence.2)
  tab_x = tab_x[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "Confidence.1", "ID_transcript.2", "Start.2", 
                    "End.2", "Strand.2", "Confidence.2", "Overlap.Perc", "Overlap", "Type")]
  TAB_x_int = rbind(TAB_x_int, tab_x)
  
  rm(list = c("tab_x"))
  
  rm(list = c("spes", "spel", "db_lnc", "db_lnc_1", "db_lnc_2"))
}

rm(list = c("i"))

################################################################################
## Classify the NAT-lncRNAs
cat(paste0("\t-Classify the NAT-lncRNAs...\n"))

TAB_x_int$Class_x = NA
TAB_x_int$Class_x = ifelse(TAB_x_int$Overlap.Perc > 0 & 
                                  (TAB_x_int$Start.1 < TAB_x_int$Start.2) & (TAB_x_int$End.1 < TAB_x_int$End.2), "HEAD-TO-HEAD", TAB_x_int$Class_x)
TAB_x_int$Class_x = ifelse(TAB_x_int$Overlap.Perc > 0 & 
                                  (TAB_x_int$Start.1 > TAB_x_int$Start.2) & (TAB_x_int$End.1 > TAB_x_int$End.2), "TAIL-TO-TAIL", TAB_x_int$Class_x)
TAB_x_int$Class_x = ifelse(TAB_x_int$Overlap.Perc > 0 & 
                                  ((TAB_x_int$Start.1 == TAB_x_int$Start.2) & (TAB_x_int$End.1 == TAB_x_int$End.2) |
                                  (TAB_x_int$Start.1 > TAB_x_int$Start.2) & (TAB_x_int$End.1 < TAB_x_int$End.2) |
                                  (TAB_x_int$Start.1 == TAB_x_int$Start.2) & (TAB_x_int$End.1 < TAB_x_int$End.2) |
                                  (TAB_x_int$Start.1 > TAB_x_int$Start.2) & (TAB_x_int$End.1 == TAB_x_int$End.2)), "EMBEDDED", TAB_x_int$Class_x)
TAB_x_int$Class_x = ifelse(is.na(TAB_x_int$Class_x), "OTHER", TAB_x_int$Class_x)

################################################################################
## Filter
cat(paste0("\t-Filter...\n"))

TAB_x_int_filt = TAB_x_int[TAB_x_int$Confidence.1 == "High",]

################################################################################
## Add correlation info
cat(paste0("\t-Add correlation info...\n"))

X = merge(TAB_x_corr[TAB_x_corr$Confidence.1 == "High",], TAB_x_int_filt[, c("Specie", "ID_transcript.1", "ID_transcript.2", "Class_x")], by = c("Specie", "ID_transcript.1", "ID_transcript.2"), all.x = T, all.y = F)

################################################################################
## Global correlation BOXPLOT.
cat(paste0("\t-Figure (BOXPLOT)...\n"))

ggX = ggplot(X, aes(x = Class_x, y = Cor.Deseq.Vst, fill = Class_x)) +
  geom_boxplot() +
  #stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#A6D5A2", "#75B3C3", "#e1ad60", "#C1BFBE")) + 
  facet_grid(Specie~., scales = "free_y") +
  ylab("correlation coefficient (NAT-lncRNAs)") + xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top",
        legend.title = element_blank())

ggsave(paste0(WD_corr, "/ALL/Global_correlation-BOXPLOT-ANTISENSE_CLASSES.png"), height = 10, width = 12, dpi = 600)
ggsave(paste0(WD_corr, "/ALL/Global_correlation-BOXPLOT-ANTISENSE_CLASSES.pdf"), height = 10, width = 12, dpi = 600)

rm(list = c("ggX"))

