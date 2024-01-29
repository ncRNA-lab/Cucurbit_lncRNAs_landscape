


################################################################################
################################################################################
################################### MODULES ####################################
################################################################################
################################################################################


suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(DESeq2))
suppressMessages(library(plotly))
suppressMessages(library(htmltools))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tibble))
# devtools::install_github("azhu513/apeglm")
suppressMessages(library(scales))
suppressMessages(library(ggExtra))
suppressMessages(library(data.table))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(limma))
suppressMessages(library(ggvenn))
suppressMessages(library(ggpubr))
library(palmerpenguins)
library(ggtext)
library(colorspace)
library(ragg)
options(bitmapType='cairo')










################################################################################
################################################################################
################################## VARIABLES ###################################
################################################################################
################################################################################


path_DEA = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/16-DEA"
path_PosCon = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/nr/05-Figures/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"
path_TS = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity/approach_1/ALL/nr/STEP3/mean-TAU.tsv"
alpha_value = 0.05
species = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
species_long = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")










################################################################################
################################################################################
################################## FIGURES #####################################
################################################################################
################################################################################

cat(paste0("Figures..."))

if (!dir.exists(paste0(path_DEA, "/ALL"))){
  dir.create(paste0(path_DEA, "/ALL"))
}










###################################

if (!dir.exists(paste0(path_DEA, "/ALL/Events"))){
  dir.create(paste0(path_DEA, "/ALL/Events"))
}

TAB = data.frame()
for (spe in species) {
  if (spe %in% list.files(path_DEA)) {
    path = paste0(path_DEA, "/", spe, "/04-DEA")
    tab = read.table(paste0(path, "/SUMMARY.tsv"), sep = "\t", header = T, quote = "\"")
    tab$"ID" = paste0(tab$Experiment, ".", tab$n_contrast)
    tab$Specie = spe
    tab[tab == "Abiotic stress"] = "Stress"
    tab[tab == "Biotic stress"] = "Stress"
    TAB = rbind(TAB, tab)
  }
}

TAB$Specie = factor(TAB$Specie, levels = species)
TAB$Class = factor(TAB$Class, levels = c("Development", "Stress"))

TAB_collap = TAB %>% 
  group_by(Specie, Class, .drop = F) %>% 
  summarise(Number_Exp_cont = n_distinct(ID), 
            Number_Projects = n_distinct(SRA_study))

TAB_collap = as.data.frame(TAB_collap)

write.table(TAB_collap, paste0(path_DEA, "/ALL/Events/Events.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

TAB_collap_filt = TAB_collap[TAB_collap$Number_Exp_cont != 0,]
species_mod = as.character(TAB_collap_filt[duplicated(TAB_collap_filt$Specie), "Specie"])

rm(list = c("TAB", "TAB_collap", "TAB_collap_filt", "tab", "path", "spe"))

## Additional table (PAPER)
TAB_add = data.frame()
for (spe in species) {
  if (spe %in% list.files(path_DEA)) {
    path = paste0(path_DEA, "/", spe, "/04-DEA")
    tab = read.table(paste0(path, "/SUMMARY.tsv"), sep = "\t", header = T, quote = "\"")
    tab$Specie = species_long[grepl(spe, species, fixed = T)]
    tab$"ID" = paste0(tab$Experiment, ".", tab$n_contrast)
    tab[tab == "Abiotic stress"] = "Stress"
    tab[tab == "Biotic stress"] = "Stress"
    tab = tab[, c("Specie", "Class", "SRA_study", "Stress", "ID")]
    TAB_add = rbind(TAB_add, tab)
  }
}

TAB_add_dev = TAB_add[TAB_add$Class == "Development",]
TAB_add_st = TAB_add[TAB_add$Class == "Stress",]

TAB_add_dev$Specie = factor(TAB_add_dev$Specie, levels = species_long)
TAB_add_st$Specie = factor(TAB_add_st$Specie, levels = species_long)
TAB_add_dev$Class = factor(TAB_add_dev$Class, levels = c("Development"))
TAB_add_st$Class = factor(TAB_add_st$Class, levels = c("Stress"))

TAB_add_dev_collap = TAB_add_dev %>%
  group_by(Specie, Class, SRA_study, .drop = F) %>%
  summarise(Number.Events = n_distinct(ID))
TAB_add_dev_collap = as.data.frame(TAB_add_dev_collap)
TAB_add_dev_collap = TAB_add_dev_collap[TAB_add_dev_collap$Number.Events != 0,]

write.table(TAB_add_dev_collap, paste0(path_DEA, "/ALL/Events/Events_dev.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

TAB_add_st_collap = TAB_add_st %>%
  group_by(Specie, Class, SRA_study, Stress, .drop = F) %>%
  summarise(Number.Events = n_distinct(ID))
TAB_add_st_collap = as.data.frame(TAB_add_st_collap)
TAB_add_st_collap = TAB_add_st_collap[TAB_add_st_collap$Number.Events != 0,]

write.table(TAB_add_st_collap, paste0(path_DEA, "/ALL/Events/Events_st.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("TAB_add", "TAB_add_dev", "TAB_add_st", "TAB_add_dev_collap", "TAB_add_st_collap", "tab", "path", "spe"))


















###################################

if (!dir.exists(paste0(path_DEA, "/ALL/Percentage_DE"))){
  dir.create(paste0(path_DEA, "/ALL/Percentage_DE"))
}

TAB_1 = data.frame()
TAB_2 = data.frame()

for (i in 1:length(species)) {
  
  spes = species[i]
  spel = species_long[i]
  
  if (spes %in% list.files(path_DEA)) {
    path = paste0(path_DEA, "/", spes, "/05-Tables_and_Figures/ALL")
    
    tab_1 = read.table(paste0(path, "/Percentage_DE/TAB-Violin-Percentage_DE-1.tsv"), sep = "\t", header = T, quote = "\"")
    tab_1$Specie = spes
    tab_1$Specie_long = spel
    tab_2 = read.table(paste0(path, "/Percentage_DE/TAB-Violin-Percentage_DE-2.tsv"), sep = "\t", header = T, quote = "\"")
    tab_2$Specie = spes
    tab_2$Specie_long = spel
    
    TAB_1 = rbind(TAB_1, tab_1)
    TAB_2 = rbind(TAB_2, tab_2)
  }
}

rm(list = c("spes", "spel", "i", "tab_1", "tab_2", "path"))

TAB_1[TAB_1 == "Abiotic stress"] = "Stress"
TAB_1[TAB_1 == "Biotic stress"] = "Stress"
TAB_2[TAB_2 == "Abiotic stress"] = "Stress"
TAB_2[TAB_2 == "Biotic stress"] = "Stress"

TAB_1 = TAB_1[TAB_1$Significance_level == "PC genes" | TAB_1$Significance_level == "HC-lncRNAs",]
rownames(TAB_1) = NULL
TAB_2 = TAB_2[TAB_2$Significance_level == "PC genes" | TAB_2$Significance_level == "HC-lncRNAs",]
rownames(TAB_2) = NULL

TAB_1$Class_code = factor(TAB_1$Class_code, levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
TAB_1$Class = factor(TAB_1$Class, levels = c("Development", "Stress"))
TAB_1$Specie = factor(TAB_1$Specie, levels = species)
TAB_2$Significance_level = factor(TAB_2$Significance_level, levels = c("PC genes", "HC-lncRNAs"))
TAB_2$Class = factor(TAB_2$Class, levels = c("Development", "Stress"))
TAB_2$Specie = factor(TAB_2$Specie, levels = species)

TAB_1 = TAB_1[order(TAB_1$Specie, TAB_1$ID, TAB_1$Class),]
rownames(TAB_1) = NULL
TAB_2 = TAB_2[order(TAB_2$Specie, TAB_2$ID, TAB_2$Class),]
rownames(TAB_2) = NULL

TAB_1_save = TAB_1[, c("Specie_long", "ID", "Class_code", "Number", "Total", "Percentage", "Class")]
colnames(TAB_1_save) = c("Specie", "Experiment", "Type", "Number", "Total", "Percentage", "Class")
TAB_1_save_dev = TAB_1_save[TAB_1_save$Class == "Development", c("Specie", "Experiment", "Type", "Number", "Total", "Percentage")]
write.table(TAB_1_save_dev, paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-1-Dev.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
TAB_1_save_str = TAB_1_save[TAB_1_save$Class == "Stress", c("Specie", "Experiment", "Type", "Number", "Total", "Percentage")]
write.table(TAB_1_save_str, paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-1-Str.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

TAB_2_save = TAB_2[, c("Specie_long", "ID", "Significance_level", "Number", "Total", "Percentage", "Class")]
colnames(TAB_2_save) = c("Specie", "Experiment", "Type", "Number", "Total", "Percentage", "Class")
TAB_2_save_dev = TAB_2_save[TAB_2_save$Class == "Development", c("Specie", "Experiment", "Type", "Number", "Total", "Percentage")]
write.table(TAB_2_save_dev, paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-2-Dev.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
TAB_2_save_str = TAB_2_save[TAB_2_save$Class == "Stress", c("Specie", "Experiment", "Type", "Number", "Total", "Percentage")]
write.table(TAB_2_save_str, paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-2-Str.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

TAB_1_collap = TAB_1 %>%
  group_by(Class_code, Class) %>%
  summarise(Mean.Percentage = round(mean(Percentage), 2))
TAB_2_collap = TAB_2 %>%
  group_by(Significance_level, Class) %>%
  summarise(Mean.Percentage = round(mean(Percentage), 2))

write.table(TAB_1_collap, paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-1-collap.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(TAB_2_collap, paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-2-collap.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

combinations = as.data.frame(t(combn(c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("cl1", "cl2")

TAB_STATISTICS_FINAL_1 = data.frame()
for (class in c("Development", "Stress")) {
  for (i in 1:nrow(combinations)) {
    cl1 = combinations[i, "cl1"]
    cl2 = combinations[i, "cl2"]
    subset1 = TAB_1[TAB_1$Class == class & TAB_1$Class_code == cl1, "Percentage"]
    subset2 = TAB_1[TAB_1$Class == class & TAB_1$Class_code == cl2, "Percentage"]
    L1 = length(subset1)
    L2 = length(subset2)
    if (L1 > 0 & L2 > 0) {
      test = wilcox.test(subset1, subset2, paired = TRUE)
      Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
      row = data.frame(Class = class, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
    }
    else{
      row = data.frame(Class = class, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
    }
    TAB_STATISTICS_FINAL_1 = rbind(TAB_STATISTICS_FINAL_1, row)
  }
}

rownames(TAB_STATISTICS_FINAL_1) = NULL

write.table(TAB_STATISTICS_FINAL_1, paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-1-STAT.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("class", "i", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))

combinations = as.data.frame(t(combn(c("PC genes", "HC-lncRNAs"), 2)))
rownames(combinations) = NULL
colnames(combinations) = c("cl1", "cl2")

TAB_STATISTICS_FINAL_2 = data.frame()
for (class in c("Development", "Stress")) {
  for (i in 1:nrow(combinations)) {
    cl1 = combinations[i, "cl1"]
    cl2 = combinations[i, "cl2"]
    subset1 = TAB_2[TAB_2$Class == class & TAB_2$Significance_level == cl1, "Percentage"]
    subset2 = TAB_2[TAB_2$Class == class & TAB_2$Significance_level == cl2, "Percentage"]
    L1 = length(subset1)
    L2 = length(subset2)
    if (L1 > 0 & L2 > 0) {
      test = wilcox.test(subset1, subset2, paired = TRUE)
      Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
      row = data.frame(Class = class, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
    }
    else{
      row = data.frame(Class = class, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
    }
    TAB_STATISTICS_FINAL_2 = rbind(TAB_STATISTICS_FINAL_2, row)
  }
}

rownames(TAB_STATISTICS_FINAL_2) = NULL

write.table(TAB_STATISTICS_FINAL_2, paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-2-STAT.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("class", "i", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))


for (cl in unique(TAB_1$Class)) {
  
  temp_1 = TAB_1[TAB_1$Class == cl,]
  
  for (i in 1:length(species)) {
    spes = species[i]
    spel = species_long[i]
    
    if (!(spes %in% temp_1$Specie)) {
      df = data.frame(Specie = spes, ID = "A", Class_code = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"),
                      Significance_level = c("PC genes", "HC-lncRNAs", "HC-lncRNAs", "HC-lncRNAs", "HC-lncRNAs"), Type = NA, Number = NA,
                      Total = NA, Percentage = NA, Class = NA, Specie_long = spel)
      temp_1 = rbind(temp_1, df)
    }
  }
  
  shapes = c(15,16,17,18,0,1,2,3,4)
  #shapes = c(1:length(species_long))
  shapes_assigned = setNames(shapes, species)
  temp_1 = temp_1[order(temp_1$Specie),]
  temp_1$Shape = shapes_assigned[temp_1$Specie]
  temp_1$Shape = factor(temp_1$Shape, levels = shapes)
  
  add_sample = function(x){return(c(y = max(x) + 1, label = mean(x)))}
  
  gg1 = ggplot(temp_1, aes(x = Class_code, y = Percentage)) + 
    ggdist::stat_halfeye(
      aes(color = Class_code,
          fill = after_scale(lighten(color, .5))),
      adjust = .5, 
      width = .75, 
      .width = 0,
      justification = -.4, 
      point_color = NA) + 
    geom_boxplot(
      aes(color = Class_code,
          color = after_scale(darken(color, .1, space = "HLS")),
          fill = after_scale(desaturate(lighten(color, .8), .4))),
      width = .42, 
      outlier.shape = NA
    ) +
    geom_point(
      aes(color = Class_code,
          color = after_scale(darken(color, .1, space = "HLS")),
          shape = Shape),
      fill = "white",
      stroke = .4,
      size = 2.4,
      position = position_jitter(seed = 1, width = .2)
    ) + 
    geom_point(
      aes(fill = Class_code, shape = Shape),
      color = "transparent",
      stroke = .4,
      size = 2.4,
      alpha = .3,
      position = position_jitter(seed = 1, width = .2)
    ) + 
    stat_summary(
      geom = "text",
      fun.data = add_sample,
      aes(label = paste("Mean =", round(after_stat(label), 2))),
      color = "black",
      family = "Roboto Condensed",
      fontface = "bold",
      size = 7,
      hjust = 0.4
    ) +
    xlab("") + ylab("") +
    scale_y_continuous(limits = c(-0.05, 32), breaks = seq(0, 30, 5), expand = c(0.03, 0.03)) +
    scale_color_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"),
      values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66"), 
      guide = "none"
    ) +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"),
      values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66"),
      guide = "none"
    ) +
    scale_shape_manual(values = shapes, labels = species_long) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 20, margin = margin(r = 0.9, l = -0.1, unit = 'cm'), face = "italic"),
          legend.text.align = 0,
          legend.key.size = unit(1, 'cm')) +
    guides(shape = guide_legend(title = element_blank(), nrow = 1, override.aes = list(size = 3.4))) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 22))
  
  temp_2 = TAB_2[TAB_2$Class == cl,]
  
  for (i in 1:length(species)) {
    spes = species[i]
    spel = species_long[i]
    
    if (!(spes %in% temp_2$Specie)) {
      df = data.frame(Specie = spes, ID = "A", Significance_level = c("PC genes", "HC-lncRNAs", "HC-lncRNAs", "HC-lncRNAs", "HC-lncRNAs"), 
                      Type = NA, Number = NA, Total = NA, Percentage = NA, Class = NA, Specie_long = spel)
      temp_2 = rbind(temp_2, df)
    }
  }
  
  shapes = c(15,16,17,18,0,1,2,3,4)
  #shapes = c(1:length(species_long))
  shapes_assigned = setNames(shapes, species)
  temp_2 = temp_2[order(temp_2$Specie),]
  temp_2$Shape = shapes_assigned[temp_2$Specie]
  temp_2$Shape = factor(temp_2$Shape, levels = shapes)
  
  add_sample = function(x){return(c(y = max(x) + 1, label = mean(x)))}
  
  gg2 = ggplot(temp_2, aes(x = Significance_level, y = Percentage)) + 
    ggdist::stat_halfeye(
      aes(color = Significance_level,
          fill = after_scale(lighten(color, .5))),
      adjust = .5, 
      width = .75, 
      .width = 0,
      justification = -.4, 
      point_color = NA) + 
    geom_boxplot(
      aes(color = Significance_level,
          color = after_scale(darken(color, .1, space = "HLS")),
          fill = after_scale(desaturate(lighten(color, .8), .4))),
      width = .42, 
      outlier.shape = NA
    ) +
    geom_point(
      aes(color = Significance_level,
          color = after_scale(darken(color, .1, space = "HLS")),
          shape = Shape),
      fill = "white",
      stroke = .4,
      size = 2.4,
      position = position_jitter(seed = 1, width = .2)
    ) + 
    geom_point(
      aes(fill = Significance_level, shape = Shape),
      color = "transparent",
      stroke = .4,
      size = 2.4,
      alpha = .3,
      position = position_jitter(seed = 1, width = .2)
    ) + 
    stat_summary(
      geom = "text",
      fun.data = add_sample,
      aes(label = paste("Mean =", round(after_stat(label), 2))),
      color = "black",
      family = "Roboto Condensed",
      fontface = "bold",
      size = 7,
      hjust = 0.4
    ) +
    xlab("") + ylab("Differentially expressed transcripts (%)") +
    scale_y_continuous(limits = c(-0.05, 32), breaks = seq(0, 30, 5), expand = c(0.03, 0.03)) +
    scale_color_manual(
      labels = c("PC genes", "HC-lncRNAs"),
      values = c("#7cc1cf", "#edd0ee"), 
      guide = "none"
    ) +
    scale_fill_manual(
      labels = c("PC genes", "HC-lncRNAs"),
      values = c("#7cc1cf", "#edd0ee"),
      guide = "none"
    ) +
    scale_shape_manual(values = shapes, labels = species_long) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 20, margin = margin(r = 0.9, l = -0.1, unit = 'cm'), face = "italic"),
          legend.text.align = 0,
          legend.key.size = unit(1, 'cm')) +
    guides(shape = guide_legend(title = element_blank(), nrow = 1, override.aes = list(size = 3.4))) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 22))
  
  # Join figures
  gg_final = ggarrange(gg2, 
                       ggparagraph(text="   ", face = "italic", size = 16, color = "black"), 
                       gg1,
                       ncol = 3, 
                       widths = c(0.38, 0.02, 0.6),
                       common.legend = T,
                       legend = "top")
  
  ggsave(paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-", cl, ".png"), height = 9, width = 20, dpi = 800, bg = "white")
  #ggsave(paste0(path_DEA, "/ALL/Percentage_DE/Violin-Percentage-", cl, ".pdf"), height = 9, width = 20, dpi = 800, bg = "white")
  
  rm(list = c("temp_1", "temp_2", "gg1", "gg2", "gg_final", "shapes", "shapes_assigned"))
}

rm(list = c("cl", "TAB_1", "TAB_2", "TAB_1_collap", "TAB_2_collap", "TAB_STATISTICS_FINAL_1", "TAB_STATISTICS_FINAL_2",
            "add_sample"))









###################################

if (!dir.exists(paste0(path_DEA, "/ALL/VennDiagram"))){
  dir.create(paste0(path_DEA, "/ALL/VennDiagram"))
}

TAB = data.frame()

for (i in 1:length(species_mod)) {
  
  spe = species_mod[i]
  path = paste0(path_DEA, "/", spe, "/05-Tables_and_Figures/ALL")
  
  tab = read.table(paste0(path, "/VennDiagram/TAB-Venn_filt.tsv"), sep = "\t", header = T, quote = "\"")
  tab$Specie = spe
  tab$ID_transcript_spe = paste0(tab$ID_transcript, "-", spe)
  
  TAB = rbind(TAB, tab)
}

for (sl in c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs")) {
  x = list(
    Development = unique(TAB[TAB$Significance_level == sl & TAB$Item == "Norm.Development", "ID_transcript_spe"]),
    Stress = unique(TAB[TAB$Significance_level == sl & TAB$Item == "Norm.Stress", "ID_transcript_spe"]),
    `Tissue specificity` = unique(TAB[TAB$Significance_level == sl & TAB$Item == "Tissue.Specificity", "ID_transcript_spe"])
  )
  ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"), stroke_size = 0.5, set_name_size = 5, text_size = 4.2)
  ggsave(paste0(path_DEA, "/ALL/VennDiagram/Venn-", gsub(" ", "-", sl), ".png"), height = 6, width = 6, dpi = 800, bg = "white")
}

rm(list = c("sl", "x"))


                  