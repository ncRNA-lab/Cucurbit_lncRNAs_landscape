################################################################################
#
# ALL: TISSUE SPECIFICITY STUDY: APPROACH 1 - STEP 5
#
# Collapse TAU values by:
#   - SRA.Study and confidence (Protein-coding, Low-confidence lncRNA, 
#     Medium-confidence lncRNA and High-confidence lncRNA).
#   - SRA.Study and confidence (Protein-coding, Low-confidence lncRNA, 
#     Medium-confidence lncRNA and High-confidence lncRNA) and class_code.
# 
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggExtra))
suppressMessages(library(data.table))
suppressMessages(options(bitmapType='cairo'))

## 1. PATHS

# Own computer
path_tissue_specificity = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
flag = "nr"
species_short_name = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
species_long_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")

# # Garnatxa
# path_tissue_specificity = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
# flag = "nr"
# species_short_name = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
# species_long_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")


## 2. MEAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEAN TAU (TABLE): \n"))

# Directories.
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP5"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP5"))
}

# Load table.
TAB_mean = read.table(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP3/mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")

# Change species format name.
species_tab = data.frame(Spe_short = species_short_name, Spe_long = species_long_name)
TAB_mean = merge(TAB_mean, species_tab, by.x = "Spe", by.y = "Spe_short", all.x = T)

# Convert to factors.
TAB_mean$Confidence = factor(TAB_mean$Confidence, levels = c("Protein-coding", "Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
TAB_mean$Class_code = factor(TAB_mean$Class_code, levels = c("Intergenic lncRNA (u)", "Antisense lncRNA (x)", "Intronic lncRNA (i)", "Sense lncRNA (o/e)"))

# Collapse table (1).
TAB_mean_red_1 = TAB_mean %>% 
  group_by(Spe_long, SRA.Study, Confidence) %>% 
  summarise(
    Mean.TAU = mean(TAU)
  )
TAB_mean_red_1 = as.data.frame(TAB_mean_red_1)
TAB_mean_red_1$"SRA.Study-Specie" = paste0(TAB_mean_red_1$SRA.Study, "-", TAB_mean_red_1$Spe_long)
TAB_mean_red_1$`SRA.Study-Specie` = factor(TAB_mean_red_1$`SRA.Study-Specie`, levels = sort(unique(TAB_mean_red_1$`SRA.Study-Specie`)))

# Collapse table (2).
TAB_mean_red_2 = TAB_mean %>% 
  group_by(Spe_long, SRA.Study, Confidence, Class_code) %>% 
  summarise(
    Mean.TAU = mean(TAU)
  )
TAB_mean_red_2 = as.data.frame(TAB_mean_red_2)
TAB_mean_red_2$"SRA.Study-Specie" = paste0(TAB_mean_red_2$SRA.Study, "-", TAB_mean_red_2$Spe_long)
TAB_mean_red_2$`SRA.Study-Specie` = factor(TAB_mean_red_2$`SRA.Study-Specie`, levels = sort(unique(TAB_mean_red_2$`SRA.Study-Specie`)))

# Figure (1).
gg1 = ggplot(TAB_mean_red_1, aes(x = Confidence, y = Mean.TAU, fill = Confidence)) + 
  geom_boxplot() + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  geom_jitter(aes(shape = `SRA.Study-Specie`), width = 0.3, size = 2) +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee")) +
  scale_shape_manual(values = 1:length(levels(TAB_mean_red_1$`SRA.Study-Specie`))) +
  scale_x_discrete(labels=c("Protein-coding" = "Protein-coding", "Low-confidence lncRNA" = "LC-lncRNA", "Medium-confidence lncRNA" = "MC-lncRNA", "High-confidence lncRNA" = "HC-lncRNA")) +
  xlab("") + ylab("Mean(TAU)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "right") +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP5/mean-TAU-Boxplot.png"), height = 10, width = 10, dpi = 600)
ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP5/mean-TAU-Boxplot.pdf"), height = 10, width = 10, dpi = 600)

# Figure (2).
gg2 = ggplot(TAB_mean_red_2, aes(x = Confidence, y = Mean.TAU, fill = Confidence)) + 
  geom_boxplot() + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  geom_jitter(aes(shape = `SRA.Study-Specie`), width = 0.3, size = 2) +
  facet_grid(. ~ Class_code) +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee")) +
  scale_shape_manual(values = 1:length(levels(TAB_mean_red_1$`SRA.Study-Specie`))) +
  scale_x_discrete(labels=c("Protein-coding" = "Protein-coding", "Low-confidence lncRNA" = "LC-lncRNA", "Medium-confidence lncRNA" = "MC-lncRNA", "High-confidence lncRNA" = "HC-lncRNA")) +
  xlab("") + ylab("Mean(TAU)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "right") +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14))

ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP5/mean-TAU-Boxplot_by_class_code.png"), height = 10, width = 20, dpi = 600)
ggsave(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP5/mean-TAU-Boxplot_by_class_code.pdf"), height = 10, width = 20, dpi = 600)

