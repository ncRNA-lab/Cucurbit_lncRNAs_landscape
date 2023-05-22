################################################################################
#
# Percentage of nucleotides by transcript coming from some kind of transposon or 
# repetitive element that can be found in the genome of each specie. We take always 
# into account the class code of each transcript.
#
# STEP 3
#
################################################################################

rm(list = ls())

## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Create_Tables_and_Figures_C-STEP3.R path_res folder_new folder_RepMask Flag confidence

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggExtra))
suppressMessages(options(bitmapType='cairo'))
options(stringsAsFactors = F)

## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments must be supplied.", call.=FALSE)
} else {
  path_res = args[1]
  folder_new = args[2]
  folder_RepMask = args[3]
  Flag = args[4]
  confidence = args[5]
}

# path_res = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results"
# folder_new = "02-Comparison_Genes_LncRNAs"
# folder_RepMask = "05-RepeatMasker"
# Flag = "NR"
# confidence = "High"


## 2. PIPELINE

# Species
species = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
species_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")
species_tab = data.frame(spe = species, name = species_name)


### 2.2 CREATE SUBTABLES

cat("\nCreating the subtables...\n")

### Paths.
WD_01 = paste0(path_res, "/11-TEs_and_genomic_repeats/", folder_new)
WD_02 = paste0(path_res, "/05-predict_lncRNAs")

### Load tables.
TAB_FINAL1_1 = read.table(paste0(WD_01, "/Figures_and_Tables/C/Final_tab-Repeat-", Flag, "-", confidence, "-COLLAPSED_REPEAT_ALL.tsv"), header = T, sep = "\t", quote = "\"")
TAB_FINAL1_2 = read.table(paste0(WD_01, "/Figures_and_Tables/C/Final_tab-Repeat-", Flag, "-", confidence, "-COLLAPSED_REPEAT.tsv"), header = T, sep = "\t", quote = "\"")
TAB_FINAL2 = read.table(paste0(WD_01, "/Figures_and_Tables/C/Final_tab-Transposon-", Flag, "-", confidence, "-COLLAPSED_TRANSPOSON.tsv"), header = T, sep = "\t", quote = "\"")

### Convert to factor.
TAB_FINAL1_1$spe = factor(TAB_FINAL1_1$spe, levels = species)
TAB_FINAL1_2$spe = factor(TAB_FINAL1_2$spe, levels = species)
TAB_FINAL2$spe = factor(TAB_FINAL2$spe, levels = species)

TAB_FINAL1_1$Class_code = factor(TAB_FINAL1_1$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_FINAL1_2$Class_code = factor(TAB_FINAL1_2$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_FINAL2$Class_code = factor(TAB_FINAL2$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))

TAB_FINAL1_1$Repeat_type_2_mod = factor(TAB_FINAL1_1$Repeat_type_2_mod, levels = c("Simple repeat", "Low complexity repeat", "LINE/SINE", "LTR", "DNA transposon", "Unknown", "No repeat"))

### Add percentage masked genome.
masked_perc = read.table(paste0(path_res, "/11-TEs_and_genomic_repeats/01-Repeat_calling/", folder_RepMask, "/masked_genome_percentage.tsv"), header = T, sep = "\t", quote = "\"")
TAB_FINAL1_1 = merge(TAB_FINAL1_1, masked_perc, by = "spe", all = T)


### 2.3 FIGURE

cat("\nDrawing the TAB_FINAL1_1 table...\n")

TAB_FINAL1_1 = merge(TAB_FINAL1_1, species_tab, by = "spe", all.x = T, all.y = F)
TAB_FINAL1_1$name = factor(TAB_FINAL1_1$name, levels = species_name)
TAB_FINAL1_1 = TAB_FINAL1_1[order(TAB_FINAL1_1$name),]
TAB_FINAL1_1$"facet" = paste0(TAB_FINAL1_1$name, " (", TAB_FINAL1_1$masked_perc, "%)")
TAB_FINAL1_1$facet = factor(TAB_FINAL1_1$facet, levels = unique(TAB_FINAL1_1$facet))

gg1 = ggplot(TAB_FINAL1_1, aes(x = Repeat_type_2_mod, y = overlap_per, fill = Repeat_type_2_mod)) +
  geom_boxplot(outlier.size=0.5) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 1, color = "black", fill = "black") +
  facet_grid(facet ~ Class_code) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000")) +
  xlab("") +
  ylab("Repeat content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, angle = 80)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 13),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 9, face = "bold.italic"))

ggsave(paste0(WD_01, "/Figures_and_Tables/C/REPEAT_ALL-", Flag, "-", confidence, ".png"), height = 16, width = 12, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/REPEAT_ALL-", Flag, "-", confidence, ".pdf"), height = 16, width = 12, dpi = 600)


gg1 = ggplot(TAB_FINAL1_1, aes(x = Repeat_type_2_mod, y = Log.overlap_per, fill = Repeat_type_2_mod)) +
  geom_boxplot(outlier.size=0.5) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 1, color = "black", fill = "black") +
  facet_grid(facet ~ Class_code) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000")) +
  xlab("") +
  ylab("Log2(Repeat content (%) + 1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, angle = 80)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 13),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 9, face = "bold.italic"))

ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_REPEAT_ALL-", Flag, "-", confidence, ".png"), height = 16, width = 12, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_REPEAT_ALL-", Flag, "-", confidence, ".pdf"), height = 16, width = 12, dpi = 600)




cat("\nDrawing the TAB TAB_FINAL1_2 table...\n")

TAB_FINAL1_2 = merge(TAB_FINAL1_2, species_tab, by = "spe", all.x = T, all.y = F)
TAB_FINAL1_2$name = factor(TAB_FINAL1_2$name, levels = species_name)

gg2 = ggplot(TAB_FINAL1_2, aes(x = Class_code, y = overlap_per, fill = Class_code)) +
  geom_boxplot() + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  facet_wrap(name~., nrow = 1) +
  xlab("") +
  ylab("Repeat content (%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 17, face = "bold.italic"))

ggsave(paste0(WD_01, "/Figures_and_Tables/C/REPEAT-", Flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/REPEAT-", Flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)

gg2 = gg2 +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(paste0(WD_01, "/Figures_and_Tables/C/REPEAT_mod-", Flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/REPEAT_mod-", Flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)


gg2 = ggplot(TAB_FINAL1_2, aes(x = Class_code, y = Log.overlap_per, fill = Class_code)) +
  geom_boxplot() + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  facet_wrap(name~., nrow = 1) +
  xlab("") +
  ylab("Log2(Repeat content (%) + 1)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 17, face = "bold.italic"))

ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_REPEAT-", Flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_REPEAT-", Flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)

gg2 = gg2 +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_REPEAT_mod-", Flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_REPEAT_mod-", Flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)




cat("\nDrawing the TAB TAB_FINAL2 table...\n")

TAB_FINAL2 = merge(TAB_FINAL2, species_tab, by = "spe", all.x = T, all.y = F)
TAB_FINAL2$name = factor(TAB_FINAL2$name, levels = species_name)

gg3 = ggplot(TAB_FINAL2, aes(x = Class_code, y = overlap_per, fill = Class_code)) +
  geom_boxplot() + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  facet_wrap(name~., nrow = 1) +
  xlab("") +
  ylab("Transposon content (%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 17, face = "bold.italic"))

ggsave(paste0(WD_01, "/Figures_and_Tables/C/TRANSPOSON-", Flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/TRANSPOSON-", Flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)

gg3 = gg3 +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(paste0(WD_01, "/Figures_and_Tables/C/TRANSPOSON_mod-", Flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/TRANSPOSON_mod-", Flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)


gg3 = ggplot(TAB_FINAL2, aes(x = Class_code, y = Log.overlap_per, fill = Class_code)) +
  geom_boxplot() + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  facet_wrap(name~., nrow = 1) +
  xlab("") +
  ylab("Log2(Transposon content (%) + 1)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 17, face = "bold.italic"))

ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_TRANSPOSON-", Flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_TRANSPOSON-", Flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)

gg3 = gg3 +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_TRANSPOSON_mod-", Flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/C/LOG_TRANSPOSON_mod-", Flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)

rm(list = c("gg1", "gg2", "gg3"))

