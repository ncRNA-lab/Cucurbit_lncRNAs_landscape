################################################################################
#
# CREATE THE SATURATION CURVES
#
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("Cairo"))
suppressMessages(options(bitmapType='cairo'))
suppressMessages(library("ggplot2"))


## 1. VARIABLES

WD = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/13-saturation_plots/cme"
range = 5

## 2. FIGURES

cat("\nFIGURES...\n")

if (!dir.exists(paste0(WD, "/Figures_CS"))){
  dir.create(paste0(WD, "/Figures_CS"))
}

cat("\nFigure 1...\n")

Batch_tab_1_RED = read.table(paste0(WD, "/Tables_CS/FINAL-Batch_tab_1-", range, ".tsv"), sep = "\t", header = T, quote = "\"")

Batch_tab_1_RED$Class_code = factor(Batch_tab_1_RED$Class_code, levels = c("intergenic (u)", "antisense (x)", "intronic (i)", "sense (o/e)"))
Batch_tab_1_RED$Confidence_level = factor(Batch_tab_1_RED$Confidence_level, levels = c("Low", "Medium", "High"))

gg1 = ggplot(Batch_tab_1_RED, aes(x = Size, y = Mean.Counts, color = Class_code)) + 
  geom_point() + 
  geom_line() + 
  facet_grid(Confidence_level ~ ., scales = 'free') +
  scale_color_manual(values = c("#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b")) +
  scale_x_continuous(breaks = seq(range,max(Batch_tab_1_RED$Size),range)) +
  xlab("Number of samples") +
  ylab("Number of LncRNAs") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave(paste0(WD, "/Figures_CS/Saturation_curve_1.png"), height = 10, width = 20, dpi = 600)

cat("\nFigure 2...\n")

Batch_tab_2_RED = read.table(paste0(WD, "/Tables_CS/FINAL-Batch_tab_2-", range, ".tsv"), sep = "\t", header = T, quote = "\"")

Batch_tab_2_RED$Class_code = factor(Batch_tab_2_RED$Class_code, levels = c("intergenic (u)", "antisense (x)", "intronic (i)", "sense (o/e)"))

gg2 = ggplot(Batch_tab_2_RED, aes(x = Size, y = Mean.Counts, color = Class_code)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b")) +
  scale_x_continuous(breaks = seq(range,max(Batch_tab_2_RED$Size),range)) +
  xlab("Number of samples") +
  ylab("Number of LncRNAs") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave(paste0(WD, "/Figures_CS/Saturation_curve_2.png"), height = 6, width = 12, dpi = 600)

rm(list = c("gg1", "gg2"))

