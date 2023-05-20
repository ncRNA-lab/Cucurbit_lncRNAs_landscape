################################################################################
# GENOME-WIDE LNCRNAS AND GENES DISTRIBUTION
# Create circos plots about the genome-wide lncRNAs and genes distribution
################################################################################

#https://mran.microsoft.com/snapshot/2016-08-07/web/packages/circlize/vignettes/genomic_plot.pdf
#https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html
#https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0027121.g001


## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(dplyr))
suppressMessages(library(circlize))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggExtra))
suppressMessages(options(bitmapType='cairo'))





## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments must be supplied.", call.=FALSE)
} else {
  spe = args[1]
  lncRNAs_tab = args[2]
  genes_tab = args[3]
  WD = args[4]
  AI = args[5]
}

# spe = "cme"
# lncRNAs_tab="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme/STEP-FINAL/Database/Database_LncRNAs_NR.tsv"
# genes_tab="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv"
# WD = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/10-Distribution/Transcript_density/nr/cme"
# AI = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Additional_info"





## 2. CHROMOSOME'S SIZE

chrs = read.table(paste0(AI, "/Chromosomes/", spe, "_chrs.txt"), header = F, sep = "\t", quote = "\"")
table_conv = data.frame(chrs.old = chrs$V1, chrs.new = sprintf("%02d", 1:length(chrs$V1)))
chrs.sizes = read.table(paste0(WD, "/sizes_genome.txt"), header = F, sep = "\t", quote = "\"")
chrs.sizes = merge(chrs.sizes, table_conv, by.x = "V1", by.y = "chrs.old", all.x = F)
chrs.sizes$V1 = chrs.sizes$chrs.new
chrs.sizes$chrs.new = NULL
colnames(chrs.sizes) = c("chr", "start", "end")

rm(list = c("chrs"))





#### 3. LNCRNAS DATABASE AND GENES FILE

DB = read.table(lncRNAs_tab, header = T, sep = "\t", quote = "\"")
DB = merge(DB, table_conv, by.x = "Chr", by.y = "chrs.old", all.x = F)
DB$Chr = DB$chrs.new
DB$chrs.new = NULL
Genes = read.table(genes_tab, header = T, sep = "\t", quote = "\"")
Genes = merge(Genes, table_conv, by.x = "Chr", by.y = "chrs.old", all.x = F)
Genes$Chr = Genes$chrs.new
Genes$chrs.new = NULL





## 4. CIRCOS - ALL CHROMOSOMES

#### 4.1 BY CLASS CODE

## Create bed files.
bed_u = DB[DB$Class_code == "u", c("Chr", "Start", "End")]
colnames(bed_u) = c("chr", "start", "end")
bed_u = bed_u[order(bed_u$chr, bed_u$start),]
bed_x = DB[DB$Class_code == "x", c("Chr", "Start", "End")]
colnames(bed_x) = c("chr", "start", "end")
bed_x = bed_x[order(bed_x$chr, bed_x$start),]
bed_i = DB[DB$Class_code == "i", c("Chr", "Start", "End")]
colnames(bed_i) = c("chr", "start", "end")
bed_i = bed_i[order(bed_i$chr, bed_i$start),]
bed_o_e = DB[DB$Class_code == "o" | DB$Class_code == "e", c("Chr", "Start", "End")]
colnames(bed_o_e) = c("chr", "start", "end")
bed_o_e = bed_o_e[order(bed_o_e$chr, bed_o_e$start),]
bed_genes = Genes[,c("Chr", "Start", "End")]
colnames(bed_genes) = c("chr", "start", "end")
bed_genes = bed_genes[order(bed_genes$chr, bed_genes$start),]



## Figure (NUMBER): Number of transcripts in the window.
png(filename = paste0(WD, "/Density-ALL-Number.png"), height = 8000, width = 8000, res = 800)
circos.genomicInitialize(chrs.sizes, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.axis(h = 1, major.at = seq(0, xlim[2] + 5000000, 5000000), labels = as.character(seq(0, xlim[2], 5000000)/1000000), labels.cex = 0.8)
  circos.rect(xlim[1], 0, xlim[2], 1, col = "#E6EDF1")
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.9, col = "black", facing = "inside", niceFacing = TRUE)
}, bg.border = NA, track.height = 0.07)


if (dim(bed_u)[1] != 0) {
  circos.genomicDensity(bed_u, col = c("#ece68f"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_x)[1] != 0) {
  circos.genomicDensity(bed_x, col = c("#decfa9"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_i)[1] != 0) {
  circos.genomicDensity(bed_i, col = c("#dcb8b8"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
} 
if (dim(bed_o_e)[1] != 0) {
  circos.genomicDensity(bed_o_e, col = c("#b5cf9b"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_genes)[1] != 0) {
  circos.genomicDensity(bed_genes, col = c("#a2ded9"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
}

circos.clear()

invisible(dev.off())

## Table (NUMBER): Number of transcripts in the window.
Tab_u = genomicDensity(bed_u, window.size = 100000, count_by = "number")
colnames(Tab_u)[colnames(Tab_u) == "value"] = "u"
Tab_x = genomicDensity(bed_x, window.size = 100000, count_by = "number")
colnames(Tab_x)[colnames(Tab_x) == "value"] = "x"
Tab_i = genomicDensity(bed_i, window.size = 100000, count_by = "number")
colnames(Tab_i)[colnames(Tab_i) == "value"] = "i"
Tab_o_e = genomicDensity(bed_o_e, window.size = 100000, count_by = "number")
colnames(Tab_o_e)[colnames(Tab_o_e) == "value"] = "o/e"
Tab_genes = genomicDensity(bed_genes, window.size = 100000, count_by = "number")
colnames(Tab_genes)[colnames(Tab_genes) == "value"] = "pc"

Tab_number_wide = merge(Tab_u, Tab_x, by = c("chr", "start", "end"), all = T)
Tab_number_wide = merge(Tab_number_wide, Tab_i, by = c("chr", "start", "end"), all = T)
Tab_number_wide = merge(Tab_number_wide, Tab_o_e, by = c("chr", "start", "end"), all = T)
Tab_number_wide = merge(Tab_number_wide, Tab_genes, by = c("chr", "start", "end"), all = T)
Tab_number_wide[is.na(Tab_number_wide)] = 0

write.table(Tab_number_wide, paste0(WD, "/Density-ALL-TAB-Number.tsv"), col.names = T, sep = "\t", row.names = F, quote = F)

Tab_number_long = melt(setDT(Tab_number_wide), id.vars = c("chr", "start", "end"), variable.name = "class_code")
colnames(Tab_number_long) = c("chr", "start", "end", "class_code", "counts")
Tab_number_mean = Tab_number_long %>% group_by(class_code) %>% summarise(mean = mean(counts), threshold = ceiling(mean(counts)) * 2)
Tab_number_mean = as.data.frame(Tab_number_mean)
Tab_number_long = merge(Tab_number_long, Tab_number_mean, by = "class_code", all = T)
Tab_number_long$"high_density_window" = ifelse(Tab_number_long$counts >= Tab_number_long$threshold, 1, 0)
Tab_number_long$class_code = factor(Tab_number_long$class_code, levels = c("u", "x", "i", "o/e", "pc"))
Tab_number_long$chr = factor(Tab_number_long$chr, levels = chrs.sizes$chr)
Tab_number_summary = Tab_number_long %>% group_by(class_code, chr, threshold, .drop = F) %>% summarise(high_density_windows = sum(high_density_window))

write.table(Tab_number_summary, paste0(WD, "/Density-ALL-TAB_Summary-Number.tsv"), col.names = T, sep = "\t", row.names = F, quote = F)

gg = ggplot(Tab_number_summary, aes(x = class_code, y = high_density_windows, fill = class_code)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#a2ded9")) +
  facet_grid(.~chr) +
  xlab("") +
  ylab("Number of high density windows") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5)) +
  theme(legend.position = "none")

ggsave(paste0(WD, "/High_density_windows-ALL-Number.png"), height = 5, width = 20, dpi = 600)



## Figure (PERCENT): Percent of window covered by transcripts.
png(filename = paste0(WD, "/Density-ALL-Percent.png"), height = 8000, width = 8000, res = 800)
circos.genomicInitialize(chrs.sizes, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.axis(h = 1, major.at = seq(0, xlim[2] + 5000000, 5000000), labels = as.character(seq(0, xlim[2], 5000000)/1000000), labels.cex = 0.8)
  circos.rect(xlim[1], 0, xlim[2], 1, col = "#E6EDF1")
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.9, col = "black", facing = "inside", niceFacing = TRUE)
}, bg.border = NA, track.height = 0.07)

if (dim(bed_u)[1] != 0) {
  circos.genomicDensity(bed_u, col = c("#ece68f"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_x)[1] != 0) {
  circos.genomicDensity(bed_x, col = c("#decfa9"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_i)[1] != 0) {
  circos.genomicDensity(bed_i, col = c("#dcb8b8"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
} 
if (dim(bed_o_e)[1] != 0) {
  circos.genomicDensity(bed_o_e, col = c("#b5cf9b"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_genes)[1] != 0) {
  circos.genomicDensity(bed_genes, col = c("#a2ded9"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
} else {
  circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
}

circos.clear()

invisible(dev.off())

## Table (PERCENT): Percent of window covered by transcripts.
Tab_u = genomicDensity(bed_u, window.size = 100000, count_by = "percent")
colnames(Tab_u)[colnames(Tab_u) == "value"] = "u"
Tab_x = genomicDensity(bed_x, window.size = 100000, count_by = "percent")
colnames(Tab_x)[colnames(Tab_x) == "value"] = "x"
Tab_i = genomicDensity(bed_i, window.size = 100000, count_by = "percent")
colnames(Tab_i)[colnames(Tab_i) == "value"] = "i"
Tab_o_e = genomicDensity(bed_o_e, window.size = 100000, count_by = "percent")
colnames(Tab_o_e)[colnames(Tab_o_e) == "value"] = "o/e"
Tab_genes = genomicDensity(bed_genes, window.size = 100000, count_by = "percent")
colnames(Tab_genes)[colnames(Tab_genes) == "value"] = "pc"

Tab_percent_wide = merge(Tab_u, Tab_x, by = c("chr", "start", "end"), all = T)
Tab_percent_wide = merge(Tab_percent_wide, Tab_i, by = c("chr", "start", "end"), all = T)
Tab_percent_wide = merge(Tab_percent_wide, Tab_o_e, by = c("chr", "start", "end"), all = T)
Tab_percent_wide = merge(Tab_percent_wide, Tab_genes, by = c("chr", "start", "end"), all = T)
Tab_percent_wide[is.na(Tab_percent_wide)] = 0

write.table(Tab_percent_wide, paste0(WD, "/Density-ALL-TAB-Percent.tsv"), col.names = T, sep = "\t", row.names = F, quote = F)





#### 4.2 BY CLASS CODE AND CONFIDENCE-LEVEL

for (confidence in c("Low", "Medium", "High")) {
  
  ## Create bed files.
  bed_u = DB[DB$Class_code == "u" & DB$Significance_level == confidence, c("Chr", "Start", "End")]
  colnames(bed_u) = c("chr", "start", "end")
  bed_u = bed_u[order(bed_u$chr, bed_u$start),]
  bed_x = DB[DB$Class_code == "x" & DB$Significance_level == confidence, c("Chr", "Start", "End")]
  colnames(bed_x) = c("chr", "start", "end")
  bed_x = bed_x[order(bed_x$chr, bed_x$start),]
  bed_i = DB[DB$Class_code == "i" & DB$Significance_level == confidence, c("Chr", "Start", "End")]
  colnames(bed_i) = c("chr", "start", "end")
  bed_i = bed_i[order(bed_i$chr, bed_i$start),]
  bed_o_e = DB[(DB$Class_code == "o" | DB$Class_code == "e") & DB$Significance_level == confidence, c("Chr", "Start", "End")]
  colnames(bed_o_e) = c("chr", "start", "end")
  bed_o_e = bed_o_e[order(bed_o_e$chr, bed_o_e$start),]
  bed_genes = Genes[, c("Chr", "Start", "End")]
  colnames(bed_genes) = c("chr", "start", "end")
  bed_genes = bed_genes[order(bed_genes$chr, bed_genes$start),]
  
  ## Figure (NUMBER): Number of transcripts in the window.
  png(filename = paste0(WD, "/Density-", confidence, "-Number.png"), height = 8000, width = 8000, res = 800)
  circos.genomicInitialize(chrs.sizes, plotType = NULL)
  
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.axis(h = 1, major.at = seq(0, xlim[2] + 5000000, 5000000), labels = as.character(seq(0, xlim[2], 5000000)/1000000), labels.cex = 0.8)
    circos.rect(xlim[1], 0, xlim[2], 1, col = "#E6EDF1")
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.9, col = "black", facing = "inside", niceFacing = TRUE)
  }, bg.border = NA, track.height = 0.07)
  
  if (dim(bed_u)[1] != 0) {
    circos.genomicDensity(bed_u, col = c("#ece68f"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_x)[1] != 0) {
    circos.genomicDensity(bed_x, col = c("#decfa9"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_i)[1] != 0) {
    circos.genomicDensity(bed_i, col = c("#dcb8b8"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  } 
  if (dim(bed_o_e)[1] != 0) {
    circos.genomicDensity(bed_o_e, col = c("#b5cf9b"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_genes)[1] != 0) {
    circos.genomicDensity(bed_genes, col = c("#a2ded9"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  }
  
  circos.clear()
  
  invisible(dev.off())
  
  ## Table (NUMBER): Number of transcripts in the window.
  Tab_u = genomicDensity(bed_u, window.size = 100000, count_by = "number")
  colnames(Tab_u)[colnames(Tab_u) == "value"] = "u"
  Tab_x = genomicDensity(bed_x, window.size = 100000, count_by = "number")
  colnames(Tab_x)[colnames(Tab_x) == "value"] = "x"
  Tab_i = genomicDensity(bed_i, window.size = 100000, count_by = "number")
  colnames(Tab_i)[colnames(Tab_i) == "value"] = "i"
  Tab_o_e = genomicDensity(bed_o_e, window.size = 100000, count_by = "number")
  colnames(Tab_o_e)[colnames(Tab_o_e) == "value"] = "o/e"
  Tab_genes = genomicDensity(bed_genes, window.size = 100000, count_by = "number")
  colnames(Tab_genes)[colnames(Tab_genes) == "value"] = "pc"
  
  Tab_number_wide = merge(Tab_u, Tab_x, by = c("chr", "start", "end"), all = T)
  Tab_number_wide = merge(Tab_number_wide, Tab_i, by = c("chr", "start", "end"), all = T)
  Tab_number_wide = merge(Tab_number_wide, Tab_o_e, by = c("chr", "start", "end"), all = T)
  Tab_number_wide = merge(Tab_number_wide, Tab_genes, by = c("chr", "start", "end"), all = T)
  Tab_number_wide[is.na(Tab_number_wide)] = 0
  
  write.table(Tab_number_wide, paste0(WD, "/Density-", confidence, "-TAB-Number.tsv"), col.names = T, sep = "\t", row.names = F, quote = F)
  
  Tab_number_long = melt(setDT(Tab_number_wide), id.vars = c("chr", "start", "end"), variable.name = "class_code")
  colnames(Tab_number_long) = c("chr", "start", "end", "class_code", "counts")
  Tab_number_mean = Tab_number_long %>% group_by(class_code) %>% summarise(mean = mean(counts), threshold = ceiling(mean(counts)) * 2)
  Tab_number_mean = as.data.frame(Tab_number_mean)
  Tab_number_long = merge(Tab_number_long, Tab_number_mean, by = "class_code", all = T)
  Tab_number_long$"high_density_window" = ifelse(Tab_number_long$counts >= Tab_number_long$threshold, 1, 0)
  Tab_number_long$class_code = factor(Tab_number_long$class_code, levels = c("u", "x", "i", "o/e", "pc"))
  Tab_number_long$chr = factor(Tab_number_long$chr, levels = chrs.sizes$chr)
  Tab_number_summary = Tab_number_long %>% group_by(class_code, chr, threshold, .drop = F) %>% summarise(high_density_windows = sum(high_density_window))
  
  write.table(Tab_number_summary, paste0(WD, "/Density-", confidence, "-TAB_Summary-Number.tsv"), col.names = T, sep = "\t", row.names = F, quote = F)
  
  gg = ggplot(Tab_number_summary, aes(x = class_code, y = high_density_windows, fill = class_code)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#a2ded9")) +
    facet_grid(.~chr) +
    xlab("") +
    ylab("Number of high density windows") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5)) +
    theme(legend.position = "none")
  
  ggsave(paste0(WD, "/High_density_windows-", confidence, "-Number.png"), height = 5, width = 20, dpi = 600)
  
  
  
  ## Figure (PERCENT): Percent of window covered by transcripts.
  png(filename = paste0(WD, "/Density-", confidence, "-Percent.png"), height = 8000, width = 8000, res = 800)
  circos.genomicInitialize(chrs.sizes, plotType = NULL)
  
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.axis(h = 1, major.at = seq(0, xlim[2] + 5000000, 5000000), labels = as.character(seq(0, xlim[2], 5000000)/1000000), labels.cex = 0.8)
    circos.rect(xlim[1], 0, xlim[2], 1, col = "#E6EDF1")
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.9, col = "black", facing = "inside", niceFacing = TRUE)
  }, bg.border = NA, track.height = 0.07)
  
  if (dim(bed_u)[1] != 0) {
    circos.genomicDensity(bed_u, col = c("#ece68f"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_x)[1] != 0) {
    circos.genomicDensity(bed_x, col = c("#decfa9"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_i)[1] != 0) {
    circos.genomicDensity(bed_i, col = c("#dcb8b8"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  } 
  if (dim(bed_o_e)[1] != 0) {
    circos.genomicDensity(bed_o_e, col = c("#b5cf9b"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_genes)[1] != 0) {
    circos.genomicDensity(bed_genes, col = c("#a2ded9"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "percent")
  } else {
    circos.track(factors = chrs.sizes$chr, ylim = c(0, 1), track.height = 0.1)
  }
  
  circos.clear()
  
  invisible(dev.off())
  
  ## Table (PERCENT): Percent of window covered by transcripts.
  Tab_u = genomicDensity(bed_u, window.size = 100000, count_by = "percent")
  colnames(Tab_u)[colnames(Tab_u) == "value"] = "u"
  Tab_x = genomicDensity(bed_x, window.size = 100000, count_by = "percent")
  colnames(Tab_x)[colnames(Tab_x) == "value"] = "x"
  Tab_i = genomicDensity(bed_i, window.size = 100000, count_by = "percent")
  colnames(Tab_i)[colnames(Tab_i) == "value"] = "i"
  Tab_o_e = genomicDensity(bed_o_e, window.size = 100000, count_by = "percent")
  colnames(Tab_o_e)[colnames(Tab_o_e) == "value"] = "o/e"
  Tab_genes = genomicDensity(bed_genes, window.size = 100000, count_by = "percent")
  colnames(Tab_genes)[colnames(Tab_genes) == "value"] = "pc"
  
  Tab_percent_wide = merge(Tab_u, Tab_x, by = c("chr", "start", "end"), all = T)
  Tab_percent_wide = merge(Tab_percent_wide, Tab_i, by = c("chr", "start", "end"), all = T)
  Tab_percent_wide = merge(Tab_percent_wide, Tab_o_e, by = c("chr", "start", "end"), all = T)
  Tab_percent_wide = merge(Tab_percent_wide, Tab_genes, by = c("chr", "start", "end"), all = T)
  Tab_percent_wide[is.na(Tab_percent_wide)] = 0
  
  write.table(Tab_percent_wide, paste0(WD, "/Density-", confidence, "-TAB-Percent.tsv"), col.names = T, sep = "\t", row.names = F, quote = F)
}

