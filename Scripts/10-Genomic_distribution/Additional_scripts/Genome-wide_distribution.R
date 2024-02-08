################################################################################
#
# GENOME-WIDE LNCRNA AND PCG DISTRIBUTION
#
# Create circos plots about the genome-wide lncRNA and PCG distribution.
#
# @author: pasviber - Pascual Villalba Bermell
#
################################################################################

rm(list = ls())



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
# lncRNAs_tab = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP-FINAL/Database/Database_LncRNAs_NR.tsv"
# genes_tab = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv"
# WD = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/10-Genomic_distribution/cme/Dist/nr"
# AI = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"





## 2. CHROMOSOME'S SIZE

chrs = read.table(paste0(AI, "/Chromosomes/", spe, "_chrs.txt"), header = F, sep = "\t", quote = "\"")
table_conv = data.frame(chrs.old = chrs$V1, chrs.new = sprintf("%02d", 1:length(chrs$V1)))
table_conv$chrs.new = paste0("chr", table_conv$chrs.new)
chrs.sizes = read.table(paste0(WD, "/", spe, ".sizes_genome.txt"), header = F, sep = "\t", quote = "\"")
chrs.sizes = merge(chrs.sizes, table_conv, by.x = "V1", by.y = "chrs.old", all.x = F)
chrs.sizes$V1 = chrs.sizes$chrs.new
chrs.sizes$chrs.new = NULL
colnames(chrs.sizes) = c("Chr", "Start", "End")

rm(list = c("chrs"))





#### 3. LNCRNAS DATABASE AND GENES FILE

DB = read.table(lncRNAs_tab, header = T, sep = "\t", quote = "\"")
DB = merge(DB, table_conv, by.x = "Chr", by.y = "chrs.old", all.x = F)
DB$Chr = DB$chrs.new
DB$chrs.new = NULL
Genes = read.table(genes_tab, header = T, sep = "\t", quote = "\"")
Genes$Class_code = "pc"
Genes = merge(Genes, table_conv, by.x = "Chr", by.y = "chrs.old", all.x = F)
Genes$Chr = Genes$chrs.new
Genes$chrs.new = NULL

rm(list = c("table_conv"))





## 4. CIRCOS - ALL CHROMOSOMES

#### 4.1 BY CLASS CODE

###### 4.1.1 Create bed files.

bed_u = DB[DB$Class_code == "u", c("Chr", "Start", "End")]
bed_u = bed_u[order(bed_u$Chr, bed_u$Start),]
bed_x = DB[DB$Class_code == "x", c("Chr", "Start", "End")]
bed_x = bed_x[order(bed_x$Chr, bed_x$Start),]
bed_i = DB[DB$Class_code == "i", c("Chr", "Start", "End")]
bed_i = bed_i[order(bed_i$Chr, bed_i$Start),]
bed_o_e = DB[DB$Class_code == "o" | DB$Class_code == "e", c("Chr", "Start", "End")]
bed_o_e = bed_o_e[order(bed_o_e$Chr, bed_o_e$Start),]
bed_genes = Genes[,c("Chr", "Start", "End")]
bed_genes = bed_genes[order(bed_genes$Chr, bed_genes$Start),]

###### 4.1.2  Create a table representing the gene density as the number of transcripts found in 250Kb size windows.

Tab_u = genomicDensity(bed_u, window.size = 250000, count_by = "number")
colnames(Tab_u)[colnames(Tab_u) == "value"] = "u"
Tab_x = genomicDensity(bed_x, window.size = 250000, count_by = "number")
colnames(Tab_x)[colnames(Tab_x) == "value"] = "x"
Tab_i = genomicDensity(bed_i, window.size = 250000, count_by = "number")
colnames(Tab_i)[colnames(Tab_i) == "value"] = "i"
Tab_o_e = genomicDensity(bed_o_e, window.size = 250000, count_by = "number")
colnames(Tab_o_e)[colnames(Tab_o_e) == "value"] = "o/e"
Tab_genes = genomicDensity(bed_genes, window.size = 250000, count_by = "number")
colnames(Tab_genes)[colnames(Tab_genes) == "value"] = "pc"

Tab_number_wide = merge(Tab_u, Tab_x, by = c("chr", "start", "end"), all = T)
Tab_number_wide = merge(Tab_number_wide, Tab_i, by = c("chr", "start", "end"), all = T)
Tab_number_wide = merge(Tab_number_wide, Tab_o_e, by = c("chr", "start", "end"), all = T)
Tab_number_wide = merge(Tab_number_wide, Tab_genes, by = c("chr", "start", "end"), all = T)
Tab_number_wide[is.na(Tab_number_wide)] = 0

colnames(Tab_number_wide) = c("Chr", "Start", "End", "Number_transcripts.u", "Number_transcripts.x", "Number_transcripts.i", "Number_transcripts.o/e", "Number_transcripts.pc")

write.table(Tab_number_wide, paste0(WD, "/", spe, "-ALL-Number-250Kb.tsv"), col.names = T, sep = "\t", row.names = F, quote = F)

rm(list = c("Tab_u", "Tab_x", "Tab_i", "Tab_o_e", "Tab_genes"))

###### 4.1.3  Create a figure representing the gene density as the number of transcripts found in 250Kb size windows.

## Figure (NUMBER): Number of transcripts in the window.
png(filename = paste0(WD, "/", spe, "-ALL-Number-250Kb.png"), height = 8000, width = 8000, res = 800)
circos.genomicInitialize(chrs.sizes, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.axis(h = 1, major.at = seq(0, xlim[2] + 5000000, 5000000), labels = as.character(seq(0, xlim[2], 5000000)/1000000), labels.cex = 1.2)
  circos.rect(xlim[1], 0, xlim[2], 1, col = "#E6EDF1")
  circos.text(mean(xlim), mean(ylim), chr, cex = 1.2, col = "black", facing = "inside", niceFacing = TRUE, font = 2)
}, bg.border = NA, track.height = 0.08)


if (dim(bed_u)[1] != 0) {
  circos.genomicDensity(bed_u, col = c("#f1ea0e"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_x)[1] != 0) {
  circos.genomicDensity(bed_x, col = c("#f29e0b"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_i)[1] != 0) {
  circos.genomicDensity(bed_i, col = c("#f1280d"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
} 
if (dim(bed_o_e)[1] != 0) {
  circos.genomicDensity(bed_o_e, col = c("#177300"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_genes)[1] != 0) {
  circos.genomicDensity(bed_genes, col = c("#19bde2"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
} else {
  circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
}

circos.clear()

invisible(dev.off())

rm(list = c("bed_u", "bed_x", "bed_i", "bed_o_e", "bed_genes", "Tab_number_wide"))





#### 4.2 BY CLASS CODE AND CONFIDENCE-LEVEL

for (confidence in c("Low", "Medium", "High")) {
  
  ###### 4.2.1 Create bed files.
  
  bed_u = DB[DB$Class_code == "u" & DB$Confidence == confidence, c("Chr", "Start", "End")]
  bed_u = bed_u[order(bed_u$Chr, bed_u$Start),]
  bed_x = DB[DB$Class_code == "x" & DB$Confidence == confidence, c("Chr", "Start", "End")]
  bed_x = bed_x[order(bed_x$Chr, bed_x$Start),]
  bed_i = DB[DB$Class_code == "i" & DB$Confidence == confidence, c("Chr", "Start", "End")]
  bed_i = bed_i[order(bed_i$Chr, bed_i$Start),]
  bed_o_e = DB[(DB$Class_code == "o" | DB$Class_code == "e") & DB$Confidence == confidence, c("Chr", "Start", "End")]
  bed_o_e = bed_o_e[order(bed_o_e$Chr, bed_o_e$Start),]
  bed_genes = Genes[, c("Chr", "Start", "End")]
  bed_genes = bed_genes[order(bed_genes$Chr, bed_genes$Start),]
  
  ###### 4.1.2  Create a table representing the gene density as the number of transcripts found in 250Kb size windows.
  
  Tab_u = genomicDensity(bed_u, window.size = 250000, count_by = "number")
  colnames(Tab_u)[colnames(Tab_u) == "value"] = "u"
  Tab_x = genomicDensity(bed_x, window.size = 250000, count_by = "number")
  colnames(Tab_x)[colnames(Tab_x) == "value"] = "x"
  Tab_i = genomicDensity(bed_i, window.size = 250000, count_by = "number")
  colnames(Tab_i)[colnames(Tab_i) == "value"] = "i"
  Tab_o_e = genomicDensity(bed_o_e, window.size = 250000, count_by = "number")
  colnames(Tab_o_e)[colnames(Tab_o_e) == "value"] = "o/e"
  Tab_genes = genomicDensity(bed_genes, window.size = 250000, count_by = "number")
  colnames(Tab_genes)[colnames(Tab_genes) == "value"] = "pc"
  
  Tab_number_wide = merge(Tab_u, Tab_x, by = c("chr", "start", "end"), all = T)
  Tab_number_wide = merge(Tab_number_wide, Tab_i, by = c("chr", "start", "end"), all = T)
  Tab_number_wide = merge(Tab_number_wide, Tab_o_e, by = c("chr", "start", "end"), all = T)
  Tab_number_wide = merge(Tab_number_wide, Tab_genes, by = c("chr", "start", "end"), all = T)
  Tab_number_wide[is.na(Tab_number_wide)] = 0
  
  colnames(Tab_number_wide) = c("Chr", "Start", "End", "Number_transcripts.u", "Number_transcripts.x", "Number_transcripts.i", "Number_transcripts.o/e", "Number_transcripts.pc")
  
  write.table(Tab_number_wide, paste0(WD, "/", spe, "-", confidence, "-Number-250Kb.tsv"), col.names = T, sep = "\t", row.names = F, quote = F)
  
  rm(list = c("Tab_u", "Tab_x", "Tab_i", "Tab_o_e", "Tab_genes"))
  
  ###### 4.2.3  Create a figure representing the gene density as the number of transcripts found in 250Kb size windows.
  
  png(filename = paste0(WD, "/", spe, "-", confidence, "-Number-250Kb.png"), height = 8000, width = 8000, res = 800)
  circos.genomicInitialize(chrs.sizes, plotType = NULL)
  
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.axis(h = 1, major.at = seq(0, xlim[2] + 5000000, 5000000), labels = as.character(seq(0, xlim[2], 5000000)/1000000), labels.cex = 1.2)
    circos.rect(xlim[1], 0, xlim[2], 1, col = "#E6EDF1")
    circos.text(mean(xlim), mean(ylim), chr, cex = 1.2, col = "black", facing = "inside", niceFacing = TRUE, font = 2)
  }, bg.border = NA, track.height = 0.08)
  
  if (dim(bed_u)[1] != 0) {
    circos.genomicDensity(bed_u, col = c("#f1ea0e"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_x)[1] != 0) {
    circos.genomicDensity(bed_x, col = c("#f29e0b"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_i)[1] != 0) {
    circos.genomicDensity(bed_i, col = c("#f1280d"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
  } 
  if (dim(bed_o_e)[1] != 0) {
    circos.genomicDensity(bed_o_e, col = c("#177300"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_genes)[1] != 0) {
    circos.genomicDensity(bed_genes, col = c("#19bde2"), track.height = 0.1, baseline = 0, window.size = 250000, count_by = "number")
  } else {
    circos.track(factors = chrs.sizes$Chr, ylim = c(0, 1), track.height = 0.1)
  }
  
  circos.clear()
  
  invisible(dev.off())
  
  rm(list = c("bed_u", "bed_x", "bed_i", "bed_o_e", "bed_genes", "Tab_number_wide"))
}

rm(list = c("Genes", "DB", "chrs.sizes"))
