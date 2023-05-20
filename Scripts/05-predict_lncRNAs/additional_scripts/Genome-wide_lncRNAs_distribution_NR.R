################################################################################
# GENOME-WIDE LNCRNAS DISTRIBUTION
# Create circos plots about the genome-wide lncRNAs distribution
################################################################################

#https://mran.microsoft.com/snapshot/2016-08-07/web/packages/circlize/vignettes/genomic_plot.pdf
#https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html
#https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0027121.g001


## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(dplyr))
suppressMessages(library(circlize))
suppressMessages(options(bitmapType='cairo'))



## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied.", call.=FALSE)
} else {
  specie = args[1]
  WD = args[2]
  AI = args[3]
}

# specie = "cla"
# WD = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cla/STEP-FINAL"
# AI = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Additional_info"


## 2. CHROMOSOME'S SIZE

chrs = read.table(paste0(AI, "/Chromosomes/", specie, "_chrs.txt"), header = F, sep = "\t", quote = "\"")
chrs = chrs$V1
chrs.sizes = read.table(paste0(WD, "/Figures/Circos/sizes_genome.txt"), header = F, sep = "\t", quote = "\"")
chrs.sizes = chrs.sizes[chrs.sizes$V1 %in% chrs,]



#### 3. LNCRNAS DATABASE

DB = read.table(paste0(WD, "/Database/Database_LncRNAs_NR.tsv"), header = T, sep = "\t", quote = "\"")



## 4. CIRCOS - ALL CHROMOSOMES

#### 4.1 BY CLASS CODE

bed_u = DB[DB$Class_code == "u", c("Chr", "Start", "End")]
colnames(bed_u) = c("chr", "start", "end")
bed_x = DB[DB$Class_code == "x", c("Chr", "Start", "End")]
colnames(bed_x) = c("chr", "start", "end")
bed_i = DB[DB$Class_code == "i", c("Chr", "Start", "End")]
colnames(bed_i) = c("chr", "start", "end")
bed_o_e = DB[DB$Class_code == "o" | DB$Class_code == "e", c("Chr", "Start", "End")]
colnames(bed_o_e) = c("chr", "start", "end")

png(filename = paste0(WD, "/Figures/Circos/nr/density_by_class_code-ALL.png"), height = 8000, width = 8000, res = 800)
circos.initializeWithIdeogram(chrs.sizes, plotType = c("axis", "labels"), labels.cex = 2)

if (dim(bed_u)[1] != 0) {
  circos.genomicDensity(bed_u, col = c("#ece68f"), track.height = 0.1, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_x)[1] != 0) {
  circos.genomicDensity(bed_x, col = c("#decfa9"), track.height = 0.1, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.1)
}
if (dim(bed_i)[1] != 0) {
  circos.genomicDensity(bed_i, col = c("#dcb8b8"), track.height = 0.1, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.1)
} 
if (dim(bed_o_e)[1] != 0) {
  circos.genomicDensity(bed_o_e, col = c("#b5cf9b"), track.height = 0.1, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.1)
}

circos.clear()

invisible(dev.off())


#### 4.2 BY CLASS CODE AND CONFIDENCE-LEVEL

bed_u_Low = DB[DB$Class_code == "u" & DB$Significance_level == "Low", c("Chr", "Start", "End")]
colnames(bed_u_Low) = c("chr", "start", "end")
bed_u_Medium = DB[DB$Class_code == "u" & DB$Significance_level == "Medium", c("Chr", "Start", "End")]
colnames(bed_u_Medium) = c("chr", "start", "end")
bed_u_High = DB[DB$Class_code == "u" & DB$Significance_level == "High", c("Chr", "Start", "End")]
colnames(bed_u_High) = c("chr", "start", "end")

bed_x_Low = DB[DB$Class_code == "x" & DB$Significance_level == "Low", c("Chr", "Start", "End")]
colnames(bed_x_Low) = c("chr", "start", "end")
bed_x_Medium = DB[DB$Class_code == "x" & DB$Significance_level == "Medium", c("Chr", "Start", "End")]
colnames(bed_x_Medium) = c("chr", "start", "end")
bed_x_High = DB[DB$Class_code == "x" & DB$Significance_level == "High", c("Chr", "Start", "End")]
colnames(bed_x_High) = c("chr", "start", "end")

bed_i_Low = DB[DB$Class_code == "i" & DB$Significance_level == "Low", c("Chr", "Start", "End")]
colnames(bed_i_Low) = c("chr", "start", "end")
bed_i_Medium = DB[DB$Class_code == "i" & DB$Significance_level == "Medium", c("Chr", "Start", "End")]
colnames(bed_i_Medium) = c("chr", "start", "end")
bed_i_High = DB[DB$Class_code == "i" & DB$Significance_level == "High", c("Chr", "Start", "End")]
colnames(bed_i_High) = c("chr", "start", "end")

bed_o_e_Low = DB[(DB$Class_code == "o" | DB$Class_code == "e") & DB$Significance_level == "Low", c("Chr", "Start", "End")]
colnames(bed_o_e_Low) = c("chr", "start", "end")
bed_o_e_Medium = DB[(DB$Class_code == "o" | DB$Class_code == "e") & DB$Significance_level == "Medium", c("Chr", "Start", "End")]
colnames(bed_o_e_Medium) = c("chr", "start", "end")
bed_o_e_High = DB[(DB$Class_code == "o" | DB$Class_code == "e") & DB$Significance_level == "High", c("Chr", "Start", "End")]
colnames(bed_o_e_High) = c("chr", "start", "end")

png(filename = paste0(WD, "/Figures/Circos/nr/density_by_class_code-LOW.png"), height = 8000, width = 8000, res = 800)
circos.initializeWithIdeogram(chrs.sizes, plotType = c("axis", "labels"), labels.cex = 2)

if (dim(bed_u_Low)[1] != 0) {
  circos.genomicDensity(bed_u_Low, col = c("#ece68f"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_x_Low)[1] != 0) {
  circos.genomicDensity(bed_x_Low, col = c("#decfa9"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_i_Low)[1] != 0) {
  circos.genomicDensity(bed_i_Low, col = c("#dcb8b8"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_o_e_Low)[1] != 0) {
  circos.genomicDensity(bed_o_e_Low, col = c("#b5cf9b"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}

circos.clear()

invisible(dev.off())


png(filename = paste0(WD, "/Figures/Circos/nr/density_by_class_code-MEDIUM.png"), height = 8000, width = 8000, res = 800)
circos.initializeWithIdeogram(chrs.sizes, plotType = c("axis", "labels"), labels.cex = 2)

if (dim(bed_u_Medium)[1] != 0) {
  circos.genomicDensity(bed_u_Medium, col = c("#ece68f"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_x_Medium)[1] != 0) {
  circos.genomicDensity(bed_x_Medium, col = c("#decfa9"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_i_Medium)[1] != 0) {
  circos.genomicDensity(bed_i_Medium, col = c("#dcb8b8"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_o_e_Medium)[1] != 0) {
  circos.genomicDensity(bed_o_e_Medium, col = c("#b5cf9b"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}

circos.clear()

invisible(dev.off())


png(filename = paste0(WD, "/Figures/Circos/nr/density_by_class_code-HIGH.png"), height = 8000, width = 8000, res = 800)
circos.initializeWithIdeogram(chrs.sizes, plotType = c("axis", "labels"), labels.cex = 2)

if (dim(bed_u_High)[1] != 0) {
  circos.genomicDensity(bed_u_High, col = c("#ece68f"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_x_High)[1] != 0) {
  circos.genomicDensity(bed_x_High, col = c("#decfa9"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_i_High)[1] != 0) {
  circos.genomicDensity(bed_i_High, col = c("#dcb8b8"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}
if (dim(bed_o_e_High)[1] != 0) {
  circos.genomicDensity(bed_o_e_High, col = c("#b5cf9b"), track.height = 0.05, baseline = 0)
} else {
  circos.track(factors = chrs.sizes$V1, ylim = c(0, 1), track.height = 0.05)
}

circos.clear()

invisible(dev.off())




## 5. CIRCOS BY CHROMOSOME

for (chr in chrs) {
  
  chrs.size = chrs.sizes[chrs.sizes$V1 == chr,]
  
  bed_u = DB[DB$Class_code == "u" & DB$Chr == chr, c("Chr", "Start", "End")]
  colnames(bed_u) = c("chr", "start", "end")
  bed_x = DB[DB$Class_code == "x" & DB$Chr == chr, c("Chr", "Start", "End")]
  colnames(bed_x) = c("chr", "start", "end")
  bed_i = DB[DB$Class_code == "i" & DB$Chr == chr, c("Chr", "Start", "End")]
  colnames(bed_i) = c("chr", "start", "end")
  bed_o_e = DB[(DB$Class_code == "o" | DB$Class_code == "e")& DB$Chr == chr, c("Chr", "Start", "End")]
  colnames(bed_o_e) = c("chr", "start", "end")
  
  png(filename = paste0(WD, "/Figures/Circos/nr/density_by_class_code-", chr, ".png"), height = 8000, width = 8000, res = 800)
  circos.initializeWithIdeogram(chrs.size, plotType = c("axis", "labels"), labels.cex = 2)
  
  if (dim(bed_u)[1] != 0) {
    circos.genomicDensity(bed_u, col = c("#ece68f"), track.height = 0.1, baseline = 0)
  } else {
    circos.track(factors = chrs.size$V1, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_x)[1] != 0) {
    circos.genomicDensity(bed_x, col = c("#decfa9"), track.height = 0.1, baseline = 0)
  } else {
    circos.track(factors = chrs.size$V1, ylim = c(0, 1), track.height = 0.1)
  }
  if (dim(bed_i)[1] != 0) {
    circos.genomicDensity(bed_i, col = c("#dcb8b8"), track.height = 0.1, baseline = 0)
  } else {
    circos.track(factors = chrs.size$V1, ylim = c(0, 1), track.height = 0.1)
  } 
  if (dim(bed_o_e)[1] != 0) {
    circos.genomicDensity(bed_o_e, col = c("#b5cf9b"), track.height = 0.1, baseline = 0)
  } else {
    circos.track(factors = chrs.size$V1, ylim = c(0, 1), track.height = 0.1)
  }
  
  circos.clear()
  
  invisible(dev.off())
}
