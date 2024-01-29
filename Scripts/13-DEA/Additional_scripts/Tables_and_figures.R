


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
options(bitmapType='cairo')











################################################################################
################################################################################
################################## VARIABLES ###################################
################################################################################
################################################################################


## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments must be supplied.", call.=FALSE)
} else {
  path_DEA = args[1]
  path_db_lnc = args[2]
  path_db_PCG = args[3]
  alpha_value = as.numeric(args[4])
  spe = args[5]
}

# path_DEA = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/16-DEA/cme"
# path_db_lnc = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme/STEP-FINAL/Database/Database_LncRNAs_NR.tsv"
# path_db_PCG = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv"
# alpha_value = 0.05
# spe = "cme"










################################################################################
################################################################################
################################## FIGURES #####################################
################################################################################
################################################################################

cat(paste0("Figures..."))

if (!dir.exists(paste0(path_DEA, "/05-Tables_and_Figures/ALL"))){
  dir.create(paste0(path_DEA, "/05-Tables_and_Figures/ALL"))
}
if (!dir.exists(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Database"))){
  dir.create(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Database"))
}
if (!dir.exists(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables"))){
  dir.create(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables"))
}
if (!dir.exists(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Volcano"))){
  dir.create(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Volcano"))
}
if (!dir.exists(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Percentage_DE"))){
  dir.create(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Percentage_DE"))
}
if (!dir.exists(paste0(path_DEA, "/05-Tables_and_Figures/ALL/VennDiagram"))){
  dir.create(paste0(path_DEA, "/05-Tables_and_Figures/ALL/VennDiagram"))
}










############################################
## 0. GLOBAL DATABASE
############################################

# Por si fuera necesaria la informacion de cada uno de los transcritos acerca de
# la clase de transcrito o su nivel de confianza, procedemos a a unir la base de
# datos de los lncRNAs y de los PCGs.

cat(paste0("\n\n0-GLOBAL DATABASE..."))

# Load lncRNAs database.
db_lnc = read.table(path_db_lnc, sep = "\t", header = T, quote = "\"")
db_lnc = db_lnc[, c("ID_transcript", "Class_code", "Confidence")]

db_lnc[db_lnc == "u"] = "lincRNAs"
db_lnc[db_lnc == "x"] = "NAT-lncRNAs"
db_lnc[db_lnc == "i"] = "int-lncRNAs"
db_lnc[db_lnc == "o"] = "SOT-lncRNAs"
db_lnc[db_lnc == "e"] = "SOT-lncRNAs"

db_lnc[db_lnc == "High"] = "HC-lncRNAs"
db_lnc[db_lnc == "Medium"] = "MC-lncRNAs"
db_lnc[db_lnc == "Low"] = "LC-lncRNAs"

# Load PCGs database.
db_PCG = read.table(path_db_PCG, sep = "\t", header = T, quote = "\"")
db_PCG$"Class_code" = "PC genes"
db_PCG$"Confidence" = "PC genes"
db_PCG = db_PCG[, c("ID_transcript", "Class_code", "Confidence")]

# Create global database.
db = rbind(db_lnc, db_PCG)

# Save
write.table(db, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Database/Database.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("db_lnc", "db_PCG"))










############################################
## 1. BINARY TABLE
############################################

# Aqui calculamos una tabla binaria (1 o 0) en la que las filas son los transcritos
# (lncRNAs y genes) y las columnas son cada uno de los experimentos-contraste en los
# que se ha llevado a cabo un analisis de expresion diferencial (DEA). La primera 
# tabla contendrá todos los transcritos esten DE o Non-DE. La segunda tabla (filt) 
# contendrá aquellos transcritos que estan DE en al menos alguno de los experiemtnos-
# contraste. Paralelamente se crean dos tablas iguales que las anteriores que contendran 
# el log2FC en la posicion del 1. Por ultimo, se crean dos tablas resumen en las 
# que se colapsan los experimentos-contraste por el tipo de estudio (Desarrollo, 
# estres abiotico y estres biotico). Estas dos ultimas tablas se crean partiendo 
# de la tabla binaria filtrada.

cat(paste0("\n\n1-BINARY TABLE..."))

if (!dir.exists(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables"))){
  dir.create(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables"))
}

# Load DEA summary table (Experiment-Contrast information).
Summary_tab = read.table(paste0(path_DEA, "/04-DEA/SUMMARY.tsv"), header = T, sep = "\t", quote = "\"")
Summary_tab$"ID" = paste0(Summary_tab$Experiment, ".", Summary_tab$n_contrast)

# Create an empty binary table.
Binary_tab = db
Binary_LFC_tab = db

# Experiments.
Experiments = list.files(paste0(path_DEA, "/04-DEA/02-DEA_sig"))

# Fill the binary table by Experiment-Contrast.
for (experiment in Experiments) {
  
  cat(paste0("\n\nExperiment: ", experiment))
  
  path_sig = paste0(path_DEA, "/04-DEA/02-DEA_sig/", experiment)
  path_fig = paste0(path_DEA, "/05-Tables_and_Figures")
  
  for (i in 1:length(list.files(path_sig))) {
    
    cat(paste0("\n\tContrast: ", i))
    
    # Load DEA table with DE transcripts.
    tab_DE = read.table(paste0(path_sig, "/", experiment, "-dea_sig-", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Add a new column called as the Experiment-Contrast to the empty binary table. 
    # 1: DE and 0: Non-DE.
    Binary_tab[[paste0(experiment, ".", i)]] = ifelse(Binary_tab$ID_transcript %in% tab_DE$ID_transcript, 1, 0)
    Binary_LFC_tab = merge(Binary_LFC_tab, tab_DE[, c("ID_transcript", "Shrunkenlog2FoldChange")], by = "ID_transcript", all = T)
    Binary_LFC_tab[is.na(Binary_LFC_tab)] = 0
    colnames(Binary_LFC_tab)[colnames(Binary_LFC_tab) == "Shrunkenlog2FoldChange"] = paste0(experiment, ".", i)
    
    rm(list = c("tab_DE"))
  }
  
  rm(list = c("i", "path_sig", "path_fig"))
}

rm(list = c("experiment", "Experiments"))

# Add a new column called "Sum" which refers to the sum of rows.
if (length(Summary_tab$ID) > 1) {
  Binary_tab$Sum = rowSums(Binary_tab[, Summary_tab$ID])
} else {
  Binary_tab$Sum = Binary_tab[, Summary_tab$ID]
}

write.table(Binary_tab, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_tab.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(Binary_LFC_tab, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_LFC_tab.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

# SUMMARY TABLE 1. Collapse the Experiment-Contrast columns according to the class 
# of experiment (Development, Abiotic stress and biotic stress).

Binary_tab_SUMMARY_1 = Binary_tab[, c("ID_transcript", "Class_code", "Confidence")]
for (class in c("Development", "Abiotic stress", "Biotic stress")) {
  Experiments = Summary_tab[Summary_tab$Class == class, "ID"]
  n_exp = length(Experiments)
  if (n_exp != 0) {
    if (n_exp == 1) {
      Binary_tab_SUMMARY_1[, gsub(" ", ".", class)] = Binary_tab[, Experiments]
    } else {
      Binary_tab_SUMMARY_1[, gsub(" ", ".", class)] = rowSums(Binary_tab[, Experiments])
    }
    Binary_tab_SUMMARY_1[, paste0("Total.", gsub(" ", ".", class))] = n_exp
    Binary_tab_SUMMARY_1[, paste0("Norm.", gsub(" ", ".", class))] = Binary_tab_SUMMARY_1[, gsub(" ", ".", class)]/n_exp
  }
  else {
    Binary_tab_SUMMARY_1[, gsub(" ", ".", class)] = 0
    Binary_tab_SUMMARY_1[, paste0("Total.", gsub(" ", ".", class))] = 0
    Binary_tab_SUMMARY_1[, paste0("Norm.", gsub(" ", ".", class))] = 0
  }
}

write.table(Binary_tab_SUMMARY_1, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_tab_filt-SUMMARY_1.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

# SUMMARY TABLE 2. Collapse the experiment-contrast columns according to the class 
# of experiment (Development and Stress).

Binary_tab_SUMMARY_2 = Binary_tab[, c("ID_transcript", "Class_code", "Confidence")]
for (class in c("Development", "Stress")) {
  
  if (class == "Development") {
    Experiments = Summary_tab[Summary_tab$Class == class, "ID"]
    n_exp = length(Experiments)
    if (n_exp != 0) {
      if (n_exp == 1) {
        Binary_tab_SUMMARY_2[, gsub(" ", ".", class)] = Binary_tab[, Experiments]
      } else {
        Binary_tab_SUMMARY_2[, gsub(" ", ".", class)] = rowSums(Binary_tab[, Experiments])
      }
      Binary_tab_SUMMARY_2[, paste0("Total.", gsub(" ", ".", class))] = n_exp
      Binary_tab_SUMMARY_2[, paste0("Norm.", gsub(" ", ".", class))] = Binary_tab_SUMMARY_2[, gsub(" ", ".", class)]/n_exp
    } else {
      Binary_tab_SUMMARY_2[, gsub(" ", ".", class)] = 0
      Binary_tab_SUMMARY_2[, paste0("Total.", gsub(" ", ".", class))] = 0
      Binary_tab_SUMMARY_2[, paste0("Norm.", gsub(" ", ".", class))] = 0
    }
  } else {
    Experiments = Summary_tab[Summary_tab$Class == "Biotic stress" | Summary_tab$Class == "Abiotic stress", "ID"]
    n_exp = length(Experiments)
    if (n_exp != 0) {
      if (n_exp == 1) {
        Binary_tab_SUMMARY_2[, "Stress"] = Binary_tab[, Experiments]
      } else {
        Binary_tab_SUMMARY_2[, "Stress"] = rowSums(Binary_tab[, Experiments])
      }
      Binary_tab_SUMMARY_2[, "Total.Stress"] = n_exp
      Binary_tab_SUMMARY_2[, "Norm.Stress"] = Binary_tab_SUMMARY_2[, "Stress"]/n_exp
    } else {
      Binary_tab_SUMMARY_2[, "Stress"] = 0
      Binary_tab_SUMMARY_2[, "Total.Stress"] = 0
      Binary_tab_SUMMARY_2[, "Norm.Stress"] = 0
    }
  }
}

write.table(Binary_tab_SUMMARY_2, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_tab_filt-SUMMARY_2.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

# Filter the tables. Keep only transcripts DE in at least one Experiment-Contrast.
Binary_tab_filt = Binary_tab[Binary_tab$Sum != 0,]
Binary_LFC_tab_filt = Binary_LFC_tab[Binary_LFC_tab$ID_transcript %in% Binary_tab_filt$ID_transcript,]
Binary_tab_SUMMARY_1_filt = Binary_tab_SUMMARY_1[Binary_tab_SUMMARY_1$ID_transcript %in% Binary_tab_filt$ID_transcript,]
Binary_tab_SUMMARY_2_filt = Binary_tab_SUMMARY_2[Binary_tab_SUMMARY_2$ID_transcript %in% Binary_tab_filt$ID_transcript,]

write.table(Binary_tab_filt, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_tab_filt.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(Binary_LFC_tab_filt, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_LFC_tab_filt.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(Binary_tab_SUMMARY_1_filt, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_tab-SUMMARY_1_filt.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(Binary_tab_SUMMARY_2_filt, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_tab-SUMMARY_2_filt.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("Binary_tab", "Binary_tab_filt", "Binary_LFC_tab", "Binary_LFC_tab_filt", 
            "Binary_tab_SUMMARY_1", "Binary_tab_SUMMARY_1_filt", "Binary_tab_SUMMARY_2", 
            "Binary_tab_SUMMARY_2_filt", "Experiments", "Summary_tab", "class", "n_exp"))










############################################
## 2. VOLCANO PLOT
############################################

# Usando los transcritos tanto DE como Non-DE se representa el log10BaseMean vs el
# log2FC.

cat(paste0("\n\n2-VOLCANO PLOT..."))

# Load DEA summary table (Experiment-Contrast information).
Summary_tab = read.table(paste0(path_DEA, "/04-DEA/SUMMARY.tsv"), header = T, sep = "\t", quote = "\"")
Summary_tab$"ID" = paste0(Summary_tab$Experiment, ".", Summary_tab$n_contrast)

# Experiments.
Experiments = list.files(paste0(path_DEA, "/04-DEA/01-DEA_raw"))

######## VOLCANO BY EXPERIMENT-CONTRAST.

TAB = data.frame()
TAB_DE = data.frame()

for (experiment in Experiments) {
  
  cat(paste0("\n\nExperiment: ", experiment))
  
  # Paths.
  path_raw = paste0(path_DEA, "/04-DEA/01-DEA_raw/", experiment)
  path_fig = paste0(path_DEA, "/05-Tables_and_Figures/", experiment)
  
  # Create output directory.
  if (!dir.exists(path_fig)) {
    dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
  }
  
  for (i in 1:length(list.files(path_raw))) {
    
    cat(paste0("\n\tContrast: ", i))
    
    # Load DEA table with DE and Non-DE transcripts.
    tab_DEA = read.table(paste0(path_raw, "/", experiment, "-dea_raw-", i, ".tsv"), sep = "\t", header = T, quote = "\"")
    
    # Merge DEA table and global db. ALL.y = F
    tab_DEA = merge(tab_DEA, db, by = "ID_transcript", all.x = T, all.y = F)
    
    # Add parameters.
    tab_DEA$"Size" = ifelse(tab_DEA$padj > alpha_value, 0.8, ifelse(tab_DEA$Class_code == "PC genes", 0.8, 1.5))
    tab_DEA$"Legend.1" = ifelse(tab_DEA$padj > alpha_value, "Non-significant", tab_DEA$Class_code)
    tab_DEA$"Legend.2" = ifelse(tab_DEA$padj > alpha_value, "Non-significant", tab_DEA$Confidence)
    
    tab_DEA$Legend.1 = factor(tab_DEA$Legend.1, levels = c("Non-significant", "PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
    tab_DEA$Legend.2 = factor(tab_DEA$Legend.2, levels = c("Non-significant", "PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))
    
    tab_DEA$ID = paste0(experiment, ".", i)
    
    TAB = rbind(TAB, tab_DEA)
    
    ## All the transcripts (DE and non-DE).
    
    # Volcano by class code.
    gg1 = ggplot(tab_DEA, aes(x = baseMean, y = Shrunkenlog2FoldChange, colour= Legend.1)) + 
      geom_point(size = tab_DEA$Size) +
      scale_color_manual(values=c("grey", "#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66")) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                    labels = trans_format("log10", math_format(10^.x))) +
      xlab("BaseMean") + ylab("Log2FC") + theme_bw() + removeGrid(x = T, y = T) +
      geom_hline(yintercept=0, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=1, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=-1, color = "black", linetype="dashed", linewidth=0.2) +
      theme(legend.position = "top", 
            legend.title = element_blank(), 
            axis.text = element_text(size = 15), 
            axis.title = element_text(size = 19), 
            legend.text = element_text(size = 13))
    
    ggsave(paste0(path_fig, "/", experiment, "-", i, "-Volcano_class_code-ALL.png"), height = 8, width = 12, dpi = 800)
    #ggsave(paste0(path_fig, "/", experiment, "-", i, "-Volcano_class_code-ALL.pdf"), height = 8, width = 12, dpi = 800)
    
    # Volcano by confidence level.
    gg2 = ggplot(tab_DEA, aes(x = baseMean, y = Shrunkenlog2FoldChange, colour= Legend.2)) + 
      geom_point(size = tab_DEA$Size) +
      scale_color_manual(values=c("grey", "#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee")) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                    labels = trans_format("log10", math_format(10^.x))) +
      xlab("BaseMean") + ylab("Log2FC") + theme_bw() + removeGrid(x = T, y = T) +
      geom_hline(yintercept=0, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=1, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=-1, color = "black", linetype="dashed", linewidth=0.2) +
      theme(legend.position = "top", 
            legend.title = element_blank(), 
            axis.text = element_text(size = 15), 
            axis.title = element_text(size = 19), 
            legend.text = element_text(size = 13))
    
    ggsave(paste0(path_fig, "/", experiment, "-", i, "-Volcano_confidence_level-ALL.png"), height = 8, width = 12, dpi = 800)
    #ggsave(paste0(path_fig, "/", experiment, "-", i, "-Volcano_confidence_level-ALL.pdf"), height = 8, width = 12, dpi = 800)
    
    ## Only DE transcripts.
    
    tab_DE = tab_DEA[tab_DEA$padj <= 0.05,]
    
    tab_DE$Legend.1 = factor(tab_DE$Legend.1, levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
    tab_DE$Legend.2 = factor(tab_DE$Legend.2, levels = c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))
    
    TAB_DE = rbind(TAB_DE, tab_DE)
    
    # Volcano by class code.
    gg3 = ggplot(tab_DE, aes(x = baseMean, y = Shrunkenlog2FoldChange, colour= Legend.1)) + 
      geom_point(size = tab_DE$Size) +
      scale_color_manual(values=c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66")) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                    labels = trans_format("log10", math_format(10^.x))) +
      xlab("BaseMean") + ylab("Log2FC") + theme_bw() + removeGrid(x = T, y = T) +
      geom_hline(yintercept=0, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=1, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=-1, color = "black", linetype="dashed", linewidth=0.2) +
      theme(legend.position = "top", 
            legend.title = element_blank(), 
            axis.text = element_text(size = 15), 
            axis.title = element_text(size = 19), 
            legend.text = element_text(size = 13))
    
    ggsave(paste0(path_fig, "/", experiment, "-", i, "-Volcano_class_code.png"), height = 8, width = 12, dpi = 800)
    #ggsave(paste0(path_fig, "/", experiment, "-", i, "-Volcano_class_code.pdf"), height = 8, width = 12, dpi = 800)
    
    # Volcano by confidence level.
    gg4 = ggplot(tab_DE, aes(x = baseMean, y = Shrunkenlog2FoldChange, colour= Legend.2)) + 
      geom_point(size = tab_DE$Size) +
      scale_color_manual(values=c("#7cc1cf", "#ce4fd3", "#eb92ef", "#edd0ee")) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                    labels = trans_format("log10", math_format(10^.x))) +
      xlab("BaseMean") + ylab("Log2FC") + theme_bw() + removeGrid(x = T, y = T) +
      geom_hline(yintercept=0, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=1, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=-1, color = "black", linetype="dashed", linewidth=0.2) +
      theme(legend.position = "top", 
            legend.title = element_blank(), 
            axis.text = element_text(size = 15), 
            axis.title = element_text(size = 19), 
            legend.text = element_text(size = 13))
    
    ggsave(paste0(path_fig, "/", experiment, "-", i, "-Volcano_confidence_level.png"), height = 8, width = 12, dpi = 800)
    #ggsave(paste0(path_fig, "/", experiment, "-", i, "-Volcano_confidence_level.pdf"), height = 8, width = 12, dpi = 800)
    
    rm(list = c("tab_DEA", "tab_DE", "gg1", "gg2", "gg3", "gg4"))
  }
  
  rm(list = c("i", "path_raw", "path_fig"))
}

rm(list = c("experiment", "Experiments"))

# Modify final volcano table.
TAB_mod = merge(TAB, Summary_tab[, c("ID", "Class")], by = "ID", all = T)
TAB_DE_mod = merge(TAB_DE, Summary_tab[, c("ID", "Class")], by = "ID", all = T)

write.table(TAB_mod, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Volcano/TAB-volcano.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(TAB_DE_mod, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Volcano/TAB_DE-volcano.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

######## VOLCANO BY CLASS OF EXPERIMENT (DEVELOPMENT, ABIOTIC STRESS AND BIOTIC STRESS).

# Volcano by confidence-level and class (Development, Abiotic stress and biotic stress).
for (co in unique(c("LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))) {
  for (class in unique(c("Development", "Abiotic stress", "Biotic stress"))) {
    
    # LncRNAs.
    tab1 = TAB_DE_mod[TAB_DE_mod$Class == class & TAB_DE_mod$Legend.2 == co,]
    tab1$Legend.1 = factor(as.character(tab1$Legend.1), levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
    tab1 = tab1[order(tab1$Legend.1),]
    
    gg1 = ggplot(tab1, aes(x = baseMean, y = Shrunkenlog2FoldChange, colour= Legend.1)) + 
      geom_point(size = tab1$Size) +
      scale_color_manual(values=c("#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66")) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                    labels = trans_format("log10", math_format(10^.x))) +
      xlab("BaseMean") + ylab("Log2FC") + theme_bw() + removeGrid(x = T, y = T) +
      geom_hline(yintercept=0, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=1, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=-1, color = "black", linetype="dashed", linewidth=0.2) +
      theme(legend.position = "top", 
            legend.title = element_blank(), 
            axis.text = element_text(size = 15), 
            axis.title = element_text(size = 19), 
            legend.text = element_text(size = 13))
    
    ggsave(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Volcano/Volcano_", class, "_", co, "-1.png"), height = 8, width = 12, dpi = 800)
    
    # PCGs and LncRNAs.
    tab2 = TAB_DE_mod[TAB_DE_mod$Class == class & (TAB_DE_mod$Legend.2 == co | TAB_DE_mod$Legend.2 == "PC genes"),]
    tab2$Legend.1 = factor(as.character(tab2$Legend.1), levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
    tab2 = tab2[order(tab2$Legend.1),]
    
    gg2 = ggplot(tab2, aes(x = baseMean, y = Shrunkenlog2FoldChange, colour = Legend.1)) + 
      geom_point(size = tab2$Size) +
      scale_color_manual(values=c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66")) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                    labels = trans_format("log10", math_format(10^.x))) +
      xlab("BaseMean") + ylab("Log2FC") + theme_bw() + removeGrid(x = T, y = T) +
      geom_hline(yintercept=0, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=1, color = "black", linetype="dashed", linewidth=0.2) +
      geom_hline(yintercept=-1, color = "black", linetype="dashed", linewidth=0.2) +
      theme(legend.position = "top", 
            legend.title = element_blank(), 
            axis.text = element_text(size = 15), 
            axis.title = element_text(size = 19), 
            legend.text = element_text(size = 13))
    
    ggsave(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Volcano/Volcano_", class, "_", co, "-2.png"), height = 8, width = 12, dpi = 800)
  }
}

######## VOLCANO BY CLASS OF EXPERIMENT (DEVELOPMENT AND STRESS).

# Volcano by confidence-level and Stress class.
for (co in unique(c("LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))) {
  
  # LncRNAs.
  tab1 = TAB_DE_mod[(TAB_DE_mod$Class == "Abiotic stress" | TAB_DE_mod$Class == "Biotic stress") & TAB_DE_mod$Legend.2 == co,]
  tab1$Legend.1 = factor(as.character(tab1$Legend.1), levels = c("lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  tab1 = tab1[order(tab1$Legend.1),]
  
  gg1 = ggplot(tab1, aes(x = baseMean, y = Shrunkenlog2FoldChange, colour= Legend.1)) + 
    geom_point(size = tab1$Size) +
    scale_color_manual(values=c("#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66")) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("BaseMean") + ylab("Log2FC") + theme_bw() + removeGrid(x = T, y = T) +
    geom_hline(yintercept=0, color = "black", linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept=1, color = "black", linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept=-1, color = "black", linetype="dashed", linewidth=0.2) +
    theme(legend.position = "top", 
          legend.title = element_blank(), 
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 19), 
          legend.text = element_text(size = 13))
  
  ggsave(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Volcano/Volcano_Stress_", co, "-1.png"), height = 8, width = 12, dpi = 800)
  
  # PCGs and LncRNAs.
  tab2 = TAB_DE_mod[(TAB_DE_mod$Class == "Abiotic stress" | TAB_DE_mod$Class == "Biotic stress") & (TAB_DE_mod$Legend.2 == co | TAB_DE_mod$Legend.2 == "PC genes"),]
  tab2$Legend.1 = factor(as.character(tab2$Legend.1), levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
  tab2 = tab2[order(tab2$Legend.1),]
  
  gg2 = ggplot(tab2, aes(x = baseMean, y = Shrunkenlog2FoldChange, colour = Legend.1)) + 
    geom_point(size = tab2$Size) +
    scale_color_manual(values=c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66")) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("BaseMean") + ylab("Log2FC") + theme_bw() + removeGrid(x = T, y = T) +
    geom_hline(yintercept=0, color = "black", linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept=1, color = "black", linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept=-1, color = "black", linetype="dashed", linewidth=0.2) +
    theme(legend.position = "top", 
          legend.title = element_blank(), 
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 19), 
          legend.text = element_text(size = 13))
  
  ggsave(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Volcano/Volcano_Stress_", co, "-2.png"), height = 8, width = 12, dpi = 800)
}

rm(list = c("TAB", "TAB_mod", "TAB_DE", "TAB_DE_mod", "tab1", "tab2", "gg1", "gg2", "class", "co", "Summary_tab"))










############################################
## 3. PERCENTAGE AND NUMBER OF DE TRANSCRIPTS
############################################

# Se calcula el porcentaje de transcrito DE usando como total tanto transcritos
# DE como Non-DE

cat(paste0("\n\n3-PERCENTAGE AND NUMBER OF DE TRANSCRIPTS..."))

# Load summary table.
Summary_tab = read.table(paste0(path_DEA, "/04-DEA/SUMMARY.tsv"), header = T, sep = "\t", quote = "\"")
Summary_tab$"ID" = paste0(Summary_tab$Experiment, ".", Summary_tab$n_contrast)

# Experiments.
Experiments = list.files(paste0(path_DEA, "/04-DEA/01-DEA_raw"))

TAB_1 = data.frame()
TAB_2 = data.frame()
TAB_1_mod = data.frame()
TAB_2_mod = data.frame()

# Lollipop by experiment.
for (experiment in Experiments) {

  cat(paste0("\n\nExperiment: ", experiment))

  # Paths.
  path_raw = paste0(path_DEA, "/04-DEA/01-DEA_raw/", experiment)
  path_fig = paste0(path_DEA, "/05-Tables_and_Figures/", experiment)

  # Create output directory.
  if (!dir.exists(path_fig)) {
    dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
  }

  # Volcano by contrast.
  for (i in 1:length(list.files(path_raw))) {

    cat(paste0("\n\tContrast: ", i))

    # Load DEA table with DE and Non-DE transcripts.
    tab_DEA = read.table(paste0(path_raw, "/", experiment, "-dea_raw-", i, ".tsv"), sep = "\t", header = T, quote = "\"")

    # Merge DEA table and db. all.y = T.
    tab_ALL = merge(tab_DEA, db, by = "ID_transcript", all.x = T, all.y = T)
    tab_ALL$padj = ifelse(is.na(tab_ALL$padj), 1, tab_ALL$padj)
    tab_ALL$Shrunkenlog2FoldChange = ifelse(is.na(tab_ALL$Shrunkenlog2FoldChange), 0, tab_ALL$Shrunkenlog2FoldChange)
    
    rm(list = c("tab_DEA"))
    
    # Add type info.
    tab_ALL$Type = ifelse(tab_ALL$padj <= 0.05, "Significant", "Non-significant")

    # Factors
    tab_ALL$Class_code = factor(tab_ALL$Class_code, levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
    tab_ALL$Confidence = factor(tab_ALL$Confidence, levels = c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))
    tab_ALL$Type = factor(tab_ALL$Type, levels = c("Non-significant", "Significant"))

    # Collapse 1.
    tab_ALL_collap_1 = suppressMessages(tab_ALL %>%
                                          group_by(Class_code, Confidence, Type, .drop = F) %>%
                                          summarise(Number = n_distinct(ID_transcript)) %>%
                                          group_by(Class_code, Confidence, .drop = F) %>%
                                          mutate(Total = sum(Number),
                                                 Percentage = round((Number*100)/sum(Number), 2)))
    tab_ALL_collap_1 = tab_ALL_collap_1[tab_ALL_collap_1$Total != 0,]
    tab_ALL_collap_1_S = tab_ALL_collap_1[tab_ALL_collap_1$Type == "Significant",]
    
    tab_ALL_collap_1_S$ID = paste0(experiment, ".", i)
    
    TAB_1 = rbind(TAB_1, tab_ALL_collap_1_S)
    
    temp = tab_ALL_collap_1_S[tab_ALL_collap_1_S$Confidence == "PC genes",]
    temp_LC = temp
    temp_LC$Confidence = "LC-lncRNAs"
    temp_MC = temp
    temp_MC$Confidence = "MC-lncRNAs"
    temp_HC = temp
    temp_HC$Confidence = "HC-lncRNAs"

    tab_ALL_collap_1_S_mod = tab_ALL_collap_1_S[tab_ALL_collap_1_S$Confidence != "PC genes", ]
    tab_ALL_collap_1_S_mod = rbind(tab_ALL_collap_1_S_mod, temp_LC, temp_MC, temp_HC)
    tab_ALL_collap_1_S_mod$Confidence = factor(tab_ALL_collap_1_S_mod$Confidence, levels = c("LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))
    
    TAB_1_mod = rbind(TAB_1_mod, tab_ALL_collap_1_S_mod)
    
    # Lollipop
    gg1 = ggplot(tab_ALL_collap_1_S_mod, aes(x = Class_code, y = Percentage, fill = Class_code)) +
      geom_segment(aes(x=Class_code, xend=Class_code, y=0, yend=Percentage)) +
      geom_point( size=5, alpha=0.7, shape=21, stroke=2) +
      scale_fill_manual(values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")) +
      scale_y_continuous(limits = c(0, max(tab_ALL_collap_1_S_mod$Percentage) + max(tab_ALL_collap_1_S_mod$Percentage)*0.1), breaks = seq(0, max(tab_ALL_collap_1_S_mod$Percentage), 10), expand = c(0, 1)) +
      facet_wrap(Confidence~., nrow = 1) +
      xlab("") +
      ylab("Percentage") +
      geom_text(aes(label = round(Percentage, 2)), vjust = -1.3) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 14, angle = 50),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 18),
            strip.text.x = element_text(size = 17, face = "bold"))

    ggsave(paste0(path_fig, "/", experiment, "-", i, "-Lollipop-Percentage.png"), height = 6, width = 12, dpi = 800)
    #ggsave(paste0(path_fig, "/", experiment, "-", i, "-Lollipop-Percentage.pdf"), height = 6, width = 12, dpi = 800)

    gg2 = ggplot(tab_ALL_collap_1_S_mod, aes(x = Class_code, y = Number, fill = Class_code)) +
      geom_segment(aes(x=Class_code, xend=Class_code, y=0, yend=Number)) +
      geom_point( size=5, alpha=0.7, shape=21, stroke=2) +
      scale_fill_manual(values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66", "#5d65b4")) +
      scale_y_continuous(limits = c(0, max(tab_ALL_collap_1_S_mod$Number) + max(tab_ALL_collap_1_S_mod$Number)*0.1), breaks = seq(0, max(tab_ALL_collap_1_S_mod$Number), 500), expand = c(0, 1)) +
      facet_wrap(Confidence~., nrow = 1) +
      xlab("") +
      ylab("Number") +
      geom_text(aes(label = Number), vjust = -1.3) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 14, angle = 50),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 18),
            strip.text.x = element_text(size = 17, face = "bold"))

    ggsave(paste0(path_fig, "/", experiment, "-", i, "-Lollipop-Number.png"), height = 6, width = 12, dpi = 800)
    #ggsave(paste0(path_fig, "/", experiment, "-", i, "-Lollipop-Number.pdf"), height = 6, width = 12, dpi = 800)
    
    rm(list = c("tab_ALL_collap_1", "tab_ALL_collap_1_S", "tab_ALL_collap_1_S_mod",
                "temp", "temp_LC", "temp_MC", "temp_HC", "gg1", "gg2"))
    
    # Collapse 2.
    tab_ALL_collap_2 = suppressMessages(tab_ALL %>%
                                          group_by(Confidence, Type, .drop = F) %>%
                                          summarise(Number = n_distinct(ID_transcript)) %>%
                                          group_by(Confidence, .drop = F) %>%
                                          mutate(Total = sum(Number),
                                                 Percentage = round((Number*100)/sum(Number), 2)))
    tab_ALL_collap_2 = tab_ALL_collap_2[tab_ALL_collap_2$Total != 0,]
    tab_ALL_collap_2_S = tab_ALL_collap_2[tab_ALL_collap_2$Type == "Significant",]
    
    tab_ALL_collap_2_S$ID = paste0(experiment, ".", i)
    
    TAB_2 = rbind(TAB_2, tab_ALL_collap_2_S)
    
    temp = tab_ALL_collap_2_S[tab_ALL_collap_2_S$Confidence == "PC genes",]
    temp_LC = temp
    temp_LC$Confidence = "LC-lncRNAs"
    temp_MC = temp
    temp_MC$Confidence = "MC-lncRNAs"
    temp_HC = temp
    temp_HC$Confidence = "HC-lncRNAs"
    
    tab_ALL_collap_2_S_mod = tab_ALL_collap_2_S[tab_ALL_collap_2_S$Confidence != "PC genes", ]
    tab_ALL_collap_2_S_mod = rbind(tab_ALL_collap_2_S_mod, temp_LC, temp_MC, temp_HC)
    tab_ALL_collap_2_S_mod$Confidence = factor(tab_ALL_collap_2_S_mod$Confidence, levels = c("LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))
    
    TAB_2_mod = rbind(TAB_2_mod, tab_ALL_collap_2_S_mod)
    
    rm(list = c("tab_ALL_collap_2", "tab_ALL_collap_2_S", "tab_ALL_collap_2_S_mod",
                "temp", "temp_LC", "temp_MC", "temp_HC"))
    
    rm(list = c("tab_ALL"))
  }

  rm(list = c("i", "path_raw", "path_fig"))
}

rm(list = c("experiment", "Experiments"))

TAB_1 = merge(TAB_1, Summary_tab[, c("ID", "Class")], by = "ID", all.x = T, all.y = F)
TAB_1_mod = merge(TAB_1_mod, Summary_tab[, c("ID", "Class")], by = "ID", all.x = T, all.y = F)
TAB_2 = merge(TAB_2, Summary_tab[, c("ID", "Class")], by = "ID", all.x = T, all.y = F)
TAB_2_mod = merge(TAB_2_mod, Summary_tab[, c("ID", "Class")], by = "ID", all.x = T, all.y = F)

write.table(TAB_1, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Percentage_DE/TAB-Violin-Percentage_DE-1.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(TAB_1_mod, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Percentage_DE/TAB-Violin-Percentage_DE_mod-1.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(TAB_2, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Percentage_DE/TAB-Violin-Percentage_DE-2.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(TAB_2_mod, paste0(path_DEA, "/05-Tables_and_Figures/ALL/Percentage_DE/TAB-Violin-Percentage_DE_mod-2.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

TAB_1_mod$Class_code = factor(TAB_1_mod$Class_code, levels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs"))
TAB_1_mod$Confidence = factor(TAB_1_mod$Confidence, levels = c("LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs"))
TAB_1_mod$Class = factor(TAB_1_mod$Class, levels = c("Development", "Abiotic stress", "Biotic stress"))

gg = ggplot(TAB_1_mod, aes(x = Class_code, y = Percentage, fill = Class_code)) +
  geom_violin() +
  geom_jitter(size = 0.5) +
  scale_fill_manual(values = c("#7cc1cf", "#e5dd6c", "#e1ad60", "#da6d6d", "#89ab66")) +
  facet_grid(Confidence~Class) +
  xlab("") +
  ylab("Percentage") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 14, angle = 50),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 17, face = "bold"))

ggsave(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Percentage_DE/Violin-Percentage-1.png"), height = 12, width = 12, dpi = 800)
#ggsave(paste0(path_fig, "/05-Tables_and_Figures/ALL/Percentage_DE/Violin-Percentage-1.pdf"), height = 12, width = 12, dpi = 800)

rm(list = c("TAB_1", "TAB_1_mod", "TAB_2", "TAB_2_mod", "gg", "Summary_tab"))











############################################
## 4. VENNDIAGRAM
############################################

# Utilizando como total tanto transcritos DE como non-DE generamos un venn diagram.

cat(paste0("\n\n4. VENNDIAGRAM..."))

# Load binary table summary.
Binary_tab_SUMMARY_2_filt = read.table(paste0(path_DEA, "/05-Tables_and_Figures/ALL/Binary_tables/Binary_tab-SUMMARY_2_filt.tsv"), sep = "\t", header = T, quote = "\"")

# Load tissue-specificity info.
tab_TS = read.table(path_TS, sep = "\t", header = T, quote = "\"")
tab_TS = tab_TS[tab_TS$Spe == spe,]
rownames(tab_TS) = NULL
IDs_TS = unique(tab_TS[tab_TS$TAU >= 0.8, "ID_transcript"])

# Join tables.
tab_joined = merge(Binary_tab_SUMMARY_2_filt[,c("ID_transcript", "Class_code", "Confidence", "Norm.Development", "Norm.Stress")], db, by = c("ID_transcript", "Class_code", "Confidence"), all = T)
tab_joined[is.na(tab_joined)] = 0
tab_joined$Tissue.Specificity = ifelse(tab_joined$ID_transcript %in% IDs_TS, 1, 0)
tab_joined_long = melt(setDT(tab_joined), id.vars = c("ID_transcript", "Class_code", "Confidence"), variable.name = "Item")
tab_joined_long = as.data.frame(tab_joined_long)

write.table(tab_joined_long, paste0(path_DEA, "/05-Tables_and_Figures/ALL/VennDiagram/TAB-Venn.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

tab_joined_long = tab_joined_long[tab_joined_long$value != 0,]

write.table(tab_joined_long, paste0(path_DEA, "/05-Tables_and_Figures/ALL/VennDiagram/TAB-Venn_filt.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("Binary_tab_SUMMARY_2_filt", "tab_TS", "IDs_TS", "tab_joined"))

for (sl in c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs")) {
  x = list(
    Development = unique(tab_joined_long[tab_joined_long$Confidence == sl & tab_joined_long$Item == "Norm.Development", "ID_transcript"]),
    Stress = unique(tab_joined_long[tab_joined_long$Confidence == sl & tab_joined_long$Item == "Norm.Stress", "ID_transcript"]),
    `Tissue specificity` = unique(tab_joined_long[tab_joined_long$Confidence == sl & tab_joined_long$Item == "Tissue.Specificity", "ID_transcript"])
  )
  ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"), stroke_size = 0.5, set_name_size = 5, text_size = 4.2)
  ggsave(paste0(path_DEA, "/05-Tables_and_Figures/ALL/VennDiagram/Venn-", gsub(" ", "-", sl), ".png"), height = 6, width = 6, dpi = 800, bg = "white")
}

rm(list = c("sl", "x", "db", "tab_joined_long"))

