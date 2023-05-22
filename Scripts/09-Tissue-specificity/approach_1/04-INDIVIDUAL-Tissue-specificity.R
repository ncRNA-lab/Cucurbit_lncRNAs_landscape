################################################################################
#
# INDIVIDUAL: TISSUE SPECIFICITY STUDY: APPROACH 1 - STEP 4
#
# Draw heatmaps using lncRNAs with TAU > 0.8.
# 
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggExtra))
suppressMessages(library(data.table))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(circlize))
suppressMessages(options(bitmapType='cairo'))

## 1. PATHS

# Own computer
path_tissue_specificity = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
path_quant = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/06-quantification"
flag = "nr"
species = c("car", "cla", "cma", "cme", "cmo", "cpe", "csa", "lsi", "mch")

# # Garnatxa
# path_tissue_specificity = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
# path_quant = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-quantification"
# flag = "nr"
# species = c("car", "cla", "cma", "cme", "cmo", "cpe", "csa", "lsi", "mch")


## 2. MEAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEAN TAU: \n"))

TAB_mean_TAU_filt_FINAL = data.frame()
for (spe in species) {
  cat(paste0("------ Specie: ", spe, "\n"))
  
  # Directories.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4"))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4"))
  }
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe))
  }
  
  # Get files with TAU metric values and log-transformed expression values.
  files = list.files(path = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe), pattern = ".tsv")
  
  if (length(files) > 0) {
    SRA.Studies = unique(unlist(lapply(strsplit(unique(unlist(lapply(strsplit(files, "_"), `[[`, 1))), "-"), `[[`, 2)))
    
    for (SRA.Study in SRA.Studies) {
      cat(paste0("--------- SRA.Study: ", SRA.Study, "\n"))
      
      # Load table with the tissue specificity results and the log-tranformed expression values. 
      tab_mean_TAU_G = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/GENES-", SRA.Study, "_mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU_L = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/LNCRNAS-", SRA.Study, "_mean-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_mean_TAU = rbind(tab_mean_TAU_G, tab_mean_TAU_L)
      
      # Convert to factors.
      tab_mean_TAU$Confidence = factor(tab_mean_TAU$Confidence, levels = c("Protein-coding", "Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
      tab_mean_TAU$Class_code = factor(tab_mean_TAU$Class_code, levels = c("pc", "u", "x", "i", "o/e"))
      
      # Heatmap plot with TAU >= 0.8 lncRNAs and without genes.
      tab_mean_TAU_filt = tab_mean_TAU[tab_mean_TAU$TAU >= 0.8 & tab_mean_TAU$Class_code != "pc", ]
      rownames(tab_mean_TAU_filt) = tab_mean_TAU_filt$ID_transcript
      tab_mean_TAU_filt$ID_transcript = NULL
      mat = tab_mean_TAU_filt[, c(colnames(tab_mean_TAU)[grepl("EXP.", colnames(tab_mean_TAU), fixed = T)])]
      colnames(mat) = gsub("EXP.", "", colnames(mat))
      heat = t(scale(t(mat)))
      
      myCol = colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
      myBreaks = seq(-3, 3, length.out = 100)
      zscore_col_fun = colorRamp2(myBreaks, myCol)
      
      hmap = Heatmap(heat, 
                     
                     # Rows
                     cluster_rows = TRUE,
                     show_row_dend = FALSE,
                     row_title = 'Transcripts TAU > 0.8',
                     row_title_side = 'right',
                     row_title_gp = gpar(fontsize = 15,  fontface = 'bold'),
                     row_title_rot = 90,
                     show_row_names = FALSE,
                     row_dend_width = unit(25,'mm'),
                     
                     # Columns
                     cluster_columns = FALSE,
                     show_column_dend = FALSE,
                     column_title = 'Tissues',
                     column_title_side = 'bottom',
                     column_title_gp = gpar(fontsize = 15, fontface = 'bold'),
                     column_title_rot = 0,
                     show_column_names = TRUE,
                     column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                     column_names_max_height = unit(10, 'cm'),
                     column_dend_height = unit(25,'mm'),
                     
                     name = 'Z-score',
                     col = zscore_col_fun,
                     
                     heatmap_legend_param = list(
                       color_bar = 'continuous',
                       legend_direction = 'horizontal',
                       legend_width = unit(8, 'cm'),
                       legend_height = unit(5.0, 'cm'),
                       title_position = 'topcenter',
                       title_gp=gpar(fontsize = 10, fontface = 'bold'),
                       labels_gp=gpar(fontsize = 9))
                     )
      
      png(filename = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe, "/", SRA.Study, "_mean-TAU-Heatmap_without_pc.png"), height = 7000, width = 4000, res = 800)
      draw(hmap,
           heatmap_legend_side = 'top',
           annotation_legend_side = 'right',
           row_sub_title_side = 'left')
      invisible(dev.off())
      
      rm(list = c("mat", "heat", "myCol", "myBreaks", "hmap", "zscore_col_fun", "tab_mean_TAU_filt"))
      
      # Heatmap plot with TAU >= 0.8 lncRNAs and with TAU >= 0.8 genes.
      tab_mean_TAU_filt = tab_mean_TAU[tab_mean_TAU$TAU >= 0.8, ]
      rownames(tab_mean_TAU_filt) = tab_mean_TAU_filt$ID_transcript
      tab_mean_TAU_filt$ID_transcript = NULL
      mat = tab_mean_TAU_filt[, c(colnames(tab_mean_TAU)[grepl("EXP.", colnames(tab_mean_TAU), fixed = T)])]
      colnames(mat) = gsub("EXP.", "", colnames(mat))
      heat = t(scale(t(mat)))
      
      myCol = colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
      myBreaks = seq(-3, 3, length.out = 100)
      zscore_col_fun = colorRamp2(myBreaks, myCol)
      
      hmap = Heatmap(heat,
                     
                     # Rows
                     cluster_rows = TRUE,
                     show_row_dend = FALSE,
                     row_title = 'Transcripts TAU > 0.8',
                     row_title_side = 'right',
                     row_title_gp = gpar(fontsize = 15,  fontface = 'bold'),
                     row_title_rot = 90,
                     show_row_names = FALSE,
                     row_dend_width = unit(25,'mm'),
                     
                     # Columns
                     cluster_columns = FALSE,
                     show_column_dend = FALSE,
                     column_title = 'Tissues',
                     column_title_side = 'bottom',
                     column_title_gp = gpar(fontsize = 15, fontface = 'bold'),
                     column_title_rot = 0,
                     show_column_names = TRUE,
                     column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                     column_names_max_height = unit(10, 'cm'),
                     column_dend_height = unit(25,'mm'),
                     
                     name = 'Z-score',
                     col = zscore_col_fun,
                     
                     heatmap_legend_param = list(
                       color_bar = 'continuous',
                       legend_direction = 'horizontal',
                       legend_width = unit(8, 'cm'),
                       legend_height = unit(5.0, 'cm'),
                       title_position = 'topcenter',
                       title_gp=gpar(fontsize = 10, fontface = 'bold'),
                       labels_gp=gpar(fontsize = 9))
                     )
      
      png(filename = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe, "/", SRA.Study, "_mean-TAU-Heatmap_with_pc.png"), height = 7000, width = 4000, res = 800)
      draw(hmap,
           heatmap_legend_side = 'top',
           annotation_legend_side = 'right',
           row_sub_title_side = 'left')
      invisible(dev.off())
      
      rm(list = c("mat", "heat", "myCol", "myBreaks", "hmap", "zscore_col_fun", "tab_mean_TAU_filt"))
      
      # Table TAU > 0.8
      tab_mean_TAU_filt = tab_mean_TAU[, c("ID_transcript", "Confidence", "Class_code", "TAU")]
      tab_mean_TAU_filt$"TAU>0.8" = ifelse(tab_mean_TAU_filt$TAU >= 0.8, "YES", "NO")
      tab_mean_TAU_filt$`TAU>0.8` = factor(tab_mean_TAU_filt$`TAU>0.8`, levels = c("YES", "NO"))
      
      tab_mean_TAU_filt_FINAL = tab_mean_TAU_filt %>% 
        group_by(Confidence, Class_code, `TAU>0.8`, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(ID_transcript)) %>%
        mutate(Total = sum(Counts),
               Perc = round(Counts/sum(Counts) * 100, 2))
      tab_mean_TAU_filt_FINAL = as.data.frame(tab_mean_TAU_filt_FINAL)
      
      tab_mean_TAU_filt_FINAL = tab_mean_TAU_filt_FINAL[!is.na(tab_mean_TAU_filt_FINAL$Perc),]
      
      write.table(tab_mean_TAU_filt_FINAL, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe, "/", SRA.Study, "_mean-0.8.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      
      tab_mean_TAU_filt_FINAL$"Spe" = spe
      tab_mean_TAU_filt_FINAL$"SRA.Study" = SRA.Study
      tab_mean_TAU_filt_FINAL = tab_mean_TAU_filt_FINAL[,c("Spe", "SRA.Study", "Confidence", "Class_code", "TAU>0.8", "Total", "Counts", "Perc")]
      TAB_mean_TAU_filt_FINAL = rbind(TAB_mean_TAU_filt_FINAL, tab_mean_TAU_filt_FINAL)
      
      rm(list = c("tab_mean_TAU_filt", "tab_mean_TAU_filt_FINAL"))
    }
  }
}

write.table(TAB_mean_TAU_filt_FINAL, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/mean-0.8.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

rm(list = c("files", "spe", "SRA.Study", "SRA.Studies", "TAB_mean_TAU_filt_FINAL"))


## 3. MEDIAN TAU - TABLES AND FIGURES

cat(paste0("\n\n--- MEDIAN TAU: \n"))

TAB_median_TAU_filt_FINAL = data.frame()
for (spe in species) {
  cat(paste0("------ Specie: ", spe, "\n"))
  
  # Directories.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4"))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4"))
  }
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe))
  }
  
  # Get files with TAU metric values and log-transformed expression values.
  files = list.files(path = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe), pattern = ".tsv")
  
  if (length(files) > 0) {
    SRA.Studies = unique(unlist(lapply(strsplit(unique(unlist(lapply(strsplit(files, "_"), `[[`, 1))), "-"), `[[`, 2)))
    
    for (SRA.Study in SRA.Studies) {
      cat(paste0("--------- SRA.Study: ", SRA.Study, "\n"))
      
      # Load table with the tissue specificity results and the log-tranformed expression values. 
      tab_median_TAU_G = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/GENES-", SRA.Study, "_median-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_median_TAU_L = read.table(paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP2/", spe, "/LNCRNAS-", SRA.Study, "_median-TAU.tsv"), sep = "\t", header = T, quote = "\"")
      tab_median_TAU = rbind(tab_median_TAU_G, tab_median_TAU_L)
      
      # Convert to factors.
      tab_median_TAU$Confidence = factor(tab_median_TAU$Confidence, levels = c("Protein-coding", "Low-confidence lncRNA", "Medium-confidence lncRNA", "High-confidence lncRNA"))
      tab_median_TAU$Class_code = factor(tab_median_TAU$Class_code, levels = c("pc", "u", "x", "i", "o/e"))
      
      # Heatmap plot with TAU >= 0.8 lncRNAs and without genes.
      tab_median_TAU_filt = tab_median_TAU[tab_median_TAU$TAU >= 0.8 & tab_median_TAU$Class_code != "pc", ]
      rownames(tab_median_TAU_filt) = tab_median_TAU_filt$ID_transcript
      tab_median_TAU_filt$ID_transcript = NULL
      mat = tab_median_TAU_filt[, c(colnames(tab_median_TAU)[grepl("EXP.", colnames(tab_median_TAU), fixed = T)])]
      colnames(mat) = gsub("EXP.", "", colnames(mat))
      heat = t(scale(t(mat)))
      
      myCol = colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
      myBreaks = seq(-3, 3, length.out = 100)
      zscore_col_fun = colorRamp2(myBreaks, myCol)
      
      hmap = Heatmap(heat, 
                     
                     # Rows
                     cluster_rows = TRUE,
                     show_row_dend = FALSE,
                     row_title = 'Transcripts TAU > 0.8',
                     row_title_side = 'right',
                     row_title_gp = gpar(fontsize = 15,  fontface = 'bold'),
                     row_title_rot = 90,
                     show_row_names = FALSE,
                     row_dend_width = unit(25,'mm'),
                     
                     # Columns
                     cluster_columns = FALSE,
                     show_column_dend = FALSE,
                     column_title = 'Tissues',
                     column_title_side = 'bottom',
                     column_title_gp = gpar(fontsize = 15, fontface = 'bold'),
                     column_title_rot = 0,
                     show_column_names = TRUE,
                     column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                     column_names_max_height = unit(10, 'cm'),
                     column_dend_height = unit(25,'mm'),
                     
                     name = 'Z-score',
                     col = zscore_col_fun,
                     
                     heatmap_legend_param = list(
                       color_bar = 'continuous',
                       legend_direction = 'horizontal',
                       legend_width = unit(8, 'cm'),
                       legend_height = unit(5.0, 'cm'),
                       title_position = 'topcenter',
                       title_gp=gpar(fontsize = 10, fontface = 'bold'),
                       labels_gp=gpar(fontsize = 9))
      )
      
      png(filename = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe, "/", SRA.Study, "_median-TAU-Heatmap_without_pc.png"), height = 7000, width = 4000, res = 800)
      draw(hmap,
           heatmap_legend_side = 'top',
           annotation_legend_side = 'right',
           row_sub_title_side = 'left')
      invisible(dev.off())
      
      rm(list = c("mat", "heat", "myCol", "myBreaks", "hmap", "zscore_col_fun", "tab_median_TAU_filt"))
      
      # Heatmap plot with TAU >= 0.8 lncRNAs and with TAU >= 0.8 genes.
      tab_median_TAU_filt = tab_median_TAU[tab_median_TAU$TAU >= 0.8, ]
      rownames(tab_median_TAU_filt) = tab_median_TAU_filt$ID_transcript
      tab_median_TAU_filt$ID_transcript = NULL
      mat = tab_median_TAU_filt[, c(colnames(tab_median_TAU)[grepl("EXP.", colnames(tab_median_TAU), fixed = T)])]
      colnames(mat) = gsub("EXP.", "", colnames(mat))
      heat = t(scale(t(mat)))
      
      myCol = colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
      myBreaks = seq(-3, 3, length.out = 100)
      zscore_col_fun = colorRamp2(myBreaks, myCol)
      
      hmap = Heatmap(heat, 
                     
                     # Rows
                     cluster_rows = TRUE,
                     show_row_dend = FALSE,
                     row_title = 'Transcripts TAU > 0.8',
                     row_title_side = 'right',
                     row_title_gp = gpar(fontsize = 15,  fontface = 'bold'),
                     row_title_rot = 90,
                     show_row_names = FALSE,
                     row_dend_width = unit(25,'mm'),
                     
                     # Columns
                     cluster_columns = FALSE,
                     show_column_dend = FALSE,
                     column_title = 'Tissues',
                     column_title_side = 'bottom',
                     column_title_gp = gpar(fontsize = 15, fontface = 'bold'),
                     column_title_rot = 0,
                     show_column_names = TRUE,
                     column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                     column_names_max_height = unit(10, 'cm'),
                     column_dend_height = unit(25,'mm'),
                     
                     name = 'Z-score',
                     col = zscore_col_fun,
                     
                     heatmap_legend_param = list(
                       color_bar = 'continuous',
                       legend_direction = 'horizontal',
                       legend_width = unit(8, 'cm'),
                       legend_height = unit(5.0, 'cm'),
                       title_position = 'topcenter',
                       title_gp=gpar(fontsize = 10, fontface = 'bold'),
                       labels_gp=gpar(fontsize = 10, fontface = 'bold'))
      )
      
      png(filename = paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe, "/", SRA.Study, "_median-TAU-Heatmap_with_pc.png"), height = 7000, width = 4000, res = 800)
      draw(hmap,
           heatmap_legend_side = 'top',
           annotation_legend_side = 'right',
           row_sub_title_side = 'left')
      invisible(dev.off())
      
      rm(list = c("mat", "heat", "myCol", "myBreaks", "hmap", "zscore_col_fun", "tab_median_TAU_filt"))
      
      # Table TAU > 0.8
      tab_median_TAU_filt = tab_median_TAU[, c("ID_transcript", "Confidence", "Class_code", "TAU")]
      tab_median_TAU_filt$"TAU>0.8" = ifelse(tab_median_TAU_filt$TAU >= 0.8, "YES", "NO")
      tab_median_TAU_filt$`TAU>0.8` = factor(tab_median_TAU_filt$`TAU>0.8`, levels = c("YES", "NO"))
      
      tab_median_TAU_filt_FINAL = tab_median_TAU_filt %>% 
        group_by(Confidence, Class_code, `TAU>0.8`, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(ID_transcript)) %>%
        mutate(Total = sum(Counts),
               Perc = round(Counts/sum(Counts) * 100, 2))
      tab_median_TAU_filt_FINAL = as.data.frame(tab_median_TAU_filt_FINAL)
      
      tab_median_TAU_filt_FINAL = tab_median_TAU_filt_FINAL[!is.na(tab_median_TAU_filt_FINAL$Perc),]
      
      write.table(tab_median_TAU_filt_FINAL, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/", spe, "/", SRA.Study, "_median-0.8.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      
      tab_median_TAU_filt_FINAL$"Spe" = spe
      tab_median_TAU_filt_FINAL$"SRA.Study" = SRA.Study
      tab_median_TAU_filt_FINAL = tab_median_TAU_filt_FINAL[,c("Spe", "SRA.Study", "Confidence", "Class_code", "TAU>0.8", "Total", "Counts", "Perc")]
      TAB_median_TAU_filt_FINAL = rbind(TAB_median_TAU_filt_FINAL, tab_median_TAU_filt_FINAL)
      
      rm(list = c("tab_median_TAU_filt", "tab_median_TAU_filt_FINAL"))
    }
  }
}

write.table(TAB_median_TAU_filt_FINAL, paste0(path_tissue_specificity, "/approach_1/INDIVIDUAL/", flag, "/STEP4/median-0.8.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

rm(list = c("files", "spe", "SRA.Study", "SRA.Studies", "TAB_median_TAU_filt_FINAL"))

