################################################################################
#
# EXPLORATORY ANALYSIS
#
# Considering each of the metadata tables previously as an experiment, this script 
# is used to perform an exploratory analysis of the samples present in each experiment.
#
# @author: pasviber - Pascual Villalba Bermell
# 
################################################################################


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
options(bitmapType='cairo')










################################################################################
################################################################################
################################## VARIABLES ###################################
################################################################################
################################################################################


## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied.", call.=FALSE)
} else {
  path_DEA = args[1]
  path_annot = args[2]
  path_quant = args[3]
}

# path_DEA = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/13-DEA/cme"
# path_annot = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP-FINAL/Files/ALL/nr"
# path_quant = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification/cme/ALL/nr/03-Quant"










################################################################################
################################################################################
############################ EXPLORATORY ANALYSIS ##############################
################################################################################
################################################################################


if (!dir.exists(paste0(path_DEA, "/02-EA"))){
  dir.create(paste0(path_DEA, "/02-EA"))
}
if (!dir.exists(paste0(path_DEA, "/03-Metadata_DEA"))){
  dir.create(paste0(path_DEA, "/03-Metadata_DEA"))
}

# Experiments.
Experiments = gsub(".tsv", "", list.files(paste0(path_DEA, "/01-Metadata_EA"), pattern = "[0-9].tsv"))

DF_test_MWN = data.frame()

for (experiment in Experiments) {
  
  cat(paste0("\n\nExperiment: ", experiment))
  
  ############################################
  ## 1. LOAD COUNTS
  ############################################
  
  cat(paste0("\n\tLoad counts..."))
  
  ## Load counts tables coming from salmon using tximport.
  
  txdb = suppressMessages(makeTxDbFromGFF(paste0(path_annot, "/ALL.gtf"), format = "gtf"))
  k = keys(txdb, keytype = "GENEID")
  df = suppressMessages(AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = columns(txdb)))
  tx2gene = df[, c("TXNAME", "GENEID")]
  tx2gene = tx2gene[!duplicated(tx2gene),]
  
  metadata = read.table(paste0(path_DEA, "/01-Metadata_EA/", experiment, ".tsv"), header = T, sep = "\t")
  files = file.path(path_quant, metadata$Sample, "quant.sf")
  names(files) = paste0(metadata$Sample)
  txi.salmon = suppressMessages(tximport(files, type = "salmon", tx2gene = tx2gene, txIn = TRUE, txOut = TRUE, countsFromAbundance = "no"))
  
  rm(list = c("txdb", "k", "df", "tx2gene", "files"))
  
  ############################################
  ## 2. CREATE THE DESEQ2 OBJECT
  ############################################
  
  cat(paste0("\n\tCreate the DESeq2 object..."))
  
  ## Create the DESeq2 object taking into account the class (DEVELOPMENT and, 
  ## ABIOTIC or BIOTIC STRESS). Use the column Name.2 as design (Variable with
  ## at least two conditions).
  
  if (unique(metadata$Class) == "Development") {
    ctrl = unique(metadata[metadata$Condition == "Control", "Name.2"])
    conds = unique(metadata[metadata$Condition == "Stage", "Name.2"])
    metadata$Name.2 = factor(metadata$Name.2, levels = c(ctrl, conds))
    dds = suppressMessages(DESeqDataSetFromTximport(txi = txi.salmon, colData = metadata, design = ~Name.2))
    dds$Name.2 = relevel(dds$Name.2, ref = ctrl)
  }
  if (unique(metadata$Class) == "Abiotic stress" | unique(metadata$Class) == "Biotic stress") {
    ctrl = unique(metadata[metadata$Condition == "Control", "Name.2"])
    conds = unique(metadata[metadata$Condition == "Treated", "Name.2"])
    metadata$Name.2 = factor(metadata$Name.2, levels = c(ctrl, conds))
    dds = suppressMessages(DESeqDataSetFromTximport(txi = txi.salmon, colData = metadata, design = ~Name.2))
    dds$Name.2 = relevel(dds$Name.2, ref = ctrl)
  }
  
  # Note: levels of factors in the design contain characters other than letters, 
  # numbers, '_' and '.'. It is recommended (but not required) to use only letters, 
  # numbers, and delimiters '_' or '.', as these are safe characters for column 
  # names in R. [This is a message, not a warning or an error] using counts and 
  # average transcript lengths from tximport.
  
  rm(list = c("txi.salmon", "ctrl", "conds"))
  
  ############################################
  ## 3. FILTER COUNTS
  ############################################
  
  cat(paste0("\n\tFilter counts..."))
  
  ## Filter the counts using the estimated counts loaded using tximport R package.
  ## Keep only PC genes or lncRNAs with at least 1 count in at least two samples
  ## by condition. We have decided to use this soft filter due to 
  ##     - The low expression of lncRNAs.
  ##     - The independent filter: You typically don't need to pre-filter because 
  ##       independent filtering occurs within results() to save you from multiple 
  ##       test correction on transcripts with no power.
  
  cat(paste0("\n\t\tNumber of transcripts pre-filtering: ", nrow(dds)))
  
  KEEP = c()
  for (name in unique(metadata$Name.2)) {
    cols = metadata[metadata$Name.2 == name, "Sample"]
    dt = counts(dds)[, cols]
    keep = rowSums(dt > 0) >= 2
    keep_name = rownames(dt)[keep]
    KEEP = unique(c(KEEP, keep))
  }
  
  dds = dds[KEEP,]
  
  cat(paste0("\n\t\tNumber of transcripts post-filtering: ", nrow(dds)))
  
  rm(list = c("cols", "dt", "keep", "keep_name", "KEEP", "name"))
  
  ############################################
  ## 4. CALCULATE SIZE FACTORS 
  ############################################
  
  cat(paste0("\n\tCalculate size factors..."))
  
  ## We see that the larger size factors correspond to the samples with higher 
  ## sequencing depth, which makes sense, because to generate our normalized 
  ## counts we need to divide the counts by the size factors. This accounts for 
  ## the differences in sequencing depth between samples.
  
  dds = suppressMessages(DESeq2::estimateSizeFactors(dds))
  
  ############################################
  ## 5. VARIANCE STABILIZING TRANSFORMATION (VST)
  ############################################
  
  cat(paste0("\n\tVariance stabilizing transformation..."))
  
  ## RNA-Seq workflow: gene-level exploratory analysis and differential expression
  ## (F1000Research)
  
  ## This transformation produces transformed data on the log2 scale which has 
  ## been normalized with respect to library size or other normalization factors.
  
  ## We specify blind=FALSE, which means that differences across donors and treatment 
  ## should not add to the variance-mean profile of the experiment. However, the 
  ## experimental design is not used directly in the transformation, only in estimating 
  ## the global amount of variability in the counts.
  
  assay(dds, 'counts.norm.VST') = suppressMessages(as.data.frame(assay(varianceStabilizingTransformation(dds, blind = FALSE))))
  
  ############################################
  ## 6. MEAN VS VARIANCE PLOT
  ############################################
  
  cat(paste0("\n\tMean vs variance plot..."))
  
  ## Many common statistical methods for exploratory analysis of multidimensional 
  ## data work best for data that generally has the same range of variance at 
  ## different ranges of the mean values. For RNA-seq raw counts, however, the 
  ## variance grows with the mean.
  
  ## We want to check if the variance stabilizing transformations worked well.
  
  # This function calculates the variance of each row in a matrix or data frame
  # "a" using the R function "apply".
  rowVar = function(a){apply(a,1,var)}
  
  # Plot.
  png(file = paste(path_DEA, "/02-EA/", experiment, '_meanvsvar.png', sep = ''), 
      width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)
  par(mfrow=c(1,2))
  
  plot(log2(rowMeans(assay(dds,'counts'))+1), log2(rowVar(assay(dds, 'counts'))+1),
       xlab=expression('Log'[2]~'Mean count'), ylab=expression('Log'[2]~'Variance'), main='Raw Counts')
  
  plot(rowMeans(assay(dds,'counts.norm.VST')), rowVar(assay(dds, 'counts.norm.VST')),
       xlab='Mean count', ylab='Variance', main='VST')
  
  invisible(dev.off())
  
  rm(list = c("rowVar"))
  
  ############################################
  ## 7. HEATMAP 
  ############################################
  
  cat(paste0("\n\tHeatmap..."))
  
  ## A useful step in an RNA-seq analysis is often to assess overall similarity 
  ## between samples: Which samples are similar to each other, which are different? 
  ## Does this fit to the expectation from the experiment’s design?
  
  # Creates a matrix in dataframe format with the Euclidean distances between samples.
  sampleDists = dist(t(assay(dds,'counts.norm.VST')))
  sampleDistDF = as.data.frame(as.matrix(sampleDists))
  
  # Modify colnames and rownames.
  conversion = metadata[, c("Sample", "Name.1")]
  sampleDistDF$"Sample" = rownames(sampleDistDF)
  sampleDistDF = merge(sampleDistDF, conversion, by = "Sample", all = T)
  sampleDistDF$Sample = NULL
  rownames(sampleDistDF) = sampleDistDF$Name.1
  sampleDistDF$Name.1 = NULL
  colnames(sampleDistDF) = NULL
  
  # Colors.
  colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  # Plot.
  png(file = paste(path_DEA, "/02-EA/", experiment, '_heatmap.png', sep = ''), 
      width = 7, height = 4, units = "in", res = 1200)
  
  pheatmap(sampleDistDF, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, 
           col=colors, fontsize_row = 3)
  
  invisible(dev.off())
  
  rm(list = c("sampleDists", "sampleDistDF", "conversion", "colors"))
  
  ############################################
  ## 8. PRINCIPAL COMPONENT ANALYSIS (PCA)
  ############################################
  
  cat(paste0("\n\tPCA..."))
  
  ## Another way to visualize sample-to-sample distances is a principal components 
  ## analysis (PCA). In this ordination method, the data points (here, the samples) 
  ## are projected onto the 2D plane such that they spread out in the two directions 
  ## that capture the most of the variance across samples.
  
  #############
  ## 8.1 PLOT
  
  cat(paste0("\n\t\tPlot..."))
  
  # Create dataframe.
  pca_res = prcomp(x=t(assay(dds,'counts.norm.VST')), rank. = 6)
  pca_df = as.data.frame(pca_res$x)
  
  # Create the conditions column for coloring the plot.
  conversion = metadata[, c("Sample", "Name.2")]
  pca_df$"Sample" = rownames(pca_df)
  pca_df = merge(pca_df, conversion, by = "Sample", all = T)
  Conditions = pca_df$Name.2
  pca_df$Sample = NULL
  
  # Create interactive plot.
  plot = plot_ly(pca_df, x = pca_df[,1], y = pca_df[,2], z = pca_df[,3], type='scatter3d', mode='markers', color=~pca_df$Name.2, colors = 'Paired') %>% 
    layout(scene = list(xaxis = list(title = 'PC1'), yaxis = list(title = 'PC2'), zaxis = list(title = 'PC3')), showlegend = TRUE, legend = list(font = list(size = 20)))
  
  # Create html file with interactive plot.
  htmlwidgets::saveWidget(widget = plot, file = paste(path_DEA, "/02-EA/", experiment, '_PCA.html', sep=''), selfcontained = FALSE)
  
  #############
  ## 8.2 MANN-WHITNEY-WILCOXON TEST
  
  cat(paste0("\n\t\tMann-whitney wilcoxon test..."))
  
  # Create Null matrix.
  pca_comp = pca_res$x
  distance_matrix = matrix(data = rep(0, nrow(pca_comp) * nrow(pca_comp)), nrow = nrow(pca_comp), ncol = nrow(pca_comp))
  rownames(distance_matrix) = rownames(pca_comp)
  colnames(distance_matrix) = rownames(pca_comp)
  
  # INTRA- and INTER-group distance vectors.
  intra_group = c()
  inter_group = c()
  
  # This function calculates the Euclidean distance between two points.
  euclideanDist = function(a, b){return(sqrt(sum((a - b) ^ 2)))}
  
  # Complete matrix and vectors with the corresponding distances.
  for (i in 1:nrow(pca_comp)) {
    for (j in 1:nrow(pca_comp)){
      # Save distances in matrix
      distance_matrix[i,j] = euclideanDist(pca_comp[i,][1:3], pca_comp[j,][1:3]);
      # Prevent the presence of repeated distances in the vectors.
      if (j > i){
        # INTRA-group distances
        if (Conditions[i] == Conditions[j]) {
          intra_group = c(intra_group, euclideanDist(pca_comp[i,][1:3], pca_comp[j,][1:3]))
        }
        # INTER-group distances
        else{
          inter_group = c(inter_group, euclideanDist(pca_comp[i,][1:3], pca_comp[j,][1:3]))
        }
      }
    }
  }
  
  # Execute wilcoxon test.
  mww_res = wilcox.test(inter_group, intra_group, paired = FALSE)
  p_value = mww_res[3][[1]]
  
  #############
  ## 8.3 PROPORTION OF VARIANCE EXPLAINED BY EACH COMPONENT 
  
  cat(paste0("\n\t\tProportion of variance explained by each component..."))
  
  # Plot.
  png(file = paste(path_DEA, "/02-EA/", experiment, '_PCA_variance_by_component.png', sep=''), 
      width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)
  
  data = head(round(pca_res$sdev ^ 2 / sum(pca_res$sdev ^ 2) * 100,2), n=6)
  plot = barplot(data, ylim = c(0, max(data) * 1.2), las=2, names.arg=colnames(pca_res$x), ylab='% variance explained')
  text(x = plot, y = data, labels = data, pos = 3)
  
  invisible(dev.off())
  
  n_ceros = 6 - length(data)
  ceros = rep(0.00, n_ceros)
  data = c(data, ceros)
  
  # Create a table with the results.
  df_test_MWN = data.frame(Specie = unique(metadata$Specie), Project = experiment, PC1 = data[1], 
                           PC2 = data[2], PC3 = data[3], PC4 = data[4],  PC5 = data[5], 
                           PC6 = data[6], 'P.value' = p_value)
  DF_test_MWN = rbind(DF_test_MWN, df_test_MWN)
  
  ## Filter of MWN test.
  if (p_value < 0.05) {
    cat("\n\tThis experiment passes the filter. It has a p.value < 0.05.")
  }else {
    cat("\n\tThis experiment doesn't pass the filter. It has a p.value > 0.05.")
  }
  
  ## Save the metadata in the 03-Metadata_DEA folder.
  write.table(metadata, paste0(path_DEA, "/03-Metadata_DEA/", experiment, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(list = c("dds", "pca_res", "pca_df", "conversion", "Conditions", "plot", "pca_comp", "distance_matrix", 
              "intra_group", "inter_group", "euclideanDist", "i", "j", "mww_res", "p_value", "data", "n_ceros",
              "ceros", "df_test_MWN"))
}

write.table(DF_test_MWN, paste0(path_DEA, "/02-EA/Summary_table_MWN.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("metadata", "experiment", "Experiments"))

cat("\n\n\n--------------------------------------------------------------------")

