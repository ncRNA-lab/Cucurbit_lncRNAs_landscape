
######### MODULES

suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(DESeq2))


######### VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop("At least 7 arguments must be supplied.", call.=FALSE)
} else {
  experiment = args[1]
  WD_corr_S2 = args[2]
  WD_corr_S3 = args[3]
  WD_pred = args[4]
  WD_quant = args[5]
  WD_DEA = args[6]
  flag = args[7]
}


######### PIPELINE

## STEP 3: Create table with pearson correlation.

## Load TAB_CIS.
TAB_CIS = read.table(paste0(WD_corr_S2, "/TAB_CIS.tsv"), sep = "\t", header = T, quote = "\"")
## Load summary table with all the information about the experiments and contrasts.
metadata = read.table(paste0(WD_DEA, "/03-Metadata_DEA/", experiment, ".tsv"), sep = "\t", header = T, quote = "\"")
## Create subset table.
subset = TAB_CIS
subset$Experiment = experiment

if (dim(subset)[1] != 0) {
  
  ## Load counts tables coming from salmon using tximport.
  txdb = suppressMessages(makeTxDbFromGFF(paste0(WD_pred, "/STEP-FINAL/Files/Joined/ALL/", flag, "/ALL.gtf"), format = "gtf"))
  k = keys(txdb, keytype = "GENEID")
  df = suppressMessages(AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = columns(txdb)))
  tx2gene = df[, c("TXNAME", "GENEID")]
  tx2gene = tx2gene[!duplicated(tx2gene),]
  
  files = file.path(paste0(WD_quant, "/Salmon/ALL/", flag, "/03-Quant/"), metadata$Sample, "quant.sf")
  names(files) = paste0(metadata$Sample)
  txi.salmon = suppressMessages(tximport(files, type = "salmon", tx2gene = tx2gene, txIn = TRUE, txOut = TRUE, countsFromAbundance = "no"))
  
  rm(list = c("txdb", "k", "df", "tx2gene", "files"))
  
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
  
  rm(list = c("ctrl", "conds", "txi.salmon"))
  
  ## Filter the counts using the estimated counts loaded using tximport R package.
  dds_filt = dds[apply(counts(dds), 1, mad) > quantile(apply(counts(dds), 1, mad), 0.75),]
  
  rm(list = c("dds"))
  
  ## Calculate size factors.
  dds = suppressMessages(DESeq2::estimateSizeFactors(dds_filt))
  
  rm(list = c("dds_filt"))
  
  ## Variance stabilizing transformation. This transformation produces transformed 
  ## data on the log2 scale which has been normalized with respect to library 
  ## size or other normalization factors.
  assay(dds, 'counts.norm.VST') = suppressMessages(as.data.frame(assay(varianceStabilizingTransformation(dds, blind = FALSE))))
  
  ## Correlation method
  vst = as.data.frame(assay(dds, 'counts.norm.VST'))
  vst_filt = vst[rownames(vst) %in% unique(c(subset$ID_transcript.1, subset$ID_transcript.2)),]
  tab_cor = suppressWarnings(cor(t(vst_filt), method = "pearson"))
  tab_cor = melt(tab_cor)
  colnames(tab_cor) = c("ID_transcript.1", "ID_transcript.2", "Cor.Deseq.Vst")
  
  rm(list = c("metadata", "dds", "vst", "vst_filt"))
  
  ## Obtain Correlation index data.
  subset = merge(subset, tab_cor, by = c("ID_transcript.1", "ID_transcript.2"), all.x = F, all.y = F)
  
  rm(list = c("tab_cor"))
  
  ## Generate the final table.
  TAB_CIS_corr =  subset[, c("Specie", "Chr", "ID_transcript.1", "Start.1", "End.1", "Strand.1", "ID_transcript.2", 
                             "Start.2", "End.2", "Strand.2", "Overlap.Perc", "Overlap", "Type", "Experiment", 
                             "Cor.Deseq.Vst")]
  rownames(TAB_CIS_corr) = NULL
  
  ## Save.
  write.table(TAB_CIS_corr, paste0(WD_corr_S3, "/TAB_CIS-", experiment, "-PEARSON.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
  
}else {
  cat(paste0("\tThere is no pairs\n"))
}

rm(list = c("TAB_CIS", "experiment", "subset"))

