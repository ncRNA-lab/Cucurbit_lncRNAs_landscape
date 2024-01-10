


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
if (length(args) < 4) {
  stop("At least 4 arguments must be supplied.", call.=FALSE)
} else {
  path_DEA = args[1]
  path_annot = args[2]
  path_quant = args[3]
  alpha_value = as.numeric(args[4])
}

# path_DEA = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/16-DEA/cla"
# path_annot = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cla/STEP-FINAL/Files/Joined/ALL/nr"
# path_quant = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/06-quantification/cla/Salmon/ALL/nr/03-Quant"
# alpha_value = 0.05









################################################################################
################################################################################
##################################### DEA ######################################
################################################################################
################################################################################


cat(paste0("DEA..."))

if (!dir.exists(paste0(path_DEA, "/04-DEA"))){
  dir.create(paste0(path_DEA, "/04-DEA"))
}

# Experiments.
Experiments = gsub(".tsv", "", list.files(paste0(path_DEA, "/03-Metadata_DEA"), pattern = "[0-9].tsv"))

SUMMARY = data.frame()

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
  
  metadata = read.table(paste0(path_DEA, "/03-Metadata_DEA/", experiment, ".tsv"), header = T, sep = "\t")
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
  ## 4. DEA
  ############################################
  
  cat(paste0("\n\tDEA..."))
  
  # Paths.
  path_raw_out = paste0(path_DEA, "/04-DEA/01-DEA_raw/", experiment)
  path_sig_out = paste0(path_DEA, "/04-DEA/02-DEA_sig/", experiment)
  path_nosig_out = paste0(path_DEA, "/04-DEA/03-DEA_nosig/", experiment)
  
  # Create output directories.
  if (!dir.exists(path_raw_out)) {
    dir.create(path_raw_out, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path_sig_out)) {
    dir.create(path_sig_out, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path_nosig_out)) {
    dir.create(path_nosig_out, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Differential expression analysis.
  dds = suppressMessages(DESeq(dds))
  
  # Save results.
  results_names = resultsNames(dds)
  i = 1
  summary = data.frame()
  for (contrast in results_names){
    if (contrast != 'Intercept'){
      
      cat(paste0("\n\tContrast: ", i))
      
      # Get results
      res_deseq = results(dds, name = contrast, alpha = alpha_value)
      resLFC = suppressMessages(lfcShrink(dds, coef = contrast, res = res_deseq))
      
      # Add Shrunken LFC to results
      final_res = na.omit(cbind(res_deseq, Shrunkenlog2FoldChange = resLFC$log2FoldChange, ShrunkenlfcSE = resLFC$lfcSE))
      final_res = final_res %>% 
        data.frame() %>%
        rownames_to_column(var = "ID_transcript")
      
      # Extract significant differentially expressed transcripts.
      final_res_sig = final_res %>% 
        as_tibble() %>%
        dplyr::filter(padj <= alpha_value) %>%
        data.frame()
      
      # Extract non-significant differentially expressed transcripts.
      final_res_nosig = final_res %>% 
        as_tibble() %>%
        dplyr::filter(padj > alpha_value) %>%
        data.frame()
      
      # Save results.
      write.table(final_res, paste0(path_raw_out, "/", experiment, "-dea_raw-", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
      write.table(final_res_sig, paste0(path_sig_out, "/", experiment, "-dea_sig-", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
      write.table(final_res_nosig, paste0(path_nosig_out, "/", experiment, "-dea_nosig-", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
      
      # Add data to dataframe. When adding the columns of the table with _info there 
      # may be characteristics (Tissue, Cultivar...) that we did not want to compare 
      # in the differential expression analysis but that appear in the columns with 
      # _info. As we want to keep this information in the SUMMARY table and we do not 
      # want to duplicate or triplicate rows, paste is used to join the values of 
      # these columns. But normally it is not necessary.
      sum = data.frame(Specie = unique(metadata$Specie), 
                       Class = unique(metadata$Class), 
                       SRA_study = unlist(strsplit(experiment, "_"))[1],
                       Experiment = experiment, 
                       n_contrast = i, 
                       contrast = gsub("Name.2_", "", contrast), 
                       Stress = gsub("_vs_Control", "", unlist(strsplit(gsub("Name.2_", "", contrast), "__"))[5]), 
                       Tissue = paste(gsub(" ", "_", unique(metadata$Tissue_info)), collapse = "__"), # It's not necessary 
                       Technique = paste(unique(metadata$Technique_info), collapse = "__"),  
                       Cultivar = paste(unique(metadata$Cultivar_info), collapse = "__"), 
                       Genotype = paste(unique(metadata$Genotype_info), collapse = "__"),
                       n_sig = nrow(final_res_sig), 
                       n_nosig = nrow(final_res_nosig), 
                       n_total = nrow(final_res))
      summary = rbind(summary, sum)
      
      i = i + 1
    }
  }
  
  SUMMARY = rbind(SUMMARY, summary)
    
  rm(list = c("dds", "final_res", "final_res_sig", "final_res_nosig", "metadata", "res_deseq", "resLFC", "summary", "contrast", "i",
              "path_raw_out", "path_sig_out", "path_nosig_out", "results_names", "sum"))
}

write.table(SUMMARY, paste0(path_DEA, "/04-DEA/SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("experiment", "Experiments"))

cat("\n\n\n--------------------------------------------------------------------")

