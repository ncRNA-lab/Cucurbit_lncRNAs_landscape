################################################################################
#
# ALL: TISSUE SPECIFICITY STUDY: APPROACH 1 - STEP 1
#
# We're going to create different tables with this R script:
# - TAB (TAB_TPMs.tsv): TPMs by ID_transcript, Sample and Specie.
# - TAB_red (TAB_TPMs_summary.tsv): Mean, median, max and min TPMs by ID_transcript, 
#   SRA.Study, Tissue and Specie.
# - TAB_red_filt (TAB_TPMs_summary_filtered.tsv): Mean, median, max and min TPMs 
#   by ID_transcript, SRA.Study, Tissue and Specie from SRA.studies which have 3 
#   or more than 3 tissues.
# - TAB_experiments (TAB_experiments.tsv): SRA.Studies classified as YES (pass the 
#   tissues filter) or NO (doesn't pass the tissues filter).
# - TAB_experiments_summary (TAB_experiments_summary.tsv): SRA.Studies classified 
#   as YES (pass the tissues filter) or NO (doesn't pass the tissues filter) by 
#   specie.
# - Create an expression table (ID_transcript + Confidence + Class_code ~ Tissue)
#   by SRA.study filtered (pass the tissues filter).
# 
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

## 1. PATHS

# Own computer.
path_tissue_specificity = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
path_pred = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs"
path_quant = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/06-quantification"
path_AI = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Additional_info"
flag = "nr"
species = c("car", "cla", "cma", "cme", "cmo", "cpe", "csa", "lsi", "mch")

# # Garnatxa
# path_tissue_specificity = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/09-Tissue-specificity"
# path_pred = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs"
# path_quant = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-quantification"
# path_AI = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info"
# flag = "nr"
# species = c("car", "cla", "cma", "cme", "cmo", "cpe", "csa", "lsi", "mch")

## 2. CREATE THE DIRECTORIES

if (!dir.exists(path_tissue_specificity)){
    dir.create(path_tissue_specificity)
}
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1"))
}
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/ALL"))
}
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag))
}
if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1"))){
  dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1"))
}

## 3. CREATE TABLE

# PIPELINE

TAB = data.frame()
TAB_red = data.frame()
TAB_red_filt = data.frame()
TAB_experiments = data.frame()

for (spe in species) {
  cat(paste0("\n\n\n\n\n----- Specie: ", spe, "\n\n"))
  
  # Transcripts info.
  gtf_file = paste0(path_pred, "/", spe, "/STEP-FINAL/Files/Joined/ALL/", flag, "/ALL.gtf")
  txdb = makeTxDbFromGFF(gtf_file, format = "gtf")
  k = keys(txdb, keytype = "GENEID")
  df = AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = columns(txdb))
  tx2gene = df[, c("TXNAME", "GENEID")]
  tx2gene = tx2gene[!duplicated(tx2gene),]
  
  # Metadata info.
  metadata_file = paste0(path_AI, "/Metadata_modified/", spe, "-SraRunTablestandard.csv")
  metadata = read.csv(metadata_file, header = T)
  metadata = metadata[, c("Run", "Tissue.Standard", "SRA.Study")]
  colnames(metadata) = c("Sample", "Tissue", "SRA.Study")
  
  # Create txi.salmon object.
  quant_dir = paste0(path_quant, "/", spe, "/Salmon/ALL/", flag, "/03-Quant/")
  files = file.path(quant_dir, metadata$Sample, "quant.sf")
  names(files) = paste0(metadata$Sample)
  txi.salmon = tximport(files, type = "salmon", tx2gene = tx2gene, txIn = TRUE, txOut = TRUE, countsFromAbundance = "no")
  
  # Extract TPMs table.
  TPMs = as.data.frame(txi.salmon$abundance)
  TPMs$"ID_transcript" = rownames(TPMs)
  TPMs = TPMs[, c("ID_transcript", metadata$Sample)]
  rownames(TPMs) = NULL
  
  # Load the LncRNA database and gene info.
  if (flag == "nr") {
    LncRNAs_info = read.table(paste0(path_pred, "/", spe, "/STEP-FINAL/Database/Database_LncRNAs_NR.tsv"), sep = "\t", header = T, quote = "\"")
    LncRNAs_info = LncRNAs_info[, c("ID_transcript", "Class_code", "Significance_level")]
    colnames(LncRNAs_info) = c("ID_transcript", "Class_code", "Confidence")
  } 
  if (flag == "r") {
    LncRNAs_info = read.table(paste0(path_pred, "/", spe, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), sep = "\t", header = T, quote = "\"")
    LncRNAs_info = LncRNAs_info[, c("ID_transcript", "Class_code", "Significance_level")]
    colnames(LncRNAs_info) = c("ID_transcript", "Class_code", "Confidence")
  }
  LncRNAs_info[LncRNAs_info == 'o'] = 'o/e'
  LncRNAs_info[LncRNAs_info == 'e'] = 'o/e'
  LncRNAs_info[LncRNAs_info == 'Low'] = 'Low-confidence lncRNA'
  LncRNAs_info[LncRNAs_info == 'Medium'] = 'Medium-confidence lncRNA'
  LncRNAs_info[LncRNAs_info == 'High'] = 'High-confidence lncRNA'
  Genes_info = read.table(paste0(path_pred, "/", spe, "/STEP-FINAL/Files/Genes/ORIGINAL_GENES.tsv"), sep = "\t", header = T, quote = "\"")
  Genes_info$"Class_code" = "pc"
  Genes_info$"Confidence" = "Protein-coding"
  Genes_info = Genes_info[, c("ID_transcript", "Class_code", "Confidence")]
  info = rbind(LncRNAs_info, Genes_info)
  
  # Create tab table.
  cat(paste0("\n\nCreating tab table\n"))
  TPMs = merge(TPMs, info, by = "ID_transcript", all = T)
  TPMs = TPMs[, c("ID_transcript", "Confidence", "Class_code", metadata$Sample)]
  
  tab = melt(setDT(TPMs), id.vars = c("ID_transcript", "Confidence", "Class_code"), variable.name = "Sample", value.name = "TPMs")
  tab$Sample = as.character(tab$Sample)
  
  tab = merge(tab, metadata, by = "Sample", all.x = T)
  tab = tab[order(tab$ID_transcript), c("ID_transcript", "Confidence", "Class_code", "Sample", "SRA.Study", "Tissue", "TPMs")]
  
  # Create tab_red table.
  cat(paste0("\n\nCreating tab_red table\n"))
  tab_red = tab %>% 
    group_by(ID_transcript, Confidence, Class_code, SRA.Study, Tissue) %>%
    summarise(
      Number_samples_by_tissue = n_distinct(Sample),
      Mean.TPMs = mean(TPMs),
      Median.TPMs = median(TPMs),
      Max.TPMs = max(TPMs),
      Min.TPMs = min(TPMs)) %>%
    mutate(Number_tissues_by_SRA.Study = n_distinct(Tissue))
  tab_red = as.data.frame(tab_red)
  tab_red = tab_red[, c("ID_transcript", "Confidence", "Class_code", "SRA.Study", "Number_samples_by_tissue", "Number_tissues_by_SRA.Study", "Tissue", "Mean.TPMs", "Median.TPMs", "Max.TPMs", "Min.TPMs")]
  
  # Create tab_red_filt table.
  cat(paste0("\n\nCreating tab_red_filt table\n"))
  tab_red_filt = tab_red[tab_red$Number_tissues_by_SRA.Study >= 3,]
  
  # Create tab_experiments table.
  cat(paste0("\n\nCreating tab_experiments table\n"))
  tab_experiments_temp = tab_red[, c("SRA.Study", "Tissue", "Number_tissues_by_SRA.Study", "Number_samples_by_tissue")]
  tab_experiments_temp = tab_experiments_temp[!duplicated(tab_experiments_temp), c("SRA.Study", "Number_tissues_by_SRA.Study", "Number_samples_by_tissue")]
  tab_experiments = tab_experiments_temp %>% group_by(SRA.Study, Number_tissues_by_SRA.Study) %>% summarise(Number_samples = sum(Number_samples_by_tissue))
  tab_experiments = as.data.frame(tab_experiments)
  tab_experiments$"KEPT" = ifelse(tab_experiments$Number_tissues_by_SRA.Study >= 3, "YES", "NO")
  
  # Add specie column to each table.
  tab$"Specie" = spe
  tab_red$"Specie" = spe
  if (dim(tab_red_filt)[1] > 0) {
    tab_red_filt$"Specie" = spe
  }
  tab_experiments$"Specie" = spe
  
  # Sort the columns.
  tab = tab[,c("Specie", "ID_transcript", "Confidence", "Class_code", "Sample", "SRA.Study", "Tissue", "TPMs")]
  tab_red = tab_red[,c("Specie", "ID_transcript", "Confidence", "Class_code", "SRA.Study", "Number_samples_by_tissue", "Number_tissues_by_SRA.Study", "Tissue", "Mean.TPMs", "Median.TPMs", "Max.TPMs", "Min.TPMs")]
  if (dim(tab_red_filt)[1] > 0) {
    tab_red_filt = tab_red_filt[,c("Specie", "ID_transcript", "Confidence", "Class_code", "SRA.Study", "Number_samples_by_tissue", "Number_tissues_by_SRA.Study", "Tissue", "Mean.TPMs", "Median.TPMs", "Max.TPMs", "Min.TPMs")]
  }
  tab_experiments = tab_experiments[,c("Specie", "SRA.Study", "Number_tissues_by_SRA.Study", "Number_samples", "KEPT")]
  
  # Join the tabs to the TABS.
  TAB = rbind(TAB, tab)
  TAB_red = rbind(TAB_red, tab_red)
  if (dim(tab_red_filt)[1] > 0) {
    TAB_red_filt = rbind(TAB_red_filt, tab_red_filt)
  }
  TAB_experiments = rbind(TAB_experiments, tab_experiments)
  
  # Create a new folder called as specie's name.
  if (!dir.exists(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/", spe))){
    dir.create(paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/", spe))
  }
  
  # Create a expression table (transcripts x tissues) by SRA Study filtered (>= 3 Tissues)
  cat(paste0("\n\nCreating table by SRA.study filtered\n"))
  SRA.studies = unique(tab_red_filt$SRA.Study)
  
  if (length(SRA.studies) > 0) {
    for (SRA.study in SRA.studies) {
      cat(paste0("\n", SRA.study))
      tab_SRA_mean_long = tab_red_filt[tab_red_filt$SRA.Study == SRA.study, c("ID_transcript", "Confidence", "Class_code", "Tissue", "Mean.TPMs")]
      tab_SRA_median_long = tab_red_filt[tab_red_filt$SRA.Study == SRA.study, c("ID_transcript", "Confidence", "Class_code", "Tissue", "Median.TPMs")]
      
      tab_SRA_mean_wide = dcast(setDT(tab_SRA_mean_long), ID_transcript + Confidence + Class_code ~ Tissue, value.var = "Mean.TPMs")
      tab_SRA_median_wide = dcast(setDT(tab_SRA_median_long), ID_transcript + Confidence + Class_code ~ Tissue, value.var = "Median.TPMs")
      
      tab_SRA_mean_wide = as.data.frame(tab_SRA_mean_wide)
      tab_SRA_median_wide = as.data.frame(tab_SRA_median_wide)
      
      write.table(tab_SRA_mean_wide, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/", spe, "/", SRA.study, "_mean.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
      write.table(tab_SRA_median_wide, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/", spe, "/", SRA.study, "_median.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    }
    
    rm(list = c("tab_SRA_mean_long", "tab_SRA_median_long", "tab_SRA_mean_wide", "tab_SRA_median_wide"))
  }
  
  rm(list = c("gtf_file", "txdb", "k", "df", "tx2gene", "metadata_file", "files", "txi.salmon", "LncRNAs_info", 
              "metadata", "TPMs", "tab_experiments_temp", "tab", "tab_red", "tab_red_filt", "tab_experiments", 
              "Genes_info", "info"))
}

# Convert to factors.
TAB_experiments$Specie = factor(TAB_experiments$Specie, levels = species)
TAB_experiments$KEPT = factor(TAB_experiments$KEPT, levels = c("NO", "YES"))

# Create TAB_experiments_summary table.
cat(paste0("\n\nCreating TAB_experiments_summary table\n"))
TAB_experiments_summary = TAB_experiments %>% group_by(Specie, KEPT, .drop = F) %>% summarise(Number_experiments = n_distinct(SRA.Study))
TAB_experiments_summary = as.data.frame(TAB_experiments_summary)

# Save tables.
write.table(TAB, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/TAB_TPMs.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_red, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/TAB_TPMs_summary.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_red_filt, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/TAB_TPMs_summary_filtered.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_experiments, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/TAB_experiments.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_experiments_summary, paste0(path_tissue_specificity, "/approach_1/ALL/", flag, "/STEP1/TAB_experiemnts_summary.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

