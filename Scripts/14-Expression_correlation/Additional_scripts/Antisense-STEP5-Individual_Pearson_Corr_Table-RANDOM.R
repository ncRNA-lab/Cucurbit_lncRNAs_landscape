
######### MODULES

suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(DESeq2))


######### VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
  stop("At least 8 arguments must be supplied.", call.=FALSE)
} else {
  spel = args[1]
  experiment = args[2]
  WD_corr_S5 = args[3]
  WD_pred = args[4]
  WD_quant = args[5]
  WD_DEA = args[6]
  n_random_pairs = as.numeric(args[7])
  flag = args[8]
}

# spel = "C. melo"
# experiment = "SRP198631_1_1"
# WD_corr_S5 = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/17-Correlation/cme/intergenic/STEP5"
# WD_pred = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme"
# WD_quant = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/06-quantification/cme"
# WD_DEA = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/16-DEA/cme"
# n_random_pairs = 1000
# flag = "nr"

######### PIPELINE

## STEP 5: Create table with pearson correlation.

## Load summary table with all the information about the experiments and contrasts.
metadata = read.table(paste0(WD_DEA, "/03-Metadata_DEA/", experiment, ".tsv"), sep = "\t", header = T, quote = "\"")
  
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

## Select random pairs.
set.seed(42)
vst = as.data.frame(assay(dds, 'counts.norm.VST'))
random_pairs = data.frame(ID_transcript.1 = row.names(slice_sample(vst, n = n_random_pairs, replace = FALSE)),
                          ID_transcript.2 = row.names(slice_sample(vst, n = n_random_pairs, replace = FALSE)))
random_pairs_filt = random_pairs[!duplicated(random_pairs),]

rm(list = c("dds", "random_pairs"))

## Correlation method.
vst_filt = vst[rownames(vst) %in% unique(c(random_pairs_filt$ID_transcript.1, random_pairs_filt$ID_transcript.2)),]
tab_cor = suppressWarnings(cor(t(vst_filt), method = "pearson"))
tab_cor = melt(tab_cor)
colnames(tab_cor) = c("ID_transcript.1", "ID_transcript.2", "Cor.Deseq.Vst")
tab_cor_random = merge(random_pairs_filt, tab_cor, by = c("ID_transcript.1", "ID_transcript.2"), all.x = F, all.y = F)

rm(list = c("metadata", "vst", "vst_filt", "tab_cor", "random_pairs_filt"))

## Generate the final tables.
TAB_CIS_closest_corr =  tab_cor_random
TAB_CIS_closest_corr$Type = "Random pairs"
TAB_CIS_closest_corr$Experiment = experiment
TAB_CIS_closest_corr$Specie = spel
TAB_CIS_closest_corr =  TAB_CIS_closest_corr[, c("Specie", "ID_transcript.1", "ID_transcript.2", "Type", "Experiment", "Cor.Deseq.Vst")]
rownames(TAB_CIS_closest_corr) = NULL

rm(list = c("tab_cor_random"))

## Save.
write.table(TAB_CIS_corr_1, paste0(WD_corr_S5, "/TAB_CIS-", experiment, "-PEARSON-random.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

