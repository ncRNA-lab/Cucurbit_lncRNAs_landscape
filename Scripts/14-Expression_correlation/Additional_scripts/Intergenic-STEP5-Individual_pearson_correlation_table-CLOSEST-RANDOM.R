################################################################################
#
# STEP5: CREATE CIS CORRELATION TABLE (CLOSEST-RANDOM)
#
# For each subexperiment, we loaded the samples (obtained with Salmon), filtered 
# the expression using median absolute deviation (Palos et al., 2022 in The Plant 
# cell) and applied variance stabilizing transformation. Then, we generated a 
# Pearson correlation matrix with all PCGs and lincRNAs present in the cis-interactions 
# table. Finally, we extract from this table random pairs and we generate the 
# cis-correlation table for each subexperiment. The table shows lncRNA-PCG and 
# PCG-PCG interactions. In this case, we use the strategy "the closest gene".
#
# @author: pasviber - Pascual Villalba Bermell
# 
################################################################################


######### MODULES

suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(DESeq2))


######### VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 9) {
  stop("At least 9 arguments must be supplied.", call.=FALSE)
} else {
  spel = args[1]
  experiment = args[2]
  WD_corr_S2 = args[3]
  WD_corr_S5 = args[4]
  WD_pred = args[5]
  WD_quant = args[6]
  WD_DEA = args[7]
  n_random_pairs = as.numeric(args[8])
  flag = args[9]
}

# spel = "C. melo"
# experiment = "SRP139690_1_1"
# WD_corr_S2 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/intergenic/nr/STEP2"
# WD_corr_S5 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/14-Expression_correlation/cme/intergenic/nr/STEP5"
# WD_pred = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme"
# WD_quant = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification/cme"
# WD_DEA = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/13-DEA/cme"
# n_random_pairs = 1500
# flag = "nr"


######### PIPELINE

cat(paste0("-Closest...\n"))

## Load TAB_CIS_closest.
TAB_CIS_closest = read.table(paste0(WD_corr_S2, "/TAB_CIS_closest.tsv"), sep = "\t", header = T, quote = "\"")
## Load summary table with all the information about the experiments and contrasts.
metadata = read.table(paste0(WD_DEA, "/03-Metadata_DEA/", experiment, ".tsv"), sep = "\t", header = T, quote = "\"")
## Create subset table.
subset = TAB_CIS_closest
  
## Load counts tables coming from salmon using tximport.
txdb = suppressMessages(makeTxDbFromGFF(paste0(WD_pred, "/STEP-FINAL/Files/ALL/", flag, "/ALL.gtf"), format = "gtf"))
k = keys(txdb, keytype = "GENEID")
df = suppressMessages(AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = columns(txdb)))
tx2gene = df[, c("TXNAME", "GENEID")]
tx2gene = tx2gene[!duplicated(tx2gene),]

files = file.path(paste0(WD_quant, "/ALL/", flag, "/03-Quant/"), metadata$Sample, "quant.sf")
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
vst_filt_1 = vst[rownames(vst) %in% unique(c(subset$ID_transcript.1, subset$ID_transcript.2)),]
random_pairs = data.frame(ID_transcript.1 = row.names(slice_sample(vst_filt_1, n = n_random_pairs, replace = FALSE)),
                          ID_transcript.2 = row.names(slice_sample(vst_filt_1, n = n_random_pairs, replace = FALSE)))
random_pairs_filt = random_pairs[!duplicated(random_pairs),]

rm(list = c("dds", "vst_filt_1", "random_pairs"))

## Correlation method.
vst_filt_2 = vst[rownames(vst) %in% unique(c(random_pairs_filt$ID_transcript.1, random_pairs_filt$ID_transcript.2)),]
tab_cor = suppressWarnings(cor(t(vst_filt_2), method = "pearson"))
tab_cor = melt(tab_cor)
colnames(tab_cor) = c("ID_transcript.1", "ID_transcript.2", "Cor.Deseq.Vst")
tab_cor_random = merge(random_pairs_filt, tab_cor, by = c("ID_transcript.1", "ID_transcript.2"), all = F)

rm(list = c("metadata", "vst", "vst_filt_2", "tab_cor", "random_pairs_filt"))

## Generate the final tables.
TAB_CIS_closest_corr =  tab_cor_random
TAB_CIS_closest_corr$Type = "Random pairs"
TAB_CIS_closest_corr$Experiment = experiment
TAB_CIS_closest_corr$Specie = spel
TAB_CIS_closest_corr =  TAB_CIS_closest_corr[, c("Specie", "ID_transcript.1", "ID_transcript.2", "Type", "Experiment", "Cor.Deseq.Vst")]
rownames(TAB_CIS_closest_corr) = NULL

rm(list = c("tab_cor_random"))

## Save.
write.table(TAB_CIS_closest_corr, paste0(WD_corr_S5, "/TAB_CIS-", experiment, "-PEARSON_closest-random.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = c("TAB_CIS_closest", "subset"))

