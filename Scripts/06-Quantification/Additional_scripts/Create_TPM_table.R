################################################################################
#
# CREATE TPM TABLE (TRANSCRIPT-LEVEL)
#
# Load the quantification files coming from salmon and create a TPM table to
# transcript-level.
#
# @author: pasviber - Pascual Villalba Bermell
#
################################################################################

rm(list = ls())



## 0. INSTALL AND LOAD LIBRARIES

# BiocManager::install("tximport")
# BiocManager::install("GenomicFeatures")
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures))



## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("At least 4 arguments must be supplied.", call.=FALSE)
} else {
  WD = args[1]
  quants = args[2]
  gtf = args[3]
  samples_file = args[4]
}

# WD = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification/cme/ALL/nr"
# quants = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-Quantification/cme/ALL/nr/03-Quant"
# gtf = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/05-LncRNAs_prediction/cme/STEP-FINAL/Files/ALL/nr/ALL.gtf"
# samples_file = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Additional_info/sra-info/accession_list/cme-SRR_Acc_List-Filter_2.txt"

if (!dir.exists(paste0(WD, "/04-Table"))){
  dir.create(paste0(WD, "/04-Table"))
}



## 2. LOAD COUNTS TABLES COMING FROM SALMON USING TXIMPORT.

txdb = makeTxDbFromGFF(gtf, format = "gtf")
k = keys(txdb, keytype = "GENEID")
df = select(txdb, keys = k, keytype = "GENEID", columns = columns(txdb))
tx2gene = df[, c("TXNAME", "GENEID")]
tx2gene = tx2gene[!duplicated(tx2gene),]

samples = read.table(samples_file, header = F)
samples = samples$V1
files = file.path(quants, samples, "quant.sf")
names(files) = paste0(samples)
txi.salmon.no = tximport(files, type = "salmon", tx2gene = tx2gene, txIn = TRUE, txOut = TRUE, countsFromAbundance = "no")



## 3. CREATE THE TPM TABLE.

abund.no = txi.salmon.no$abundance
count.no = txi.salmon.no$counts
length.no = txi.salmon.no$length

abund.no.transID = cbind(data.frame(ID_transcript = rownames(abund.no)), abund.no)
rownames(abund.no.transID) = NULL
write.table(abund.no.transID, paste0(WD, "/04-Table/TPMs.tsv"), col.names = T, row.names = F, sep = "\t", quote = F)



## 4. CREATE THE SUMMARY TPM TABLE.

sum_TPMs = data.frame(
  ID_transcript = rownames(abund.no), 
  TPMs.mean = rowSums(abund.no)/ncol(abund.no), 
  log2.TPMs.1.mean = log((rowSums(abund.no)/ncol(abund.no)) + 1, 2), 
  log2.TPMs.mean = log(rowSums(abund.no)/ncol(abund.no), 2)
  )
#sum_TPMs = sum_TPMs[sum_TPMs$TPMs.mean >= 1,]
write.table(sum_TPMs, paste0(WD, "/04-Table/TPMs_summary.tsv"), col.names = T, row.names = F, sep = "\t", quote = F)

