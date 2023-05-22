################################################################################
# CHOOSE THE BEST X HITS FROM THE BLASTN TABLES
################################################################################


## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(dplyr))

## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments must be supplied.", call.=FALSE)
} else {
  spe1 = args[1]
  spe2 = args[2]
  tab_in = args[3]
  tab_out = args[4]
  n = as.integer(args[5])
}

# spe1 = "0"
# spe2 = "1"
# tab_in="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sequence_level/OrthoFinder/nr/Prueba_max_5//03-Blastn/High/ALL/Blast0_1_temp1.txt"
# tab_out="/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sequence_level/OrthoFinder/nr/Prueba_max_5//03-Blastn/High/ALL/Blast0_1_temp2.txt"
# n = 5

## 2. PIPELINE

tabin = tryCatch({read.table(tab_in, header = F, sep = "\t", quote = "\"")}, error = function(e) {NULL})

if (!is.null(tabin)) {
  colnames(tabin) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                      "qend", "sstart", "send", "evalue", "bitscore")
  
  tabin_filt_1 = tabin %>% 
    group_by(qseqid, sseqid) %>% 
    arrange(qseqid, sseqid, evalue, desc(pident), desc(length)) %>% 
    slice(1)
  
  if (spe1 == spe2) {
    # In this study, we only want to study conservation between species and not within species. 
    # Therefore, the logical thing to do is to run Blast with all possible paired combinations 
    # without necessarily combining a species against itself. However, programs such as Orthofinder 
    # and OrthoMCL require a blast of each species against itself. To avoid the introduction of 
    # paralogs, we will limit the acceptance of only hits that are of a transcript against itself.
    tabin_filt_2 = tabin_filt_1[tabin_filt_1$qseqid == tabin_filt_1$sseqid,]
  }
  
  if (spe1 != spe2) {
    tabin_filt_2 = tabin_filt_1 %>% 
      group_by(qseqid) %>% 
      arrange(qseqid, evalue, desc(pident), desc(length)) %>%
      slice_head(n = n)
  }
  
  write.table(tabin_filt_2, tab_out, col.names = F, sep = "\t", row.names = F, quote = F)
} 

if (is.null(tabin)) {
  log = file.copy(from = tab_in, to = tab_out)
}

