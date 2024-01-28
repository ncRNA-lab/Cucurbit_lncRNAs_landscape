################################################################################
#
# CREATE THE REDUNDANT LNCRNA DATABASE
#
# Create a lncRNA database using the information obtained in the lncRNA prediction
# pipeline.
#
################################################################################

rm(list = ls())



## 0. INSTALL AND LOAD LIBRARIES.

# remotes::install_github("jokergoo/ComplexHeatmap")
# install.packages("UpSetR")
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("limma"))
suppressMessages(library("Cairo"))
suppressMessages(options(bitmapType='cairo'))
suppressMessages(library("UpSetR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))



## 1. VARIABLES.

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("At least 4 arguments must be supplied.", call.=FALSE)
} else {
  specie = args[1]
  WD = args[2]
  AI = args[3]
  SP = args[4]
}

# specie = "cme"
# WD = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/cme"
# AI = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Additional_info"
# SP = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Scripts/Pascual/05-predict_lncRNAs/Softwares_prediction"



## 2. CREATE THE LNCRNAS DATABASE.

### 2.1 STEP 1: Load initial lncRNA database (genomic location and some molecular features).
cat("LOAD STEP 1 DATA...\n")

DB = read.table(paste0(WD, "/STEP1/Potential_lncRNAs/POTENTIAL_LNCRNAS.tsv"), header = T, sep = "\t", quote = "\"")


### 2.2 STEP 2: Add coding potential information (CPC2, FEELnc and CPAT).
cat("LOAD STEP 2 DATA...\n")

## CPC2.
cat("\t- CPC2...\n")
CPC2 = tryCatch({read.table(paste0(WD, "/STEP2/CPC2/noncoding_ids_CPC2.txt"), quote = "\"")}, error = function(e) {NULL})
if (!is.null(CPC2)) {
  DB$"CPC2" = ifelse(DB$ID_transcript %in% CPC2$V1, "NC", "C")
  rm(list = c("CPC2"))
}else {
  DB$"CPC2" = "C"
}

## FEELnc.
cat("\t- FEELnc...\n")
FEELnc = tryCatch({read.table(paste0(WD, "/STEP2/FEELnc/noncoding_and_no_ORF_ids_FEELnc.txt"), quote = "\"")}, error = function(e) {NULL})
if (!is.null(FEELnc)) {
  DB$"FEELnc" = ifelse(DB$ID_transcript %in% FEELnc$V1, "NC", "C")
  rm(list = c("FEELnc"))
}else {
  DB$"FEELnc" = "C"
}

## CPAT.
cat("\t- CPAT...\n")
CPAT = tryCatch({read.table(paste0(WD, "/STEP2/CPAT/noncoding_and_no_ORF_ids_CPAT.txt"), quote = "\"")}, error = function(e) {NULL})
if (!is.null(CPAT)) {
  DB$"CPAT" = ifelse(DB$ID_transcript %in% CPAT$V1, "NC", "C")
  rm(list = c("CPAT"))
}else {
  DB$"CPAT" = "C"
}


### 2.3 STEP 3: Add coding potential information (SwissProt and Pfam).
cat("LOAD STEP 3 DATA...\n")

## SwissProt.
cat("\t- SwissProt...\n")
SwissProt = tryCatch({read.table(paste0(WD, "/STEP3/SwissProt/diamond_output.tsv"), header = T, sep = "\t", quote = "\"")}, error = function(e) {NULL})
if (!is.null(SwissProt)) {
  DB$"SwissProt" = ifelse(DB$ID_transcript %in% SwissProt$qseqid, "C", "NC")
  write.table(SwissProt, paste0(WD, "/STEP-FINAL/Database/Extra/SwissProt.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  rm(list = c("SwissProt"))
}else {
  DB$"SwissProt" = "NC"
}

## Pfam.
cat("\t- Pfam...\n")
Pfam = tryCatch({read.table(paste0(WD, "/STEP3/Pfam/Hmmer/Results_PFAM_domtblout.tsv"), header = T, sep = "\t", quote = "\"")}, error = function(e) {NULL})
if (!is.null(Pfam)) {
  Pfam_dat = read.table(paste0(SP, "/PFAM/Pfam-A.hmm.dat.tsv"), header = T, sep = "\t", quote = "\"")
  Pfam_add = merge(Pfam, Pfam_dat, by.x = "query_name", by.y = "ID", all.x = T, all.y = F)
  write.table(Pfam_add, paste0(WD, "/STEP-FINAL/Database/Extra/Pfam.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  #
  Pfam_add_mod = Pfam_add
  Pfam_add_mod$target_name = gsub(".p[0-9]+", "", Pfam_add_mod$target_name)
  #
  DB$"Pfam" = ifelse(DB$ID_transcript %in% Pfam_add_mod$target_name, "C", "NC")
  rm(list = c("Pfam_dat", "Pfam_add", "Pfam_add_mod", "Pfam"))
}else {
  DB$"Pfam" = "NC"
}


### 2.4 STEP 4: Add information about potential peptides (length >80aa, >100aa and >120aa).
cat("LOAD STEP 4 DATA...\n")

## Transdecoder (> 80aa).
cat("\t- Transdecoder 80aa...\n")
Transdecoder_80 = tryCatch({read.table(paste0(WD, "/STEP4/ORF/Transdecoder_80/ids_more_than_80_aa.txt"), quote = "\"")}, error = function(e) {NULL})
if (!is.null(Transdecoder_80)) {
  DB$"ORF_80" = ifelse(DB$ID_transcript %in% Transdecoder_80$V1, "YES", "NO")
  rm(list = c("Transdecoder_80"))
}else {
  DB$"ORF_80" = "NO"
}

## Transdecoder (> 100aa).
cat("\t- Transdecoder 100aa...\n")
Transdecoder_100 = tryCatch({read.table(paste0(WD, "/STEP4/ORF/Transdecoder_100/ids_more_than_100_aa.txt"), quote = "\"")}, error = function(e) {NULL})
if (!is.null(Transdecoder_100)) {
  DB$"ORF_100" = ifelse(DB$ID_transcript %in% Transdecoder_100$V1, "YES", "NO")
  rm(list = c("Transdecoder_100"))
}else {
  DB$"ORF_100" = "NO"
}

## Transdecoder (> 120aa).
cat("\t- Transdecoder 120aa...\n")
Transdecoder_120 = tryCatch({read.table(paste0(WD, "/STEP4/ORF/Transdecoder_120/ids_more_than_120_aa.txt"), quote = "\"")}, error = function(e) {NULL})
if (!is.null(Transdecoder_120)) {
  DB$"ORF_120" = ifelse(DB$ID_transcript %in% Transdecoder_120$V1, "YES", "NO")
  rm(list = c("Transdecoder_120"))
}else {
  DB$"ORF_120" = "NO"
}


### 2.5 STEP 5: Add information about housekeeping RNAs and miRNA precursors (Databases: RNAcentral, miRBase, PmiREN).
cat("LOAD STEP 5 DATA...\n")

## RNAcentral.
cat("\t- RNAcentral (rRNA, tRNA, snRNA, snoRNA)...\n")
RNAcentral = tryCatch({read.table(paste0(WD, "/STEP5/RNAcentral/output_blastn_no_bacterial.tsv"), header = T, sep = "\t", quote = "\"")}, error = function(e) {NULL})
if (!is.null(RNAcentral)) {
  RNAcentral_conv = read.table(paste0(AI, "/RNAcentral/ALL/cucurbitaceae_RNAcentral_conversion.tsv"), header = F, sep = "\t", quote = "\"")
  colnames(RNAcentral_conv) = c("Term_long", "Annotation")
  lst = strsplit(RNAcentral_conv$Term_long, "\\s+")
  RNAcentral_conv$"sseqid" = sapply(lst ,`[`, 1)
  RNAcentral = left_join(RNAcentral, RNAcentral_conv[,c("sseqid", "Annotation")], by = "sseqid")
  write.table(RNAcentral, paste0(WD, "/STEP-FINAL/Database/Extra/RNAcentral.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  #
  RNAcentral_red = RNAcentral[,c("qseqid", "Annotation")]
  colnames(RNAcentral_red) = c("ID_transcript", "RNAcentral")
  RNAcentral_red = RNAcentral_red[!duplicated(RNAcentral_red),]
  #
  RNAcentral_red_rRNA = RNAcentral_red[RNAcentral_red$RNAcentral == "rRNA", "ID_transcript"]
  RNAcentral_red_tRNA = RNAcentral_red[RNAcentral_red$RNAcentral == "tRNA", "ID_transcript"]
  RNAcentral_red_snRNA = RNAcentral_red[RNAcentral_red$RNAcentral == "snRNA", "ID_transcript"]
  RNAcentral_red_snoRNA = RNAcentral_red[RNAcentral_red$RNAcentral == "snoRNA", "ID_transcript"]
  #
  DB$"RNAcentral_rRNA" = ifelse(DB$ID_transcript %in% RNAcentral_red_rRNA, "rRNA", NA)
  DB$"RNAcentral_tRNA" = ifelse(DB$ID_transcript %in% RNAcentral_red_tRNA, "tRNA", NA)
  DB$"RNAcentral_snRNA" = ifelse(DB$ID_transcript %in% RNAcentral_red_snRNA, "snRNA", NA)
  DB$"RNAcentral_snoRNA" = ifelse(DB$ID_transcript %in% RNAcentral_red_snoRNA, "snoRNA", NA)
  rm(list = c("RNAcentral_conv", "RNAcentral", "RNAcentral_red", "lst", "RNAcentral_red_rRNA", 
              "RNAcentral_red_tRNA", "RNAcentral_red_snRNA", "RNAcentral_red_snoRNA"))
}else {
  DB$"RNAcentral_rRNA" = NA
  DB$"RNAcentral_tRNA" = NA
  DB$"RNAcentral_snRNA" = NA
  DB$"RNAcentral_snoRNA" = NA
}

## miRBase.
cat("\t- miRBase (precursor-miRNAs)...\n")
miRBase_prec = tryCatch({read.table(paste0(WD, "/STEP5/miRBase/prec_ID.txt"), quote = "\"")}, error = function(e) {NULL})
if (!is.null(miRBase_prec)) {
  DB$"miRBase" = ifelse(DB$ID_transcript %in% miRBase_prec$V1, "precursor-miRNA", NA)
  rm(list = c("miRBase_prec"))
}else {
  DB$"miRBase" = NA
}

## PmiREN.
cat("\t- PmiREN (precursor-miRNAs)...\n")
PmiREN_prec = tryCatch({read.table(paste0(WD, "/STEP5/PmiREN/prec_ID.txt"), quote = "\"")}, error = function(e) {NULL})
if (!is.null(PmiREN_prec)) {
  DB$"PmiREN" = ifelse(DB$ID_transcript %in% PmiREN_prec$V1, "precursor-miRNA", NA)
  rm(list = c("PmiREN_prec"))
}else {
  DB$"PmiREN" = NA
}


### 2.6 STEP 6: Add information about known potential lncRNAs (Databases: Cantatadb, Plncdb and Greenc).
cat("LOAD STEP 6 DATA...\n")

## Cantatadb.
cat("\t- Cantatadb...\n")
Cantatadb = tryCatch({read.table(paste0(WD, "/STEP6/CANTATAdb/output_blastn.tsv"), header = T, sep = "\t", quote = "\"")}, error = function(e) {NULL})
if (!is.null(Cantatadb)) {
  DB$"CANTATAdb" = ifelse(DB$ID_transcript %in% Cantatadb$qseqid, "YES", "NO")
  rm(list = c("Cantatadb"))
}else {
  DB$"CANTATAdb" = "NO"
}

## Plncdb.
cat("\t- Plncdb...\n")
Plncdb = tryCatch({read.table(paste0(WD, "/STEP6/PLncDB/output_blastn.tsv"), header = T, sep = "\t", quote = "\"")}, error = function(e) {NULL})
if (!is.null(Plncdb)) {
  DB$"PLncDB" = ifelse(DB$ID_transcript %in% Plncdb$qseqid, "YES", "NO")
  rm(list = c("Plncdb"))
}else {
  DB$"PLncDB" = "NO"
}

## Greenc.
cat("\t- Greenc...\n")
Greenc = tryCatch({read.table(paste0(WD, "/STEP6/GreeNC/output_blastn.tsv"), header = T, sep = "\t", quote = "\"")}, error = function(e) {NULL})
if (!is.null(Greenc)) {
  DB$"GreeNC" = ifelse(DB$ID_transcript %in% Greenc$qseqid, "YES", "NO")
  rm(list = c("Greenc"))
}else {
  DB$"GreeNC" = "NO"
}


### 2.7 N's.
cat("CALCULATE N'S PERCENTAGE...\n")
seq = read.table(paste0(WD, "/STEP1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta"), sep = "\t", header = F, quote = "\"")
SEQ = data.frame(ID = gsub(">", "", seq$V1[grepl(">", seq$V1)]), seq = seq$V1[!grepl(">", seq$V1)])
SEQ$"Ns" = (str_count(toupper(SEQ$seq), "N")*100)/nchar(toupper(SEQ$seq))
DB = merge(DB, SEQ[, c("ID", "Ns")], by.x = "ID_transcript", by.y = "ID")

rm(list = c("seq", "SEQ"))


### 2.8 ASSIGN CONFIDENCE LEVEL.
cat("ASSIGN CONFIDENCE LEVELS TO THE POTENTIAL LNCRNAS...\n")
P = rowSums(DB[,c("CPC2", "FEELnc", "CPAT")] == "NC")
B = rowSums(DB[,c("SwissProt", "Pfam")] == "NC")

DB$"Confidence" = ifelse(P == 3 & B == 2,
                                 "High",
                                 ifelse((P == 3 & B == 1) | (P == 2 & B == 2),
                                        "Medium",
                                        ifelse((P == 1 & B == 1) | (P == 2 & B == 0) | (P == 2 & B == 1),
                                               "Low",
                                               "Unclassified")))

rm(list = c("P", "B"))


### 2.9 CREATE THE DIFFERENT FINAL DATABASES.
cat("CREATE THE DIFFERENT FINAL DATABASES...\n")

# All transcripts predicted as non-coding by some of the prediction software.
write.table(DB, paste0(WD, "/STEP-FINAL/Database/Database_ALL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

# All transcripts with assigned confidence level (Unclassified transcripts are removed).
DB_filt = DB[DB$Confidence != "Unclassified",]
rownames(DB_filt) = NULL

write.table(DB_filt, paste0(WD, "/STEP-FINAL/Database/Database_Classified.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

# All transcripts with assigned confidence level and not classified as rRNA, tRNA, 
# snRNA, snoRNA or miRNA precursor.
DB_filt_filt = DB_filt
DB_filt_filt = DB_filt_filt[is.na(DB_filt_filt$RNAcentral_rRNA),]
DB_filt_filt = DB_filt_filt[is.na(DB_filt_filt$RNAcentral_tRNA),]
DB_filt_filt = DB_filt_filt[is.na(DB_filt_filt$RNAcentral_snRNA),]
DB_filt_filt = DB_filt_filt[is.na(DB_filt_filt$RNAcentral_snoRNA),]
DB_filt_filt = DB_filt_filt[is.na(DB_filt_filt$PmiREN),]
DB_filt_filt = DB_filt_filt[is.na(DB_filt_filt$miRBase),]

write.table(DB_filt_filt, paste0(WD, "/STEP-FINAL/Database/Database_LncRNAs_R.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)



## 3. FIGURES.

### 3.1 VennDiagram: Non-coding transcripts according to each prediction software.
cat("VENNDIAGRAM: NON-CODING TRANSCRIPTS ACCORDING TO EACH PREDICTION SOFTWARE...\n")

x = list(
  CPC2 = unique(DB[DB$CPC2 == "NC", "ID_transcript"]),
  FEELnc = unique(DB[DB$FEELnc == "NC", "ID_transcript"]),
  CPAT = unique(DB[DB$CPAT == "NC", "ID_transcript"]),
  PFAM = unique(DB[DB$Pfam == "NC", "ID_transcript"]),
  SwissProt = unique(DB[DB$SwissProt == "NC", "ID_transcript"])
)

A = list_to_matrix(x)
Z = vennCounts(A, include = "both")

png(filename=paste0(WD, '/STEP-FINAL/Figures/Venn_Diagram_NC_Prediction.png'), height = 5000, width = 5000, res=600)
vennDiagram(Z, 
            include="both", 
            names=NULL, 
            mar=rep(1,4), 
            cex=c(1.5,1,0.7), 
            lwd=1,
            circle.col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), 
            counts.col=NULL, 
            show.include=NULL
            )
invisible(dev.off())

rm(list = c("A", "Z", "x"))


### 3.2 Upset plot: Classified non-coding transcripts (assigned confidence level)
### by type of non-coding sequence (rRNA, tRNA, snRNA, snoRNA, miRNA precursor,
### knonwn lncRNA or lncRNA).
cat("UPSET PLOT: CLASSIFIED NON-CODING TRANSCRIPTS BY TYPE OF NON-CODING SEQUENCE...\n")

x = data.frame(
  matrix(vector(), 
         nrow = dim(DB_filt)[1], 
         ncol = 8, 
         dimnames=list(c(), 
                       c("ID_transcript", "rRNA", "tRNA", "snRNA", "snoRNA", "precursor.miRNA", "known.lncRNA", "lncRNA"))),
  stringsAsFactors=F)

x$ID_transcript = DB_filt$ID_transcript
x$rRNA = ifelse(!is.na(DB_filt$RNAcentral_rRNA), ifelse(DB_filt$RNAcentral_rRNA == "rRNA", 1, 0), 0)
x$tRNA = ifelse(!is.na(DB_filt$RNAcentral_tRNA), ifelse(DB_filt$RNAcentral_tRNA == "tRNA", 1, 0), 0)
x$snRNA = ifelse(!is.na(DB_filt$RNAcentral_snRNA), ifelse(DB_filt$RNAcentral_snRNA == "snRNA", 1, 0), 0)
x$snoRNA = ifelse(!is.na(DB_filt$RNAcentral_snoRNA), ifelse(DB_filt$RNAcentral_snoRNA == "snoRNA", 1, 0), 0)
x$precursor.miRNA = ifelse(!is.na(DB_filt$PmiREN), ifelse(DB_filt$PmiREN == "precursor-miRNA" , 1, 0), 0)
x$precursor.miRNA = ifelse(!is.na(DB_filt$miRBase), ifelse(DB_filt$miRBase == "precursor-miRNA", 1, x$precursor.miRNA), x$precursor.miRNA)
x$known.lncRNA = ifelse(DB_filt$CANTATAdb == "YES" | DB_filt$PLncDB == "YES" | DB_filt$GreeNC == "YES", 1, 0)
x$lncRNA = ifelse(rowSums(x[,2:7]) == 0, 1, 0)

names(x)[names(x) == "precursor.miRNA"] = "precursor-miRNA"
names(x)[names(x) == "known.lncRNA"] = "Homologous to \nknown lncRNAs"
names(x)[names(x) == "lncRNA"] = "Novel potential \nlncRNAs"

png(filename=paste0(WD, '/STEP-FINAL/Figures/Upset_Plot_NC_Classification.png'), height = 5500, width = 7000, res=600)
upset(x, 
      sets = c("rRNA", "tRNA", "snRNA", "snoRNA", "precursor-miRNA", "Homologous to \nknown lncRNAs", "Novel potential \nlncRNAs"), 
      sets.bar.color = c(rep("#d16262", 2), rep("#86c3e7", 5)),
      order.by = "freq", 
      #empty.intersections = "on",
      text.scale = 1.5,
      mainbar.y.label = "Number of transcripts (Intersection Size)",
      sets.x.label = "Number of transcripts (Set Size)"
      )
invisible(dev.off())

rm(list = c("x"))


### 3.3 Bar plot: LncRNA classification by class code.
cat("BAR PLOT: LNCRNA CLASSIFICATION BY CLASS CODE...\n")

Tab = DB_filt_filt[,c("ID_transcript", "Class_code")]
Tab$"Counts" = 1
Tab[Tab == "u"] = "intergenic (u)"
Tab[Tab == "x"] = "antisense (x)"
Tab[Tab == "i"] = "intronic (i)"
Tab[Tab == "o"] = "sense (o/e)"
Tab[Tab == "e"] = "sense (o/e)"
Tab$Class_code = factor(Tab$Class_code, levels = c("intergenic (u)", "antisense (x)", "intronic (i)", "sense (o/e)"))
Tab_FINAL = Tab %>% 
  group_by(Class_code, .drop=FALSE) %>%
  summarise(
    sum.Counts = sum(Counts))
Tab_FINAL = as.data.frame(Tab_FINAL)

gg = ggplot(Tab_FINAL, aes(x = Class_code, y = sum.Counts, fill = Class_code)) + 
  geom_bar(position=position_dodge(), aes(y=sum.Counts), stat="identity") +
  scale_fill_manual(values = c("#80c0e5", "#e2e89d", "#baf4b4", "#eeb8a1")) +
  xlab("") +
  ylab("Number of LncRNAs") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5)) + 
  theme(legend.position = "none") + 
  geom_text(aes(label = sum.Counts), nudge_y = 100)

ggsave(paste0(WD, "/STEP-FINAL/Figures/Bar_Plot_LncRNA_Classification.png"), height = 7, width = 8, dpi = 600)

rm(list = c("Tab", "Tab_FINAL", "gg"))


