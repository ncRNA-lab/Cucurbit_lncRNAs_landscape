################################################################################
#
# CREATE THE LNCRNA DATABASE
#
# Create a lncRNA database from the information obtained in the lncRNA prediction
# pipeline.
#
################################################################################


## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Create_customed_LncRNA_db.R $specie $WD2/04-Predict_lncRNAs "Batch-1" $AI $SP 

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

# if (!require(remotes)) {
#   install.packages("remotes")
# }
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


## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments must be supplied.", call.=FALSE)
} else {
  specie = args[1]
  WD = args[2]
  Batch = args[3]
  AI = args[4]
  SP = args[5]
}

# specie = "cme"
# WD = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/13-saturation_plots/cme/04-Predict_lncRNAs"
# Batch = "Batch-1"
# AI = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Additional_info"
# SP = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Scripts/Pascual/05-predict_lncRNAs/Softwares_prediction"



## 2. CREATE THE LNCRNAS DATABASE

### 2.1 STEP 1
cat("LOAD STEP 1 DATA...\n")
LncRNAs = read.table(paste0(WD, "/STEP1/", Batch, "/Potential_lncRNAs/POTENTIAL_LNCRNAS.tsv"), header = T, sep = "\t", quote = "\"")
DB = LncRNAs


### 2.2 STEP 2
cat("LOAD STEP 2 DATA...\n")
CPC2 = read.table(paste0(WD, "/STEP2/", Batch, "/CPC2/noncoding_ids_CPC2.txt"), quote = "\"")
CPC2 = CPC2$V1
CPAT = read.table(paste0(WD, "/STEP2/", Batch, "/CPAT/noncoding_and_no_ORF_ids_CPAT.txt"), quote = "\"")
CPAT = CPAT$V1
FEELnc = read.table(paste0(WD, "/STEP2/", Batch, "/FEELnc/noncoding_and_no_ORF_ids_FEELnc.txt"), quote = "\"")
FEELnc = FEELnc$V1

DB$"CPC2" = ifelse(DB$ID_transcript %in% CPC2, "NC", "C")
DB$"CPAT" = ifelse(DB$ID_transcript %in% CPAT, "NC", "C")
DB$"FEELnc" = ifelse(DB$ID_transcript %in% FEELnc, "NC", "C")

rm(list = c("CPC2", "CPAT", "FEELnc"))


### 2.3 STEP 3
cat("LOAD STEP 3 DATA...\n")
SwissProt = read.table(paste0(WD, "/STEP3/", Batch, "/SwissProt/diamond_output.tsv"), header = T, sep = "\t", quote = "\"")
Pfam = read.table(paste0(WD, "/STEP3/", Batch, "/Pfam/Hmmer/Results_PFAM_domtblout.tsv"), header = T, sep = "\t", quote = "\"")
Pfam_dat = read.table(paste0(SP, "/PFAM/Pfam-A.hmm.dat.tsv"), header = T, sep = "\t", quote = "\"")

Pfam_add = merge(Pfam, Pfam_dat, by.x = "query_name", by.y = "ID", all.x = T, all.y = F)
Pfam_add_mod = Pfam_add
Pfam_add_mod$target_name = gsub(".p[0-9]+", "", Pfam_add_mod$target_name)

DB$"SwissProt" = ifelse(DB$ID_transcript %in% SwissProt$qseqid, "C", "NC")
DB$"Pfam" = ifelse(DB$ID_transcript %in% Pfam_add_mod$target_name, "C", "NC")

write.table(SwissProt, paste0(WD, "/STEP-FINAL/", Batch, "/Database/Extra/SwissProt.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(Pfam_add, paste0(WD, "/STEP-FINAL/", Batch, "/Database/Extra/Pfam.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("Pfam_dat", "SwissProt", "Pfam_add", "Pfam_add_mod", "Pfam"))


### 2.4 STEP 4
cat("LOAD STEP 4 DATA...\n")
Transdecoder_80 = read.table(paste0(WD, "/STEP4/", Batch, "/ORF/Transdecoder_80/ids_more_than_80_aa.txt"), quote = "\"")
Transdecoder_80 = Transdecoder_80$V1
Transdecoder_100 = read.table(paste0(WD, "/STEP4/", Batch, "/ORF/Transdecoder_100/ids_more_than_100_aa.txt"), quote = "\"")
Transdecoder_100 = Transdecoder_100$V1
Transdecoder_120 = read.table(paste0(WD, "/STEP4/", Batch, "/ORF/Transdecoder_120/ids_more_than_120_aa.txt"), quote = "\"")
Transdecoder_120 = Transdecoder_120$V1

DB$"ORF>80" = ifelse(DB$ID_transcript %in% Transdecoder_80, "YES", "NO")
DB$"ORF>100" = ifelse(DB$ID_transcript %in% Transdecoder_100, "YES", "NO")
DB$"ORF>120" = ifelse(DB$ID_transcript %in% Transdecoder_120, "YES", "NO")

rm(list = c("Transdecoder_80", "Transdecoder_100", "Transdecoder_120"))


### 2.5 STEP 5
cat("LOAD STEP 5 DATA...\n")
# RNAcentral
cat("\t- RNAcentral (rRNA, tRNA, snRNA, snoRNA)...\n")
RNAcentral_conv = read.table(paste0(AI, "/RNAcentral/ALL/cucurbitaceae_RNAcentral_conversion.tsv"), header = F, sep = "\t", quote = "\"")
colnames(RNAcentral_conv) = c("Term_long", "Annotation")
lst = strsplit(RNAcentral_conv$Term_long, "\\s+")
RNAcentral_conv$"sseqid" = sapply(lst ,`[`, 1)

RNAcentral = read.table(paste0(WD, "/STEP5/", Batch, "/RNAcentral/output_blastn_no_bacterial.tsv"), header = T, sep = "\t", quote = "\"")
#RNAcentral = RNAcentral[RNAcentral$evalue < 1e-10,]

RNAcentral = left_join(RNAcentral, RNAcentral_conv[,c("sseqid", "Annotation")], by = "sseqid")
RNAcentral_red = RNAcentral[,c("qseqid", "Annotation")]
colnames(RNAcentral_red) = c("ID_transcript", "RNAcentral")
RNAcentral_red = RNAcentral_red[!duplicated(RNAcentral_red),]
RNAcentral_red_rRNA = RNAcentral_red[RNAcentral_red$RNAcentral == "rRNA", "ID_transcript"]
RNAcentral_red_tRNA = RNAcentral_red[RNAcentral_red$RNAcentral == "tRNA", "ID_transcript"]
RNAcentral_red_snRNA = RNAcentral_red[RNAcentral_red$RNAcentral == "snRNA", "ID_transcript"]
RNAcentral_red_snoRNA = RNAcentral_red[RNAcentral_red$RNAcentral == "snoRNA", "ID_transcript"]

write.table(RNAcentral, paste0(WD, "/STEP-FINAL/", Batch, "/Database/Extra/RNAcentral.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

# Add RNAcentral info to the database
DB$"RNAcentral_rRNA" = ifelse(DB$ID_transcript %in% RNAcentral_red_rRNA, "rRNA", NA)
DB$"RNAcentral_tRNA" = ifelse(DB$ID_transcript %in% RNAcentral_red_tRNA, "tRNA", NA)
DB$"RNAcentral_snRNA" = ifelse(DB$ID_transcript %in% RNAcentral_red_snRNA, "snRNA", NA)
DB$"RNAcentral_snoRNA" = ifelse(DB$ID_transcript %in% RNAcentral_red_snoRNA, "snoRNA", NA)

# miRBase
cat("\t- miRBase (precursor-miRNAs)...\n")
miRBase_prec = read.table(paste0(WD, "/STEP5/", Batch, "/miRBase/prec_ID.txt"), quote = "\"")
miRBase_prec = miRBase_prec$V1

# Add miRBase info to the database
DB$"miRBase" = ifelse(DB$ID_transcript %in% miRBase_prec, "precursor-miRNA", NA)

# PmiREN
cat("\t- PmiREN (precursor-miRNAs)...\n")
PmiREN_prec = read.table(paste0(WD, "/STEP5/", Batch, "/PmiREN/prec_ID.txt"), quote = "\"")
PmiREN_prec = PmiREN_prec$V1

# Add PmiREN info to the database
DB$"PmiREN" = ifelse(DB$ID_transcript %in% PmiREN_prec, "precursor-miRNA", NA)

rm(list = c("RNAcentral_conv", "RNAcentral", "PmiREN_prec", "miRBase_prec", "RNAcentral_red", 
            "lst", "RNAcentral_red_rRNA", "RNAcentral_red_tRNA", "RNAcentral_red_snRNA", 
            "RNAcentral_red_snoRNA"))


### 2.6 STEP 6
cat("LOAD STEP 6 DATA...\n")
Cantatadb = read.table(paste0(WD, "/STEP6/", Batch, "/CANTATAdb/output_blastn.tsv"), header = T, sep = "\t", quote = "\"")
Plncdb = read.table(paste0(WD, "/STEP6/", Batch, "/PLncDB/output_blastn.tsv"), header = T, sep = "\t", quote = "\"")
Greenc = read.table(paste0(WD, "/STEP6/", Batch, "/GreeNC/output_blastn.tsv"), header = T, sep = "\t", quote = "\"")

DB$"CANTATAdb" = ifelse(DB$ID_transcript %in% Cantatadb$qseqid, "YES", "NO")
DB$"PLncDB" = ifelse(DB$ID_transcript %in% Plncdb$qseqid, "YES", "NO")
DB$"GreeNC" = ifelse(DB$ID_transcript %in% Greenc$qseqid, "YES", "NO")

rm(list = c("Cantatadb", "Plncdb", "Greenc"))


### 2.7 N's
cat("CALCULATE N'S PERCENTAGE...\n")
seq = read.table(paste0(WD, "/STEP1/", Batch, "/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta"), sep = "\t", header = F, quote = "\"")
SEQ = data.frame(ID = gsub(">", "", seq$V1[grepl(">", seq$V1)]), seq = seq$V1[!grepl(">", seq$V1)])
SEQ$"Ns" = (str_count(toupper(SEQ$seq), "N")*100)/nchar(toupper(SEQ$seq))
DB = merge(DB, SEQ[, c("ID", "Ns")], by.x = "ID_transcript", by.y = "ID")

rm(list = c("seq", "SEQ"))


### 2.8 ASSIGN CONFIDENCE LEVEL
cat("ASSIGN CONFIDENCE LEVELS TO THE POTENTIAL LNCRNAS...\n")
P = rowSums(DB[,c("CPC2", "FEELnc", "CPAT")] == "NC")
B = rowSums(DB[,c("SwissProt", "Pfam")] == "NC")

DB$"Significance_level" = ifelse(P == 3 & B == 2,
                                 "High",
                                 ifelse((P == 3 & B == 1) | (P == 2 & B == 2),
                                        "Medium",
                                        ifelse((P == 1 & B == 1) | (P == 2 & B == 0) | (P == 2 & B == 1),
                                               "Low",
                                               "Unclassified")))

rm(list = c("P", "B"))


### 2.9 ADD A CODE TO IDENTIFY THE SEQUENCES AS LNCRNAS
cat("ASSIGN AN OWN ID TO THE POTENTIAL LNCRNAS...\n")
L = DB[DB$Significance_level != "Unclassified", c("ID_transcript", "Significance_level", "RNAcentral_rRNA", "RNAcentral_tRNA", "RNAcentral_snRNA", "RNAcentral_snoRNA", "miRBase", "PmiREN")]
L = L[is.na(L$RNAcentral_rRNA),]
L = L[is.na(L$RNAcentral_tRNA),]
L = L[is.na(L$RNAcentral_snRNA),]
L = L[is.na(L$RNAcentral_snoRNA),]
L = L[is.na(L$miRBase),]
L = L[is.na(L$PmiREN),]
L$"ID" = paste0(toupper(specie), "-LNC", 1:length(L$ID_transcript))
DB = merge(DB, L[,c("ID_transcript","ID")], by = "ID_transcript", all.x = T, all.y = F)
DB_filt = DB[DB$Significance_level != "Unclassified",]
rownames(DB_filt) = NULL

write.table(DB, paste0(WD, "/STEP-FINAL/", Batch, "/Database/Database_ALL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(DB_filt, paste0(WD, "/STEP-FINAL/", Batch, "/Database/Database_NC.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

rm(list = c("L"))



## 3. FIGURES

### 3.1 VennDiagram: Non-coding transcripts according to predictor software.
# https://rdrr.io/bioc/limma/man/venn.html
# https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/venn
cat("CREATE NON-CODING PREDICTION VENNDIAGRAM...\n")

#DB = read.table(paste0(WD, "/STEP-FINAL/", Batch, "/Database/Database_ALL.tsv"), sep = "\t", header = T, quote = "\"")

x = list(
  CPC2 = unique(DB[DB$CPC2 == "NC", "ID_transcript"]),
  FEELnc = unique(DB[DB$FEELnc == "NC", "ID_transcript"]),
  CPAT = unique(DB[DB$CPAT == "NC", "ID_transcript"]),
  PFAM = unique(DB[DB$Pfam == "NC", "ID_transcript"]),
  SwissProt = unique(DB[DB$SwissProt == "NC", "ID_transcript"])
)

A = list_to_matrix(x)
Z = vennCounts(A, include = "both")

png(filename=paste0(WD, '/STEP-FINAL/', Batch, '/Figures/venn_diagram_NC_prediction.png'), height = 5000, width = 5000, res=600)
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


### 3.2 Non-Coding RNA Classification
# https://github.com/hms-dbmi/UpSetR
cat("CREATE THE NON-CODING CLASSIFICATION FIGURE...\n")

#DB_filt = read.table(paste0(WD, "/STEP-FINAL/", Batch, "/Database/Database_NC.tsv"), sep = "\t", header = T, quote = "\"")

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

png(filename=paste0(WD, '/STEP-FINAL/', Batch, '/Figures/Upset_Classification.png'), height = 5500, width = 7000, res=600)
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


### 3.3 Class code LncRNA Classification
cat("CREATE THE CLASS CODE LNCRNA CLASSIFICATION FIGURE...\n")

#DB_filt = read.table(paste0(WD, "/STEP-FINAL/", Batch, "/Database/Database_NC.tsv"), sep = "\t", header = T, quote = "\"")

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

id = x[x$lncRNA == 1 | (x$known.lncRNA == 1 & rowSums(x[,2:6]) == 0), "ID_transcript"]
DB_filt_filt = DB_filt[DB_filt$ID_transcript %in% id,]

write.table(DB_filt_filt, paste0(WD, "/STEP-FINAL/", Batch, "/Database/Database_LncRNAs.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

#DB_filt_filt = read.table(paste0(WD, "/STEP-FINAL/", Batch, "/Database/Database_LncRNAs.tsv"), sep = "\t", header = T, quote = "\"")

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

ggsave(paste0(WD, "/STEP-FINAL/", Batch, "/Figures/Class_code_lncRNAs_distribution.png"), height = 7, width = 8, dpi = 600)

rm(list = c("x", "id", "Tab", "Tab_FINAL", "gg"))


