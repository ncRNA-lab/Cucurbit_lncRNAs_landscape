################################################################################
#
# CREATE THE NON-REDUNDANT LNCRNA DATABASE
#
# Create a lncRNA database from the information obtained in the lncRNA prediction
# pipeline.
#
################################################################################


## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Create_customed_LncRNA_db_NR.R $WD2

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("limma"))
suppressMessages(library("Cairo"))
suppressMessages(options(bitmapType='cairo'))
suppressMessages(library("UpSetR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("ggvenn"))


## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("At least 1 argument must be supplied.", call.=FALSE)
} else {
  WD = args[1]
}

#WD = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/05-predict_lncRNAs/car"

## 2. CREATE THE NON-REDUNDANT LNCRNAS DATABASE

### Load redundant LncRNAs database
DB = read.table(paste0(WD, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), sep = "\t", header = T, quote = "\"")

### Load results of redundancy analysis
AGAT = read.table(paste0(WD, "/STEP-FINAL/Redundancy_analysis/AGAT/POTENTIAL_LNCRNAS_filt_sort_ids.txt"), header = F, quote = "\"")
AGAT = AGAT$V1
CGAT = read.table(paste0(WD, "/STEP-FINAL/Redundancy_analysis/CGAT/POTENTIAL_LNCRNAS_sort_filt_sort_ids.txt"), header = F, quote = "\"")
CGAT = CGAT$V1
CDHIT90 = read.table(paste0(WD, "/STEP-FINAL/Redundancy_analysis/CDHIT/POTENTIAL_LNCRNAS_filt_90_ids.txt"), header = F, quote = "\"")
CDHIT90 = CDHIT90$V1
CDHIT95 = read.table(paste0(WD, "/STEP-FINAL/Redundancy_analysis/CDHIT/POTENTIAL_LNCRNAS_filt_95_ids.txt"), header = F, quote = "\"")
CDHIT95 = CDHIT95$V1

### Add CGAT info to the database.
if ("CGAT" %in% colnames(DB)) {
  DB$CGAT = NULL
} 

DB$"CGAT" = ifelse(DB$ID_transcript %in% CGAT, "non-redundant", "redundant")
DB = DB[, c("ID_transcript", "Chr", "Origin", "Start", "End", "Strand", "ID_Gene", "Class_code", "Exons", 
            "Length", "GC", "CPC2", "CPAT", "FEELnc", "SwissProt", "Pfam", "ORF.80", "ORF.100", "ORF.120", 
            "RNAcentral_rRNA", "RNAcentral_tRNA", "RNAcentral_snRNA","RNAcentral_snoRNA", "miRBase", 
            "PmiREN", "CANTATAdb", "PLncDB", "GreeNC", "Ns", "CGAT", "Significance_level", "ID")]


write.table(DB, paste0(WD, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

### Select non-redundant lncRNAs.
DB_nr = DB[DB$CGAT == "non-redundant",]
write.table(DB_nr, paste0(WD, "/STEP-FINAL/Database/Database_LncRNAs_NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


## 3. FIGURES

if (!dir.exists(paste0(WD, "/STEP-FINAL/Redundancy_analysis/Figures"))){
  dir.create(paste0(WD, "/STEP-FINAL/Redundancy_analysis/Figures"))
}

### 3.1 Prediction tools
#### 3.1.1 VennDiagram: ALL LncRNAs according to the different tools.
# https://rdrr.io/bioc/limma/man/venn.html
# https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/venn

cat("VENNDIAGRAM: ALL LNCRNAS ACCORDING TO THE PREDICTION TOOLS...\n")

x = list(
  CPC2 = unique(DB[DB$CPC2 == "NC", "ID_transcript"]),
  FEELnc = unique(DB[DB$FEELnc == "NC", "ID_transcript"]),
  CPAT = unique(DB[DB$CPAT == "NC", "ID_transcript"]),
  PFAM = unique(DB[DB$Pfam == "NC", "ID_transcript"]),
  SwissProt = unique(DB[DB$SwissProt == "NC", "ID_transcript"])
)

A = list_to_matrix(x)
Z = vennCounts(A, include = "both")

png(filename=paste0(WD, '/STEP-FINAL/Redundancy_analysis/Figures/LncRNAs_R_and_NR_prediction_tools.png'), height = 5000, width = 5000, res=600)
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

#### 3.1.2 VennDiagram: LncRNAs considered as Redundant according to the different tools.
# https://rdrr.io/bioc/limma/man/venn.html
# https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/venn

cat("VENNDIAGRAM: LNCRNAS CONSIDERED AS NON-REDUNDANT ACCORDING TO THE PREDICTION TOOLS...\n")

x = list(
  CPC2 = unique(DB_nr[DB_nr$CPC2 == "NC", "ID_transcript"]),
  FEELnc = unique(DB_nr[DB_nr$FEELnc == "NC", "ID_transcript"]),
  CPAT = unique(DB_nr[DB_nr$CPAT == "NC", "ID_transcript"]),
  PFAM = unique(DB_nr[DB_nr$Pfam == "NC", "ID_transcript"]),
  SwissProt = unique(DB_nr[DB_nr$SwissProt == "NC", "ID_transcript"])
)

A = list_to_matrix(x)
Z = vennCounts(A, include = "both")

png(filename=paste0(WD, '/STEP-FINAL/Redundancy_analysis/Figures/LncRNAs_NR_prediction_tools.png'), height = 5000, width = 5000, res=600)
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


### 3.2 Redundancy analysis tools

cat("VENNDIAGRAM: NON-REDUNDANT LNCRNAS ACCORDING TO THE TOOLS...\n")

x = list(
  AGAT = AGAT,
  CGAT = CGAT,
  CDHIT90 = CDHIT90,
  CDHIT95 = CDHIT95
)

png(filename = paste0(WD, '/STEP-FINAL/Redundancy_analysis/Figures/LncRNAs_NR_redundancy_tools_4.png'), width = 3000, height = 3000, res = 500)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
invisible(dev.off())

x = list(
  AGAT = AGAT,
  CGAT = CGAT
)

png(filename = paste0(WD, '/STEP-FINAL/Redundancy_analysis/Figures/LncRNAs_NR_redundancy_tools_2.png'), width = 3000, height = 3000, res = 500)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)
invisible(dev.off())


### 3.3 Class code LncRNA Classification
cat("CREATE THE CLASS CODE LNCRNA CLASSIFICATION FIGURE...\n")

Tab = DB_nr[, c("ID_transcript", "Class_code")]
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

ggsave(paste0(WD, "/STEP-FINAL/Redundancy_analysis/Figures/Class_code_lncRNAs_distribution_NR.png"), height = 7, width = 8, dpi = 600)

rm(list = c("x", "Tab", "Tab_FINAL", "gg"))


