################################################################################
#
# Percentage of transcripts which cover at least the 50 percent (this value can 
# be changed in the script) of some kind of transposon or repetitive element that 
# can be found in the genome of each specie. We take always into account the class
# code of each transcript.
#
################################################################################

rm(list = ls())

## EXECUTION IN COMMAND LINE:
#EXAMPLE --> Create_Tables_and_Figures_A.R path_res folder_new folder_RepMask flag confidence

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggExtra))
suppressMessages(options(bitmapType='cairo'))
options(stringsAsFactors = F)

## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments must be supplied.", call.=FALSE)
} else {
  path_res = args[1]
  folder_new = args[2]
  folder_RepMask = args[3]
  flag = args[4]
  confidence = args[5]
}

# path_res = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results"
# folder_new = "02-Comparison_Genes_LncRNAs"
# folder_RepMask = "05-RepeatMasker"
# flag = "NR"
# confidence = "High"


## 2. PIPELINE

# Species
species = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
species_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")
species_tab = data.frame(spe = species, name = species_name)


### 2.1 CREATE FINAL TABLE

cat("\nCreating the FINAL table...\n")

TAB_FINAL = data.frame()
for (spe in species) {
  cat(paste0("\t", spe, "...\n"))
  
  ### Paths.
  WD_01 = paste0(path_res, "/11-TEs_and_genomic_repeats/", folder_new)
  WD_02 = paste0(path_res, "/05-predict_lncRNAs")
  
  ### Directory.
  if (!dir.exists(paste0(WD_01, "/Figures_and_Tables"))){
    dir.create(paste0(WD_01, "/Figures_and_Tables"))
  }
  if (!dir.exists(paste0(WD_01, "/Figures_and_Tables/A"))){
    dir.create(paste0(WD_01, "/Figures_and_Tables/A"))
  }
  
  ### Genes, LncRNAs and IR intersected with Repeat regions.
  res_genes = read.table(paste0(WD_01, "/", spe, "/ORIGINAL_GENES_intersect_Rep.tsv"), header = F, sep = "\t", quote = "\"")
  colnames(res_genes) = c("chr", "start", "end", "transcript_id", "strand", "Repeat_id", "Repeat_type_1", "Repeat_type_2")
  if (flag == "NR") {
    res_lncRNAs = read.table(paste0(WD_01, "/", spe, "/POTENTIAL_LNCRNAS_intersect_Rep_NR.tsv"), header = F, sep = "\t", quote = "\"")
    colnames(res_lncRNAs) = c("chr", "start", "end", "transcript_id", "strand", "Repeat_id", "Repeat_type_1", "Repeat_type_2")
  }
  if (flag == "R") {
    res_lncRNAs = read.table(paste0(WD_01, "/", spe, "/POTENTIAL_LNCRNAS_intersect_Rep_R.tsv"), header = F, sep = "\t", quote = "\"")
    colnames(res_lncRNAs) = c("chr", "start", "end", "transcript_id", "strand", "Repeat_id", "Repeat_type_1", "Repeat_type_2")
  }
  res_IR = read.table(paste0(WD_01, "/", spe, "/Random_IR_intersect_Rep.tsv"), header = F, sep = "\t", quote = "\"")
  colnames(res_IR) = c("chr", "start", "end", "transcript_id", "strand", "Repeat_id", "Repeat_type_1", "Repeat_type_2")
  res_IR$transcript_id = paste0(res_IR$chr, ":", res_IR$start, "-", res_IR$end, "(", res_IR$strand, ")")
  RES = rbind(res_genes, res_lncRNAs, res_IR)
  
  ### Remove intersections with small RNAs.
  RES = RES[!(RES$Repeat_type_2 %in% c("rRNA", "tRNA", "snRNA", "snoRNA")),]
  
  ### All the genes, lncRNAs and IR.
  All_genes = read.table(paste0(WD_01, "/", spe, "/ORIGINAL_GENES.bed"), header = F, sep = "\t", quote = "\"")
  All_genes = All_genes[,c(1:4,6)]
  colnames(All_genes) = c("chr", "start", "end", "transcript_id", "strand")
  All_genes$"Type" = "Protein Coding"
  if (flag == "NR") {
    All_lncRNAs = read.table(paste0(WD_01, "/", spe, "/POTENTIAL_LNCRNAS_pred_NR.bed"), header = F, sep = "\t", quote = "\"")
    All_lncRNAs = All_lncRNAs[,c(1:4,6)]
    colnames(All_lncRNAs) = c("chr", "start", "end", "transcript_id", "strand")
    All_lncRNAs$"Type" = "LncRNA"
    LncRNAs_db = read.table(paste0(WD_02, "/", spe, "/STEP-FINAL/Database/Database_LncRNAs_NR.tsv"), header = T, sep = "\t", quote = "\"")
    LncRNAs_db = LncRNAs_db[LncRNAs_db$Significance_level == confidence,]
    All_lncRNAs = All_lncRNAs[All_lncRNAs$transcript_id %in% LncRNAs_db$ID_transcript,]
  }
  if (flag == "R") {
    All_lncRNAs = read.table(paste0(WD_01, "/", spe, "/POTENTIAL_LNCRNAS_pred_R.bed"), header = F, sep = "\t", quote = "\"")
    All_lncRNAs = All_lncRNAs[,c(1:4,6)]
    colnames(All_lncRNAs) = c("chr", "start", "end", "transcript_id", "strand")
    All_lncRNAs$"Type" = "LncRNA"
    LncRNAs_db = read.table(paste0(WD_02, "/", spe, "/STEP-FINAL/Database/Database_LncRNAs.tsv"), header = T, sep = "\t", quote = "\"")
    LncRNAs_db = LncRNAs_db[LncRNAs_db$Significance_level == confidence,]
    All_lncRNAs = All_lncRNAs[All_lncRNAs$transcript_id %in% LncRNAs_db$ID_transcript,]
  }
  All_IR = read.table(paste0(WD_01, "/", spe, "/Random_IR.bed"), header = F, sep = "\t", quote = "\"")
  All_IR = All_IR[,c(1:4,6)]
  colnames(All_IR) = c("chr", "start", "end", "transcript_id", "strand")
  All_IR$transcript_id = paste0(All_IR$chr, ":", All_IR$start, "-", All_IR$end, "(", All_IR$strand, ")")
  All_IR$"Type" = "Intergenic Region"
  ALL = rbind(All_genes, All_lncRNAs, All_IR)
  
  ### Join all the info ALL and RES. 
  TAB = merge(ALL, RES[,c("transcript_id", "Repeat_id", "Repeat_type_1", "Repeat_type_2")], by = "transcript_id", all.x = T, all.y = F)
  
  ### Convert all NA values to No_overlap. They are the genes, lncRNAs and IR which haven't overlap on some repeat region.
  TAB[is.na(TAB)] = "No repeat"
  
  ### Modify the Repeat_type_2 values to make easier the downstream analysis.
  TAB$"Repeat_type_2_mod" = ifelse(grepl("LINE", TAB$Repeat_type_2, fixed = T), "LINE/SINE",
                                   ifelse(grepl("SINE", TAB$Repeat_type_2, fixed = T), "LINE/SINE",
                                          ifelse(grepl("LTR", TAB$Repeat_type_2, fixed = T), "LTR",
                                                 ifelse(grepl("DNA/", TAB$Repeat_type_2, fixed = T), "DNA transposon",
                                                        ifelse(grepl("RC/", TAB$Repeat_type_2, fixed = T), "DNA transposon",
                                                               ifelse(TAB$Repeat_type_2 == "Simple_repeat", "Simple repeat",
                                                                      ifelse(TAB$Repeat_type_2 == "Low_complexity", "Low complexity repeat",
                                                                             ifelse(TAB$Repeat_type_2 == "No repeat", "No repeat",
                                                                                    "Unknown"))))))))
  
  ### Add class code info to the table. u (intergenic), x (antisense), i (intronic), o/e (sense), pc (protein coding) and ir (intergenic region).
  TAB = merge(TAB, LncRNAs_db[,c("ID_transcript", "Class_code")], by.x = "transcript_id", by.y = "ID_transcript", all.x = T, all.y = F)
  TAB$Class_code = ifelse(TAB$Type == "Protein Coding", "pc", ifelse(TAB$Type == "Intergenic Region", "ir", TAB$Class_code))
  TAB[TAB == "o"] = "o/e"
  TAB[TAB == "e"] = "o/e"
  
  ### Create factor classes.
  TAB$Class_code = factor(TAB$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  TAB$Repeat_type_2_mod = factor(TAB$Repeat_type_2_mod, levels = c("Simple repeat", "Low complexity repeat", "LINE/SINE", 
                                                                   "LTR", "DNA transposon", "Unknown", "No repeat"))
  TAB = TAB[order(TAB$Class_code),]
  
  ### Add specie info.
  TAB$"spe" = spe
  
  ## Join TAB to the TAB_FINAL table.
  TAB_FINAL = rbind(TAB_FINAL, TAB)
}

### Convert spe to factor.
TAB_FINAL$spe = factor(TAB_FINAL$spe, levels = species)

### Save TAB_FINAL.
write.table(TAB_FINAL, paste0(WD_01, "/Figures_and_Tables/A/Final_tab-", flag, "-", confidence, ".tsv"), col.names = T, row.names = F, sep = "\t", quote = F)

### 2.2 CREATE SUMMARY TABLE

cat("\nCreating the SUMMARY table...\n")

### Create summary table. Use n_distinct() to remove the redundancy. I mean, a transcript which overlap with
### the same kind of repeat several times because you can find this repeat sequence in several genome positions.
TAB_SUMMARY_individual = TAB_FINAL %>% 
  group_by(spe, Class_code, Repeat_type_2_mod, .drop=FALSE) %>% 
  summarise(sum.Counts = n_distinct(transcript_id))
TAB_SUMMARY_individual = as.data.frame(TAB_SUMMARY_individual)

TAB_SUMMARY_total = TAB_FINAL %>% 
  group_by(spe, Class_code, .drop=FALSE) %>% 
  summarise(Total = n_distinct(transcript_id))
TAB_SUMMARY_total = as.data.frame(TAB_SUMMARY_total)

TAB_SUMMARY = merge(TAB_SUMMARY_individual, TAB_SUMMARY_total, by = c("spe", "Class_code"), all = T)
TAB_SUMMARY = TAB_SUMMARY[order(TAB_SUMMARY$spe, TAB_SUMMARY$Class_code, TAB_SUMMARY$Repeat_type_2_mod),]

### Get percentage.
TAB_SUMMARY$"Percentage" = round((TAB_SUMMARY$sum.Counts*100)/TAB_SUMMARY$Total, 1)

### Add percentage masked genome.
masked_perc = read.table(paste0(path_res, "/11-TEs_and_genomic_repeats/01-Repeat_calling/", folder_RepMask, "/masked_genome_percentage.tsv"), header = T, sep = "\t", quote = "\"")
TAB_SUMMARY = merge(TAB_SUMMARY, masked_perc, by = "spe", all = T)

### Subdivide in two tables.
TAB_SUMMARY_1 = TAB_SUMMARY
TAB_SUMMARY_2 = TAB_SUMMARY[TAB_SUMMARY$Repeat_type_2_mod == "No repeat", c("spe", "Class_code", "Percentage")]
colnames(TAB_SUMMARY_2) = c("spe", "Class_code", "No_repeat")
TAB_SUMMARY_2$Repeat = 100 - TAB_SUMMARY_2$No_repeat
TAB_SUMMARY_2$No_repeat = NULL

### Save TAB_SUMMARY.
write.table(TAB_SUMMARY_1, paste0(WD_01, "/Figures_and_Tables/A/Summary_tab_1-", flag, "-", confidence, ".tsv"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(TAB_SUMMARY_2, paste0(WD_01, "/Figures_and_Tables/A/Summary_tab_2-", flag, "-", confidence, ".tsv"), col.names = T, row.names = F, sep = "\t", quote = F)

### Remove variables.
rm(list = c("ALL", "All_genes", "All_IR", "All_lncRNAs", "LncRNAs_db", "RES", "res_genes", 
            "res_IR", "res_lncRNAs", "spe", "TAB", "WD_02", "TAB_SUMMARY_total", 
            "TAB_SUMMARY_individual", "masked_perc", "TAB_SUMMARY"))

### 2.3 FIGURE

cat("\nDrawing the TAB SUMMARY 1 table...\n")

TAB_SUMMARY_1 = merge(TAB_SUMMARY_1, species_tab, by = "spe", all.x = T, all.y = F)
TAB_SUMMARY_1$name = factor(TAB_SUMMARY_1$name, levels = species_name)
TAB_SUMMARY_1 = TAB_SUMMARY_1[order(TAB_SUMMARY_1$name),]
TAB_SUMMARY_1$"facet" = paste0(TAB_SUMMARY_1$name, " (", TAB_SUMMARY_1$masked_perc, "%)")
TAB_SUMMARY_1$facet = factor(TAB_SUMMARY_1$facet, levels = unique(TAB_SUMMARY_1$facet))

gg1 = ggplot(TAB_SUMMARY_1, aes(x = Repeat_type_2_mod, y = Percentage, fill = Repeat_type_2_mod)) +
  geom_bar(position=position_dodge(), aes(y=Percentage), stat="identity") +
  facet_grid(facet ~ Class_code) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000")) +
  scale_y_continuous(limits = c(0, max(TAB_SUMMARY_1$Percentage) + 10), breaks = seq(0, max(TAB_SUMMARY_1$Percentage) + 10, by = 10), expand = c(0.02, 5)) +
  xlab("") +
  ylab("Transcripts with repeat content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, angle = 80)) +
  theme(legend.position = "none") +
  geom_text(aes(label=paste0(Percentage, "\n(", sum.Counts,")")), vjust=-0.2, size = 2)

ggsave(paste0(WD_01, "/Figures_and_Tables/A/Percentage_repeat_1-", flag, "-", confidence, ".png"), height = 16, width = 12, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/A/Percentage_repeat_1-", flag, "-", confidence, ".pdf"), height = 16, width = 12, dpi = 600)

cat("\nDrawing the TAB SUMMARY 2 table...\n")

TAB_SUMMARY_2 = merge(TAB_SUMMARY_2, species_tab, by = "spe", all.x = T, all.y = F)
TAB_SUMMARY_2$name = factor(TAB_SUMMARY_2$name, levels = species_name)

gg2 = ggplot(TAB_SUMMARY_2, aes(x = Class_code, y = Repeat, fill = Class_code)) +
  geom_bar(position=position_dodge(), aes(y=Repeat), colour="black", stat="identity") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10), expand = c(0.02, 0.02)) +
  facet_wrap(name~., nrow = 1) +
  xlab("") +
  ylab("Transcripts with repeat content (%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 14),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 17, face = "bold.italic"))

ggsave(paste0(WD_01, "/Figures_and_Tables/A/Percentage_repeat_2-", flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/A/Percentage_repeat_2-", flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)

gg2 = gg2 +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(paste0(WD_01, "/Figures_and_Tables/A/Percentage_repeat_2_mod-", flag, "-", confidence, ".png"), height = 5, width = 25, dpi = 600)
ggsave(paste0(WD_01, "/Figures_and_Tables/A/Percentage_repeat_2_mod-", flag, "-", confidence, ".pdf"), height = 5, width = 25, dpi = 600)

rm(list = c("gg1", "gg2"))

