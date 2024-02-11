################################################################################
#
# MOTIF LEVEL: CREATE FIGURES AND TABLES (ADDITIONAL)
#
# This script is used to create some summary tables and figures of motif 
# conservation at both the lncRNA and family level.
#
# @author: pasviber - Pascual Villalba Bermell
#
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggridges"))
suppressMessages(library("pRoloc"))

options(bitmapType='cairo')
options(stringsAsFactors = F)

## 1. VARIABLES

## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 10) {
  stop("At least 10 arguments must be supplied.", call.=FALSE)
} else {
  WD_in_1 = args[1]
  WD_in_2 = args[2]
  WD_out = args[3]
  classes = unlist(strsplit(args[4], " "))
  confidences = unlist(strsplit(args[5], " "))
  strictnesses = unlist(strsplit(args[6], " "))
  nonmatches = unlist(strsplit(args[7], " "))
  widths = unlist(strsplit(args[8], " "))
  modes = unlist(strsplit(args[9], " "))
  n_sim = as.numeric(args[10])
}

# WD_in_1 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Motif_level/nr/05-Summary"
# WD_in_2 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Positional_level/nr/05-Figures_and_tables"
# WD_out = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/11-Comparative_genomics/Motif_level/nr/06-Figures_and_tables"
# classes = c("intergenic", "antisense", "intronic", "sense")
# confidences = c("Low", "Medium", "High")
# strictnesses = c("ORIGINAL")
# nonmatches = c("no")
# widths = c("6-15", "6-50")
# modes = c("oops")
# n_sim = 50





################################################################################
################################################################################
## 2. CREATE THE DIRECTORY.

cat(paste0("\n\nBUILD THE DIRECTORY\n\n"))

if (!dir.exists(WD_out)){
  dir.create(WD_out)
}
if (!dir.exists(paste0(WD_out, "/Tables"))){
  dir.create(paste0(WD_out, "/Tables"))
}
if (!dir.exists(paste0(WD_out, "/Figures"))){
  dir.create(paste0(WD_out, "/Figures"))
}
for (st in strictnesses) {
  for (no in nonmatches) {
    
    ## Tables.
    if (!dir.exists(paste0(WD_out, "/Tables/", st))){
      dir.create(paste0(WD_out, "/Tables/", st))
    }
    if (!dir.exists(paste0(WD_out, "/Tables/", st, "/", no))){
      dir.create(paste0(WD_out, "/Tables/", st, "/", no))
    }
    
    ## Figures.
    if (!dir.exists(paste0(WD_out, "/Figures/", st))){
      dir.create(paste0(WD_out, "/Figures/", st))
    }
    if (!dir.exists(paste0(WD_out, "/Figures/", st, "/", no))){
      dir.create(paste0(WD_out, "/Figures/", st, "/", no))
    }
  }
}



################################################################################
################################################################################
## 4. NUMBER/PERCENTAGE OF SYNTENIC FAMILIES WITH MOTIFS

#### 4.1 Create the tables.

cat(paste0("\n\nNUMBER/PERCENTAGE OF SYNTENIC FAMILIES WITH MOTIFS --> TABLES\n\n"))

for (st in strictnesses) {
  for (no in nonmatches) {
    
    cat(paste0("\nStrictness: ", st))
    cat(paste0("\nNonmatch: ", no))
    
    ## Positional level conservation table.
    Poslev = read.table(paste0(WD_in_2, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
    Poslev = Poslev[Poslev$Strictness == st & Poslev$NonMatch == no,]
    
    ## Create a global table
    TAB_RED = data.frame()
    for (co in confidences) {
      cat(paste0("\n-", co, "\n"))
      for (cl in classes) {
        cat(paste0("\t-", cl, "\n"))
        for (mo in modes) {
          cat(paste0("\t\t-", mo, "\n"))
          for (wi in widths) {
            cat(paste0("\t\t\t-", wi, "\n"))
            
            #### LOAD TABLES.
            ## Meme.
            meme = read.table(paste0(WD_in_1, "/MEME/", st, "/", no, "/", co, "/", cl, "/", mo, "/", wi, "/MEME-SUMMARY.tsv"), sep = "\t", header = T, quote = "\"")
            meme$"Mode" = mo
            meme$"Width" = wi
            meme$"Confidence" = co
            meme$"Comparison" = cl
            
            ## Real families.
            Families_REAL = meme[meme$Type == "REAL", c("Family", "LncRNA")]
            Families_REAL = Families_REAL[!duplicated(Families_REAL),]
            
            ## Positional level conservation info.
            Poslev_mod = Poslev[Poslev$Confidence == co & Poslev$Class == cl, c("Member", "Type", "Conserved_level", "Number_species_by_family")]
            colnames(Poslev_mod) = c("LncRNA", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3")
            Poslev_mod = Poslev_mod[Poslev_mod$Type.Conservation.3 != 1,]
            Poslev_mod = merge(Poslev_mod, Families_REAL, by = "LncRNA", all = T)
            Poslev_mod = Poslev_mod[, c("Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3")]
            Poslev_mod = Poslev_mod[!duplicated(Poslev_mod),]
            rownames(Poslev_mod) = NULL
            
            #### JOIN TABLES
            ## Join the tables: Meme and Poslev_mod.
            tab = merge(meme, Poslev_mod, by = "Family", all = T)
            
            ## Convert to factors.
            tab$Type = factor(tab$Type, levels = c("REAL", paste0("SIMULATION_", 1:n_sim)))
            tab$Family = factor(tab$Family, levels = c(paste0("Fam", 1:length(unique(tab$Family)))))
            tab$Type.Conservation.1 = factor(tab$Type.Conservation.1, levels = c("Conserved"))
            tab$Type.Conservation.2 = factor(tab$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
            tab$Type.Conservation.3 = factor(tab$Type.Conservation.3, levels = 2:9)
            
            ## Get number of motifs.
            tab = tab[, c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3", "Meme_Motif.Identifier")]
            tab[tab == ""] = NA
            tab = tab[!duplicated(tab),]
            tab$"Count" = ifelse(!is.na(tab$Meme_Motif.Identifier), 1, 0)
            tab_red = tab %>% group_by(Confidence, Comparison, Mode, Width, Type, Family, Type.Conservation.1, Type.Conservation.2, Type.Conservation.3) %>% summarise(Motifs.Number = sum(Count))
            tab_red = as.data.frame(tab_red)
            tab_red$"Motifs.Number.Binary" = ifelse(tab_red$Motifs.Number > 0, 1, 0)
            
            ### FINAL JOIN.
            TAB_RED = rbind(TAB_RED, tab_red)
          }
        }
      }
    }
    
    rm(list = c("meme", "Families_REAL", "Poslev_mod", "Poslev", "tab", "tab_red", "co", "cl", "mo", "wi"))
    
    ## Convert to factors.
    TAB_RED$Confidence = factor(TAB_RED$Confidence, levels = confidences)
    TAB_RED$Comparison = factor(TAB_RED$Comparison, levels = classes)
    TAB_RED$Width = factor(TAB_RED$Width, levels = widths)
    TAB_RED$Type = factor(TAB_RED$Type, levels = c("REAL", paste0("SIMULATION_", 1:n_sim)))
    TAB_RED$Type.Conservation.1 = factor(TAB_RED$Type.Conservation.1, levels = c("Conserved"))
    TAB_RED$Type.Conservation.2 = factor(TAB_RED$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
    TAB_RED$Type.Conservation.3 = factor(TAB_RED$Type.Conservation.3, levels = 2:9)
    
    ## Collapse the table 1 (Taking into account Comparison: ALL, intergenic, antisense, intronic, sense).
    TAB_1_temp = TAB_RED %>% 
      group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.1, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = n_distinct(Family), 
        Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
        Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
        Mean.Accum.Motifs = mean(Motifs.Number)
      )
    TAB_1_temp = as.data.frame(TAB_1_temp)
    TAB_1_temp$"Type_mod" = ifelse(TAB_1_temp$Type == "REAL", "REAL", "SIMULATION")
    TAB_1 = TAB_1_temp %>% 
      group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.1, Families.Number.Total, .drop=FALSE) %>% 
      summarise(
        Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
        Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
        Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
      )
    TAB_1 = as.data.frame(TAB_1)
    
    rm(list = c("TAB_1_temp"))
    
    ## Collapse the table 1 (Taking into account Comparison_mod: ALL, CLASSES).
    TAB_1_temp = TAB_RED %>% 
      group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.1, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = n_distinct(Family), 
        Families.Number.With.Motifs = sum(Motifs.Number.Binary)
      )
    TAB_1_temp = as.data.frame(TAB_1_temp)
    TAB_1_temp$"Comparison_mod" = ifelse(TAB_1_temp$Comparison == "ALL", "ALL", "CLASSES")
    TAB_1_temp_temp = TAB_1_temp %>% 
      group_by(Confidence, Comparison_mod, Mode, Width, Type, Type.Conservation.1, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = sum(Families.Number.Total), 
        Families.Number.With.Motifs = sum(Families.Number.With.Motifs),
        Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2)
      )
    TAB_1_temp_temp = as.data.frame(TAB_1_temp_temp)
    TAB_1_temp_temp$"Type_mod" = ifelse(TAB_1_temp_temp$Type == "REAL", "REAL", "SIMULATION")
    TAB_1_mod = TAB_1_temp_temp %>% 
      group_by(Confidence, Comparison_mod, Mode, Width, Type_mod, Type.Conservation.1, Families.Number.Total, .drop=FALSE) %>% 
      summarise(
        Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
        Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2)
      )
    TAB_1_mod = as.data.frame(TAB_1_mod)
    
    rm(list = c("TAB_1_temp", "TAB_1_temp_temp"))
    
    ## Collapse the table 2 (Taking into account Comparison: ALL, intergenic, antisense, intronic, sense).
    TAB_2_temp = TAB_RED %>% 
      group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.2, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = n_distinct(Family), 
        Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
        Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
        Mean.Accum.Motifs = mean(Motifs.Number)
      )
    TAB_2_temp = as.data.frame(TAB_2_temp)
    TAB_2_temp$"Type_mod" = ifelse(TAB_2_temp$Type == "REAL", "REAL", "SIMULATION")
    TAB_2 = TAB_2_temp %>% 
      group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.2, Families.Number.Total, .drop=FALSE) %>% 
      summarise(
        Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
        Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
        Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
      )
    TAB_2 = as.data.frame(TAB_2)
    
    rm(list = c("TAB_2_temp"))
    
    ## Collapse the table 2 (Taking into account Comparison_mod: ALL, CLASSES).
    TAB_2_temp = TAB_RED %>% 
      group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.2, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = n_distinct(Family), 
        Families.Number.With.Motifs = sum(Motifs.Number.Binary)
      )
    TAB_2_temp = as.data.frame(TAB_2_temp)
    TAB_2_temp$"Comparison_mod" = ifelse(TAB_2_temp$Comparison == "ALL", "ALL", "CLASSES")
    TAB_2_temp_temp = TAB_2_temp %>% 
      group_by(Confidence, Comparison_mod, Mode, Width, Type, Type.Conservation.2, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = sum(Families.Number.Total), 
        Families.Number.With.Motifs = sum(Families.Number.With.Motifs),
        Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2)
      )
    TAB_2_temp_temp = as.data.frame(TAB_2_temp_temp)
    TAB_2_temp_temp$"Type_mod" = ifelse(TAB_2_temp_temp$Type == "REAL", "REAL", "SIMULATION")
    TAB_2_mod = TAB_2_temp_temp %>% 
      group_by(Confidence, Comparison_mod, Mode, Width, Type_mod, Type.Conservation.2, Families.Number.Total, .drop=FALSE) %>% 
      summarise(
        Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
        Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2)
      )
    TAB_2_mod = as.data.frame(TAB_2_mod)
    
    rm(list = c("TAB_2_temp", "TAB_2_temp_temp"))
    
    ## Collapse the table 3 (Taking into account Comparison: ALL, intergenic, antisense, intronic, sense).
    TAB_3_temp = TAB_RED %>% 
      group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.3, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = n_distinct(Family), 
        Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
        Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
        Mean.Accum.Motifs = mean(Motifs.Number)
      )
    TAB_3_temp = as.data.frame(TAB_3_temp)
    TAB_3_temp$"Type_mod" = ifelse(TAB_3_temp$Type == "REAL", "REAL", "SIMULATION")
    TAB_3 = TAB_3_temp %>% 
      group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.3, Families.Number.Total, .drop=FALSE) %>% 
      summarise(
        Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
        Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
        Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
      )
    TAB_3 = as.data.frame(TAB_3)
    
    rm(list = c("TAB_3_temp"))
    
    ## Collapse the table 3 (Taking into account Comparison_mod: ALL, CLASSES).
    TAB_3_temp = TAB_RED %>% 
      group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.3, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = n_distinct(Family), 
        Families.Number.With.Motifs = sum(Motifs.Number.Binary)
      )
    TAB_3_temp = as.data.frame(TAB_3_temp)
    TAB_3_temp$"Comparison_mod" = ifelse(TAB_3_temp$Comparison == "ALL", "ALL", "CLASSES")
    TAB_3_temp_temp = TAB_3_temp %>% 
      group_by(Confidence, Comparison_mod, Mode, Width, Type, Type.Conservation.3, .drop=FALSE) %>% 
      summarise(
        Families.Number.Total = sum(Families.Number.Total), 
        Families.Number.With.Motifs = sum(Families.Number.With.Motifs),
        Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2)
      )
    TAB_3_temp_temp = as.data.frame(TAB_3_temp_temp)
    TAB_3_temp_temp$"Type_mod" = ifelse(TAB_3_temp_temp$Type == "REAL", "REAL", "SIMULATION")
    TAB_3_mod = TAB_3_temp_temp %>% 
      group_by(Confidence, Comparison_mod, Mode, Width, Type_mod, Type.Conservation.3, Families.Number.Total, .drop=FALSE) %>% 
      summarise(
        Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
        Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2)
      )
    TAB_3_mod = as.data.frame(TAB_3_mod)
    
    rm(list = c("TAB_3_temp", "TAB_3_temp_temp"))
    
    ## Change NA by 0.
    TAB_1[is.na(TAB_1)] = 0
    TAB_2[is.na(TAB_2)] = 0
    TAB_3[is.na(TAB_3)] = 0
    TAB_1_mod[is.na(TAB_1_mod)] = 0
    TAB_2_mod[is.na(TAB_2_mod)] = 0
    TAB_3_mod[is.na(TAB_3_mod)] = 0
    
    ## Save the tables.
    write.table(TAB_1, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_1.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_2, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_2.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_3, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_3.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_1_mod, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_1_mod.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_2_mod, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_2_mod.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_3_mod, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_3_mod.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
}

rm(list = c("st", "no", "TAB_1", "TAB_2", "TAB_3", "TAB_1_mod", "TAB_2_mod", "TAB_3_mod", "TAB_RED"))



################################################################################
################################################################################
## 6. GC CONTENT, LENGTH AND E-VALUE MOTIFS ANALYSIS

#### 6.1 Create the tables.

cat(paste0("\n\nGC CONTENT, LENGTH AND E-VALUE MOTIFS ANALYSIS --> TABLES\n\n"))

for (st in strictnesses) {
  for (no in nonmatches) {
    
    cat(paste0("\nStrictness: ", st))
    cat(paste0("\nNonmatch: ", no))
    
    ## Positional level conservation table.
    Poslev = read.table(paste0(WD_in_2, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
    Poslev = Poslev[Poslev$Strictness == st & Poslev$NonMatch == no,]
    
    ## Create a global table
    TAB_GC = data.frame()
    TAB_Length = data.frame()
    TAB_E_value = data.frame()
    for (co in confidences) {
      cat(paste0("\n-", co, "\n"))
      for (cl in classes) {
        cat(paste0("\t-", cl, "\n"))
        for (mo in modes) {
          cat(paste0("\t\t-", mo, "\n"))
          for (wi in widths) {
            cat(paste0("\t\t\t-", wi, "\n"))
            
            #### LOAD TABLES.
            ## Meme.
            meme = read.table(paste0(WD_in_1, "/MEME/", st, "/", no, "/", co, "/", cl, "/", mo, "/", wi, "/MEME-SUMMARY.tsv"), sep = "\t", header = T, quote = "\"")
            meme$"Mode" = mo
            meme$"Width" = wi
            meme$"Confidence" = co
            meme$"Comparison" = cl
            
            ## Real families.
            Families_REAL = meme[meme$Type == "REAL", c("Family", "LncRNA")]
            Families_REAL = Families_REAL[!duplicated(Families_REAL),]
            
            ## Positional level conservation info.
            Poslev_mod = Poslev[Poslev$Confidence == co & Poslev$Class == cl, c("Member", "Type", "Conserved_level", "Number_species_by_family")]
            colnames(Poslev_mod) = c("LncRNA", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3")
            Poslev_mod = Poslev_mod[Poslev_mod$Type.Conservation.3 != 1,]
            Poslev_mod = merge(Poslev_mod, Families_REAL, by = "LncRNA", all = T)
            Poslev_mod = Poslev_mod[, c("Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3")]
            Poslev_mod = Poslev_mod[!duplicated(Poslev_mod),]
            rownames(Poslev_mod) = NULL
            
            #### JOIN TABLES
            ## Join the tables: Meme and Poslev_mod.
            tab = merge(meme, Poslev_mod, by = "Family", all = T)
            
            ## Convert to factors.
            tab$Type = factor(tab$Type, levels = c("REAL", paste0("SIMULATION_", 1:n_sim)))
            tab$Family = factor(tab$Family, levels = c(paste0("Fam", 1:length(unique(tab$Family)))))
            tab$Type.Conservation.1 = factor(tab$Type.Conservation.1, levels = c("Conserved"))
            tab$Type.Conservation.2 = factor(tab$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
            tab$Type.Conservation.3 = factor(tab$Type.Conservation.3, levels = 2:9)
            
            ## Get number of motifs.
            tab = tab[, c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "LncRNA", "Type.Conservation.1", "Type.Conservation.2", 
                          "Type.Conservation.3", "Meme_Motif.Identifier", "Meme_Motif.LncRNA", "Meme_Motif.Width", "Meme_Motif.E.value", "Meme_Motif.Sites")]
            tab[tab == ""] = NA
            tab = tab[!is.na(tab$Meme_Motif.Identifier),]
            tab = tab[!duplicated(tab),]
            
            ## GC
            tab$"Meme_Motif.GC" = round(((nchar(gsub("[AT]", "", tab$Meme_Motif.LncRNA)) * 100)/tab$Meme_Motif.Width), 2)
            tab_GC = tab %>% 
              group_by(Confidence, Comparison, Mode, Width, Type, Family, Type.Conservation.1, Type.Conservation.2, Type.Conservation.3, Meme_Motif.Identifier) %>% 
              summarise(Mean.Meme_Motif.GC = mean(Meme_Motif.GC))
            tab_GC = as.data.frame(tab_GC)
            tab_GC$"Type_mod" = ifelse(tab_GC$Type == "REAL", "REAL", "SIMULATION")
            
            ## Length
            tab_Length = tab[,c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3", "Meme_Motif.Identifier", "Meme_Motif.Width")]
            tab_Length = tab_Length[!duplicated(tab_Length),]
            tab_Length$"Log2.Meme_Motif.Width" = log(tab_Length$Meme_Motif.Width, 2)
            tab_Length$"Type_mod" = ifelse(tab_Length$Type == "REAL", "REAL", "SIMULATION")
            
            ## E-value
            tab_E_value = tab[,c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3", "Meme_Motif.Identifier", "Meme_Motif.E.value")]
            tab_E_value = tab_E_value[!duplicated(tab_E_value),]
            tab_E_value$"-Log10.Meme_Motif.E.value" = -log(tab_E_value$Meme_Motif.E.value, 10)
            tab_E_value$"-Log2.Meme_Motif.E.value" = -log(tab_E_value$Meme_Motif.E.value, 2)
            tab_E_value$"Type_mod" = ifelse(tab_E_value$Type == "REAL", "REAL", "SIMULATION")

            ### FINAL JOIN.
            TAB_GC = rbind(TAB_GC, tab_GC)
            TAB_Length = rbind(TAB_Length, tab_Length)
            TAB_E_value = rbind(TAB_E_value, tab_E_value)
          }
        }
      }
    }
    
    rm(list = c("meme", "Families_REAL", "Poslev_mod", "Poslev", "tab", "tab_GC", "tab_Length", "tab_E_value", "co", "cl", "mo", "wi"))
    
    TAB_GC_mod = TAB_GC
    TAB_Length_mod = TAB_Length
    TAB_E_value_mod = TAB_E_value
    
    TAB_GC_mod$"Comparison_mod" = ifelse(TAB_GC_mod$Comparison == "ALL", "ALL", "CLASSES")
    TAB_Length_mod$"Comparison_mod" = ifelse(TAB_Length_mod$Comparison == "ALL", "ALL", "CLASSES")
    TAB_E_value_mod$"Comparison_mod" = ifelse(TAB_E_value_mod$Comparison == "ALL", "ALL", "CLASSES")
    
    ## Save the tables.
    write.table(TAB_GC, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_GC.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_Length, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_LENGTH.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_E_value, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_E_VALUE.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_GC_mod, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_GC_mod.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_Length_mod, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_LENGTH_mod.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_E_value_mod, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_E_VALUE_mod.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
}

rm(list = c("st", "no", "TAB_GC", "TAB_Length", "TAB_E_value", "TAB_GC_mod", "TAB_Length_mod", "TAB_E_value_mod"))


