################################################################################
#
# TABLES AND FIGURES
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

flag = "nr"
WD_in_1 = paste0("/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Motif_level/", flag, "/Positional_conserved/05-Summary")
WD_in_2 = paste0("/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/07-comparison_lncRNAs_vs_coding_genes/ALL")
WD_in_3 = paste0("/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Positional_level/Approach_2/", flag, "/05-Figures")
WD_out = paste0("/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Motif_level/", flag, "/Positional_conserved/06-Figures_and_tables")
classes = c("ALL", "intergenic", "antisense", "intronic", "sense")
confidences = c("Low", "Medium", "High")
strictnesses = c("ORIGINAL")
nonmatches = c("no")
widths = c("6-15", "6-50")
modes = c("oops")
n_sim = 50





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
## 3. GLOBAL TABLEs

cat(paste0("\n\nGLOBAL TABLE\n\n"))

#### 3.1 MEME, GOMO and Positional level conservation info.

for (st in strictnesses) {
  for (no in nonmatches) {
    
    cat(paste0("\nStrictness: ", st))
    cat(paste0("\nNonmatch: ", no))
    
    ## Positional level conservation table.
    Poslev = read.table(paste0(WD_in_3, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
    Poslev = Poslev[Poslev$Strictness == st & Poslev$NonMatch == no,]
    
    ## Create a global table
    TAB_all = data.frame()
    TAB_only_motifs = data.frame()
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
            
            ## Gomo.
            gomo = tryCatch(read.table(paste0(WD_in_1, "/GOMO/", st, "/", no, "/", co, "/", cl, "/", mo, "/", wi, "/GOMO-SUMMARY.tsv"), sep = "\t", header = T, quote = "\""), error=function(e) NULL)
            
            ## Positional level conservation info.
            Poslev_mod = Poslev[Poslev$Confidence == co & Poslev$Class == cl, c("Member", "Type", "Conserved_level", "Number_species_by_family")]
            colnames(Poslev_mod) = c("LncRNA", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3")
            Poslev_mod = Poslev_mod[Poslev_mod$Type.Conservation.3 != 1,]
            Poslev_mod = merge(Poslev_mod, Families_REAL, by = "LncRNA", all = T)
            Poslev_mod = Poslev_mod[, c("Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3")]
            Poslev_mod = Poslev_mod[!duplicated(Poslev_mod),]
            rownames(Poslev_mod) = NULL
            
            #### JOIN TABLES
            ## Join the tables: Meme, Gomo and Poslev_mod.
            tab = merge(meme, Poslev_mod, by = "Family", all = T)
            if (!is.null(gomo)){
              tab = merge(tab, gomo, by.x = c("Type", "Family", "Meme_Motif.Identifier"), by.y = c("Type", "Family", "Gomo_Motif.Identifier"), all = T)
            }
            if (is.null(gomo)){
              tab$"Gomo_GO.Term.Identifier" = NA
              tab$"Gomo_Score" = NA
              tab$"Gomo_P.value" = NA
              tab$"Gomo_Q.value" = NA
            }
            tab = tab[,c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "LncRNA", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3", 
                         "Meme_Motif.Identifier", "Meme_Motif.Regex", "Meme_Motif.Conserved", "Meme_Motif.E.value", "Meme_Motif.LncRNA", "Meme_Motif.Width", "Meme_Motif.Start", 
                         "Meme_Motif.P.value", "Meme_Motif.Sites", "Gomo_GO.Term.Identifier", "Gomo_Score", "Gomo_P.value", "Gomo_Q.value")]
            
            ## Change "" by NA
            tab[tab == ""] = NA
            
            ## Filter: Keep only families with at least one conserved DNA motif.
            tab_all = tab
            tab_only_motifs = tab[!is.na(tab$Meme_Motif.Identifier),]
            
            TAB_all = rbind(TAB_all, tab_all)
            TAB_only_motifs = rbind(TAB_only_motifs, tab_only_motifs)
          }
        }
      }
    }
    
    rm(list = c("meme", "gomo", "Families_REAL", "Poslev_mod", "Poslev", "tab", "tab_all", "tab_only_motifs", "co", "cl", "mo", "wi"))
    
    ## Save the table.
    write.table(TAB_all, paste0(WD_out, "/Tables/", st, "/", no, "/GLOBAL_TABLE_MEME-GOMO-POSLEV.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_only_motifs, paste0(WD_out, "/Tables/", st, "/", no, "/GLOBAL_TABLE_MEME-GOMO-POSLEV_Filtered.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
}

rm(list = c("st", "no", "TAB_all", "TAB_only_motifs"))





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
    Poslev = read.table(paste0(WD_in_3, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
    Poslev = Poslev[Poslev$Strictness == st & Poslev$NonMatch == no,]
    
    ## Create a global table
    TAB_1 = data.frame()
    TAB_2 = data.frame()
    TAB_3 = data.frame()
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
            
            ## Collapse the table 1.
            tab_red_1 = tab_red %>% 
              group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.1, .drop=FALSE) %>% 
              summarise(
                Families.Number.Total = n_distinct(Family), 
                Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
                Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
                Mean.Accum.Motifs = mean(Motifs.Number)
                )
            tab_red_1 = as.data.frame(tab_red_1)
            tab_red_1$"Type_mod" = ifelse(tab_red_1$Type == "REAL", "REAL", "SIMULATION")
            tab_red_1_red = tab_red_1 %>% 
              group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.1, Families.Number.Total, .drop=FALSE) %>% 
              summarise(
                Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
                Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
                Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
                )
            tab_red_1_red = as.data.frame(tab_red_1_red)
            
            ## Collapse the table 2.
            tab_red_2 = tab_red %>% 
              group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.2, .drop=FALSE) %>% 
              summarise(
                Families.Number.Total = n_distinct(Family), 
                Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
                Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
                Mean.Accum.Motifs = mean(Motifs.Number)
                )
            tab_red_2 = as.data.frame(tab_red_2)
            tab_red_2$"Type_mod" = ifelse(tab_red_2$Type == "REAL", "REAL", "SIMULATION")
            tab_red_2_red = tab_red_2 %>% 
              group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.2, Families.Number.Total, .drop=FALSE) %>% 
              summarise(
                Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
                Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
                Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
                )
            tab_red_2_red = as.data.frame(tab_red_2_red)
            
            ## Collapse the table 3.
            tab_red_3 = tab_red %>% 
              group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.3, .drop=FALSE) %>% 
              summarise(
                Families.Number.Total = n_distinct(Family), 
                Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
                Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
                Mean.Accum.Motifs = mean(Motifs.Number)
                )
            tab_red_3 = as.data.frame(tab_red_3)
            tab_red_3$"Type_mod" = ifelse(tab_red_3$Type == "REAL", "REAL", "SIMULATION")
            tab_red_3_red = tab_red_3 %>% 
              group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.3, Families.Number.Total, .drop=FALSE) %>% 
              summarise(
                Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
                Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
                Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
                )
            tab_red_3_red = as.data.frame(tab_red_3_red)
            
            ### FINAL JOIN.
            TAB_1 = rbind(TAB_1, tab_red_1_red)
            TAB_2 = rbind(TAB_2, tab_red_2_red)
            TAB_3 = rbind(TAB_3, tab_red_3_red)
          }
        }
      }
    }
    
    rm(list = c("meme", "Families_REAL", "Poslev_mod", "Poslev", "tab", "tab_red", "tab_red_1", "tab_red_2", "tab_red_3", "tab_red_1_red", "tab_red_2_red", "tab_red_3_red", "co", "cl", "mo", "wi"))
    
    ## Convert to factors.
    TAB_1$Confidence = factor(TAB_1$Confidence, levels = confidences)
    TAB_1$Comparison = factor(TAB_1$Comparison, levels = classes)
    TAB_1$Width = factor(TAB_1$Width, levels = widths)
    TAB_1$Type_mod = factor(TAB_1$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_2$Confidence = factor(TAB_2$Confidence, levels = confidences)
    TAB_2$Comparison = factor(TAB_2$Comparison, levels = classes)
    TAB_2$Width = factor(TAB_2$Width, levels = widths)
    TAB_2$Type_mod = factor(TAB_2$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_2$Type.Conservation.2 = factor(TAB_2$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
    TAB_3$Confidence = factor(TAB_3$Confidence, levels = confidences)
    TAB_3$Comparison = factor(TAB_3$Comparison, levels = classes)
    TAB_3$Width = factor(TAB_3$Width, levels = widths)
    TAB_3$Type_mod = factor(TAB_3$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_3$Type.Conservation.3 = factor(TAB_3$Type.Conservation.3, levels = 2:9)
    
    ## Change NA by 0.
    TAB_1[is.na(TAB_1)] = 0
    TAB_2[is.na(TAB_2)] = 0
    TAB_3[is.na(TAB_3)] = 0
    
    ## Save the tables.
    write.table(TAB_1, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_1.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_2, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_2.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_3, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_3.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
}

rm(list = c("st", "no", "TAB_1", "TAB_2", "TAB_3"))


#### 4.2 Create the figures.

cat(paste0("\n\nNUMBER/PERCENTAGE OF SYNTENIC FAMILIES WITH MOTIFS --> FIGURES\n\n"))

for (st in strictnesses) {
  cat(paste0("-", st, "\n"))
  for (no in nonmatches) {
    cat(paste0("\t-", no, "\n"))
    for (mo in modes) {
      cat(paste0("\t\t-", mo, "\n"))
      
      ## Load tables
      TAB_1 = read.table(paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_1.tsv"), sep = "\t", header = T, quote = "\"", check.names = F)
      TAB_1 = TAB_1[TAB_1$Mode == mo,]
      TAB_2 = read.table(paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_2.tsv"), sep = "\t", header = T, quote = "\"", check.names = F)
      TAB_2 = TAB_2[TAB_2$Mode == mo,]
      
      ## Convert to factors.
      TAB_1$Confidence = factor(TAB_1$Confidence, levels = confidences)
      TAB_1$Comparison = factor(TAB_1$Comparison, levels = classes)
      TAB_1$Width = factor(TAB_1$Width, levels = widths)
      TAB_1$Type_mod = factor(TAB_1$Type_mod, levels = c("REAL", "SIMULATION"))
      TAB_2$Confidence = factor(TAB_2$Confidence, levels = confidences)
      TAB_2$Comparison = factor(TAB_2$Comparison, levels = classes)
      TAB_2$Width = factor(TAB_2$Width, levels = widths)
      TAB_2$Type_mod = factor(TAB_2$Type_mod, levels = c("REAL", "SIMULATION"))
      TAB_2$Type.Conservation.2 = factor(TAB_2$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
      
      ## Number (TAB_1)
      gg1 = ggplot(TAB_1, aes(x = Width, y = Mean.Families.Number.With.Motifs)) +
        geom_point(aes(color = Type_mod), size=4, alpha=0.6) +
        scale_color_manual(values = c("#0072b5", "#bc3c29")) +
        scale_x_discrete(limits = rev) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width") + ylab("Number of syntenic families with at least one DNA motif") +
        coord_flip() +
        guides(color = guide_legend(nrow = 1, title = "Dataset")) +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_number_1.png"), height = 6, width = 12, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_number_1.pdf"), height = 6, width = 12, dpi = 600)
      
      ## Number (TAB_2)
      gg2 = with(TAB_2, TAB_2[order(Width, Type.Conservation.2),]) %>% 
        mutate(x = paste(Width, Type.Conservation.2, sep = " ")) %>%
        mutate(x = factor(x, levels = paste(rep(widths, each = length(levels(TAB_2$Type.Conservation.2))), rep(levels(TAB_2$Type.Conservation.2), length(widths), sep = " ")))) %>%
        ggplot(aes(x = x, y = Mean.Families.Number.With.Motifs)) +
        geom_point(aes(color = Type.Conservation.2, alpha = Type_mod), size = 3) +
        scale_color_manual(values = c("#6f99ad", "#efc000", "#D87758")) +
        scale_x_discrete(limits = rev) +
        scale_alpha_discrete(range = c(1, 0.5)) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width and Positional conservation level") + ylab("Number of syntenic families with at least one DNA motif") +
        coord_flip() +
        guides(color = guide_legend(nrow = 1)) +
        labs(color = "Conservation level", alpha = "Dataset") +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_number_2.png"), height = 8, width = 15, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_number_2.pdf"), height = 8, width = 15, dpi = 600)
      
      ## Percentage (TAB_1)
      gg1 = ggplot(TAB_1, aes(x = Width, y = Mean.Families.Percentage.With.Motifs)) +
        geom_point(aes(color = Type_mod), size=4, alpha=0.6) +
        scale_color_manual(values = c("#0072b5", "#bc3c29")) +
        scale_x_discrete(limits = rev) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width") + ylab("Percentage of syntenic families with at least one DNA motif") +
        coord_flip() +
        guides(color = guide_legend(nrow = 1, title = "Dataset")) +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_percentage_1.png"), height = 6, width = 12, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_percentage_1.pdf"), height = 6, width = 12, dpi = 600)
      
      ## Percentage (TAB_2)
      gg2 = with(TAB_2, TAB_2[order(Width, Type.Conservation.2),]) %>% 
        mutate(x = paste(Width, Type.Conservation.2, sep = " ")) %>%
        mutate(x = factor(x, levels = paste(rep(widths, each = length(levels(TAB_2$Type.Conservation.2))), rep(levels(TAB_2$Type.Conservation.2), length(widths), sep = " ")))) %>%
        ggplot(aes(x = x, y = Mean.Families.Percentage.With.Motifs)) +
        geom_point(aes(color = Type.Conservation.2, alpha = Type_mod), size = 3) +
        scale_color_manual(values = c("#6f99ad", "#efc000", "#D87758")) +
        scale_x_discrete(limits = rev) +
        scale_alpha_discrete(range = c(1, 0.5)) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width and Positional conservation level") + ylab("Percentage of syntenic families with at least one DNA motif") +
        coord_flip() +
        guides(color = guide_legend(nrow = 1)) +
        labs(color = "Conservation level", alpha = "Dataset") +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_percentage_2.png"), height = 8, width = 15, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_percentage_2.pdf"), height = 8, width = 15, dpi = 600)
    }
  }
}

rm(list = c("st", "no", "mo", "gg1", "gg2", "TAB_1", "TAB_2"))





################################################################################
################################################################################
## 5. NUMBER/PERCENTAGE OF SYNTENIC FAMILIES WITH PUTATIVE FUNCTIONAL MOTIFS

#### 5.1 Create the tables.

cat(paste0("\n\nNUMBER/PERCENTAGE OF SYNTENIC FAMILIES WITH PUTATIVE FUNCTIONAL MOTIFS --> TABLES\n\n"))

for (st in strictnesses) {
  for (no in nonmatches) {
    
    cat(paste0("\nStrictness: ", st))
    cat(paste0("\nNonmatch: ", no))
    
    ## Positional level conservation table.
    Poslev = read.table(paste0(WD_in_3, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
    Poslev = Poslev[Poslev$Strictness == st & Poslev$NonMatch == no,]
    
    ## Create a global table
    TAB_F1 = data.frame()
    TAB_F2 = data.frame()
    TAB_F3 = data.frame()
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
            
            ## Gomo.
            gomo = tryCatch(read.table(paste0(WD_in_1, "/GOMO/", st, "/", no, "/", co, "/", cl, "/", mo, "/", wi, "/GOMO-SUMMARY.tsv"), sep = "\t", header = T, quote = "\""), error=function(e) NULL)
            
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
            
            if (!is.null(gomo)){
              tab = merge(tab, gomo, by.x = c("Type", "Family", "Meme_Motif.Identifier"), by.y = c("Type", "Family", "Gomo_Motif.Identifier"), all = T)
            }
            if (is.null(gomo)){
              tab$"Gomo_GO.Term.Identifier" = NA
            }
            
            ## Convert to factors.
            tab$Type = factor(tab$Type, levels = c("REAL", paste0("SIMULATION_", 1:n_sim)))
            tab$Family = factor(tab$Family, levels = c(paste0("Fam", 1:length(unique(tab$Family)))))
            tab$Type.Conservation.1 = factor(tab$Type.Conservation.1, levels = c("Conserved"))
            tab$Type.Conservation.2 = factor(tab$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
            tab$Type.Conservation.3 = factor(tab$Type.Conservation.3, levels = 2:9)
            
            ## Get number of motifs.
            tab = tab[,c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3", "Meme_Motif.Identifier", "Gomo_GO.Term.Identifier")]
            tab[tab == ""] = NA
            tab = tab[!duplicated(tab),]
            tab$Meme_Motif.Identifier = ifelse(is.na(tab$Gomo_GO.Term.Identifier), NA, tab$Meme_Motif.Identifier)
            tab = tab[,c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3", "Meme_Motif.Identifier")]
            tab$"Count" = ifelse(!is.na(tab$Meme_Motif.Identifier), 1, 0)
            tab_red = tab %>% group_by(Confidence, Comparison, Mode, Width, Type, Family, Type.Conservation.1, Type.Conservation.2, Type.Conservation.3) %>% summarise(Motifs.Number = sum(Count))
            tab_red = as.data.frame(tab_red)
            tab_red$"Motifs.Number.Binary" = ifelse(tab_red$Motifs.Number > 0, 1, 0)
            
            ## Collapse the table 1.
            tab_red_1 = tab_red %>% 
              group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.1, .drop=FALSE) %>% 
              summarise(
                Families.Number.Total = n_distinct(Family), 
                Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
                Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
                Mean.Accum.Motifs = mean(Motifs.Number)
              )
            tab_red_1 = as.data.frame(tab_red_1)
            tab_red_1$"Type_mod" = ifelse(tab_red_1$Type == "REAL", "REAL", "SIMULATION")
            tab_red_1_red = tab_red_1 %>% 
              group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.1, Families.Number.Total, .drop=FALSE) %>% 
              summarise(
                Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
                Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
                Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
              )
            tab_red_1_red = as.data.frame(tab_red_1_red)
            
            ## Collapse the table 2.
            tab_red_2 = tab_red %>% 
              group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.2, .drop=FALSE) %>% 
              summarise(
                Families.Number.Total = n_distinct(Family), 
                Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
                Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
                Mean.Accum.Motifs = mean(Motifs.Number)
              )
            tab_red_2 = as.data.frame(tab_red_2)
            tab_red_2$"Type_mod" = ifelse(tab_red_2$Type == "REAL", "REAL", "SIMULATION")
            tab_red_2_red = tab_red_2 %>% 
              group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.2, Families.Number.Total, .drop=FALSE) %>% 
              summarise(
                Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
                Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
                Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
              )
            tab_red_2_red = as.data.frame(tab_red_2_red)
            
            ## Collapse the table 3.
            tab_red_3 = tab_red %>% 
              group_by(Confidence, Comparison, Mode, Width, Type, Type.Conservation.3, .drop=FALSE) %>% 
              summarise(
                Families.Number.Total = n_distinct(Family), 
                Families.Number.With.Motifs = sum(Motifs.Number.Binary), 
                Families.Percentage.With.Motifs = round(Families.Number.With.Motifs/Families.Number.Total * 100, 2),
                Mean.Accum.Motifs = mean(Motifs.Number)
              )
            tab_red_3 = as.data.frame(tab_red_3)
            tab_red_3$"Type_mod" = ifelse(tab_red_3$Type == "REAL", "REAL", "SIMULATION")
            tab_red_3_red = tab_red_3 %>% 
              group_by(Confidence, Comparison, Mode, Width, Type_mod, Type.Conservation.3, Families.Number.Total, .drop=FALSE) %>% 
              summarise(
                Mean.Families.Number.With.Motifs = round(mean(Families.Number.With.Motifs), 0), 
                Mean.Families.Percentage.With.Motifs = round(mean(Families.Percentage.With.Motifs), 2),
                Mean.Mean.Accum.Motifs = round(mean(Mean.Accum.Motifs), 2)
              )
            tab_red_3_red = as.data.frame(tab_red_3_red)
            
            ### FINAL JOIN.
            TAB_F1 = rbind(TAB_F1, tab_red_1_red)
            TAB_F2 = rbind(TAB_F2, tab_red_2_red)
            TAB_F3 = rbind(TAB_F3, tab_red_3_red)
          }
        }
      }
    }
    
    rm(list = c("meme", "gomo", "Families_REAL", "Poslev_mod", "Poslev", "tab", "tab_red", "tab_red_1", "tab_red_2", "tab_red_3", "tab_red_1_red", "tab_red_2_red", "tab_red_3_red", "co", "cl", "mo", "wi"))
    
    ## Convert to factors.
    TAB_F1$Confidence = factor(TAB_F1$Confidence, levels = confidences)
    TAB_F1$Comparison = factor(TAB_F1$Comparison, levels = classes)
    TAB_F1$Width = factor(TAB_F1$Width, levels = widths)
    TAB_F1$Type_mod = factor(TAB_F1$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_F2$Confidence = factor(TAB_F2$Confidence, levels = confidences)
    TAB_F2$Comparison = factor(TAB_F2$Comparison, levels = classes)
    TAB_F2$Width = factor(TAB_F2$Width, levels = widths)
    TAB_F2$Type_mod = factor(TAB_F2$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_F2$Type.Conservation.2 = factor(TAB_F2$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
    TAB_F3$Confidence = factor(TAB_F3$Confidence, levels = confidences)
    TAB_F3$Comparison = factor(TAB_F3$Comparison, levels = classes)
    TAB_F3$Width = factor(TAB_F3$Width, levels = widths)
    TAB_F3$Type_mod = factor(TAB_F3$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_F3$Type.Conservation.3 = factor(TAB_F3$Type.Conservation.3, levels = 2:9)
    
    ## Change NA by 0.
    TAB_F1[is.na(TAB_F1)] = 0
    TAB_F2[is.na(TAB_F2)] = 0
    TAB_F3[is.na(TAB_F3)] = 0
    
    ## Save the tables.
    write.table(TAB_F1, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_1_FUNCTIONAL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_F2, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_2_FUNCTIONAL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_F3, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_3_FUNCTIONAL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
}

rm(list = c("st", "no", "TAB_F1", "TAB_F2", "TAB_F3"))


#### 5.2 Create the figures.

cat(paste0("\n\nNUMBER/PERCENTAGE OF SYNTENIC FAMILIES WITH PUTATIVE FUNCTIONAL MOTIFS --> FIGURES\n\n"))

for (st in strictnesses) {
  cat(paste0("-", st, "\n"))
  for (no in nonmatches) {
    cat(paste0("\t-", no, "\n"))
    for (mo in modes) {
      cat(paste0("\t\t-", mo, "\n"))
      
      ## Load tables
      TAB_F1 = read.table(paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_1_FUNCTIONAL.tsv"), sep = "\t", header = T, quote = "\"", check.names = F)
      TAB_F1 = TAB_F1[TAB_F1$Mode == mo,]
      TAB_F2 = read.table(paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_COLLAPSED_2_FUNCTIONAL.tsv"), sep = "\t", header = T, quote = "\"", check.names = F)
      TAB_F2 = TAB_F2[TAB_F2$Mode == mo,]
      
      ## Convert to factors.
      TAB_F1$Confidence = factor(TAB_F1$Confidence, levels = confidences)
      TAB_F1$Comparison = factor(TAB_F1$Comparison, levels = classes)
      TAB_F1$Width = factor(TAB_F1$Width, levels = widths)
      TAB_F1$Type_mod = factor(TAB_F1$Type_mod, levels = c("REAL", "SIMULATION"))
      TAB_F2$Confidence = factor(TAB_F2$Confidence, levels = confidences)
      TAB_F2$Comparison = factor(TAB_F2$Comparison, levels = classes)
      TAB_F2$Width = factor(TAB_F2$Width, levels = widths)
      TAB_F2$Type_mod = factor(TAB_F2$Type_mod, levels = c("REAL", "SIMULATION"))
      TAB_F2$Type.Conservation.2 = factor(TAB_F2$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
      
      ## Number (TAB_F1)
      gg1 = ggplot(TAB_F1, aes(x = Width, y = Mean.Families.Number.With.Motifs)) +
        geom_point(aes(color = Type_mod), size=4, alpha=0.6) +
        scale_color_manual(values = c("#0072b5", "#bc3c29")) +
        scale_x_discrete(limits = rev) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width") + ylab("Number of syntenic families with at least one functional DNA motif") +
        coord_flip() +
        guides(color = guide_legend(nrow = 1, title = "Dataset")) +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_number_1_FUNCTIONAL.png"), height = 6, width = 12, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_number_1_FUNCTIONAL.pdf"), height = 6, width = 12, dpi = 600)
      
      ## Number (TAB_F2)
      gg2 = with(TAB_F2, TAB_F2[order(Width, Type.Conservation.2),]) %>% 
        mutate(x = paste(Width, Type.Conservation.2, sep = " ")) %>%
        mutate(x = factor(x, levels = paste(rep(widths, each = length(levels(TAB_F2$Type.Conservation.2))), rep(levels(TAB_F2$Type.Conservation.2), length(widths), sep = " ")))) %>%
        ggplot(aes(x = x, y = Mean.Families.Number.With.Motifs)) +
        geom_point(aes(color = Type.Conservation.2, alpha = Type_mod), size = 3) +
        scale_color_manual(values = c("#6f99ad", "#efc000", "#D87758")) +
        scale_x_discrete(limits = rev) +
        scale_alpha_discrete(range = c(1, 0.5)) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width and Positional conservation level") + ylab("Number of syntenic families with at least one functional DNA motif") +
        coord_flip() +
        guides(color = guide_legend(nrow = 1)) +
        labs(color = "Conservation level", alpha = "Dataset") +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_number_2_FUNCTIONAL.png"), height = 8, width = 15, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_number_2_FUNCTIONAL.pdf"), height = 8, width = 15, dpi = 600)
      
      ## Percentage (TAB_F1)
      gg1 = ggplot(TAB_F1, aes(x = Width, y = Mean.Families.Percentage.With.Motifs)) +
        geom_point(aes(color = Type_mod), size=4, alpha=0.6) +
        scale_color_manual(values = c("#0072b5", "#bc3c29")) +
        scale_x_discrete(limits = rev) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width") + ylab("Percentage of syntenic families with at least one functional DNA motif") +
        coord_flip() +
        guides(color = guide_legend(nrow = 1, title = "Dataset")) +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_percentage_1_FUNCTIONAL.png"), height = 6, width = 12, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_percentage_1_FUNCTIONAL.pdf"), height = 6, width = 12, dpi = 600)
      
      ## Percentage (TAB_F2)
      gg2 = with(TAB_F2, TAB_F2[order(Width, Type.Conservation.2),]) %>% 
        mutate(x = paste(Width, Type.Conservation.2, sep = " ")) %>%
        mutate(x = factor(x, levels = paste(rep(widths, each = length(levels(TAB_F2$Type.Conservation.2))), rep(levels(TAB_F2$Type.Conservation.2), length(widths), sep = " ")))) %>%
        ggplot(aes(x = x, y = Mean.Families.Percentage.With.Motifs)) +
        geom_point(aes(color = Type.Conservation.2, alpha = Type_mod), size = 3) +
        scale_color_manual(values = c("#6f99ad", "#efc000", "#D87758")) +
        scale_x_discrete(limits = rev) +
        scale_alpha_discrete(range = c(1, 0.5)) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width and Positional conservation level") + ylab("Percentage of syntenic families with at least one functional DNA motif") +
        coord_flip() +
        guides(color = guide_legend(nrow = 1)) +
        labs(color = "Conservation level", alpha = "Dataset") +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_percentage_2_FUNCTIONAL.png"), height = 8, width = 15, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_percentage_2_FUNCTIONAL.pdf"), height = 8, width = 15, dpi = 600)
    }
  }
}

rm(list = c("st", "no", "mo", "gg1", "gg2", "TAB_F1", "TAB_F2"))





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
    Poslev = read.table(paste0(WD_in_3, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
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
    
    ## Convert to factors.
    TAB_GC$Confidence = factor(TAB_GC$Confidence, levels = confidences)
    TAB_GC$Comparison = factor(TAB_GC$Comparison, levels = classes)
    TAB_GC$Width = factor(TAB_GC$Width, levels = widths)
    TAB_GC$Type_mod = factor(TAB_GC$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_GC$Type.Conservation.2 = factor(TAB_GC$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
    TAB_GC$Type.Conservation.3 = factor(TAB_GC$Type.Conservation.3, levels = 2:9)
    TAB_Length$Confidence = factor(TAB_Length$Confidence, levels = confidences)
    TAB_Length$Comparison = factor(TAB_Length$Comparison, levels = classes)
    TAB_Length$Width = factor(TAB_Length$Width, levels = widths)
    TAB_Length$Type_mod = factor(TAB_Length$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_Length$Type.Conservation.2 = factor(TAB_Length$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
    TAB_Length$Type.Conservation.3 = factor(TAB_Length$Type.Conservation.3, levels = 2:9)
    TAB_E_value$Confidence = factor(TAB_E_value$Confidence, levels = confidences)
    TAB_E_value$Comparison = factor(TAB_E_value$Comparison, levels = classes)
    TAB_E_value$Width = factor(TAB_E_value$Width, levels = widths)
    TAB_E_value$Type_mod = factor(TAB_E_value$Type_mod, levels = c("REAL", "SIMULATION"))
    TAB_E_value$Type.Conservation.2 = factor(TAB_E_value$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
    TAB_E_value$Type.Conservation.3 = factor(TAB_E_value$Type.Conservation.3, levels = 2:9)
    
    ## Save the tables.
    write.table(TAB_GC, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_GC.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_Length, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_LENGTH.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(TAB_E_value, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_E_VALUE.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
}

rm(list = c("st", "no", "TAB_GC", "TAB_Length", "TAB_E_value"))


#### 6.2 Create the figures.

cat(paste0("\n\nGC CONTENT, LENGTH AND E-VALUE MOTIFS ANALYSIS --> FIGURES\n\n"))

for (st in strictnesses) {
  cat(paste0("-", st, "\n"))
  for (no in nonmatches) {
    cat(paste0("\t-", no, "\n"))
    for (mo in modes) {
      cat(paste0("\t\t-", mo, "\n"))
      
      ## Load tables
      TAB_GC = read.table(paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_GC.tsv"), sep = "\t", header = T, quote = "\"", check.names = F)
      TAB_GC = TAB_GC[TAB_GC$Mode == mo,]
      TAB_Length = read.table(paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_LENGTH.tsv"), sep = "\t", header = T, quote = "\"", check.names = F)
      TAB_Length = TAB_Length[TAB_Length$Mode == mo,]
      TAB_E_value = read.table(paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_E_VALUE.tsv"), sep = "\t", header = T, quote = "\"", check.names = F)
      TAB_E_value = TAB_E_value[TAB_E_value$Mode == mo,]
      
      ## Convert to factors.
      TAB_GC$Confidence = factor(TAB_GC$Confidence, levels = confidences)
      TAB_GC$Comparison = factor(TAB_GC$Comparison, levels = classes)
      TAB_GC$Width = factor(TAB_GC$Width, levels = widths)
      TAB_GC$Type_mod = factor(TAB_GC$Type_mod, levels = c("REAL", "SIMULATION"))
      TAB_GC$Type.Conservation.2 = factor(TAB_GC$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
      TAB_GC$Type.Conservation.3 = factor(TAB_GC$Type.Conservation.3, levels = 2:9)
      TAB_Length$Confidence = factor(TAB_Length$Confidence, levels = confidences)
      TAB_Length$Comparison = factor(TAB_Length$Comparison, levels = classes)
      TAB_Length$Width = factor(TAB_Length$Width, levels = widths)
      TAB_Length$Type_mod = factor(TAB_Length$Type_mod, levels = c("REAL", "SIMULATION"))
      TAB_Length$Type.Conservation.2 = factor(TAB_Length$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
      TAB_Length$Type.Conservation.3 = factor(TAB_Length$Type.Conservation.3, levels = 2:9)
      TAB_E_value$Confidence = factor(TAB_E_value$Confidence, levels = confidences)
      TAB_E_value$Comparison = factor(TAB_E_value$Comparison, levels = classes)
      TAB_E_value$Width = factor(TAB_E_value$Width, levels = widths)
      TAB_E_value$Type_mod = factor(TAB_E_value$Type_mod, levels = c("REAL", "SIMULATION"))
      TAB_E_value$Type.Conservation.2 = factor(TAB_E_value$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
      TAB_E_value$Type.Conservation.3 = factor(TAB_E_value$Type.Conservation.3, levels = 2:9)
      
      ## GC (1)
      cat(paste0("\t\t\t-GC (1)\n"))
      ggA1 = ggplot(TAB_GC, aes(x = Mean.Meme_Motif.GC, y = Width)) +
        geom_density_ridges(aes(fill = Type_mod), alpha=.4) +
        scale_fill_manual(values = c("#0072b5", "#bc3c29")) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("GC content (%)") + ylab("Density") +
        guides(color = guide_legend(nrow = 1, title = "Dataset")) +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_GC_1.png"), height = 6, width = 12, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_GC_1.pdf"), height = 6, width = 12, dpi = 600)
      
      ## GC (2)
      cat(paste0("\t\t\t-GC (2)\n"))
      ggA2 = with(TAB_GC, TAB_GC[order(Width, Type.Conservation.2),]) %>% 
        mutate(x = paste(Width, Type.Conservation.2, sep = " ")) %>%
        mutate(x = factor(x, levels = paste(rep(widths, each = length(levels(TAB_GC$Type.Conservation.2))), rep(levels(TAB_GC$Type.Conservation.2), length(widths), sep = " ")))) %>% 
        ggplot(aes(x = Mean.Meme_Motif.GC, y = x)) +
        geom_density_ridges(aes(fill = Type_mod), alpha=.4) +
        scale_fill_manual(values = c("#0072b5", "#bc3c29")) +
        facet_grid(Confidence ~ Comparison, scales = "free") +
        theme_bw() +
        xlab("Motif Width and Positional conservation level") + ylab("Density") +
        guides(color = guide_legend(nrow = 1)) +
        labs(fill = "Dataset") +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_GC_2.png"), height = 6, width = 12, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_GC_2.pdf"), height = 6, width = 12, dpi = 600)
      
      ## Length (1)
      cat(paste0("\t\t\t-Length (1)\n"))
      ggB1 = ggplot(TAB_Length, aes(x = Width, y = Meme_Motif.Width, fill = Type_mod)) + 
        geom_boxplot(aes(fill = Type_mod), size = 0.2, colour = "black", outlier.size = 0.25, position = position_dodge2(width = 0.9, preserve = "single")) + 
        stat_summary(fun = mean, geom = "point", aes(group = Type_mod), shape = 20, size = 1.2, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
        scale_fill_manual(values = c("#0072b5", "#bc3c29")) +
        scale_x_discrete(limits = rev) +
        facet_grid(Confidence ~ Comparison) +
        theme_bw() +
        xlab("Motif Width") + ylab("Length") +
        coord_flip() +
        guides(fill = guide_legend(nrow = 1, title = "Dataset")) +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_Length_1.png"), height = 6, width = 12, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_Length_1.pdf"), height = 6, width = 12, dpi = 600)
      
      ## Length (2)
      cat(paste0("\t\t\t-Length (2)\n"))
      ggB2 = with(TAB_Length, TAB_Length[order(Width, Type.Conservation.2),]) %>% 
        mutate(x = paste(Width, Type.Conservation.2, sep = " ")) %>%
        mutate(x = factor(x, levels = paste(rep(widths, each = length(levels(TAB_Length$Type.Conservation.2))), rep(levels(TAB_Length$Type.Conservation.2), length(widths), sep = " ")))) %>% 
        ggplot(aes(x = x, y = Meme_Motif.Width, fill = Type.Conservation.2, alpha = Type_mod)) + 
        geom_boxplot(aes(fill = Type.Conservation.2, alpha = Type_mod), size = 0.2, colour = "black", outlier.size = 0.25, position = position_dodge2(width = 0.9, preserve = "single")) + 
        stat_summary(fun = mean, geom = "point", aes(fill = Type.Conservation.2, alpha = Type_mod), shape = 20, size = 1.2, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
        scale_fill_manual(values = c("#6f99ad", "#efc000", "#D87758")) +
        scale_alpha_discrete(range = c(1, 0.5)) +
        scale_x_discrete(limits = rev) +
        facet_grid(Confidence ~ Comparison) +
        theme_bw() +
        xlab("Motif Width and Positional conservation level") + ylab("Length") +
        coord_flip() +
        guides(fill = guide_legend(nrow = 1)) +
        labs(fill = "Conservation level", alpha = "Dataset") +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_Length_2.png"), height = 8, width = 15, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_Length_2.pdf"), height = 8, width = 15, dpi = 600)
      
      ## E-value (1)
      cat(paste0("\t\t\t-E_value (1)\n"))
      ggC1 = ggplot(TAB_E_value, aes(x = Width, y = `-Log10.Meme_Motif.E.value`, fill = Type_mod)) + 
        geom_boxplot(aes(fill = Type_mod), size = 0.2, colour = "black", outlier.size = 0.25, position = position_dodge2(width = 0.9, preserve = "single")) + 
        stat_summary(fun = mean, geom = "point", aes(group = Type_mod), shape = 20, size = 1.2, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
        scale_fill_manual(values = c("#0072b5", "#bc3c29")) +
        scale_x_discrete(limits = rev) +
        scale_y_continuous(limits = c(0, 80)) +
        facet_grid(Confidence ~ Comparison) +
        theme_bw() +
        xlab("Motif Width") + ylab("-Log10(E-value)") +
        coord_flip() +
        guides(fill = guide_legend(nrow = 1, title = "Dataset")) +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_E_value_1.png"), height = 6, width = 12, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_E_value_1.pdf"), height = 6, width = 12, dpi = 600)
      
      ## E-value (2)
      cat(paste0("\t\t\t-E_value (2)\n"))
      ggC2 = with(TAB_E_value, TAB_E_value[order(Width, Type.Conservation.2),]) %>% 
        mutate(x = paste(Width, Type.Conservation.2, sep = " ")) %>%
        mutate(x = factor(x, levels = paste(rep(widths, each = length(levels(TAB_E_value$Type.Conservation.2))), rep(levels(TAB_E_value$Type.Conservation.2), length(widths), sep = " ")))) %>% 
        ggplot(aes(x = x, y = `-Log10.Meme_Motif.E.value`, fill = Type.Conservation.2, alpha = Type_mod)) + 
        geom_boxplot(aes(fill = Type.Conservation.2, alpha = Type_mod), size = 0.2, colour = "black", outlier.size = 0.25, position = position_dodge2(width = 0.9, preserve = "single")) + 
        stat_summary(fun = mean, geom = "point", aes(fill = Type.Conservation.2, alpha = Type_mod), shape = 20, size = 1.2, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
        scale_fill_manual(values = c("#6f99ad", "#efc000", "#D87758")) +
        scale_alpha_discrete(range = c(1, 0.5)) +
        scale_x_discrete(limits = rev) +
        scale_y_continuous(limits = c(0, 80)) +
        facet_grid(Confidence ~ Comparison) +
        theme_bw() +
        xlab("Motif Width and Positional conservation level") + ylab("-Log10(E-value)") +
        coord_flip() +
        guides(fill = guide_legend(nrow = 1)) +
        labs(fill = "Conservation level", alpha = "Dataset") +
        theme(panel.background = element_blank()) +
        theme(legend.position = "top")
      
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_E_value_2.png"), height = 8, width = 15, dpi = 600)
      ggsave(paste0(WD_out, "/Figures/", st, "/", no, "/FIGURE_E_value_2.pdf"), height = 8, width = 15, dpi = 600)
    }
  }
}

rm(list = c("st", "no", "mo", "ggA1", "ggA2", "ggB1", "ggB2", "ggC1", "ggC2", "TAB_GC", "TAB_Length", "TAB_E_value"))





################################################################################
################################################################################
## 7. FUNCTIONAL MOTIFS

cat(paste0("\n\nFUNCTIONAL MOTIFS --> TABLES\n\n"))

#### 7.1 Create the tables.

for (st in strictnesses) {
  for (no in nonmatches) {
    
    cat(paste0("\nStrictness: ", st))
    cat(paste0("\nNonmatch: ", no))
    
    ## Positional level conservation table.
    Poslev = read.table(paste0(WD_in_3, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", header = T, quote = "\"")
    Poslev = Poslev[Poslev$Strictness == st & Poslev$NonMatch == no,]
    
    ## Create a global table
    TAB_F = data.frame()
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
            
            ## Gomo.
            gomo = tryCatch(read.table(paste0(WD_in_1, "/GOMO/", st, "/", no, "/", co, "/", cl, "/", mo, "/", wi, "/GOMO-SUMMARY.tsv"), sep = "\t", header = T, quote = "\""), error=function(e) NULL)
            
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
            
            if (!is.null(gomo)){
              tab = merge(tab, gomo, by.x = c("Type", "Family", "Meme_Motif.Identifier"), by.y = c("Type", "Family", "Gomo_Motif.Identifier"), all = T)
            }
            if (is.null(gomo)){
              tab$"Gomo_GO.Term.Identifier" = NA
              tab$"Gomo_P.value" = NA
              tab$"Gomo_Q.value" = NA
            }
            
            ## Convert to factors.
            tab$Type = factor(tab$Type, levels = c("REAL", paste0("SIMULATION_", 1:n_sim)))
            tab$Family = factor(tab$Family, levels = c(paste0("Fam", 1:length(unique(tab$Family)))))
            tab$Type.Conservation.1 = factor(tab$Type.Conservation.1, levels = c("Conserved"))
            tab$Type.Conservation.2 = factor(tab$Type.Conservation.2, levels = c("Low-conserved", "Medium-conserved", "High-conserved"))
            tab$Type.Conservation.3 = factor(tab$Type.Conservation.3, levels = 2:9)
            
            ## Get number of motifs.
            tab = tab[,c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3", "Meme_Motif.Identifier", "Meme_Motif.Width", "Meme_Motif.E.value", "Gomo_GO.Term.Identifier", "Gomo_P.value", "Gomo_Q.value")]
            tab = tab[tab$Type == "REAL",]
            tab[tab == ""] = NA
            tab = tab[!is.na(tab$Gomo_GO.Term.Identifier), ]
            tab = tab[!duplicated(tab),]
            tab = tab[,c("Confidence", "Comparison", "Mode", "Width", "Type", "Family", "Type.Conservation.1", "Type.Conservation.2", "Type.Conservation.3", "Meme_Motif.Identifier", "Meme_Motif.Width", "Meme_Motif.E.value", "Gomo_GO.Term.Identifier", "Gomo_P.value", "Gomo_Q.value")]
            tab$"Gomo_GO.Term.Identifier.Long" = goIdToTerm(tab$Gomo_GO.Term.Identifier, names = FALSE)
            
            ### FINAL JOIN.
            TAB_F = rbind(TAB_F, tab)
          }
        }
      }
    }
    
    rm(list = c("meme", "gomo", "Families_REAL", "Poslev_mod", "Poslev", "tab", "co", "cl", "mo", "wi"))
    
    ## Convert to factors.
    TAB_F$Confidence = factor(TAB_F$Confidence, levels = confidences)
    TAB_F$Comparison = factor(TAB_F$Comparison, levels = classes)
    TAB_F$Width = factor(TAB_F$Width, levels = widths)
    
    ## Save the tables.
    write.table(TAB_F, paste0(WD_out, "/Tables/", st, "/", no, "/TABLE_FUNCTIONAL_MOTIFS.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
}

rm(list = c("st", "no", "TAB_F"))

