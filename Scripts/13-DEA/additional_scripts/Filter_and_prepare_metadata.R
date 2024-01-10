


################################################################################
################################################################################
################################### MODULES ####################################
################################################################################
################################################################################


suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(DESeq2))
suppressMessages(library(plotly))
suppressMessages(library(htmltools))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tibble))
# devtools::install_github("azhu513/apeglm")
suppressMessages(library(scales))
suppressMessages(library(ggExtra))
options(bitmapType='cairo')










################################################################################
################################################################################
################################## VARIABLES ###################################
################################################################################
################################################################################


## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied.", call.=FALSE)
} else {
  path_meta = args[1]
  path_DEA = args[2]
  path_quant = args[3]
}

# path_meta = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Additional_info/sra-info/metadata/Search_studies/csa/Studies"
# path_DEA = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/16-DEA/csa"
# path_quant = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/06-quantification/csa/Salmon/ALL/nr/03-Quant"










################################################################################
################################################################################
####################### METADATA: EXPLORATORY ANALYSIS #########################
################################################################################
################################################################################


cat(paste0("METADATA: EXPLORATORY ANALYSIS..."))

if (!dir.exists(paste0(path_DEA, "/01-Metadata_EA"))){
  dir.create(paste0(path_DEA, "/01-Metadata_EA"))
}

# SRA studies.
SRA_studies = gsub(".tsv", "", list.files(path_meta, pattern = "[0-9].tsv"))

############################################
## 1. SRA STUDIES
############################################

for (SRA_study in SRA_studies) {
  
  cat(paste0("\n\nSRA study: ", SRA_study))
  
  # Load the metadata and metadata_info tables for a particular SRA study.
  META = read.table(paste0(path_meta, "/", SRA_study, ".tsv"), sep = "\t", header = T, quote = "\"")
  META_info = read.table(paste0(path_meta, "/", SRA_study, "_info.tsv"), sep = "\t", header = T, quote = "\"")
  colnames(META_info) = paste0(colnames(META_info), "_info")
  # Merge both tables to join all the information.
  META = merge(META, 
               META_info, 
               by.x = c("Specie", "SRA_study", "Sample", "Comparison", "Class", "Condition", "Replicate"), 
               by.y = c("Specie_info", "SRA_study_info", "Sample_info", "Comparison_info", "Class_info", "Condition_info", "Replicate_info"),
               all = T)
  
  rm(list = c("META_info"))
  
  
  
  
  ############################################
  ## 2. COMPARISONS
  ############################################
  
  # Subdivide the main table (META) into separate tables using the "Comparison" 
  # column ("Comparison 1", "Comparison 2"...). Most of the tables have only one 
  # "Comparison" value ("Comparison 1").
 
  for (comparison in unique(META$Comparison)) {
     
    cat(paste0("\n\t-Comparison: ", comparison))
    
    i = 1
    n = gsub("Comparison ", "", comparison)
    
    # Get the table. For example, "Comparison 1".
    meta = META[META$Comparison == comparison,]
    # According to the "Class" column, this table can analyze "Development", 
    # "Abiotic stress" or "Biotic stress". Each comparison has only one class. 
    # If there is more than one class in META table, it is because there is more 
    # than one comparison in the table.
    class = unique(meta$Class)
    
    cat(paste0("\n\t\tClass: ", class))
    
    
    
    
    ############################################
    ## 3.1 DEVELOPMENT
    ############################################
    
    if (class == "Development") {

      ## First, we subdivide the table meta into separate tables by "Dev_stage",
      ## "Tissue", "Age", "Technique", "Cultivar" or "Genotype", if we can do it.
      ## We can do this if some of these columns have more than one value.
      last = meta[, c("Tissue", "Age", "Technique", "Cultivar", "Genotype")]
      items = colnames(last)[sapply(last, function(col) length(unique(col)) > 1)]
      tabs = meta %>% group_by_at(items) %>% group_split()

      cat(paste0("\n\t\tColumns with variability: ", ifelse(length(items) == 0, "None", paste(items, collapse = ", "))))
      cat(paste0("\n\t\t", ifelse(length(tabs) == 1, "1 subdivision", paste0(length(tabs), " subdivisions"))))
      
      ## Second, we iterate each possible and independent table, and save it.
      for (tab in tabs) {

        cat(paste0("\n\t\t\tCreate table ", i, "..."))

        tab = as.data.frame(tab)
        
        # Check if there are information on the "Dev_stage" column.
        if ("DS.0" %in% tab$Dev_stage) {
          cat(paste0("\n\t\t\tNO (Info: There is no information on the \"Dev_stage\" column. Give the information. Check it)."))
        }
        else {
          tab = tab[order(tab$Condition, tab$Dev_stage_info, tab$Replicate),]
          tab$"Name.1" = paste0(tab$Condition, "__", gsub("-", "_", gsub(" ", "_", tab$Dev_stage_info)), "__", tab$Treatment_time_exp, "__", tab$Exp_level, "__", tab$Stress, "__", tab$Replicate)
          tab$"Name.2" = paste0(tab$Condition, "__", gsub("-", "_", gsub(" ", "_", tab$Dev_stage_info)), "__", tab$Treatment_time_exp, "__", tab$Exp_level, "__", tab$Stress)
          
          # Check if all samples are quantified.
          cat(paste0("\n\t\t\tFilter 1..."))
          
          tab_mod_1 = data.frame()
          for (a in 1:nrow(tab)) {
            file = tab[a, "Sample"]
            if (file.exists(paste0(path_quant, "/", file))) {
              tab_mod_1 = rbind(tab_mod_1, tab[a,])
            }else {
              cat(paste0("\n\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
            }
          }

          # Check if there are more than one sample by condition.
          cat(paste0("\n\t\t\tFilter 2..."))
          
          tab_mod_2 = data.frame()
          for (name in unique(tab_mod_1$Name.2)) {
            subset = tab_mod_1[tab_mod_1$Name.2 == name,]
            if (dim(subset)[1] > 1) {
              tab_mod_2 = rbind(tab_mod_2, subset)
            }
          }
          
          rm(list = c("name", "subset", "a", "file"))
          
          if (("Control" %in% tab_mod_2$Condition) & ("Stage" %in% tab_mod_2$Condition)) {
            cat(paste0("\n\t\t\tYES"))
            write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
            i = i + 1
          }
          else {
            cat(paste0("\n\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
          }
        }
        
      } # End loop tabs (Development)

      rm(list = c("i", "n", "last", "items", "tabs", "tab", "tab_mod_1", "tab_mod_2"))
      
    } # End conditional DEVELOPMENT
    
    
    
    
    
    ############################################
    ## 3.2 ABIOTIC AND BIOTIC STRESS
    ############################################
    
    if (class == "Biotic stress" | class == "Abiotic stress") {
      
      
      
      
      ############################################
      ## 4.1 ONE STRESS
      ############################################
      
      if (length(unique(meta$Stress)) == 1) {
        
        ## First, we subdivide the table meta into separate tables by "Dev_stage", 
        ## "Tissue", "Age", "Technique", "Cultivar" or "Genotype", if we can do it. 
        ## We can do this if some of these columns have more than one value.
        last = meta[, c("Tissue", "Age", "Technique", "Cultivar", "Genotype")]
        items = colnames(last)[sapply(last, function(col) length(unique(col)) > 1)]
        tabs = meta %>% group_by_at(items) %>% group_split()
        
        cat(paste0("\n\t\tColumns with variability: ", ifelse(length(items) == 0, "None", paste(items, collapse = ", "))))
        cat(paste0("\n\t\t", ifelse(length(tabs) == 1, "1 subdivision", paste0(length(tabs), " subdivisions"))))
        
        ## Second, we iterate each possible and independent table, and save it.
        for (tab in tabs) {
          
          tab = as.data.frame(tab)
          
          
          ############################################
          ## 5.1 TREATMENT TIME EXPOSITION VALUES
          ############################################
          
          if (!("TTE.0" %in% tab$Treatment_time_exp) & ("EL.0" %in% tab$Exp_level)) {
            
            tab = tab[order(tab$Condition, tab$Treatment_time_exp, tab$Replicate),]
            tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", gsub("-", "_", gsub(" ", "_", tab$Treatment_time_exp_info)), "__", tab$Exp_level, "__", tab$Stress, "__", tab$Replicate)
            tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", gsub("-", "_", gsub(" ", "_", tab$Treatment_time_exp_info)), "__", tab$Exp_level, "__", tab$Stress)
            
            
            ############################################
            ## 6.1 ONE CONTROL
            ############################################
            
            if (length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1) {
              
              cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])), " subsubdivisions"))))
              
              cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
              
              # Check if all samples are quantified.
              cat(paste0("\n\t\t\t\tFilter 1..."))
              
              tab_mod_1 = data.frame()
              for (a in 1:nrow(tab)) {
                file = tab[a, "Sample"]
                if (file.exists(paste0(path_quant, "/", file))) {
                  tab_mod_1 = rbind(tab_mod_1, tab[a,])
                }else {
                  cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                }
              }
              
              # Check if there are more than one sample by condition.
              cat(paste0("\n\t\t\t\tFilter 2..."))
              
              tab_mod_2 = data.frame()
              for (name in unique(tab_mod_1$Name.2)) {
                subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                if (dim(subset)[1] > 1) {
                  tab_mod_2 = rbind(tab_mod_2, subset)
                }
              }
              
              rm(list = c("name", "subset", "a", "file"))
              
              if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                cat(paste0("\n\t\t\t\tYES"))
                write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                i = i + 1
              }
              else {
                cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
              }
            }
            
            
            ############################################
            ## 6.2 SEVERAL CONTROLS
            ############################################
            
            if (length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) > 1) {
              
              cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])), " subsubdivisions"))))
              
              # Check if all samples are quantified.
              cat(paste0("\n\t\t\t\tFilter 1..."))
              
              tab_mod_1 = data.frame()
              for (a in 1:nrow(tab)) {
                file = tab[a, "Sample"]
                if (file.exists(paste0(path_quant, "/", file))) {
                  tab_mod_1 = rbind(tab_mod_1, tab[a,])
                }else {
                  cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                }
              }
              
              rm(list = c("a", "file"))
              
              for (TTE in unique(tab_mod_1$Treatment_time_exp)) {
                
                cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                
                # Check if there are more than one sample by condition.
                cat(paste0("\n\t\t\t\tFilter 2..."))
                
                tab_mod_2 = data.frame()
                for (name in unique(tab_mod_1[tab_mod_1$Treatment_time_exp == TTE, "Name.2"])) {
                  subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                  if (dim(subset)[1] > 1) {
                    tab_mod_2 = rbind(tab_mod_2, subset)
                  }
                }
                
                rm(list = c("name", "subset"))
                
                if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                  cat(paste0("\n\t\t\t\tYES"))
                  write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                  i = i + 1
                }
                else {
                  cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                }
              }
              
              rm(list = c("TTE"))
              
            }
            
            rm(list = c("tab_mod_1", "tab_mod_2"))
            
          } # End conditional TREATMENT TIME EXPOSITION VALUES
          
          
          ############################################
          ## 5.2 EXPOSITION LEVEL VALUES
          ############################################
          
          else if (("TTE.0" %in% tab$Treatment_time_exp) & !("EL.0" %in% tab$Exp_level)) {
            
            tab = tab[order(tab$Condition, tab$Exp_level, tab$Replicate),]
            tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", gsub("-", "_", gsub(" ", "_", tab$Exp_level_info)), "__", tab$Stress, "__", tab$Replicate)
            tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", gsub("-", "_", gsub(" ", "_", tab$Exp_level_info)), "__", tab$Stress)
            
            
            ############################################
            ## 6.1 ONE CONTROL
            ############################################
            
            if (length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1) {
              
              cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Exp_level"])), " subsubdivisions"))))
              
              cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
              
              # Check if all samples are quantified.
              cat(paste0("\n\t\t\t\tFilter 1..."))
              
              tab_mod_1 = data.frame()
              for (a in 1:nrow(tab)) {
                file = tab[a, "Sample"]
                if (file.exists(paste0(path_quant, "/", file))) {
                  tab_mod_1 = rbind(tab_mod_1, tab[a,])
                }else {
                  cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                }
              }
              
              # Check if there are more than one sample by condition.
              cat(paste0("\n\t\t\t\tFilter 2..."))
              
              tab_mod_2 = data.frame()
              for (name in unique(tab_mod_1$Name.2)) {
                subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                if (dim(subset)[1] > 1) {
                  tab_mod_2 = rbind(tab_mod_2, subset)
                }
              }
              
              rm(list = c("name", "subset", "a", "file"))
              
              if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                cat(paste0("\n\t\t\t\tYES"))
                write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                i = i + 1
              }
              else {
                cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
              }
            }
            
            
            ############################################
            ## 6.2 SEVERAL CONTROLS
            ############################################
            
            if (length(unique(tab[tab$Condition == "Control", "Exp_level"])) > 1) {
              
              cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Exp_level"])), " subsubdivisions"))))
              
              # Check if all samples are quantified.
              cat(paste0("\n\t\t\t\tFilter 1..."))
              
              tab_mod_1 = data.frame()
              for (a in 1:nrow(tab)) {
                file = tab[a, "Sample"]
                if (file.exists(paste0(path_quant, "/", file))) {
                  tab_mod_1 = rbind(tab_mod_1, tab[a,])
                }else {
                  cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                }
              }
              
              rm(list = c("a", "file"))
              
              for (EL in unique(tab_mod_1$Exp_level)) {
                
                cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                
                # Check if there are more than one sample by condition.
                cat(paste0("\n\t\t\t\tFilter 2..."))
                
                tab_mod_2 = data.frame()
                for (name in unique(tab_mod_1[tab_mod_1$Exp_level == EL, "Name.2"])) {
                  subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                  if (dim(subset)[1] > 1) {
                    tab_mod_2 = rbind(tab_mod_2, subset)
                  }
                }
                
                rm(list = c("name", "subset"))
                
                if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                  cat(paste0("\n\t\t\t\tYES"))
                  write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                  i = i + 1
                }
                else {
                  cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                }
              }
              
              rm(list = c("EL"))
              
            }
            
            rm(list = c("tab_mod_1", "tab_mod_2"))
            
          } # end conditional EXPOSITION LEVEL VALUES
          
          
          ############################################
          ## 5.3 TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
          ############################################
          
          else if (!("TTE.0" %in% tab$Treatment_time_exp) & !("EL.0" %in% tab$Exp_level)) {
            
            cat(paste0("\n\t\t\tCreate table ", i, "..."))
            
            cat(paste0("\n\t\t\tNO (Info: This experiment contains Treatment time exposition and Exposition level values. Therefore, it is too complex)."))
            
          } # End conditional TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
          
          
          ############################################
          ## 5.4 NO TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES 
          ############################################
          
          else {
            
            tab = tab[order(tab$Condition, tab$Exp_level, tab$Replicate),]
            tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", tab$Exp_level, "__", tab$Stress, "__", tab$Replicate)
            tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", tab$Exp_level, "__", tab$Stress)
            
            cat(paste0("\n\t\t\tCreate table ", i, "..."))
            
            # Check if all samples are quantified.
            cat(paste0("\n\t\t\tFilter 1..."))
            
            tab_mod_1 = data.frame()
            for (a in 1:nrow(tab)) {
              file = tab[a, "Sample"]
              if (file.exists(paste0(path_quant, "/", file))) {
                tab_mod_1 = rbind(tab_mod_1, tab[a,])
              }else {
                cat(paste0("\n\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
              }
            }
            
            # Check if there are more than one sample by condition.
            cat(paste0("\n\t\t\tFilter 2..."))
            
            tab_mod_2 = data.frame()
            for (name in unique(tab_mod_1$Name.2)) {
              subset = tab_mod_1[tab_mod_1$Name.2 == name,]
              if (dim(subset)[1] > 1) {
                tab_mod_2 = rbind(tab_mod_2, subset)
              }
            }
            
            rm(list = c("name", "subset", "a", "file"))
            
            if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
              cat(paste0("\n\t\t\tYES"))
              write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
              i = i + 1
            }
            else {
              cat(paste0("\n\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
            }
            
            rm(list = c("tab_mod_1", "tab_mod_2"))
            
          } # End conditional NO TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
          
        } # End loop tabs (Abiotic or biotic stress; One stress)
        
      } # End conditional ONE STRESS
      
      
      
      
      ############################################
      ## 4.2 TWO OR MORE STRESSES
      ############################################
      
      if (length(unique(meta$Stress)) > 1) {
        
        
        ############################################
        ## 5.1 DEPENDENT CONTROL
        ############################################
        
        if ("TRUE" %in% grepl(":", unique(meta$Stress), fixed = T)) {
          
          ## First, we subdivide the table meta into separate tables by "Dev_stage", 
          ## "Tissue", "Age", "Technique", "Cultivar" or "Genotype", if we can do it. 
          ## We can do this if some of these columns have more than one value.
          last = meta[, c("Tissue", "Age", "Technique", "Cultivar", "Genotype")]
          items = colnames(last)[sapply(last, function(col) length(unique(col)) > 1)]
          tabs = meta %>% group_by_at(items) %>% group_split()
          
          cat(paste0("\n\t\tColumns with variability: ", ifelse(length(items) == 0, "None", paste(items, collapse = ", "))))
          cat(paste0("\n\t\t", ifelse(length(tabs) == 1, "1 subdivision", paste0(length(tabs), " subdivisions"))))
          
          for (tab in tabs) {
            
            tab = as.data.frame(tab)
            
            
            ############################################
            ## 6.1 TREATMENT TIME EXPOSITION VALUES
            ############################################
            
            if (!("TTE.0" %in% tab$Treatment_time_exp) & ("EL.0" %in% tab$Exp_level)) {
              
              tab = tab[order(tab$Condition, tab$Treatment_time_exp, tab$Replicate),]
              tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", gsub("-", "_", gsub(" ", "_", tab$Treatment_time_exp_info)), "__", tab$Exp_level, "__", tab$Stress, "__", tab$Replicate)
              tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", gsub("-", "_", gsub(" ", "_", tab$Treatment_time_exp_info)), "__", tab$Exp_level, "__", tab$Stress)
              
              ############################################
              ## 7.1 ONE CONTROL
              ############################################
              
              if (length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1) {
                
                cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])), " subsubdivisions"))))
                
                cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                
                # Check if all samples are quantified.
                cat(paste0("\n\t\t\t\tFilter 1..."))
                
                tab_mod_1 = data.frame()
                for (a in 1:nrow(tab)) {
                  file = tab[a, "Sample"]
                  if (file.exists(paste0(path_quant, "/", file))) {
                    tab_mod_1 = rbind(tab_mod_1, tab[a,])
                  }else {
                    cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                  }
                }
                
                # Check if there are more than one sample by condition.
                cat(paste0("\n\t\t\t\tFilter 2..."))
                
                tab_mod_2 = data.frame()
                for (name in unique(tab_mod_1$Name.2)) {
                  subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                  if (dim(subset)[1] > 1) {
                    tab_mod_2 = rbind(tab_mod_2, subset)
                  }
                }
                
                rm(list = c("name", "subset", "a", "file"))
                
                if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                  cat(paste0("\n\t\t\t\tYES"))
                  write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                  i = i + 1
                }
                else {
                  cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                }
              }
              
              ############################################
              ## 7.2 SEVERAL CONTROLS
              ############################################
              
              if (length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) > 1) {
                
                cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])), " subsubdivisions"))))
                
                # Check if all samples are quantified.
                cat(paste0("\n\t\t\t\tFilter 1..."))
                
                tab_mod_1 = data.frame()
                for (a in 1:nrow(tab)) {
                  file = tab[a, "Sample"]
                  if (file.exists(paste0(path_quant, "/", file))) {
                    tab_mod_1 = rbind(tab_mod_1, tab[a,])
                  }else {
                    cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                  }
                }
                
                rm(list = c("a", "file"))
                
                for (TTE in unique(tab_mod_1$Treatment_time_exp)) {
                  
                  cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                  
                  # Check if there are more than one sample by condition.
                  cat(paste0("\n\t\t\t\tFilter 2..."))
                  
                  tab_mod_2 = data.frame()
                  for (name in unique(tab_mod_1[tab_mod_1$Treatment_time_exp == TTE, "Name.2"])) {
                    subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                    if (dim(subset)[1] > 1) {
                      tab_mod_2 = rbind(tab_mod_2, subset)
                    }
                  }
                  
                  rm(list = c("name", "subset"))
                  
                  if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                    cat(paste0("\n\t\t\t\tYES"))
                    write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                    i = i + 1
                  }
                  else {
                    cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                  }
                }
                
                rm(list = c("TTE"))
                
              }
              
              rm(list = c("tab_mod_1", "tab_mod_2"))
              
            } # End conditional TREATMENT TIME EXPOSITION VALUES
            
            
            ############################################
            ## 6.2 EXPOSITION LEVEL VALUES
            ############################################
            
            else if (("TTE.0" %in% tab$Treatment_time_exp) & !("EL.0" %in% tab$Exp_level)) {
              
              tab = tab[order(tab$Condition, tab$Exp_level, tab$Replicate),]
              tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", gsub("-", "_", gsub(" ", "_", tab$Exp_level_info)), "__", tab$Stress, "__", tab$Replicate)
              tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", gsub("-", "_", gsub(" ", "_", tab$Exp_level_info)), "__", tab$Stress)
              
              ############################################
              ## 7.1 ONE CONTROL
              ############################################
              
              if (length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1) {
                
                cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Exp_level"])), " subsubdivisions"))))
                
                cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                
                # Check if all samples are quantified.
                cat(paste0("\n\t\t\t\tFilter 1..."))
                
                tab_mod_1 = data.frame()
                for (a in 1:nrow(tab)) {
                  file = tab[a, "Sample"]
                  if (file.exists(paste0(path_quant, "/", file))) {
                    tab_mod_1 = rbind(tab_mod_1, tab[a,])
                  }else {
                    cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                  }
                }
                
                # Check if there are more than one sample by condition.
                cat(paste0("\n\t\t\t\tFilter 2..."))
                
                tab_mod_2 = data.frame()
                for (name in unique(tab_mod_1$Name.2)) {
                  subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                  if (dim(subset)[1] > 1) {
                    tab_mod_2 = rbind(tab_mod_2, subset)
                  }
                }
                
                rm(list = c("name", "subset", "a", "file"))
                
                if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                  cat(paste0("\n\t\t\t\tYES"))
                  write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                  i = i + 1
                }
                else {
                  cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                }
              }
              
              ############################################
              ## 7.2 SEVERAL CONTROLS
              ############################################
              
              if (length(unique(tab[tab$Condition == "Control", "Exp_level"])) > 1) {
                
                cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Exp_level"])), " subsubdivisions"))))
                
                # Check if all samples are quantified.
                cat(paste0("\n\t\t\t\tFilter 1..."))
                
                tab_mod_1 = data.frame()
                for (a in 1:nrow(tab)) {
                  file = tab[a, "Sample"]
                  if (file.exists(paste0(path_quant, "/", file))) {
                    tab_mod_1 = rbind(tab_mod_1, tab[a,])
                  }else {
                    cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                  }
                }
                
                rm(list = c("a", "file"))
                
                for (EL in unique(tab_mod_1$Exp_level)) {
                  
                  cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                  
                  # Check if there are more than one sample by condition.
                  cat(paste0("\n\t\t\t\tFilter 2..."))
                  
                  tab_mod_2 = data.frame()
                  for (name in unique(tab_mod_1[tab_mod_1$Exp_level == EL, "Name.2"])) {
                    subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                    if (dim(subset)[1] > 1) {
                      tab_mod_2 = rbind(tab_mod_2, subset)
                    }
                  }
                  
                  rm(list = c("name", "subset"))
                  
                  if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                    cat(paste0("\n\t\t\t\tYES"))
                    write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                    i = i + 1
                  }
                  else {
                    cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                  }
                }
                
                rm(list = c("EL"))
                
              }
              
              rm(list = c("tab_mod_1", "tab_mod_2"))
              
            } # End conditional EXPOSITION LEVEL VALUES
            
            
            ############################################
            ## 6.3 TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
            ############################################
            
            else if (!("TTE.0" %in% tab$Treatment_time_exp) & !("EL.0" %in% tab$Exp_level)) {
              
              cat(paste0("\n\t\t\tCreate table ", i, "..."))
              
              cat(paste0("\n\t\t\tNO (Info: This experiment contains Treatment time exposition and Exposition level values. Therefore, it is too complex)."))
              
            } # End conditional TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
            
            
            ############################################
            ## 6.4 NO TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES 
            ############################################
            
            else {
              
              tab = tab[order(tab$Condition, tab$Exp_level, tab$Replicate),]
              tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", tab$Exp_level, "__", tab$Stress, "__", tab$Replicate)
              tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", tab$Exp_level, "__", tab$Stress)
              
              cat(paste0("\n\t\t\tCreate table ", i, "..."))
              
              # Check if all samples are quantified.
              cat(paste0("\n\t\t\tFilter 1..."))
              
              tab_mod_1 = data.frame()
              for (a in 1:nrow(tab)) {
                file = tab[a, "Sample"]
                if (file.exists(paste0(path_quant, "/", file))) {
                  tab_mod_1 = rbind(tab_mod_1, tab[a,])
                }else {
                  cat(paste0("\n\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                }
              }
              
              # Check if there are more than one sample by condition.
              cat(paste0("\n\t\t\tFilter 2..."))
              
              tab_mod_2 = data.frame()
              for (name in unique(tab_mod_1$Name.2)) {
                subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                if (dim(subset)[1] > 1) {
                  tab_mod_2 = rbind(tab_mod_2, subset)
                }
              }
              
              rm(list = c("name", "subset", "a", "file"))
              
              if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                cat(paste0("\n\t\t\tYES"))
                write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                i = i + 1
              }
              else {
                cat(paste0("\n\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
              }
              
              rm(list = c("tab_mod_1", "tab_mod_2"))
              
            } # End conditional NO TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
            
          } # End loop tabs (Abiotic or biotic stress; Two or more stresses; Dependent control)
          
        } # End conditional DEPENDENT CONTROL
        
        
        ############################################
        ## 5.2 INDEPENDENT CONTROL
        ############################################
        
        # Don't contain two points
        else if (!("TRUE" %in% grepl(":", unique(meta$Stress), fixed = T))) {
          
          ## First, we subdivide the table meta into separate tables by "Dev_stage", 
          ## "Tissue", "Age", "Technique", "Cultivar" or "Genotype", if we can do it. 
          ## We can do this if some of these columns have more than one value.
          last = meta[, c("Stress", "Tissue", "Age", "Technique", "Cultivar", "Genotype")]
          items = colnames(last)[sapply(last, function(col) length(unique(col)) > 1)]
          tabs = meta %>% group_by_at(items) %>% group_split()
          
          cat(paste0("\n\t\tColumns with variability: ", ifelse(length(items) == 0, "None", paste(items, collapse = ", "))))
          cat(paste0("\n\t\t", ifelse(length(tabs) == 1, "1 subdivision", paste0(length(tabs), " subdivisions"))))
          
          for (tab in tabs) {
            
            tab = as.data.frame(tab)
            
            
            ############################################
            ## 6.1 TREATMENT TIME EXPOSITION VALUES
            ############################################
            
            if (!("TTE.0" %in% tab$Treatment_time_exp) & ("EL.0" %in% tab$Exp_level)) {
              
              tab = tab[order(tab$Condition, tab$Treatment_time_exp, tab$Replicate),]
              tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", gsub("-", "_", gsub(" ", "_", tab$Treatment_time_exp_info)), "__", tab$Exp_level, "__", tab$Stress, "__", tab$Replicate)
              tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", gsub("-", "_", gsub(" ", "_", tab$Treatment_time_exp_info)), "__", tab$Exp_level, "__", tab$Stress)
              
              ############################################
              ## 7.1 ONE CONTROL
              ############################################
              
              if (length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1) {
                
                cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])), " subsubdivisions"))))
                
                cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                
                # Check if all samples are quantified.
                cat(paste0("\n\t\t\t\tFilter 1..."))
                
                tab_mod_1 = data.frame()
                for (a in 1:nrow(tab)) {
                  file = tab[a, "Sample"]
                  if (file.exists(paste0(path_quant, "/", file))) {
                    tab_mod_1 = rbind(tab_mod_1, tab[a,])
                  }else {
                    cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                  }
                }
                
                # Check if there are more than one sample by condition.
                cat(paste0("\n\t\t\t\tFilter 2..."))
                
                tab_mod_2 = data.frame()
                for (name in unique(tab_mod_1$Name.2)) {
                  subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                  if (dim(subset)[1] > 1) {
                    tab_mod_2 = rbind(tab_mod_2, subset)
                  }
                }
                
                rm(list = c("name", "subset", "a", "file"))
                
                if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                  cat(paste0("\n\t\t\t\tYES"))
                  write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                  i = i + 1
                }
                else {
                  cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                }
              }
              
              ############################################
              ## 7.2 SEVERAL CONTROLS
              ############################################
              
              if (length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) > 1) {
                
                cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Treatment_time_exp"])), " subsubdivisions"))))
                
                # Check if all samples are quantified.
                cat(paste0("\n\t\t\t\tFilter 1..."))
                
                tab_mod_1 = data.frame()
                for (a in 1:nrow(tab)) {
                  file = tab[a, "Sample"]
                  if (file.exists(paste0(path_quant, "/", file))) {
                    tab_mod_1 = rbind(tab_mod_1, tab[a,])
                  }else {
                    cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                  }
                }
                
                rm(list = c("a", "file"))
                
                for (TTE in unique(tab_mod_1$Treatment_time_exp)) {
                  
                  cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                  
                  # Check if there are more than one sample by condition.
                  cat(paste0("\n\t\t\t\tFilter 2..."))
                  
                  tab_mod_2 = data.frame()
                  for (name in unique(tab_mod_1[tab_mod_1$Treatment_time_exp == TTE, "Name.2"])) {
                    subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                    if (dim(subset)[1] > 1) {
                      tab_mod_2 = rbind(tab_mod_2, subset)
                    }
                  }
                  
                  rm(list = c("name", "subset"))
                  
                  if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                    cat(paste0("\n\t\t\t\tYES"))
                    write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                    i = i + 1
                  }
                  else {
                    cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                  }
                }
                
                rm(list = c("TTE"))
                
              }
              
              rm(list = c("tab_mod_1", "tab_mod_2"))
              
            } # End conditional TREATMENT TIME EXPOSITION VALUES
            
            
            ############################################
            ## 6.2 EXPOSITION LEVEL VALUES
            ############################################
            
            else if (("TTE.0" %in% tab$Treatment_time_exp) & !("EL.0" %in% tab$Exp_level)) {
              
              tab = tab[order(tab$Condition, tab$Exp_level, tab$Replicate),]
              tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", gsub("-", "_", gsub(" ", "_", tab$Exp_level_info)), "__", tab$Stress, "__", tab$Replicate)
              tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", gsub("-", "_", gsub(" ", "_", tab$Exp_level_info)), "__", tab$Stress)
              
              ############################################
              ## 7.1 ONE CONTROL
              ############################################
              
              if (length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1) {
                
                cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Exp_level"])), " subsubdivisions"))))
                
                cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                
                # Check if all samples are quantified.
                cat(paste0("\n\t\t\t\tFilter 1..."))
                
                tab_mod_1 = data.frame()
                for (a in 1:nrow(tab)) {
                  file = tab[a, "Sample"]
                  if (file.exists(paste0(path_quant, "/", file))) {
                    tab_mod_1 = rbind(tab_mod_1, tab[a,])
                  }else {
                    cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                  }
                }
                
                # Check if there are more than one sample by condition.
                cat(paste0("\n\t\t\t\tFilter 2..."))
                
                tab_mod_2 = data.frame()
                for (name in unique(tab_mod_1$Name.2)) {
                  subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                  if (dim(subset)[1] > 1) {
                    tab_mod_2 = rbind(tab_mod_2, subset)
                  }
                }
                
                rm(list = c("name", "subset", "a", "file"))
                
                if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                  cat(paste0("\n\t\t\t\tYES"))
                  write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                  i = i + 1
                }
                else {
                  cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                }
              }
              
              ############################################
              ## 7.2 SEVERAL CONTROLS
              ############################################
              
              if (length(unique(tab[tab$Condition == "Control", "Exp_level"])) > 1) {
                
                cat(paste0("\n\t\t\t", ifelse(length(unique(tab[tab$Condition == "Control", "Exp_level"])) == 1, "1 subsubdivision", paste0(length(unique(tab[tab$Condition == "Control", "Exp_level"])), " subsubdivisions"))))
                
                # Check if all samples are quantified.
                cat(paste0("\n\t\t\t\tFilter 1..."))
                
                tab_mod_1 = data.frame()
                for (a in 1:nrow(tab)) {
                  file = tab[a, "Sample"]
                  if (file.exists(paste0(path_quant, "/", file))) {
                    tab_mod_1 = rbind(tab_mod_1, tab[a,])
                  }else {
                    cat(paste0("\n\t\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                  }
                }
                
                rm(list = c("a", "file"))
                
                for (EL in unique(tab_mod_1$Exp_level)) {
                  
                  cat(paste0("\n\t\t\t\tCreate table ", i, "..."))
                  
                  # Check if there are more than one sample by condition.
                  cat(paste0("\n\t\t\t\tFilter 2..."))
                  
                  tab_mod_2 = data.frame()
                  for (name in unique(tab_mod_1[tab_mod_1$Exp_level == EL, "Name.2"])) {
                    subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                    if (dim(subset)[1] > 1) {
                      tab_mod_2 = rbind(tab_mod_2, subset)
                    }
                  }
                  
                  rm(list = c("name", "subset"))
                  
                  if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                    cat(paste0("\n\t\t\t\tYES"))
                    write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                    i = i + 1
                  }
                  else {
                    cat(paste0("\n\t\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
                  }
                }
                
                rm(list = c("EL"))
                
              }
              
              rm(list = c("tab_mod_1", "tab_mod_2"))
              
            } # End conditional EXPOSITION LEVEL VALUES
            
            
            ############################################
            ## 6.3 TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
            ############################################
            
            else if (!("TTE.0" %in% tab$Treatment_time_exp) & !("EL.0" %in% tab$Exp_level)) {
              
              cat(paste0("\n\t\t\tCreate table ", i, "..."))
              
              cat(paste0("\n\t\t\tNO (Info: This experiment contains Treatment time exposition and Exposition level values. Therefore, it is too complex)."))
              
            } # End conditional TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
            
            
            ############################################
            ## 6.4 NO TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES 
            ############################################
            
            else {
              
              tab = tab[order(tab$Condition, tab$Exp_level, tab$Replicate),]
              tab$"Name.1" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", tab$Exp_level, "__", tab$Stress, "__", tab$Replicate)
              tab$"Name.2" = paste0(tab$Condition, "__", tab$Dev_stage, "__", tab$Treatment_time_exp, "__", tab$Exp_level, "__", tab$Stress)
              
              cat(paste0("\n\t\t\tCreate table ", i, "..."))
              
              # Check if all samples are quantified.
              cat(paste0("\n\t\t\tFilter 1..."))
              
              tab_mod_1 = data.frame()
              for (a in 1:nrow(tab)) {
                file = tab[a, "Sample"]
                if (file.exists(paste0(path_quant, "/", file))) {
                  tab_mod_1 = rbind(tab_mod_1, tab[a,])
                }else {
                  cat(paste0("\n\t\t\t", file, " folder doesn't exist. This sample probably did not pass the library depth or strand specificity filter."))
                }
              }
              
              # Check if there are more than one sample by condition.
              cat(paste0("\n\t\t\tFilter 2..."))
              
              tab_mod_2 = data.frame()
              for (name in unique(tab_mod_1$Name.2)) {
                subset = tab_mod_1[tab_mod_1$Name.2 == name,]
                if (dim(subset)[1] > 1) {
                  tab_mod_2 = rbind(tab_mod_2, subset)
                }
              }
              
              rm(list = c("name", "subset", "a", "file"))
              
              if (("Control" %in% tab_mod_2$Condition) & ("Treated" %in% tab_mod_2$Condition)) {
                cat(paste0("\n\t\t\tYES"))
                write.table(tab_mod_2, paste0(path_DEA, "/01-Metadata_EA/", SRA_study, "_", n, "_", i, ".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
                i = i + 1
              }
              else {
                cat(paste0("\n\t\t\tNO (Info: There are some conditions with less than two replicates. Check it)."))
              }
              
              rm(list = c("tab_mod_1", "tab_mod_2"))
              
            } # End conditional NO TREATMENT TIME EXPOSITION AND EXPOSITION LEVEL VALUES
            
          } # End loop tabs (Abiotic or biotic stress; Two or more stresses; Independent control)
          
        } # End conditional INDEPENDENT CONTROL
        
      } # End conditional TWO OR MORE STRESSES
      
      rm(list = c("i", "n", "last", "items", "tabs", "tab"))
      
    } # End conditional ABIOTIC AND BIOTIC STRESS
    
  rm(list = c("class", "meta"))
  
  } # End loop COMPARISONS
  
  rm(list = c("META", "comparison"))
  
} # End loop SRA STUDIES

rm(list = c("SRA_study", "SRA_studies"))

cat("\n\n\n--------------------------------------------------------------------")

