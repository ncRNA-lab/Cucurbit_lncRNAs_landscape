################################################################################
#
# FIGURES
#
################################################################################

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library("UpSetR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggbreak"))

options(bitmapType='cairo')
options(stringsAsFactors = F)

## 1. VARIABLES

flag = "nr"
WD1 = paste0("/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sequence_level/Blastn/", flag, "/Prueba_filters_identity_percentage_and_alignment_length/05-Families")
WD2 = paste0("/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/08-comparative_genomics/Sequence_level/Blastn/", flag, "/Prueba_filters_identity_percentage_and_alignment_length/06-Figures")
WD3 = "/storage/ncRNA/Projects/lncRNAs/Cucurbitaceae/Results/06-quantification"
species = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
classes = c("ALL", "intergenic", "antisense", "intronic", "sense")
confidences = c("Low", "Medium", "High")

if (!dir.exists(WD2)){
  dir.create(WD2)
}

# In OrthoMCL, it's important to check if the file exists. If there are some 
# blast table without hits OrthoMCL fails and doesn't generate the all_orthomcl.out 
# table. This table is necessary to get the families with the script get_families_from_OrthoMCL.py.
# So it's important to use file.exists() with Famtab and Gentab

# In OrthoFinder, this error doesn't appear because it always generates the Orthogroups.txt
# table and then the script get_families_from_OrthoFinder.py doesn't fail. However, if there are 
# few hits we can get some error with upsetR plots. So we use tryCatch to capture the error.

################################################################################
## 2. UPSETR PLOT

# Esta figura no se puede utilizar en el paper pero se puede utilizar a modo de consulta.
# El grafico de la izquierda muestra para cada especie el numero de lncRNAs o familias que
# se estan representando. Por otro lado, el grafico de arriba muestra cual es el numero de 
# familias o lncRNAs dados por la interseccion de especies mostrada bajo. Por ejemplo, si
# se trata de una familia con lncRNAs de cme y cmo la familia sera sumada a la columna de
# la interseccion cme-cmo. En el caso de lncRNAs, si se trata de un lncRNA de cme que se
# encuentra en una familia que alberga lncRNAs de cme y cmo pues aparecerá en la columna de
# la interaccion cme-cmo. Cada familia y lncRNA pueden sumar solo una vez, es decir, no hay
# posibilidad de que esten representados varias veces.

# Hablamos en numeros absolutos.

# Si el fichero presenta _ALL es que contiene tambien las familias no conservadas. Es decir,
# las que solo albergan una especie y por tanto un unico lncRNA ya que el blastn solo es entre
# especies. La unica manera de que hayan dos lncRNAs de la misma especie en la misma familia
# es que hay al menos un lncRNA de otra especie en la misma familia y que los conecte entre
# ellos, pero en este caso estariamos hablando de una familia conservada.

# Se representa una figura para cada caso confidence level y class.

cat("\nPlotting upsetr plots...\n")

for (confidence in confidences) {
  for (class in classes) {
    
    if (!dir.exists(paste0(WD2, "/", confidence))){
      dir.create(paste0(WD2, "/", confidence))
    }
    if (!dir.exists(paste0(WD2, "/", confidence, "/", class))){
      dir.create(paste0(WD2, "/", confidence, "/", class))
    }
    
    ## FAMILIES
    Fampath = paste0(WD1, "/", confidence, "/", class, "/fam.tsv")
    if (file.exists(Fampath)) {
      Famtab_all = read.table(Fampath, header = T, sep = "\t", quote = "\"")
      pdf(file = paste0(WD2, "/", confidence, "/", class, "/Upset_Families_all.pdf"), height = 10, width = 140)
      tryCatch (
        print(
          upset(Famtab_all, 
                sets = species, 
                nintersects = NA,
                sets.bar.color = c(rep("#86c3e7", length(species))), 
                order.by = "freq", 
                #empty.intersections = "on", 
                text.scale = 1.7, 
                mainbar.y.label = "Number of families (Intersection Size)", 
                sets.x.label = "Number of families (Set Size)")
        ),
        error = function(e){message(paste0("ERROR: Families, ", confidence, "-", class))}
      )
      dev.off()
      
      # When we remove the UNCONSERVED families (families composed of a single lncRNA or several lncRNAs of a single specie 
      # if we have allowed paralogs. The last one are considered unconserved families because we are studying conservation
      # across the species), it's important to sort the columns (Species). UpsetR give an error if the first specie which 
      # appears in the table has 0 in all the rows. I mean, if the specie doesn't participate in any CONSERVED family, 
      # UpsetR will give us an error.
      Famtab = Famtab_all[rowSums(Famtab_all[,2:dim(Famtab_all)[2]]) > 1, ]
      check = as.data.frame(colSums(Famtab[,2:dim(Famtab)[2]]))
      check$"Specie" = rownames(check)
      colnames(check) = c("Size", "Specie")
      check = check[order(check$Size, decreasing = TRUE),]
      New_order = check$Specie
      Famtab = Famtab[, c("Family", New_order)]
      pdf(file = paste0(WD2, "/", confidence, "/", class, "/Upset_Families.pdf"), height = 10, width = 140)
      tryCatch (
        print(
          upset(Famtab, 
                sets = species, 
                nintersects = NA,
                sets.bar.color = c(rep("#86c3e7", length(species))), 
                order.by = "freq", 
                #empty.intersections = "on", 
                text.scale = 1.7, 
                mainbar.y.label = "Number of families (Intersection Size)", 
                sets.x.label = "Number of families (Set Size)")
        ),
        error = function(e){message(paste0("ERROR: Families, ", confidence, "-", class))}
      )
      dev.off()
    }
    else {
      cat(paste0("ERROR: File Fampath doesn't exist for families, ", confidence, "-", class, "\n"))
    }
    
    ## LNCRNAS
    Genpath = paste0(WD1, "/", confidence, "/", class, "/gen.tsv")
    if (file.exists(Genpath)) {
      Gentab_all = read.table(Genpath, header = T, sep = "\t", quote = "\"")
      pdf(file = paste0(WD2, "/", confidence, "/", class, "/Upset_LncRNAs_all.pdf"), height = 10, width = 140)
      tryCatch (
        print(
          upset(Gentab_all, 
                sets = species, 
                nintersects = NA,
                sets.bar.color = c(rep("#86c3e7", length(species))), 
                order.by = "freq", 
                #empty.intersections = "on", 
                text.scale = 1.7, 
                mainbar.y.label = "Number of LncRNAs (Intersection Size)", 
                sets.x.label = "Number of LncRNAs (Set Size)")
        ),
        error = function(e){message(paste0("ERROR: LncRNAs, ", confidence, "-", class))}
      )
      dev.off()
      
      # When we remove the UNCONSERVED lncRNAs (lncRNAs of families composed of a single lncRNA or several lncRNAs of a single 
      # specie if we have allowed paralogs. The last one are considered unconserved lncRNAs because we are studying conservation
      # across the species), it's important to sort the columns (Species). UpsetR give an error if the first specie which appears 
      # in the table has 0 in all the rows. I mean, if the specie doesn't participate in any CONSERVED family, UpsetR will give us 
      # an error.
      Gentab = Gentab_all[rowSums(Gentab_all[,3:dim(Gentab_all)[2]]) > 1, ]
      check = as.data.frame(colSums(Gentab[,3:dim(Gentab)[2]]))
      check$"Specie" = rownames(check)
      colnames(check) = c("Size", "Specie")
      check = check[order(check$Size, decreasing = TRUE),]
      New_order = check$Specie
      Gentab = Gentab[, c("Family", "Member", New_order)]
      pdf(file = paste0(WD2, "/", confidence, "/", class, "/Upset_LncRNAs.pdf"), height = 10, width = 140)
      tryCatch (
        print(
          upset(Gentab, 
                sets = species, 
                nintersects = NA,
                sets.bar.color = c(rep("#86c3e7", length(species))), 
                order.by = "freq", 
                #empty.intersections = "on", 
                text.scale = 1.7, 
                mainbar.y.label = "Number of LncRNAs (Intersection Size)", 
                sets.x.label = "Number of LncRNAs (Set Size)")
        ),
        error = function(e){message(paste0("ERROR: LncRNAs, ", confidence, "-", class))}
      )
      dev.off()
    }
    else {
      cat(paste0("ERROR: File Genpath doesn't exist for lncRNAs, ", confidence, "-", class, "\n"))
    }
  }
}

################################################################################
## 3. LOLLIPOP PLOT

# Esta figura si podria servir mas para el paper. Si no le hubiesemos dado la vuelta
# el eje 'x' representa el numero de especies que podemos encontrar dentro de la familia
# representada o de la familia donde se encuentra el lncRNA representado y en el eje 'y'
# encontramos el numero de familias o lncRNAs cumpliendo esta regla.

# Hablamos en numeros absolutos.

# Si el fichero presenta _ALL es que contiene tambien las familias no conservadas. Es decir,
# las que solo albergan una especie y por tanto un unico lncRNA ya que el blastn solo es entre
# especies. La unica manera de que hayan dos lncRNAs de la misma especie en la misma familia
# es que hay al menos un lncRNA de otra especie en la misma familia y que los conecte entre
# ellos, pero en este caso estariamos hablando de una familia conservada.

# Se representa una figura para cada caso confidence level y class.

cat("\nPlotting lollipop plots...\n")

for (confidence in confidences) {
  for (class in classes) {
    
    if (!dir.exists(paste0(WD2, "/", confidence))){
      dir.create(paste0(WD2, "/", confidence))
    }
    if (!dir.exists(paste0(WD2, "/", confidence, "/", class))){
      dir.create(paste0(WD2, "/", confidence, "/", class))
    }
    
    ## FAMILIES
    Fampath = paste0(WD1, "/", confidence, "/", class, "/fam.tsv")
    if (file.exists(Fampath)) {
      Famtab_all = read.table(Fampath, header = T, sep = "\t", quote = "\"")
      Famtab_all$"Sum" = rowSums(Famtab_all[, 2:dim(Famtab_all)[2]])
      Famtab_all$"Counts" = 1
      
      Famtab_all_temp = Famtab_all
      Famtab_all_temp$Sum = factor(Famtab_all_temp$Sum, levels = c(1:length(species)))
      
      Famtab_all_red = Famtab_all_temp %>% 
        group_by(Sum, .drop=FALSE) %>%
        summarise(
          sum.Counts = sum(Counts))
      Famtab_all_red = as.data.frame(Famtab_all_red)
      
      ggplot(Famtab_all_red, aes(x = Sum, y=sum.Counts)) + 
        geom_segment(aes(x = Sum, xend=Sum, y=0, yend = sum.Counts), color = "skyblue", linewidth = 2) +
        geom_point(color="blue", size=4, alpha=0.6) +
        scale_x_discrete(expand = c(0.03, 0.03)) + 
        theme_bw() +
        xlab("Number of species") +
        ylab("Number of families") +
        #scale_y_break(breaks=c(1250, 3500), scales = 0.15, expand = c(0.01, 0.03)) +
        scale_y_continuous(limits = c(0, max(Famtab_all_red$sum.Counts) + max(Famtab_all_red$sum.Counts)*0.05), breaks = seq(0, max(Famtab_all_red$sum.Counts), 200), expand = c(0, 1)) +
        coord_flip() +
        geom_text(aes(label = sum.Counts), hjust = -1) +
        theme(
          panel.grid.minor.y = element_blank(),
          #panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90), 
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 18)
      )
      
      ggsave(paste0(WD2, "/", confidence, "/", class, "/Lollipop_Families_all.png"), height = 5, width = 20, dpi = 600)
      ggsave(paste0(WD2, "/", confidence, "/", class, "/Lollipop_Families_all.pdf"), height = 5, width = 20, dpi = 600)
      
      Famtab = Famtab_all[Famtab_all$Sum > 1,]
      Famtab$Sum = factor(Famtab$Sum, levels = c(2:length(species)))
      
      Famtab_red = Famtab %>% 
        group_by(Sum, .drop=FALSE) %>%
        summarise(
          sum.Counts = sum(Counts))
      Famtab_red = as.data.frame(Famtab_red)
      
      ggplot(Famtab_red, aes(x = Sum, y=sum.Counts)) + 
        geom_segment(aes(x = Sum, xend=Sum, y=0, yend = sum.Counts), color = "skyblue", linewidth = 2) +
        geom_point(color="blue", size=4, alpha=0.6) +
        scale_x_discrete(expand = c(0.03, 0.03)) + 
        theme_bw() +
        xlab("Number of species") +
        ylab("Number of families") +
        #scale_y_break(breaks=c(1250, 3500), scales = 0.15, expand = c(0.01, 0.03)) +
        scale_y_continuous(limits = c(0, max(Famtab_red$sum.Counts) + max(Famtab_red$sum.Counts)*0.05), breaks = seq(0, max(Famtab_red$sum.Counts), 50), expand = c(0, 1)) +
        coord_flip() +
        geom_text(aes(label = sum.Counts), hjust = -1) +
        theme(
          panel.grid.minor.y = element_blank(),
          #panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90), 
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 18)
        )
      
      ggsave(paste0(WD2, "/", confidence, "/", class, "/Lollipop_Families.png"), height = 5, width = 20, dpi = 600)
      ggsave(paste0(WD2, "/", confidence, "/", class, "/Lollipop_Families.pdf"), height = 5, width = 20, dpi = 600)
    } 
    else {
      cat(paste0("ERROR: File Fampath doesn't exist for families, ", confidence, "-", class, "\n"))
    }
    
    ## LNCRNAS
    Genpath = paste0(WD1, "/", confidence, "/", class, "/gen.tsv")
    if (file.exists(Genpath)) {
      Gentab_all = read.table(Genpath, header = T, sep = "\t", quote = "\"")
      Gentab_all$"Sum" = rowSums(Gentab_all[, 3:dim(Gentab_all)[2]])
      Gentab_all$"Counts" = 1
      
      Gentab_all_temp = Gentab_all
      Gentab_all_temp$Sum = factor(Gentab_all_temp$Sum, levels = c(1:length(species)))
      
      Gentab_all_red = Gentab_all_temp %>% 
        group_by(Sum, .drop=FALSE) %>%
        summarise(
          sum.Counts = sum(Counts))
      Gentab_all_red = as.data.frame(Gentab_all_red)
      
      ggplot(Gentab_all_red, aes(x = Sum, y=sum.Counts)) + 
        geom_segment(aes(x = Sum, xend=Sum, y=0, yend = sum.Counts), color = "skyblue", linewidth = 2) +
        geom_point(color="blue", size=4, alpha=0.6) +
        scale_x_discrete(expand = c(0.03, 0.03)) + 
        theme_bw() +
        xlab("Number of species") +
        ylab("Number of lncRNAs") +
        #scale_y_break(breaks=c(1250, 3500), scales = 0.15, expand = c(0.01, 0.03)) +
        scale_y_continuous(limits = c(0, max(Gentab_all_red$sum.Counts) + max(Gentab_all_red$sum.Counts)*0.05), breaks = seq(0, max(Gentab_all_red$sum.Counts), 500), expand = c(0, 1)) +
        coord_flip() +
        geom_text(aes(label = sum.Counts), hjust = -1) +
        theme(
          panel.grid.minor.y = element_blank(),
          #panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90), 
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 18)
        )
      
      ggsave(paste0(WD2, "/", confidence, "/", class, "/Lollipop_LncRNAs_all.png"), height = 5, width = 20, dpi = 600)
      ggsave(paste0(WD2, "/", confidence, "/", class, "/Lollipop_LncRNAs_all.pdf"), height = 5, width = 20, dpi = 600)
      
      Gentab = Gentab_all[Gentab_all$Sum > 1,]
      Gentab$Sum = factor(Gentab$Sum, levels = c(2:length(species)))
      
      Gentab_red = Gentab %>% 
        group_by(Sum, .drop=FALSE) %>%
        summarise(
          sum.Counts = sum(Counts))
      Gentab_red = as.data.frame(Gentab_red)
      
      ggplot(Gentab_red, aes(x = Sum, y=sum.Counts)) + 
        geom_segment(aes(x = Sum, xend=Sum, y=0, yend = sum.Counts), color = "skyblue", linewidth = 2) +
        geom_point(color="blue", size=4, alpha=0.6) +
        scale_x_discrete(expand = c(0.03, 0.03)) + 
        theme_bw() +
        xlab("Number of species") +
        ylab("Number of lncRNAs") +
        #scale_y_break(breaks=c(1250, 3500), scales = 0.15, expand = c(0.01, 0.03)) +
        scale_y_continuous(limits = c(0, max(Gentab_red$sum.Counts) + max(Gentab_red$sum.Counts)*0.05), breaks = seq(0, max(Gentab_red$sum.Counts), 50), expand = c(0, 1)) +
        coord_flip() +
        geom_text(aes(label = sum.Counts), hjust = -1) +
        theme(
          panel.grid.minor.y = element_blank(),
          #panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90), 
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 18)
        )
      
      ggsave(paste0(WD2, "/", confidence, "/", class, "/Lollipop_LncRNAs.png"), height = 5, width = 20, dpi = 600)
      ggsave(paste0(WD2, "/", confidence, "/", class, "/Lollipop_LncRNAs.pdf"), height = 5, width = 20, dpi = 600)
    } 
    else {
      cat(paste0("ERROR: File Genpath doesn't exist for lncRNAs, ", confidence, "-", class, "\n"))
    }
  }
}

################################################################################
## 4. FAMILIES: CONSERVATION PERCENTAGE TABLE

# Aqui mostramos tanto una tabla como dos figuras globales con el porcentaje de familias 
# conservadas (Number of species by family > 1) y no conservadas (Number of species by 
# family = 1).

# Se representan dos figuras teniendo en cuenta:
#     -Fig 1: confidence level, class, specie y Number_species_by_family.
#     -Fig 2: confidence level, class, specie y type.

cat("\nFamilies: Building the conservation percentage table...\n")

TAB_1_FAM = data.frame()
TAB_2_FAM = data.frame()
TAB_3_FAM = data.frame()
TAB_ALL_FAM = data.frame()
for (confidence in confidences) {
  for (class in classes) {
    if (file.exists(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"))) {
      gen = read.table(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"), sep = "\t", header = T, quote = "\"")
      gen$"Number_species_by_family" = rowSums(gen[,species])
      gen$"Type" = ifelse(gen$Number_species_by_family > 1, "Conserved", "Non-conserved")
      gen$"Conserved_level" = ifelse(gen$Number_species_by_family == 2 | gen$Number_species_by_family == 3, "Low-conserved",
                                     ifelse(gen$Number_species_by_family == 4 | gen$Number_species_by_family == 5 | gen$Number_species_by_family == 6, "Medium-conserved",
                                            ifelse(gen$Number_species_by_family == 7 | gen$Number_species_by_family == 8 | gen$Number_species_by_family == 9, "High-conserved", "Non-conserved")))
      gen$"Specie" = sapply(strsplit(gen$Member, "-"), "[[", 2)
      fam = gen[,c("Family", "Type", "Conserved_level", "Number_species_by_family", "Specie")]
      fam = fam[!duplicated(fam),]
      fam$"Confidence" = confidence
      fam$"Class" = class
      
      fam$Type = factor(fam$Type, levels = c("Non-conserved", "Conserved"))
      fam$Conserved_level = factor(fam$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      fam$Number_species_by_family = factor(fam$Number_species_by_family, levels = 1:length(species))
      fam$Specie = factor(fam$Specie, levels = species)
      
      TAB_ALL_FAM = rbind(TAB_ALL_FAM, fam)
      
      fam_red_1 = fam %>% 
        group_by(Confidence, Class, Specie, Number_species_by_family, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Family)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2))
      fam_red_1 = as.data.frame(fam_red_1)
      colnames(fam_red_1) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_Families", "Total_Families", "Percentage_Families")
      TAB_1_FAM = rbind(TAB_1_FAM, fam_red_1)
      
      fam_red_2 = fam %>% 
        group_by(Confidence, Class, Specie, Conserved_level, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Family)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2))
      fam_red_2 = as.data.frame(fam_red_2)
      colnames(fam_red_2) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts_Families", "Total_Families", "Percentage_Families")
      TAB_2_FAM = rbind(TAB_2_FAM, fam_red_2)
      
      fam_red_3 = fam %>% 
        group_by(Confidence, Class, Specie, Type, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Family)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2))
      fam_red_3 = as.data.frame(fam_red_3)
      colnames(fam_red_3) = c("Confidence", "Class", "Specie", "Type", "Counts_Families", "Total_Families", "Percentage_Families")
      TAB_3_FAM = rbind(TAB_3_FAM, fam_red_3)
    }
  }
}

write.table(TAB_1_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_1.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_2_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_2.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_3_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_3.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_ALL_FAM, paste0(WD2, "/TABLE_FAMILIES_PERCENTAGE_ALL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

TAB_1_FAM$Confidence = factor(TAB_1_FAM$Confidence, levels = confidences)
TAB_1_FAM$Class = factor(TAB_1_FAM$Class, levels = classes)

gg1 = ggplot(TAB_1_FAM, aes(x = Confidence, y = Percentage_Families, fill = Number_species_by_family)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_1.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_1.pdf"), height = 18, width = 16, dpi = 600)

TAB_2_FAM$Confidence = factor(TAB_2_FAM$Confidence, levels = confidences)
TAB_2_FAM$Class = factor(TAB_2_FAM$Class, levels = classes)

gg2 = ggplot(TAB_2_FAM, aes(x = Confidence, y = Percentage_Families, fill = Conserved_level)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_2.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_2.pdf"), height = 18, width = 16, dpi = 600)

TAB_3_FAM$Confidence = factor(TAB_3_FAM$Confidence, levels = confidences)
TAB_3_FAM$Class = factor(TAB_3_FAM$Class, levels = classes)

gg3 = ggplot(TAB_3_FAM, aes(x = Confidence, y = Percentage_Families, fill = Type)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_Families), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_Families_3.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_Families_3.pdf"), height = 18, width = 16, dpi = 600)

################################################################################
## 5. LNCRNAS: CONSERVATION PERCENTAGE TABLE

# Aqui mostramos tanto una tabla como dos figuras globales con el porcentaje de lncRNAs 
# conservados (Number of species by family > 1, de la familia donde se encuentra el
# lncRNA) y no conservados (Number of species by family = 1, de la familia donde se 
# encuentra el lncRNA).

# Se representan dos figuras teniendo en cuenta:
#     -Fig 1: confidence level, class, specie y Number_species_by_family.
#     -Fig 2: confidence level, class, specie y type.

cat("\nLncRNAs: Building the conservation percentage table...\n")

TAB_1_LNCRNAS = data.frame()
TAB_2_LNCRNAS = data.frame()
TAB_3_LNCRNAS = data.frame()
TAB_ALL_LNCRNAS = data.frame()
for (confidence in confidences) {
  for (class in classes) {
    if (file.exists(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"))) {
      gen = read.table(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"), sep = "\t", header = T, quote = "\"")
      gen$"Number_species_by_family" = rowSums(gen[,species])
      gen$"Type" = ifelse(gen$Number_species_by_family > 1, "Conserved", "Non-conserved")
      gen$"Conserved_level" = ifelse(gen$Number_species_by_family == 2 | gen$Number_species_by_family == 3, "Low-conserved",
                                     ifelse(gen$Number_species_by_family == 4 | gen$Number_species_by_family == 5 | gen$Number_species_by_family == 6, "Medium-conserved",
                                            ifelse(gen$Number_species_by_family == 7 | gen$Number_species_by_family == 8 | gen$Number_species_by_family == 9, "High-conserved", "Non-conserved")))
      gen$"Specie" = sapply(strsplit(gen$Member, "-"), "[[", 2)
      gen = gen[,c("Member", "Type", "Conserved_level", "Number_species_by_family", "Specie")]
      gen$"Confidence" = confidence
      gen$"Class" = class
      
      gen$Type = factor(gen$Type, levels = c("Non-conserved", "Conserved"))
      gen$Conserved_level = factor(gen$Conserved_level, levels = c("Non-conserved", "Low-conserved", "Medium-conserved", "High-conserved"))
      gen$Number_species_by_family = factor(gen$Number_species_by_family, levels = 1:length(species))
      gen$Specie = factor(gen$Specie, levels = species)
      
      TAB_ALL_LNCRNAS = rbind(TAB_ALL_LNCRNAS, gen)
      
      gen_red_1 = gen %>% 
        group_by(Confidence, Class, Specie, Number_species_by_family, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Member)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2))
      gen_red_1 = as.data.frame(gen_red_1)
      colnames(gen_red_1) = c("Confidence", "Class", "Specie", "Number_species_by_family", "Counts_LncRNAs", "Total_LncRNAs", "Percentage_LncRNAs")
      TAB_1_LNCRNAS = rbind(TAB_1_LNCRNAS, gen_red_1)
      
      gen_red_2 = gen %>% 
        group_by(Confidence, Class, Specie, Conserved_level, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Member)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2))
      gen_red_2 = as.data.frame(gen_red_2)
      colnames(gen_red_2) = c("Confidence", "Class", "Specie", "Conserved_level", "Counts_LncRNAs", "Total_LncRNAs", "Percentage_LncRNAs")
      TAB_2_LNCRNAS = rbind(TAB_2_LNCRNAS, gen_red_2)
      
      gen_red_3 = gen %>% 
        group_by(Confidence, Class, Specie, Type, .drop=FALSE) %>%
        summarise(
          Counts = n_distinct(Member)) %>%
        mutate(Total = sum(Counts),
               perc = round(Counts/sum(Counts) * 100, 2))
      gen_red_3 = as.data.frame(gen_red_3)
      colnames(gen_red_3) = c("Confidence", "Class", "Specie", "Type", "Counts_LncRNAs", "Total_LncRNAs", "Percentage_LncRNAs")
      TAB_3_LNCRNAS = rbind(TAB_3_LNCRNAS, gen_red_3)
    }
  }
}

write.table(TAB_1_LNCRNAS, paste0(WD2, "/TABLE_LNCRNAS_PERCENTAGE_1.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_2_LNCRNAS, paste0(WD2, "/TABLE_LNCRNAS_PERCENTAGE_2.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_3_LNCRNAS, paste0(WD2, "/TABLE_LNCRNAS_PERCENTAGE_3.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_ALL_LNCRNAS, paste0(WD2, "/TABLE_LNCRNAS_PERCENTAGE_ALL.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

TAB_1_LNCRNAS$Confidence = factor(TAB_1_LNCRNAS$Confidence, levels = confidences)
TAB_1_LNCRNAS$Class = factor(TAB_1_LNCRNAS$Class, levels = classes)

gg1 = ggplot(TAB_1_LNCRNAS, aes(x = Confidence, y = Percentage_LncRNAs, fill = Number_species_by_family)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1", "#6f99ad", "#efc000", "#D87758", "#DDAFD6")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_1.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_1.pdf"), height = 18, width = 16, dpi = 600)

TAB_2_LNCRNAS$Confidence = factor(TAB_2_LNCRNAS$Confidence, levels = confidences)
TAB_2_LNCRNAS$Class = factor(TAB_2_LNCRNAS$Class, levels = classes)

gg2 = ggplot(TAB_2_LNCRNAS, aes(x = Confidence, y = Percentage_LncRNAs, fill = Conserved_level)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_2.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_2.pdf"), height = 18, width = 16, dpi = 600)

TAB_3_LNCRNAS$Confidence = factor(TAB_3_LNCRNAS$Confidence, levels = confidences)
TAB_3_LNCRNAS$Class = factor(TAB_3_LNCRNAS$Class, levels = classes)

gg3 = ggplot(TAB_3_LNCRNAS, aes(x = Confidence, y = Percentage_LncRNAs, fill = Type)) +
  geom_bar(colour = "black", position="fill", stat="identity") +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(aes(y = 1, label = Total_LncRNAs), vjust = -0.3)

ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_3.png"), height = 18, width = 16, dpi = 600)
ggsave(paste0(WD2, "/Stacked_barplot_LncRNAs_3.pdf"), height = 18, width = 16, dpi = 600)


################################################################################
## 6. LNCRNAS: ACCUMULATION TABLE

# Aqui mostramos tanto una tabla como dos figuras globales con los TPMs de lncRNAs 
# conservados (Number of species by family > 1, de la familia donde se encuentra el
# lncRNA) y no conservados (Number of species by family = 1, de la familia donde se 
# encuentra el lncRNA).

# Se representan dos figuras teniendo en cuenta:
#     -Fig 1: confidence level, class, specie y Number_species_by_family.
#     -Fig 2: confidence level, class, specie y type.

cat("\nLncRNAs: Building the accumulation table...\n")

TAB_1_LNCRNAS_TPMs = data.frame()
for (confidence in confidences) {
  for (class in classes) {
    if (file.exists(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"))) {
      gen = read.table(paste0(WD1, "/", confidence, "/", class, "/gen.tsv"), sep = "\t", header = T, quote = "\"")
      gen$"LncRNA" = sapply(strsplit(gen$Member, "-"), "[[", 1)
      gen$"Specie" = sapply(strsplit(gen$Member, "-"), "[[", 2)
      gen$"Number_species_by_family" = rowSums(gen[,species])
      gen$"Type" = ifelse(gen$Number_species_by_family > 1, "Conserved", "Non-conserved")
      gen$"Conserved_level" = ifelse(gen$Number_species_by_family == 2 | gen$Number_species_by_family == 3, "Low-conserved",
                                     ifelse(gen$Number_species_by_family == 4 | gen$Number_species_by_family == 5 | gen$Number_species_by_family == 6, "Medium-conserved",
                                            ifelse(gen$Number_species_by_family == 7 | gen$Number_species_by_family == 8 | gen$Number_species_by_family == 9, "High-conserved", "Non-conserved")))
      
      gen = gen[, c("Family", "Member", "LncRNA", "Specie", "Number_species_by_family", "Type", "Conserved_level")]
      gen$"Confidence" = confidence
      gen$"Class" = class
      for (spe in species) {
        TPMs = read.table(paste0(WD3, "/", spe, "/Salmon/ALL/", flag, "/04-Table/TPMs_summary.tsv"), sep = "\t", header = T, quote = "\"")
        TPMs = TPMs[, c("ID_transcript", "log2.TPMs.1.mean")]
        gen_spe = gen[gen$Specie == spe,]
        gen_TPMs = merge(gen_spe, TPMs, by.x = "LncRNA", by.y = "ID_transcript", all.x = T, all.y = F)
        gen_TPMs = gen_TPMs[!is.na(gen_TPMs$log2.TPMs.1.mean),]
        
        TAB_1_LNCRNAS_TPMs = rbind(TAB_1_LNCRNAS_TPMs, gen_TPMs)
      }
    }
  }
}

write.table(TAB_1_LNCRNAS_TPMs, paste0(WD2, "/TABLE_LNCRNAS_TPMS.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

TAB_1_LNCRNAS_TPMs$Confidence = factor(TAB_1_LNCRNAS_TPMs$Confidence, levels = confidences)
TAB_1_LNCRNAS_TPMs$Class = factor(TAB_1_LNCRNAS_TPMs$Class, levels = classes)
TAB_1_LNCRNAS_TPMs$Specie = factor(TAB_1_LNCRNAS_TPMs$Specie, levels = species)
TAB_1_LNCRNAS_TPMs$Number_species_by_family = factor(TAB_1_LNCRNAS_TPMs$Number_species_by_family, levels = 1:length(species))
TAB_1_LNCRNAS_TPMs$Type = factor(TAB_1_LNCRNAS_TPMs$Type, levels = c("Non-conserved", "Conserved"))

gg1 = ggplot(TAB_1_LNCRNAS_TPMs, aes(x = Number_species_by_family, y = log2.TPMs.1.mean, fill = Confidence)) +
  geom_boxplot(aes(fill = Confidence), size = 0.2, colour = "black", outlier.size = 0.05, position = position_dodge2(width = 0.9, preserve = "single")) + 
  stat_summary(fun = mean, geom = "point", aes(group=Confidence), shape = 20, size = 0.5, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5", "#e18727")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("log2(TPMs + 1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1))

ggsave(paste0(WD2, "/Boxplot_LncRNAs_TPMs_1.png"), height = 8, width = 32, dpi = 600)
ggsave(paste0(WD2, "/Boxplot_LncRNAs_TPMs_1.pdf"), height = 8, width = 32, dpi = 600)

gg2 = ggplot(TAB_1_LNCRNAS_TPMs, aes(x = Confidence, y = log2.TPMs.1.mean, fill = Conserved_level)) +
  geom_boxplot(aes(fill = Conserved_level), size = 0.2, colour = "black", outlier.size = 0.05, position = position_dodge2(width = 0.9, preserve = "single")) + 
  stat_summary(fun = mean, geom = "point", aes(group=Conserved_level), shape = 20, size = 0.5, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#FFFFFF", "#CEEBF7", "#569DBD", "#18485D")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("log2(TPMs + 1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1))

ggsave(paste0(WD2, "/Boxplot_LncRNAs_TPMs_2.png"), height = 12, width = 20, dpi = 600)
ggsave(paste0(WD2, "/Boxplot_LncRNAs_TPMs_2.pdf"), height = 12, width = 20, dpi = 600)

gg3 = ggplot(TAB_1_LNCRNAS_TPMs, aes(x = Confidence, y = log2.TPMs.1.mean, fill = Type)) +
  geom_boxplot(aes(fill = Type), size = 0.2, colour = "black", outlier.size = 0.05, position = position_dodge2(width = 0.9, preserve = "single")) + 
  stat_summary(fun = mean, geom = "point", aes(group=Type), shape = 20, size = 0.5, color = "black", fill = "black", position = position_dodge2(width = 0.75, preserve = "single")) +
  facet_grid(Class ~ Specie) +
  scale_fill_manual(values = c("#bc3c29", "#0072b5")) +
  scale_y_continuous(expand = c(0.02, 0.05)) +
  xlab("") +
  ylab("log2(TPMs + 1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1))

ggsave(paste0(WD2, "/Boxplot_LncRNAs_TPMs_3.png"), height = 12, width = 20, dpi = 600)
ggsave(paste0(WD2, "/Boxplot_LncRNAs_TPMs_3.pdf"), height = 12, width = 20, dpi = 600)

rm(list = ls())
