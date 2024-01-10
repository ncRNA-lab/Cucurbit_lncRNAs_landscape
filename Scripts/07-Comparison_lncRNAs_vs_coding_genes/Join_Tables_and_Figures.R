################################################################################
#
# FIGURES NON-REDUNDANT
#
# Join plots comparing genes and lncRNAs: GC content, Length, Exon number, TPMs 
# and Repeat content.
#
################################################################################

# ATENCION 1: 
# En este script se grafican los resultados de 5 caracxteristicas: contenido en GC, Longitud, numero de exones
# expresion y contenido en repetitivo. Contenido en GC y repetitivo no han dado ningun problema porque iban en 
# porcentaje. Sin embargo, las otras tres caracteristicas (longitud, numero de exones y expresion) han dado muchos
# problemas. Centrandonos en estas 3 caracteristicas, inicialmente se transformaron los datos a logaritmo pero se
# vio que al hacer el boxplot no es lo mismo hacer el logaritmo de la media de los valores normales que la media
# de los logaritmos. En algunos casos cambiaba el orden de mayor a menor de la media para los distintos transcritos.
# Por tanto, esta opcion se descarto por adulterar los resultados. Despues se decidio poner los valores normales
# acotando con coord_cartesian (ATENCION 2) porque sino las cajas eran lineas en el 0 debido a la distancia entre
# el valor maximo y el minimo. Pero esta opcion no eran tan buena porque no conseguias ver toda la distribucion de
# los valores. Finalmente se opto por usar los valores normales pero poner en escala logaritmica el eje. De este
# modo, tenemos todos los valores representados sin meter el sesgo de transformar los valores a logaritmo.

# ATENCION 2: 
# Es importante saber que los valores que se dejan de gráficar al usar scale_y_continuous, no son 
# incluidos en el gráfico para el computo de los boxplot, es decir, de la caja del boxplot y tampoco de la media
# en stat_summary(). La solucion es usar coord_cartesian(ylim = c(0,100)), por ejemplo. De este modo, graficas 
# solo los valores de ese rango (0, 100), pero tenemos en cuenta todos los valores de la tabla para calcular la
# caja del boxplot o la media. Ocurre lo mismo con scale_x_continuous por ejemplo si haces un grafico de densidad.

rm(list = ls())

## 0. INSTALL AND LOAD LIBRARIES

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggExtra))
suppressMessages(library(data.table))
suppressMessages(library(grid))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(options(bitmapType='cairo'))


## 1. VARIABLES

WD = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/07-comparison_lncRNAs_vs_coding_genes"
repeat_content_path = "/mnt/doctorado/3-lncRNAs/Cucurbitaceae/Results/11-TEs_and_genomic_repeats/02-Comparison_Genes_LncRNAs/Figures_and_Tables/C"
species_short_name = c("csa", "cme", "cla", "lsi", "cmo", "car", "cpe", "cma", "mch")
species_long_name = c("C. sativus", "C. melo", "C. lanatus", "L. siceraria", "C. moschata", "C. argyrosperma", "C. pepo", "C. maxima", "M. charantia")
types = c("LncRNAs low", "LncRNAs medium", "LncRNAs high")

if (!dir.exists(paste0(WD, "/ALL"))){
  dir.create(paste0(WD, "/ALL"))
}


## 2. FIGURES TAKING INTO ACCOUNT SPECIES

cat(paste0("\n\n\nFIGURES TAKING INTO ACCOUNT SPECIES..."))

TAB_LONG_FINAL = data.frame()
TAB_WIDE_FINAL = data.frame()
TAB_LONG_SUMMARY_FINAL = data.frame()
TAB_WIDE_SUMMARY_FINAL = data.frame()
TAB_STATISTICS_FINAL = data.frame()

for (type in types) {
  cat(paste0("\n\n\nType: ", type, "..."))
  
  ## 2.1 CREATE LONG AND WIDE TABLE.
  
  ### 2.1.1 TAB LONG
  
  cat(paste0("\n\n\t-LONG TABLE..."))
  TAB_LONG = data.frame()
  
  #### GC CONTENT
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-GC CONTENT: Spe: ", species_short_name[i], "..."))
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "GC")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "GC")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "GC")]
    
    tab = rbind(L, G, IR)
    tab$"Feature" = "GC content"
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  #### EXON NUMBER
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-EXON NUMBER: Spe: ", species_short_name[i], "..."))
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Exons" = NA
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons")]
    
    tab = rbind(L, G, IR)
    tab$"Feature" = "Exon number"
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  #### LENGTH
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-LENGTH: Spe: ", species_short_name[i], "..."))
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Length")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Length")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Length" = NA
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Length")]
    
    tab = rbind(L, G, IR)
    tab$"Feature" = "Length"
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  #### EXPRESSION
  
  # All the species have some NA values due to the detection of duplicated transcripts in the quantification with salmon. These duplicated 
  # transcripts are removed and then they are not quantified. This happens in the same way with LncRNAs and genes. This kind of duplication 
  # isn't intra-locus, so they weren't removed in the redundancy filter. It exists other region in the genome where exists a transcript equal 
  # to this transcript. It happens few times. They will not be plotted because they are kept as NA values.
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-EXPRESSION: Spe: ", species_short_name[i], "..."))
    
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "TPMs.mean")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "TPMs.mean")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"TPMs.mean" = NA
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "TPMs.mean")]
    
    tab = rbind(L, G, IR)
    tab$"Feature" = "Expression"
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  #### REPEAT CONTENT
  
  # The table coming from repeat content analysis only contains transcripts with more than 0% of repeat content. So, when we merge
  # tables we find NA values which will be converted to 0 because they are transcripts with 0% of repeat content.
  
  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-REPEAT CONTENT: Spe: ", species_short_name[i], "..."))
    
    tab_rep = read.table(paste0(repeat_content_path, "/Final_tab-Repeat-NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-COLLAPSED_REPEAT.tsv"), sep = "\t", header = T, quote = "\"")
    tab_rep = tab_rep[tab_rep$spe == species_short_name[i], c("transcript_id", "overlap_per")]
    colnames(tab_rep) = c("ID_transcript", "overlap_per")
  
    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie")]
    
    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie")]
    
    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")
    
    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie")]
    
    tab = rbind(L, G, IR)
    tab = merge(tab, tab_rep, by = "ID_transcript", all = T)
    tab$"Feature" = "Repeat content"
    tab = tab[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code",  "Type", "Specie", "overlap_per", "Feature")]
    colnames(tab) = c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Specie", "Value", "Feature")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Value", "Feature")]
    
    tab[is.na(tab)] = 0
    
    TAB_LONG = rbind(TAB_LONG, tab)
  }
  
  ### 2.1.2 TAB WIDE
  
  cat(paste0("\n\n\t-WIDE TABLE..."))
  TAB_WIDE = data.frame()

  for (i in 1:length(species_short_name)) {
    cat(paste0("\n\t\t-Spe: ", species_short_name[i], "..."))

    tab_rep = read.table(paste0(repeat_content_path, "/Final_tab-Repeat-NR-", str_to_title(unlist(strsplit(type, " "))[2]), "-COLLAPSED_REPEAT.tsv"), sep = "\t", header = T, quote = "\"")
    tab_rep = tab_rep[tab_rep$spe == species_short_name[i], c("transcript_id", "overlap_per")]
    colnames(tab_rep) = c("ID_transcript", "overlap_per")

    tab_1 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_1_NR.tsv"), sep = "\t", header = T, quote = "\"")

    L = tab_1[tab_1$Type == type,]
    L$"Specie" = species_long_name[i]
    L = L[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons", "Length", "GC", "TPMs.mean")]

    G = tab_1[tab_1$Type == "Genes",]
    G$"Specie" = species_long_name[i]
    G = G[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons", "Length", "GC", "TPMs.mean")]

    tab_2 = read.table(paste0(WD, "/", species_short_name[i], "/TAB_FINAL_JOIN_2_NR.tsv"), sep = "\t", header = T, quote = "\"")

    IR = tab_2[tab_2$Type == "Intergenic Regions",]
    IR$"Specie" = species_long_name[i]
    IR$"Length" = NA
    IR$"Exons" = NA
    IR$"TPMs.mean" = NA
    IR$"Class_code" = "ir"
    IR = IR[, c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type", "Specie", "Exons", "Length", "GC", "TPMs.mean")]

    tab = rbind(L, G, IR)
    tab = merge(tab, tab_rep, by = "ID_transcript", all = T)
    colnames(tab) = c("ID_transcript", "Chr", "Start", "End", "Strand", "Origin", "Class_code", "Type.1", "Specie", "Exons", "Length", "GC", "TPMs.mean", "overlap_per")
    tab$"Type.2" = ifelse(tab$Type.1 %in% types, "LncRNAs", ifelse(tab$Type.1 == "Genes", "Genes", "Intergenic Regions"))
    tab = tab[,c("Chr", "Start", "End", "Strand", "Origin", "ID_transcript", "Class_code", "Type.1", "Type.2", "Specie", "Exons", "Length", "GC", "TPMs.mean", "overlap_per")]

    TAB_WIDE = rbind(TAB_WIDE, tab)
  }
  
  # Modify the class codes.
  TAB_LONG[TAB_LONG == "intergenic (u)"] = "u"
  TAB_LONG[TAB_LONG == "antisense (x)"] = "x"
  TAB_LONG[TAB_LONG == "intronic (i)"] = "i"
  TAB_LONG[TAB_LONG == "sense (o/e)"] = "o/e"
  TAB_LONG[TAB_LONG == "gene (=)"] = "pc"
  
  TAB_WIDE[TAB_WIDE == "intergenic (u)"] = "u"
  TAB_WIDE[TAB_WIDE == "antisense (x)"] = "x"
  TAB_WIDE[TAB_WIDE == "intronic (i)"] = "i"
  TAB_WIDE[TAB_WIDE == "sense (o/e)"] = "o/e"
  TAB_WIDE[TAB_WIDE == "gene (=)"] = "pc"
  
  # Create factors
  TAB_LONG$Type.1 = factor(TAB_LONG$Type.1, levels = c("Genes", type, "Intergenic Regions"))
  TAB_LONG$Type.2 = factor(TAB_LONG$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
  TAB_LONG$Class_code = factor(TAB_LONG$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  TAB_LONG$Specie = factor(TAB_LONG$Specie, levels = species_long_name)
  TAB_LONG$Feature = factor(TAB_LONG$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
  
  TAB_WIDE$Type.1 = factor(TAB_WIDE$Type.1, levels = c("Genes", type, "Intergenic Regions"))
  TAB_WIDE$Type.2 = factor(TAB_WIDE$Type.2, levels = c("Genes", "LncRNAs", "Intergenic Regions"))
  TAB_WIDE$Class_code = factor(TAB_WIDE$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  TAB_WIDE$Specie = factor(TAB_WIDE$Specie, levels = species_long_name)
  
  write.table(TAB_LONG, paste0(WD, "/ALL/TAB_LONG_ALL_SPECIES-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(TAB_WIDE, paste0(WD, "/ALL/TAB_WIDE_ALL_SPECIES-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(list = c("tab_1", "tab_2", "tab", "i", "L", "G", "IR", "tab_rep"))
  
  
  ## 2.2 FIGURES
  
  cat(paste0("\n\n\t-FIGURES..."))
  
  my_mean <- function(x) {
    log10(mean(10^x))
  }
  
  #### GC CONTENT
  
  cat(paste0("\n\t\t-GC CONTENT..."))
  
  # Individual
  GC_content = TAB_LONG[TAB_LONG$Feature == "GC content",]
  gg = ggplot(GC_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_wrap(Specie~., nrow = 1) +
    xlab("") +
    ylab("GC content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 18),
          strip.text.x = element_text(size = 19, face = "bold.italic"))

  ggsave(paste0(WD, "/ALL/GC_content-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.png"), height = 5, width = 25, dpi = 600)
  
  rm(list = c("gg", "GC_content"))
  
  # Grid
  GC_content = TAB_LONG[TAB_LONG$Feature == "GC content",]
  gg1 = ggplot(GC_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_grid(Feature~Specie) +
    xlab("") +
    ylab("GC content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_text(size = 22, face = "bold.italic"),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(5.5, 5.5, 0, 5.5))
  
  rm(list = c("GC_content"))

  #### EXON NUMBER

  cat(paste0("\n\t\t-EXON NUMBER..."))
  
  # Individual
  Exon_number = TAB_LONG[TAB_LONG$Feature == "Exon number",]
  gg = ggplot(Exon_number, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_wrap(Specie~., nrow = 1) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Exon number") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 14),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 18),
          strip.text.x = element_text(size = 17, face = "bold.italic"))

  ggsave(paste0(WD, "/ALL/Exon_number-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.png"), height = 5, width = 25, dpi = 600)

  rm(list = c("gg", "Exon_number"))
  
  # Grid
  Exon_number = TAB_LONG[TAB_LONG$Feature == "Exon number",]
  gg2 = ggplot(Exon_number, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_grid(Feature~Specie) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Exon number") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))
  
  rm(list = c("Exon_number"))
  
  #### LENGTH
  
  cat(paste0("\n\t\t-LENGTH..."))
  
  # Individual
  Length = TAB_LONG[TAB_LONG$Feature == "Length",]
  gg = ggplot(Length, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_wrap(Specie~., nrow = 1) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Length") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 14),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 18),
          strip.text.x = element_text(size = 17, face = "bold.italic"))

  ggsave(paste0(WD, "/ALL/Length-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.png"), height = 5, width = 25, dpi = 600)

  rm(list = c("gg", "Length"))
  
  # Grid
  Length = TAB_LONG[TAB_LONG$Feature == "Length",]
  gg3 = ggplot(Length, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_grid(Feature~Specie) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("Length") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))
  
  rm(list = c("Length"))
  
  #### EXPRESSION
  
  cat(paste0("\n\t\t-EXPRESSION..."))
  
  # Individual
  Expression = TAB_LONG[TAB_LONG$Feature == "Expression",]
  gg = ggplot(Expression, aes(x = Class_code, y = Value, fill = Class_code)) + 
    geom_boxplot() + 
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_wrap(Specie~., nrow = 1) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("TPM") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 18),
          strip.text.x = element_text(size = 17, face = "bold.italic")) 
  
  ggsave(paste0(WD, "/ALL/Expression-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.png"), height = 5, width = 25, dpi = 600)
  
  rm(list = c("gg", "Expression"))
  
  # Grid
  Expression = TAB_LONG[TAB_LONG$Feature == "Expression",]
  gg4 = ggplot(Expression, aes(x = Class_code, y = Value, fill = Class_code)) + 
    geom_boxplot() + 
    stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_grid(Feature~Specie) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("") +
    ylab("TPM") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))
  
  rm(list = c("Expression"))
  
  #### REPEAT CONTENT

  cat(paste0("\n\t\t-REPEAT CONTENT..."))
  
  # Individual
  Repeat_content = TAB_LONG[TAB_LONG$Feature == "Repeat content",]
  gg = ggplot(Repeat_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
    facet_wrap(Specie~., nrow = 1) +
    xlab("") +
    ylab("Repeat content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 18),
          strip.text.x = element_text(size = 17, face = "bold.italic"))

  ggsave(paste0(WD, "/ALL/Repeat_content-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.png"), height = 5, width = 25, dpi = 600)

  rm(list = c("gg", "Repeat_content"))
  
  # Grid
  Repeat_content = TAB_LONG[TAB_LONG$Feature == "Repeat content",]
  gg5 = ggplot(Repeat_content, aes(x = Class_code, y = Value, fill = Class_code)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
    scale_fill_manual(
      labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
      values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")
    ) +
    facet_grid(Feature~Specie) +
    xlab("") +
    ylab("Repeat content (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
          legend.key.size = unit(1, 'cm')
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 24),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 23),
          strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
    theme(plot.margin = margin(0, 5.5, 5.5, 5.5))
  
  rm(list = c("Repeat_content"))

  #### FIGURE FEATURES
  
  # Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
  # el siguiente error:
  # 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
  # 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
  # Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
  # -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
  # a logaritmo serian -Inf.
  # 5: Transformation introduced infinite values in continuous y-axis.
  
  cat(paste0("\n\t\t-ALL FEATURES..."))
  gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 1, heights = c(0.21, 0.19, 0.19, 0.19, 0.22))

  ggsave(paste0(WD, "/ALL/FINAL-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.png"), height = 25, width = 25, dpi = 600)
  ggsave(paste0(WD, "/ALL/FINAL-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.pdf"), height = 25, width = 25, dpi = 600)

  rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))
  
  ## 2.3 MEAN AND MEDIAN TABLES
  
  cat(paste0("\n\n\t-MEAN AND MEDIAN TABLE..."))
  TAB_LONG_SUMMARY = TAB_LONG %>% 
    drop_na() %>%
    group_by(Specie, Type.1, Class_code, Feature) %>% 
    summarise(MEAN = mean(Value),
              MEDIAN = median(Value))
  TAB_WIDE_SUMMARY <- pivot_wider(TAB_LONG_SUMMARY, 
                                  id_cols = c("Specie", "Type.1", "Class_code"),
                                  names_from = "Feature",
                                  values_from = c("MEAN", "MEDIAN"),
                                  names_glue = "{Feature}.{.value}")
  colnames(TAB_WIDE_SUMMARY) = gsub(" ", "_", colnames(TAB_WIDE_SUMMARY))
  
  TAB_LONG_SUMMARY$Specie = factor(TAB_LONG_SUMMARY$Specie, levels = species_long_name)
  TAB_LONG_SUMMARY$Type.1 = factor(TAB_LONG_SUMMARY$Type.1, levels = c("Genes", type, "Intergenic Regions"))
  TAB_LONG_SUMMARY$Class_code = factor(TAB_LONG_SUMMARY$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  TAB_LONG_SUMMARY$Feature = factor(TAB_LONG_SUMMARY$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))
  
  TAB_WIDE_SUMMARY$Specie = factor(TAB_WIDE_SUMMARY$Specie, levels = species_long_name)
  TAB_WIDE_SUMMARY$Type.1 = factor(TAB_WIDE_SUMMARY$Type.1, levels = c("Genes", type, "Intergenic Regions"))
  TAB_WIDE_SUMMARY$Class_code = factor(TAB_WIDE_SUMMARY$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
  
  write.table(TAB_LONG_SUMMARY, paste0(WD, "/ALL/TAB_LONG_ALL_SPECIES-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(TAB_WIDE_SUMMARY, paste0(WD, "/ALL/TAB_WIDE_ALL_SPECIES-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  ## 2.4 STATISTICAL ANALYSIS
  
  # Si consideramos cada transcrito como un individuo de la población, como hay más de 30 o 40 individuos por poblacion no seria
  # necesario hacer la spruebas de normalidad y homocedasticidad pudiendo directamente utilizar una prueba parametrica como 
  # el t-test. Sin embargo, son varios los papers que utilizan prueba no parametrica, es decir, Wilcoxon signed-rank test para datos
  # indipendientes. Algunos de estos papers son:
  # - The long non-coding RNA landscape of Candida yeast pathogens (Hrant Hovhannisyan and Toni Gabaldon)
  # - Identification and functional annotation of long intergenic non-coding RNAs in Brassicaceae (Palos et al.)
  
  cat(paste0("\n\n\t-STATISTICAL ANALYSIS..."))
  combinations = as.data.frame(t(combn(c("pc", "u", "x", "i", "o/e", "ir"), 2)))
  rownames(combinations) = NULL
  colnames(combinations) = c("cl1", "cl2")
  
  TAB_STATISTICS = data.frame()
  for (feature in c("GC content", "Exon number", "Length", "Expression", "Repeat content")) {
    for (spe in species_long_name) {
      for (i in 1:nrow(combinations)) {
        cl1 = combinations[i, "cl1"]
        cl2 = combinations[i, "cl2"]
        subset1 = TAB_LONG[TAB_LONG$Feature == feature & TAB_LONG$Specie == spe & TAB_LONG$Class_code == cl1 & !is.na(TAB_LONG$Value), "Value"]
        subset2 = TAB_LONG[TAB_LONG$Feature == feature & TAB_LONG$Specie == spe & TAB_LONG$Class_code == cl2 & !is.na(TAB_LONG$Value), "Value"]
        L1 = length(subset1)
        L2 = length(subset2)
        if (L1 > 0 & L2 > 0) {
          test = wilcox.test(subset1, subset2, paired = FALSE)
          Significance = ifelse(test$p.value < 0.01, "Significant", "Non-Significant")
          row = data.frame(Type.1 = type, Specie = spe, Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = test$statistic, P.value = test$p.value, Significance = Significance, Method = test$method)
        } 
        else{
          row = data.frame(Type.1 = type, Specie = spe, Feature = feature, CL1 = cl1, CL2 = cl2, N1 = L1, N2 = L2, Statistic = NA, P.value = NA, Significance = NA, Method = NA)
        }
        TAB_STATISTICS = rbind(TAB_STATISTICS, row)
      }
    }
  }
  
  write.table(TAB_STATISTICS, paste0(WD, "/ALL/TAB_STATISTICS_ALL_SPECIES-", str_to_title(unlist(strsplit(type, " "))[2]), "-NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  rm(list = c("feature", "spe", "i", "cl1", "cl2", "subset1", "subset2", "L1", "L2", "test", "Significance", "row", "combinations"))
  
  ## 2.5 FINAL TABLES
  
  cat(paste0("\n\n\t-JOIN TABLES TO FINAL TABLES..."))
  TAB_LONG_FINAL = rbind(TAB_LONG_FINAL, TAB_LONG)
  TAB_WIDE_FINAL = rbind(TAB_WIDE_FINAL, TAB_WIDE)
  TAB_LONG_SUMMARY_FINAL = rbind(TAB_LONG_SUMMARY_FINAL, TAB_LONG_SUMMARY)
  TAB_WIDE_SUMMARY_FINAL = rbind(TAB_WIDE_SUMMARY_FINAL, TAB_WIDE_SUMMARY)
  TAB_STATISTICS_FINAL = rbind(TAB_STATISTICS_FINAL, TAB_STATISTICS)
  
  rm(list = c("TAB_LONG", "TAB_WIDE", "TAB_LONG_SUMMARY", "TAB_WIDE_SUMMARY", "TAB_STATISTICS"))
}

# Remove duplicated rows. For example, genes and intergenic regions in each iteration (confidence level) are the same, so they are repeated.
# In summary tables, we don't remove duplicated rows (genes and intergenic regions) to make easier the comparison between genes, lncRNAs and 
# intergenic regions.
TAB_LONG_FINAL = TAB_LONG_FINAL[!duplicated(TAB_LONG_FINAL),]
TAB_WIDE_FINAL = TAB_WIDE_FINAL[!duplicated(TAB_WIDE_FINAL),]

TAB_LONG_FINAL$Type.1 = factor(TAB_LONG_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_WIDE_FINAL$Type.1 = factor(TAB_WIDE_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_LONG_SUMMARY_FINAL$Type.1 = factor(TAB_LONG_SUMMARY_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))
TAB_WIDE_SUMMARY_FINAL$Type.1 = factor(TAB_WIDE_SUMMARY_FINAL$Type.1, levels = c("Genes", types, "Intergenic Regions"))

write.table(TAB_LONG_FINAL, paste0(WD, "/ALL/TAB_LONG_FINAL_ALL_SPECIES-NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_FINAL, paste0(WD, "/ALL/TAB_WIDE_FINAL_ALL_SPECIES-NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_LONG_SUMMARY_FINAL, paste0(WD, "/ALL/TAB_LONG_FINAL_ALL_SPECIES-NR-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_WIDE_SUMMARY_FINAL, paste0(WD, "/ALL/TAB_WIDE_FINAL_ALL_SPECIES-NR-SUMMARY.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(TAB_STATISTICS_FINAL, paste0(WD, "/ALL/TAB_STATISTICS_FINAL_ALL_SPECIES-NR.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


## 2.6 FIGURE FINAL FEATURES

cat(paste0("\n\n\t-FIGURE: ALL FEATURES FINAL..."))

my_mean <- function(x) {
  log10(mean(10^x))
}

#### GC CONTENT

# Grid
GC_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")) +
  facet_grid(Feature~Specie) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold.italic"),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

rm(list = c("GC_content"))

#### EXON NUMBER

# Grid
Exon_number = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Exon_number"))

#### LENGTH

# Grid
Length = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Length"))

#### EXPRESSION

# Grid
Expression = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")) +
  facet_grid(Feature~Specie) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Expression"))

#### REPEAT CONTENT

# Grid
Repeat_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"), 
    values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")
  ) +
  facet_grid(Feature~Specie) +
  xlab("") +
  ylab("Repeat content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 25, margin = margin(r = 1, unit = 'cm')),
        legend.key.size = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

rm(list = c("Repeat_content"))

#### FIGURE FEATURES

gg = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 1, heights = c(0.21, 0.19, 0.19, 0.19, 0.22))

ggsave(paste0(WD, "/ALL/FINAL-NR.png"), height = 25, width = 25, dpi = 600)
ggsave(paste0(WD, "/ALL/FINAL-NR.pdf"), height = 25, width = 25, dpi = 600)

rm(list = c("gg", "gg1", "gg2", "gg3", "gg4", "gg5"))






## 3. FIGURES WITHOUT TAKING INTO ACCOUNT SPECIES.

cat(paste0("\n\n\nFIGURES WITHOUT TAKING INTO ACCOUNT SPECIES..."))

## 3.1 PREPARE TABLE

cat(paste0("\n\n\t-PREPARE TABLE..."))

Genes = TAB_LONG_FINAL[TAB_LONG_FINAL$Type.1 == "Genes", c("ID_transcript", "Class_code", "Type.1", "Feature", "Value")]
Genes_L = Genes
Genes_L[Genes_L == "Genes"] = "LncRNAs low"
Genes_M = Genes
Genes_M[Genes_M == "Genes"] = "LncRNAs medium"
Genes_H = Genes
Genes_H[Genes_H == "Genes"] = "LncRNAs high"

IR = TAB_LONG_FINAL[TAB_LONG_FINAL$Type.1 == "Intergenic Regions", c("ID_transcript", "Class_code", "Type.1", "Feature", "Value")]
IR_L = IR
IR_L[IR_L == "Intergenic Regions"] = "LncRNAs low"
IR_M = IR
IR_M[IR_M == "Intergenic Regions"] = "LncRNAs medium"
IR_H = IR
IR_H[IR_H == "Intergenic Regions"] = "LncRNAs high"

LncRNAs = TAB_LONG_FINAL[!(TAB_LONG_FINAL$Type.1 %in% c("Genes", "Intergenic Regions")), c("ID_transcript", "Class_code", "Type.1", "Feature", "Value")]

TAB_LONG_FINAL_MOD = rbind(Genes_L, Genes_M, Genes_H, IR_L, IR_M, IR_H, LncRNAs)
TAB_LONG_FINAL_MOD$Type.1 = factor(TAB_LONG_FINAL_MOD$Type.1, levels = types)
TAB_LONG_FINAL_MOD$Class_code = factor(TAB_LONG_FINAL_MOD$Class_code, levels = c("pc", "u", "x", "i", "o/e", "ir"))
TAB_LONG_FINAL_MOD$Feature = factor(TAB_LONG_FINAL_MOD$Feature, levels = c("GC content", "Exon number", "Length", "Expression", "Repeat content"))

rm(list = c("Genes", "Genes_L", "Genes_M", "Genes_H", "IR", "IR_L", "IR_M", "IR_H", "LncRNAs"))

## 3.2 FIGURES

cat(paste0("\n\n\t-FIGURES..."))

my_mean <- function(x) {
  log10(mean(10^x))
}

#### GC CONTENT A

cat(paste0("\n\t\t-GC CONTENT A..."))

# Grid
GC_content = TAB_LONG_FINAL_MOD[TAB_LONG_FINAL_MOD$Feature == "GC content",]
gg1 = ggplot(GC_content, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  facet_grid(Feature~Type.1) +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold.italic"),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

rm(list = c("GC_content"))

#### EXON NUMBER A

cat(paste0("\n\t\t-EXON NUMBER A..."))

# Grid
Exon_number = TAB_LONG_FINAL_MOD[TAB_LONG_FINAL_MOD$Feature == "Exon number",]
gg2 = ggplot(Exon_number, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  facet_grid(Feature~Type.1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Exon_number"))

#### LENGTH A

cat(paste0("\n\t\t-LENGTH A..."))

# Grid
Length = TAB_LONG_FINAL_MOD[TAB_LONG_FINAL_MOD$Feature == "Length",]
gg3 = ggplot(Length, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  facet_grid(Feature~Type.1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Length"))

#### EXPRESSION A

cat(paste0("\n\t\t-EXPRESSION A..."))

# Grid
Expression = TAB_LONG_FINAL_MOD[TAB_LONG_FINAL_MOD$Feature == "Expression",]
gg4 = ggplot(Expression, aes(x = Class_code, y = Value, fill = Class_code)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")) +
  facet_grid(Feature~Type.1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Expression"))

#### REPEAT CONTENT A

cat(paste0("\n\t\t-REPEAT CONTENT A..."))

# Grid
Repeat_content = TAB_LONG_FINAL_MOD[TAB_LONG_FINAL_MOD$Feature == "Repeat content",]
gg5 = ggplot(Repeat_content, aes(x = Class_code, y = Value, fill = Class_code)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "lincRNAs", "NAT-lncRNAs", "int-lncRNAs", "SOT-lncRNAs", "Intergenic regions"),
    values = c("#a2ded9", "#ece68f", "#decfa9", "#dcb8b8", "#b5cf9b", "#ca9bcf")
  ) +
  facet_grid(Feature~Type.1) +
  xlab("") +
  ylab("Repeat content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

rm(list = c("Repeat_content"))

#### GC CONTENT B

cat(paste0("\n\t\t-GC CONTENT B..."))

# Grid
GC_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "GC content",]
gg6 = ggplot(GC_content, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")) +
  facet_grid(Feature~"ALL") +
  xlab("") +
  ylab("GC content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 22, face = "bold.italic"),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

rm(list = c("GC_content"))

#### EXON NUMBER B

cat(paste0("\n\t\t-EXON NUMBER B..."))

# Grid
Exon_number = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Exon number",]
gg7 = ggplot(Exon_number, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")) +
  facet_grid(Feature~.) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Exon number") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Exon_number"))

#### LENGTH B

cat(paste0("\n\t\t-LENGTH B..."))

# Grid
Length = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Length",]
gg8 = ggplot(Length, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() +
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")) +
  facet_grid(Feature~.) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Length"))

#### EXPRESSION B

cat(paste0("\n\t\t-EXPRESSION B..."))

# Grid
Expression = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Expression",]
gg9 = ggplot(Expression, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() + 
  stat_summary(fun = my_mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")) +
  facet_grid(Feature~.) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  ylab("TPM") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

rm(list = c("Expression"))

#### REPEAT CONTENT B

cat(paste0("\n\t\t-REPEAT CONTENT B..."))

# Grid
Repeat_content = TAB_LONG_FINAL[TAB_LONG_FINAL$Feature == "Repeat content",]
gg10 = ggplot(Repeat_content, aes(x = Type.1, y = Value, fill = Type.1)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black", fill = "black") +
  scale_fill_manual(
    labels = c("PC genes", "LC-lncRNAs", "MC-lncRNAs", "HC-lncRNAs", "Intergenic regions"), 
    values = c("#a2ded9", "#ce4fd3", "#eb92ef", "#edd0ee", "#cf4141")
  ) +
  facet_grid(Feature~.) +
  xlab("") +
  ylab("Repeat content (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 23),
        strip.background.y = element_rect(fill = "#c7d4de", color = "black", linewidth = 1)) +
  theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

rm(list = c("Repeat_content"))


#### FIGURE FEATURES

# Si tenemos valores de NA para ciertos transcritos, como es el caso de intergenic regions en numero de exones, se generara 
# el siguiente error:
# 1: Removed 169902 rows containing non-finite values (`stat_boxplot()`). 
# 2: Removed 169902 rows containing non-finite values (`stat_summary()`).
# Tambien puede ocurrir que tengas valores iguales a cero como en el caso de TPM. el valor 0 en escala logaritmica tiende a 
# -Inf. Por tanto como la escala del eje y es en escala logaritmica simplemente te da el aviso de que hay valores que transformados
# a logaritmo serian -Inf.
# 5: Transformation introduced infinite values in continuous y-axis.

cat(paste0("\n\t\t-ALL FEATURES A AND B..."))
gg_A = ggarrange(gg1, gg2, gg3, gg4, gg5, ncol = 1, heights = c(0.21, 0.20, 0.20, 0.20, 0.19))
gg_B = ggarrange(gg6, gg7, gg8, gg9, gg10, ncol = 1, heights = c(0.21, 0.20, 0.20, 0.20, 0.19))
gg_final = ggarrange(gg_A, ggparagraph(text="   ", face = "italic", size = 16, color = "black"), gg_B, ncol = 3, widths = c(0.65, 0.08, 0.27))

ggsave(paste0(WD, "/ALL/FINAL-COLLAPSE-NR.png"), height = 25, width = 25, dpi = 600)
ggsave(paste0(WD, "/ALL/FINAL-COLLAPSE-NR.pdf"), height = 25, width = 25, dpi = 600)

rm(list = c("gg_final", "gg_A", "gg_B", "gg1", "gg2", "gg3", "gg4", "gg5", "gg6", "gg7", "gg8", "gg9", "gg10"))

