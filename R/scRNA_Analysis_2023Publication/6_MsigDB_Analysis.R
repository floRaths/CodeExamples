library(readr)
library(qusage)
library(tidyverse)

#### Functions
{
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

split_plot <- function(){
  data %>% select(CellType, Type, mods) %>% pivot_longer(mods) %>% 
    ggplot(aes(x = CellType, y = value, fill = Type)) + 
    geom_split_violin(scale = "width", 
                      alpha = 1, 
                      trim = F,
                      lwd = 1.25, 
                      color = "grey30") + 
    
    geom_boxplot(aes(color = Type), 
                 fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0) +
    
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    style +
    
    facet_wrap("name", ncol = 10, scales =  "free_y")
  
}

vln_plot <- function(){
  data %>% select(CellType, Subcluster, Type, mods) %>% pivot_longer(mods) %>% 
    ggplot(aes(x = Subcluster, y = value, fill = Subcluster)) + 
    geom_violin(scale = "width", 
                      alpha = 1, 
                      trim = F,
                      lwd = 1.25, 
                      color = "grey30") + 
    
    geom_boxplot(fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0) +
    
    scale_fill_manual(values = pal_d3("category20")(20)) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    style +
    
    facet_wrap("name", ncol = 5, scales =  "free_y")
  
}

split_plot_subs <- function(subclusters){
  data %>% 
    select(CellType, Subcluster, Type, mods) %>% 
    pivot_longer(mods) %>% 
    filter(Subcluster %in% subclusters) %>% 
    ggplot(aes(x = factor(Subcluster, levels = levels), y = value, fill = Type)) + 
    geom_split_violin(scale = "width", 
                      alpha = 1, 
                      trim = F,
                      lwd = 1.25, 
                      color = "grey30") + 
    
    geom_boxplot(aes(color = Type), 
                 fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0) +
    
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    style + theme(axis.text.y = element_text(angle = 90)) +
    
    facet_wrap("name", ncol = 5, 
               #scales = "free_y"
               )
  
}
}

setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

GID_cols  <- c("#A6499B", "#FAA42F")


TMCF_snSeq_rerun_scores <- read_csv("2021_ReRun/output/MsigDB/TMCF_snSeq_rerun_scores_v7.4.csv")
gmt                     <- read.gmt("~/Downloads/msigdb.v7.4.symbols.gmt")
Meta                    <- readRDS ("2021_ReRun/Seurat_Objects/MetaData_08.25.2021.rds")
names                   <- colnames(TMCF_snSeq_rerun_scores)[-c(1, 2)]
gene_counts             <- filter(readRDS("2021_ReRun/output/MsigDB/Module_Gene_Counts.rds"), Module %in% names) %>% rename("nGenes" = value)
Searchable_gmt          <- readRDS("2021_ReRun/output/MsigDB/Searchable_gmt.rds")
levels                  <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymph_EC", "Vasc.Acc.",  "Myeloid",  "Lymphoid")

#{
#TMCF_snSeq_rerun_scores <- read_csv("scNuclei-Libraries/Analysis/2021_ReRun/MsigDB/TMCF_snSeq_rerun_scores_v7.4.csv")
#Meta                    <- readRDS("scNuclei-Libraries/Analysis/2021_ReRun/Seurat_Objects/MetaData_07.09.2021.rds")
#names                   <- colnames(TMCF_snSeq_rerun_scores)[-c(1, 2)]
#gene_counts             <- filter(readRDS("scNuclei-Libraries/Analysis/2021_ReRun/MsigDB/Module_Gene_Counts.rds"), Module %in% names) %>% rename(nGenes = "value")
#levels                  <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymph_EC", "Vasc.Acc.",  "Myeloid",  "Lymphoid")
#} #### for HPC



##################################################
CellT <- "Lymphoid"
Celltype_Treatment_Pathways   <- readRDS(paste0("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/MsigDB/CellType_Treatment_", CellT, ".rds"))
Celltype_Treatment_Pathways %>% filter(nGenes > 10) %>% arrange(-CF) %>% rename("score_cisfemale" = "CF", "score_transmale" = "TM") %>% write_csv(paste0("2021_ReRun/Submission/Supp.Tables/Pathways_", CellT, ".csv"))
#CellT <- unique(Sobj$CellType)

data <- TMCF_snSeq_rerun_scores %>% 
  filter(CellType %in% CellT) %>% 
  left_join(rownames_to_column(Meta, "...1"))


# Subcluster Marker Pathways
##################################################
for (i in 1:length(levels)) {

  CellT <- levels[i]
  
  data <- TMCF_snSeq_rerun_scores %>% 
    filter(CellType == CellT) %>% 
    left_join(rownames_to_column(Meta, "X1"))

  cluster <- data %>% pull(Subcluster) %>% unique() %>% as.character()
  Marker_list    <- vector("list", length = length(cluster))
  
  for (i in seq_along(Marker_list)) {
    
    x <- 
      data %>% 
      select(Subcluster, CellType, names) %>% 
      as_tibble() %>% pivot_longer(names, names_to = "Module") %>% 
      mutate(Test = Subcluster, Group = ifelse(Test == cluster[i], "group1", "group2")) %>% 
      
      group_by(Group, Module) %>% mutate(group_avg = mean(value)) %>% 
      select(Module, Group, group_avg) %>% 
      distinct() %>% 
      pivot_wider(names_from = Group, values_from = group_avg) %>% 
      mutate(FC = log2(group1/group2)) %>% 
      add_column(Cluster = cluster[i])
    
    
    y <- 
      data %>% 
      select(Subcluster, CellType, names) %>% 
      as_tibble() %>% pivot_longer(names, names_to = "Module") %>% 
      mutate(Test = Subcluster, Group = ifelse(Test == cluster[i], "group1", "group2")) %>% 
      
      group_by(Module) %>% do(w = wilcox.test(value~Group, data = ., paired = F)) %>% 
      summarise(Module, Wilcox = w$p.value) %>% add_column(Subcluster = cluster[i])
    
    
    Marker_list[[i]] <- x %>% left_join(y, by = c("Module", "Cluster" = "Subcluster"))
    
  }

Bind <- do.call(rbind.data.frame, Marker_list)
Bind %>% ungroup() %>% mutate(CellType = CellT) %>% saveRDS(paste0("2021_ReRun/output/MsigDB/Subcluster_Marker_Pathways_", CellT, ".rds"))

}

# Test per subcluster
##################################################
for (i in 1:length(levels)) {
  
  CellT <- levels[i]
  
  data <- TMCF_snSeq_rerun_scores %>% 
    filter(CellType == CellT) %>% 
    left_join(rownames_to_column(Meta, "X1"))
  
  
  
  wilcox_Msig_SC <- 
    data %>% 
    select(Type, CellType, Subcluster, names) %>% 
    as_tibble() %>% pivot_longer(names, names_to = "Module") %>% 
    group_by(CellType, Subcluster, Module) %>%  
    do(w = wilcox.test(value~Type, data = ., paired = F)) %>% 
    summarise(CellType, Subcluster, Module, Wilcox = w$p.value)
  
  data_avg_SC <- 
    data %>% 
    select(Type, CellType, Subcluster, names) %>% 
    as_tibble() %>% group_by(CellType, Subcluster, Type) %>% 
    mutate_at(vars(names), mean) %>% ungroup %>% 
    distinct()
  
  data_diff_SC <- 
    data_avg_SC %>% 
    pivot_longer(names, names_to = "Module") %>% 
    pivot_wider(names_from = Type, values_from = value) %>% 
    group_by(CellType, Subcluster, Module) %>% 
    mutate(log2FC = log2(TM/CF)) %>% 
    left_join(wilcox_Msig_SC) %>% 
    left_join(gene_counts)
  
  data_diff_SC %>% saveRDS(paste0("scNuclei-Libraries/Analysis/2021_ReRun/MsigDB/Subcluster_Treatment_", CellT, ".rds"))

}

# Test per CellType
##################################################
for (i in 1:length(levels)) {
  
  CellT <- levels[i]

  data <- TMCF_snSeq_rerun_scores %>% 
    filter(CellType == CellT) %>% 
    left_join(rownames_to_column(Meta, "X1"))
  
  
  
  wilcox_Msig_CT <- 
    data %>% 
    select(Type, CellType, names) %>% 
    as_tibble() %>% pivot_longer(names, names_to = "Module") %>% 
    group_by(CellType, Module) %>%  
    do(w = wilcox.test(value~Type, data = ., paired = F)) %>% 
    summarise(CellType, Module, Wilcox = w$p.value)
  
  data_avg_CT <- 
    data %>% 
    select(Type, CellType, names) %>% 
    as_tibble() %>% group_by(CellType, Type) %>% 
    mutate_at(vars(names), mean) %>% ungroup %>% 
    distinct()
  
  data_diff_CT <- 
    data_avg_CT %>% 
    pivot_longer(names, names_to = "Module") %>% 
    pivot_wider(names_from = Type, values_from = value) %>% 
    group_by(CellType, Module) %>% 
    mutate(log2FC = log2(TM/CF)) %>% 
    left_join(wilcox_Msig_CT) %>% 
    left_join(gene_counts)
  
  data_diff_CT %>% saveRDS(paste0("scNuclei-Libraries/Analysis/2021_ReRun/MsigDB/CellType_Treatment_", CellT, ".rds"))

}



Subcluster_Marker_Pathways    <- readRDS(paste0("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/MsigDB/Subcluster_Marker_Pathways_", CellT, ".rds"))
Celltype_Treatment_Pathways   <- readRDS(paste0("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/MsigDB/CellType_Treatment_", CellT, ".rds"))
Subcluster_Treatment_Pathways <- readRDS(paste0("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/MsigDB/Subcluster_Treatment_", CellT, ".rds"))


#mods2 <- 
Subcluster_Marker_Pathways %>% 
  left_join(gene_counts) %>% 
  group_by(Cluster) %>% 
  mutate(rank = group1 * FC) %>% 
  ungroup() %>% 
  filter(Cluster == "lun_2", 
         group1 > 0.1,
         FC > 0,
         nGenes  >= 10) %>% 
  top_n(50, FC) %>% 
  arrange(-FC) #%>%  #%>% pull(Module)


Basal_matrix {
  
  query <- 
    c("REACTOME_SMOOTH_MUSCLE_CONTRACTION",
      "REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY", 
      "KEGG_ADHERENS_JUNCTION",
      
      "REACTOME_PI3K_EVENTS_IN_ERBB4_SIGNALING",
      "BIOCARTA_PDGF_PATHWAY",
      "REACTOME_RUNX3_REGULATES_NOTCH_SIGNALING",
      "BIOCARTA_HDAC_PATHWAY",
      
      "BIOCARTA_IL6_PATHWAY",
      "BIOCARTA_CDMAC_PATHWAY",
      "REACTOME_REGULATION_OF_LOCALIZATION_OF_FOXO_TRANSCRIPTION_FACTORS",
      
      "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION"
    )
  
  
  mat <- 
    data %>% select(X1, CellType, query) %>% 
    left_join(select(rownames_to_column(Meta, "X1"), X1, Subcluster)) %>% 
    pivot_longer(query, names_to = "Module") %>% 
    group_by(Subcluster, Module) %>% 
    summarise(Mean = mean(value)) %>% 
    pivot_wider(names_from = Subcluster, values_from = Mean) %>% 
    column_to_rownames("Module")
  
  pheatmap(mat, scale = "row",
           color = viridis::magma(n = 100, begin = 0, end = 0.9), 
           clustering_distance_rows = "correlation", clustering_method =  "ward.D2",
           clustering_distance_cols = "correlation",
           fontsize_col = 15,
           fontsize_row = 12, 
           treeheight_row = 10, treeheight_col = 0, 
           #annotation_col = Ann,
           main = "Pathway avg_Score")
  
}

LUM_HRneg_matrix {
  
  query <- 
    c("REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING",
      "PID_SYNDECAN_3_PATHWAY", 
      "WP_IL10_ANTIINFLAMMATORY_SIGNALING_PATHWAY",
      
      "REACTOME_SHC1_EVENTS_IN_EGFR_SIGNALING",
      #"REACTOME_ERBB2_ACTIVATES_PTK6_SIGNALING",
       "REACTOME_PI3K_EVENTS_IN_ERBB4_SIGNALING",
    
      
      #"PID_ARF6_DOWNSTREAM_PATHWAY",
      "REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS",
      "REACTOME_SMOOTH_MUSCLE_CONTRACTION",
      
      "REACTOME_CAMK_IV_MEDIATED_PHOSPHORYLATION_OF_CREB",
      "REACTOME_IL_6_TYPE_CYTOKINE_RECEPTOR_LIGAND_INTERACTIONS",
      #"REACTOME_N_GLYCAN_ANTENNAE_ELONGATION",
      
      #"WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_EMBRYONIC_DEVELOPMENT_STAGE_1_OF_4",
      #"REACTOME_ACYL_CHAIN_REMODELLING_OF_PG",
      "REACTOME_PRE_NOTCH_PROCESSING_IN_GOLGI",
      
      
      "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
      
      "REACTOME_DNA_STRAND_ELONGATION"
    )
  
  
  mat <- 
    data %>% select(...1, CellType, query) %>% 
    left_join(select(rownames_to_column(Meta, "...1"), ...1, Subcluster)) %>% 
    pivot_longer(query, names_to = "Module") %>% 
    group_by(Subcluster, Module) %>% 
    summarise(Mean = mean(value)) %>% 
    pivot_wider(names_from = Subcluster, values_from = Mean) %>% 
    mutate(Module = str_to_lower(Module), 
           Module = str_replace(Module, "biocarta", "BIOC:"),  
           Module = str_replace(Module, "hallmark", "HALLMARK:"), 
           Module = str_replace(Module, "wp", "WP:"), 
           Module = str_replace(Module, "kegg", "KEGG:"), 
           Module = str_replace(Module, "reactome", "REAC:"), 
           Module = str_replace(Module, "pid", "PID:"), 
           Module = str_replace_all(Module, "_", " ")) %>% 
    column_to_rownames("Module")
  
  p <- pheatmap(mat[,c(-7)], scale = "row",
           color = viridis::magma(n = 100, begin = 0, end = 0.9), 
           clustering_distance_rows = "correlation", clustering_method =  "ward.D2",
           clustering_distance_cols = "correlation",
           fontsize_col = 20,
           fontsize_row = 20, 
           treeheight_row = 10, treeheight_col = 0, 
           #annotation_col = Ann,
           main = "Pathway avg_Score")
  
}

LUM_HRpos_matrix {

query <- 
  c("BIOCARTA_IGF1_PATHWAY", 
    #"REACTOME_INTERLEUKIN_27_SIGNALING", 
    "PID_ERBB4_PATHWAY", 
  #"WP_ESTROGEN_SIGNALING_PATHWAY",
  "REACTOME_MAPK3_ERK1_ACTIVATION",
  #"REACTOME_BIOTIN_TRANSPORT_AND_METABOLISM",
  #"REACTOME_IL_6_TYPE_CYTOKINE_RECEPTOR_LIGAND_INTERACTIONS",
  "WP_CANONICAL_AND_NONCANONICAL_NOTCH_SIGNALING",
  "REACTOME_INTERLEUKIN_RECEPTOR_SHC_SIGNALING",
  #"WP_ESTROGEN_RECEPTOR_PATHWAY",
  #"REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  #"HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_ANDROGEN_RESPONSE",
  "REACTOME_CGMP_EFFECTS",
  "WP_FATTY_ACID_BIOSYNTHESIS",
  #"WP_LIPID_METABOLISM_PATHWAY", 
  "BIOCARTA_LEPTIN_PATHWAY", 
  #"REACTOME_TFAP2_AP_2_FAMILY_REGULATES_TRANSCRIPTION_OF_GROWTH_FACTORS_AND_THEIR_RECEPTORS",
  #"REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING",
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  #"REACTOME_NITRIC_OXIDE_STIMULATES_GUANYLATE_CYCLASE",
  #"BIOCARTA_RANKL_PATHWAY",
  #"PID_ERBB_NETWORK_PATHWAY",
  "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_INVOLUTION_STAGE_4_OF_4",
  #"WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4",
  #"WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PUBERTY_STAGE_2_OF_4",
  "WP_CHOLESTEROL_BIOSYNTHESIS_PATHWAY",
  #"REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS",
  "REACTOME_METABOLISM_OF_STEROID_HORMONES",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "REACTOME_DNA_STRAND_ELONGATION"
  )


mat <- 
  data %>% select(X1, CellType, query) %>% 
  left_join(select(rownames_to_column(Meta, "X1"), X1, Subcluster)) %>% 
  pivot_longer(query, names_to = "Module") %>% 
  group_by(Subcluster, Module) %>% 
  summarise(Mean = mean(value)) %>% 
  pivot_wider(names_from = Subcluster, values_from = Mean) %>% 
  mutate(Module = str_to_lower(Module), 
         Module = str_replace(Module, "biocarta", "BIOC:"),  
         Module = str_replace(Module, "hallmark", "HALLMARK:"), 
         Module = str_replace(Module, "wp", "WP:"), 
         Module = str_replace(Module, "kegg", "KEGG:"), 
         Module = str_replace(Module, "reactome", "REAC:"), 
         Module = str_replace(Module, "pid", "PID:"), 
         Module = str_replace_all(Module, "_", " ")) %>% 
  column_to_rownames("Module")

p2 <- 
pheatmap(mat[, c("lup_1", "lup_3", "lup_4", "lup_2", "lup_cycling", "lup_trans")], 
         scale = "row", cluster_cols = T,
         color = viridis::magma(n = 100, begin = 0, end = 0.9), 
         clustering_distance_rows = "correlation", clustering_method =  "ward.D2",
         clustering_distance_cols = "correlation",
         fontsize_col = 15,
         fontsize_row = 15,
         treeheight_row = 10, treeheight_col = 0, 
         #annotation_col = Ann,
         main = "Pathway avg_Score")

p2 %>% save_plot(name = "fig_2_modules", scale = 1, w = 16, h = 9, svg = T)

}

mods <- c("PID_INSULIN_PATHWAY",
          "REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE",
          "REACTOME_SIGNALING_BY_INSULIN_RECEPTOR",
          "WP_INSULIN_SIGNALING",
          "KEGG_INSULIN_SIGNALING_PATHWAY"
          )


mods <- "WP_INSULIN_SIGNALING"





#mods <- 
Celltype_Treatment_Pathways %>% 
  mutate(rank = CF * (log2FC*-1)) %>% 
  ungroup() %>% 
  filter(CF > 0.1, 
         nGenes  >= 10,
         Wilcox <= 0.05) %>% 
  top_n(50, CF) %>% 
  arrange(log2FC)

#mods <- 
Celltype_Treatment_Pathways %>% 
  mutate(rank = TM * log2FC) %>% 
  ungroup() %>% 
  filter(TM > 0.1, 
         nGenes  >= 10,
         Wilcox <= 0.05) %>% 
  arrange(-log2FC) #%>% top_n(10, rank) %>% pull(Module)



### up in CF
Subcluster_Treatment_Pathways %>% 
  group_by(Subcluster) %>% 
  mutate(rank = CF * (log2FC*-1)) %>% 
  ungroup() %>% 
  filter(Subcluster == "Macrophage",
         CF > 0.1,
         nGenes  >= 10) %>% 
  top_n(-50, log2FC) %>% 
  arrange(log2FC) #%>% top_n(10, rank) #%>% pull(Module)

### up in TM
Subcluster_Treatment_Pathways %>% 
  group_by(Subcluster) %>% 
  mutate(rank = TM * log2FC) %>% 
  ungroup() %>% 
  filter(Subcluster == "Macrophage", 
         TM > 0.1,
         nGenes  >= 10) %>% 
  top_n(50, log2FC) %>% 
  arrange(-log2FC) #%>% top_n(10, rank) #%>% pull(Module)



Lymphoid_cell_plots{

mods <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE")
p1 <- split_plot_subs(c("CD4_T", "CD8_T")) + ylim(c(-0.1, 0.3)) + theme(strip.text = element_text(size = 10), legend.position = "none")

mods <- c("WP_IL18_SIGNALING_PATHWAY")
p2 <- split_plot_subs(c("CD4_T", "CD8_T")) + theme(strip.text = element_text(size = 10), legend.position = "none")

p1 + p2

p3 <- VlnPlot(Sobj, group.by = "Subcluster", features = c("RUNX2"), idents = c("CD4_T", "CD8_T"), split.by = "Type", cols = GID_cols, pt.size = 0) + theme(legend.position = "none")

p2 + p3

}
Myeloid_cell_plots{

mods <- c("REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
          "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE",
          "REACTOME_CLATHRIN_MEDIATED_ENDOCYTOSIS"
)


split_plot_subs(c("Macrophage")) + theme(strip.text = element_text(size = 10), legend.position = "none")
}


BloodEC_matrix {
  
  query <- 
    c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
      "BIOCARTA_EGFR_SMRTE_PATHWAY",
      
      "REACTOME_N_GLYCAN_ANTENNAE_ELONGATION",
      "REACTOME_REGULATION_OF_IFNG_SIGNALING",
      "PID_LPA4_PATHWAY",
      
      "REACTOME_SIGNALING_BY_NODAL",
      "WP_TRANSCRIPTIONAL_CASCADE_REGULATING_ADIPOGENESIS",
      "REACTOME_REGULATION_OF_IFNG_SIGNALING",

      "REACTOME_RUNX3_REGULATES_NOTCH_SIGNALING",
      "WP_IL9_SIGNALING_PATHWAY",
      "WP_BONE_MORPHOGENIC_PROTEIN_BMP_SIGNALLING_AND_REGULATION",
      "KEGG_NOTCH_SIGNALING_PATHWAY",
      
      "PID_LPA4_PATHWAY",
      "REACTOME_CGMP_EFFECTS"
    )
  
  
  mat <- 
    data %>% select(...1, CellType, query) %>% 
    left_join(select(rownames_to_column(Meta, "...1"), ...1, Subcluster, Type)) %>% 
    pivot_longer(query, names_to = "Module") %>% 
    group_by(Type, Subcluster, Module) %>% 
    summarise(Mean = mean(value)) %>% 
    unite(Subcluster, Type, Subcluster) %>% 
    pivot_wider(names_from = Subcluster, values_from = Mean) %>% 
    mutate(Module = str_to_lower(Module), 
           Module = str_replace(Module, "biocarta", "BIOC:"),  
           Module = str_replace(Module, "hallmark", "HALLMARK:"), 
           Module = str_replace(Module, "wp", "WP:"), 
           Module = str_replace(Module, "kegg", "KEGG:"), 
           Module = str_replace(Module, "reactome", "REAC:"), 
           Module = str_replace(Module, "pid", "PID:"), 
           Module = str_replace_all(Module, "_", " ")) %>% 
    column_to_rownames("Module")
  
  p <- pheatmap(mat, scale = "row",
           color = viridis::magma(n = 100, begin = 0, end = 0.9), 
           clustering_distance_rows = "correlation", clustering_method =  "ward.D2",
           cluster_cols = F,
           fontsize_col = 15,
           fontsize_row = 12, 
           treeheight_row = 10, treeheight_col = 0, 
           #annotation_col = Ann,
           main = "Pathway avg_Score")
  
}


mods <- c("BIOCARTA_ECM_PATHWAY", "KEGG_ECM_RECEPTOR_INTERACTION", "REACTOME_LAMININ_INTERACTIONS", "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX")


mods <- mods_all[1]

#mods <- c("KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION", "HALLMARK_ANGIOGENESIS")
mods <- c("REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
          "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE",
          "REACTOME_CLATHRIN_MEDIATED_ENDOCYTOSIS"
          )

mods <- "REACTOME_VEGFR2_MEDIATED_CELL_PROLIFERATION"
mods <- c("KEGG_ADHERENS_JUNCTION", "KEGG_FOCAL_ADHESION", "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON")

p <- split_plot() + theme(strip.text = element_text(size = 10), legend.position = "none")

levels <- c("Matrix_1", "Matrix_2", "Lipo_F", "Vasc_F")

p <- split_plot_subs(c("Matrix_2", "Lipo_F", "Matrix_1", "Vasc_F")) + theme(strip.text = element_text(size = 10), legend.position = "none")

vln_plot() + theme(strip.text = element_text(size = 10))



Sobj <- AddModuleScore(Sobj, features = list(gmt[["REACTOME_PI3K_EVENTS_IN_ERBB4_SIGNALING"]]), name = "test")
FeaturePlot(Sobj, features = c("test1"), cols = viridis(n = 100, option = "A"), pt.size = 0.5, order = T)











CellType_Treatment_Basal %>% mutate(rank = CF*(log2FC*-1)) %>% arrange(-rank) %>% filter(CF > 0.1, nGenes > 9, Wilcox <= 0.05)
CellType_Treatment_Basal %>% mutate(rank = TM*log2FC)      %>% arrange(-rank) %>% filter(TM > 0.1, nGenes > 9, Wilcox <= 0.05)





recepts <- DE_Response_by_CellType %>% 
  filter(Cluster == "Basal", Gene %in% cabello_aguilar$receptor, 
         avg_log2FC < -0.25) %>% 
  pull(Gene)


ligs <- cabello_aguilar %>% filter(receptor %in% recepts, !is.na(PMIDs)) %>% pull(ligand)


DE_Response_by_CellType %>% filter(Gene %in% ligs, avg_log2FC < -0.25)






mat <- 
data %>% 
  select(Subcluster, c(mods, mods2)) %>% 
  pivot_longer(c(mods, mods2)) %>% 
  group_by(Subcluster, name) %>% 
  summarise(mean = mean(value)) %>% 
  pivot_wider(names_from = Subcluster, values_from = mean) %>% 
  column_to_rownames("name")


pheatmap(mat, scale = "row", 
         color = colorRampPalette(brewer.pal(11, "BrBG"))(100),
         fontsize = 12)

