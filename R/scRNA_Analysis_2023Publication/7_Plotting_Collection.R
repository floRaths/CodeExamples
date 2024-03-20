header {

library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(pheatmap)
library(ggplotify)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(colorspace)

setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

col <- pal_d3("category20")(20)
options(tibble.print_max = 150, tibble.print_min = 35)
`%notin%` <- Negate(`%in%`)


Sobj <- readRDS("2021_ReRun/Seurat_Objects/Sobj_Final_Scaled_MetaUpdt.rds")

save_plot     <- function(data, name, scale, w, h, svg) {
  
  svg = svg
  
  png.path <- paste0("2021_ReRun/Submission/Figures/PNG")
  svg.path <- paste0("2021_ReRun/Submission/Figures/SVG")
  
  if (svg == TRUE) {
    ggsave(plot     = data,
           path     = svg.path,
           filename = paste0(name,".svg"),
           width    = w, 
           height   = h,
           scale    = scale,
           dpi      = 300)
  }
  
  ggsave(plot     = data,
         path     = png.path,
         filename = paste0(name,".png"),
         width    = w, 
         height   = h,
         scale    = scale,
         dpi      = 300)
  
}

DE_Response <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/DE_Response_by_CellType_Additional.rds")
DE_Response_Sub <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/DE_Response_by_Subcluster.rds")
Perc_Avg    <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/Percent&Average_GeneExpression_CellType_Type.rds")
Avg_RNA_Expression_Sample_CellType <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/Avg_RNA_Expression_Sample_CellType.rds")
cabello_aguilar <- read_csv("~/Downloads/cabello_aguilar_LR_PMID.csv")
hs_hgnc_tfs <- read_csv("~/Documents/SCENIC/pySCENIC_Full/scenicdata/Resources/hs_hgnc_tfs.txt", col_names = FALSE)
HRT_data <- read_excel("2021_ReRun/HRT-data-cleaned_2021.06.23.xlsx") %>% rename("Treatment" = "Length of Treatment (months)") 
cluster.translation <- readRDS("2021_ReRun/utilities/Cluster_Translation.rds")

CT_levels <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymph_EC", "Vasc.Acc.", "Myeloid", "Lymphoid")
CT_cols   <- c("#44aa99", "#88ccee", "#117733", "#332288", "#ddcc77", "#882255", "#aa4499", "#cc6677", "#999933", "#dddddd")

SP_levels <- c("CF-3920", "CF-1380", "CF-7780", "CF-2797", "CF-318-813", "CF-0404", "CF-428-112", "CF-2099", "CF-4014", "TM-9469", "TM-1956", "TM-6544", "TM-6477", "TM-8249", "TM-7567", "TM-2768", "TM-3937", "TM-9817")
SP_cols   <- c("#F26C66", "#FFBB78", "#F57F46", "#7FC210", "#2CA02C",    "#BCBD22", "#EDBD1F",    "#FF9896", "#C49C94", "#9467BD", "#C5B0D5", "#17BECF", "#C5D4D2", "#438ABB", "#FA7AA3", "#F7B6D2", "#A8486A", "#3FB8AF")

GID_cols  <- c("#A6499B", "#FAA42F")
Mens_cols  <- c("#FAA42F", "#B966A8",  "#994597")



style <- theme(text = element_text(family = "Lato", size = 35),  
               title = element_text(size = 15, face = "bold"), 
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_rect(fill = "white", color = "grey35", size = 2))


}


query <- c("SEC14L2")

Example_Plot{
  b <- Avg_RNA_Expression_Sample_CellType %>% 
    filter(Gene %in% query, 
           #CellType %in% c("Blood_EC")
           ) %>% 
    #filter(Gene %in% query) %>% 
    left_join(cluster.translation, by = c(CellType = "Cluster")) %>% 
    
    ggplot(aes(#x = factor(Cluster_new, levels = cluster.translation$Cluster_new), 
               x = Type,
               y = avg_RNA, 
               fill = Type)) + 
    
    geom_boxplot(outlier.alpha = 0) +
    
    geom_point(position=position_jitterdodge(dodge.width = 0.25, jitter.width = 0.01),
                aes(#shape = Type, 
                    color = factor(Sample, levels = SP_levels)), 
                size = 4) + 
    
    scale_color_manual(values = SP_cols) + 
    scale_fill_manual (values = GID_cols) + 
    scale_shape_manual(values = c(16, 18)) + 
    labs(x = NULL) + #scale_y_log10() +
    
    style + theme(legend.position = "none", strip.background = element_blank(),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("CellType", scales = "fixed", ncol = 10)
}

# Figure 1 ########################################
figure_1b   {
  
  HRT_data <- read_excel("2021_ReRun/HRT-data-cleaned_2021.06.23.xlsx") %>% rename("Treatment" = "Length of Treatment (months)") 
  
  HRT_data$Treatment <- as.numeric(HRT_data$Treatment)
    
  #p1 <- 
  HRT_data %>% 
    ggplot(aes(fill = factor(Menopause, levels = c("NA", "pre", "post")), y = Age, x = Gender)) + 
    geom_boxplot() + 
    geom_point(position=position_jitterdodge(jitter.width = 0.35), 
               aes(color = factor(Sample, levels = SP_levels)), 
               size = 6) + 
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.75),
                 color = "floralwhite") +
    
    scale_color_manual(values = SP_cols) + 
    scale_fill_manual (values = Mens_cols) + 
    scale_shape_manual(values = c(16, 18)) + 
    labs(title = "Cohort Age and Menopause", 
         x     = "Menopause", 
         y     = "Age", 
         fill  = "Menopause") +
    style + 
    theme(legend.position = "none", 
                  axis.title.x = element_blank())
  
  #p %>% save_plot(name = "fig_1d_global_umaps", scale = 5, w = 3.77, h = 3, svg = F)
  
  
  ### Plotting the scores for the whole dataset as sample averages
  {
    data_avrg  <- readRDS("2021_ReRun/output/DE_ResponseScores_For_Menopausal_Impact/DE_ResponseScores_WholeData_FC_0.25.rds")
    
    my_comparisons <- list(c("NA", "pre"), c("NA", "post"), c("pre", "post"))
    
    ### plot
    #p <- 
    data_avrg %>% 
      
      ggplot(aes(x = factor(Menopause, levels = c("NA", "pre", "post")), 
                 fill = factor(Menopause, levels = c("NA", "pre", "post")), 
                 y = value)) + 
      
      geom_boxplot(color = "grey30",
                   width = 0.5, alpha = 1, 
                   lwd = 1.25, 
                   position = position_dodge(0.9), 
                   outlier.alpha = 0.5) +
      geom_point(position=position_dodge(width = 0.35), 
                 aes(color = factor(Sample, levels = SP_levels)), 
                 size = 6) + 
      stat_summary(fun = "median", 
                   geom = "point", size = 1.5,
                   position = position_dodge(0.9),
                   color = "floralwhite") +
      
      stat_compare_means(comparisons = my_comparisons, test = "wilcox", paired = F) +
      scale_fill_manual(values = Mens_cols) +
      scale_color_manual (values = SP_cols) + 
      
      labs(title = "Menstrual Status DE_Impact", 
           x     = "Gender", 
           y     = "module score (avg/Sample)", 
           fill  = "Menopause") +
      
      #facet_wrap("CellType", scales = "free_y", ncol = 10) + 
      facet_wrap("Module",  scales = "free_y", nrow = 1) + 
      style + theme(legend.position = "bottom")
    
    }
  
  p <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 2))
  
  p %>% save_plot(name = "fig_1c_cohort", scale = 1, w = 16, h = 19, svg = T)
  
}

figure_1d   {

Sobj <- AddMetaData(Sobj, as.data.frame(Sobj@reductions$umap@cell.embeddings))

RNAdata  <- Sobj@meta.data

ATACdata <- readRDS("2021_ReRun/data_resources/umap_tileMat.rds") %>% 
  rename(UMAP_1 = "UMAP1", UMAP_2 = "UMAP2") %>% 
  mutate(Type = str_replace(SampleType, "Trans male", "TM"), Type = str_replace(Type, "Cis female", "CF"))



p1 <- RNAdata[sample(nrow(RNAdata)),] %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) + 
  geom_point(size = 2, alpha = 0.5, stroke = 0) + 
  scale_color_manual(values = CT_cols) + 
  style + theme(legend.position = "none",
                panel.background = element_rect(size = 0.25),
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank())

p2 <- RNAdata[sample(nrow(RNAdata)),] %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Type)) + 
  geom_point(size = 2, alpha = 0.5, stroke = 0) + 
  scale_color_manual(values = GID_cols) + 
  style + theme(legend.position = "none",
                panel.background = element_rect(size = 0.25),
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank())


RNAdata[sample(nrow(RNAdata)),] %>% group_by(Menopause) %>% sample_n(29541) %>% 
  ungroup() %>% sample_n(3*29541) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Menopause)) + 
  geom_point(size = 2, alpha = 0.5, stroke = 0) + 
  scale_color_manual(values = c(Mens_cols[1], col[3], Mens_cols[3])) + 
  style + theme(legend.position = "bottom",
                panel.background = element_rect(size = 0.25),
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank())


ann <- RNAdata[sample(nrow(RNAdata)),] %>% select(Sample, Menopause) %>% distinct()

ATACdata[sample(nrow(ATACdata)),] %>% 
  rownames_to_column("Cell") %>% 
  separate(Cell, into = c("Sample", "Cell"), sep = "_") %>% 
  mutate(Sample = str_replace(Sample, "CF-2797-2", "CF-2797")) %>% 
  left_join(ann) %>% 
  group_by(Menopause) %>% sample_n(10027) %>% 
  ungroup() %>% sample_n(3*10027) %>%  
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Menopause)) + 
  geom_point(size = 2, alpha = 0.5, stroke = 0) + 
    scale_color_manual(values = c(Mens_cols[1], col[3], Mens_cols[3])) + 
  style + theme(legend.position = "none", 
                panel.background = element_rect(size = 0.25),
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank())


p3 <- ATACdata[sample(nrow(ATACdata)),] %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) + 
  geom_point(size = 2, alpha = 0.5, stroke = 0) + 
  scale_color_manual(values = CT_cols) + 
  style + theme(legend.position = "none",
                panel.background = element_rect(size = 0.25),
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank())

p4 <- ATACdata[sample(nrow(ATACdata)),] %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Type)) + 
  geom_point(size = 2, alpha = 0.5, stroke = 0) + 
  scale_color_manual(values = GID_cols) + 
  style + theme(legend.position = "none", 
                panel.background = element_rect(size = 0.25),
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank())


p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 2) 
p <- p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 2) 
p %>% save_plot(name = "fig_1d_global_umaps", scale = 5, w = 3.77, h = 3, svg = F)

}

figure_1e   {

Prop_RNA <- 
  prop.table(table(Sobj$CellType)) %>% 
  as.data.frame() %>% rename(CellType = "Var1")


Prop_ATAC <- 
  prop.table(table(ATACdata$CellType)) %>% 
  as.data.frame() %>% rename(CellType = "Var1")


p1 <- 
  ggplot(data = Prop_RNA, aes(x = "", y = Freq, fill = CellType)) + 
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = round(Freq*100, 1)), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y") +
  
    scale_fill_manual (values = CT_cols) +
  style +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none")
  
p2 <-   
  ggplot(data = Prop_ATAC, aes(x = "", y = Freq, fill = CellType)) + 
    geom_bar(stat = "identity", position = position_fill()) +
    geom_text(aes(label = round(Freq*100, 1)), position = position_fill(vjust = 0.5)) +
    coord_polar(theta = "y") +
    
    scale_fill_manual (values = CT_cols) +
    style +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom")



p3 <- 
  prop.table(table(Sobj$CellType, Sobj$Type), margin = 2) %>% 
  as.data.frame() %>% rename(CellType = "Var1", Type = "Var2") %>% 
  pivot_wider(names_from = Type, values_from = Freq) %>% 
  mutate(log2FC = log2(TM/CF)) %>% 
  left_join(Prop_RNA) %>% 
  arrange(-row_number()) %>% 

  ggplot(aes(x = "", y = Freq, fill = log2FC)) + 
  geom_bar(stat = "identity", position = position_fill()) +
  coord_polar(theta = "y") +
  #geom_text(aes(label = round(log2FC, 2)), position = position_fill(vjust = 0.5)) +
  scale_fill_continuous_divergingx(palette = 'BrBG', mid = 0, limits = c(-1.35, 1.35)) +
  style +
  theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text = element_blank(),
            panel.background = element_blank(),
            legend.position = "right")
  
    
p4 <- 
  prop.table(table(ATACdata$CellType, ATACdata$Type), margin = 2) %>% 
  as.data.frame() %>% rename(CellType = "Var1", Type = "Var2") %>% 
  pivot_wider(names_from = Type, values_from = Freq) %>% 
  mutate(log2FC = log2(TM/CF)) %>% 
  left_join(Prop_ATAC) %>% 
  arrange(-row_number()) %>% 
  
  ggplot(aes(x = "", y = Freq, fill = log2FC)) + 
  geom_bar(stat = "identity", position = position_fill()) +
  coord_polar(theta = "y") +
  #geom_text(aes(label = round(log2FC, 2)), position = position_fill(vjust = 0.5)) +
  scale_fill_continuous_divergingx(palette = 'BrBG', mid = 0, limits = c(-1.35, 1.35)) +
  style +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        legend.position = "right")

p1 + p3 + p2 + p4 + patchwork::plot_layout(ncol = 2)

#p <- p1 + p3 + p2 + p4 + patchwork::plot_layout(ncol = 2)
#p %>% save_plot(name = "fig_1e_proportions_pie", scale = 1, w = 16, h = 9, svg = T)

}

figure_1fgh {

  
  sig <- 
    metadf %>% 
    group_by(CellType) %>% 
    do(w = wilcox.test(nCount_RNA~Type, data = ., paired = F)) %>% 
    summarise(CellType, Wilcox = w$p.value)
  
  
  sig <- 
    metadf %>% 
    group_by(predictedCellType) %>% 
    do(w = wilcox.test(NucleosomeRatio~SampleType, data = ., paired = F)) %>% 
    summarise(predictedCellType, Wilcox = w$p.value)
  
  sig <- 
    metadf %>% 
    group_by(CellType) %>% 
    do(w = wilcox.test(Percent.spliced~Type, data = ., paired = F)) %>% 
    summarise(CellType, Wilcox = w$p.value)
  
  
  
  UMIcountPlot {
    
    metadf <- Sobj@meta.data
    my_comparisons <- list(c("CF", "TM"))
    
    lims <- c(350, 30000)
    
    umi1 <- 
      metadf %>% 
      ggplot(aes(x = Type, y = nCount_RNA, fill = Type)) + 
      
      geom_violin(aes(fill = Type),
                  scale = "width", 
                  trim = F,
                  lwd = 1.25, 
                  color = "grey30") + 
      
      geom_boxplot(aes(color = Type), 
                   fill = "grey30", 
                   width = 0.1, alpha = 1, 
                   lwd = 1.25,
                   position=position_dodge(0.7), outlier.alpha = 0) +
      
      stat_summary(fun = "median", 
                   geom = "point", size = 3,
                   position = position_dodge(0.7),
                   color = "floralwhite") +
      stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
      
      scale_fill_manual(values = GID_cols) +
      scale_color_manual(values = c("grey30", "grey30")) +
      
      scale_y_log10(limits = lims) +
      
      style + 
      coord_flip() + 
      theme(legend.position = "none", axis.title = element_blank())
    
    
    #umi2 <- 
      metadf %>% 
      ggplot(aes(x = factor(CellType, levels = rev(CT_levels)), y = nCount_RNA, fill = Type)) + 
      
      geom_boxplot(aes(color = Type), 
                   #fill = "grey30", 
                   width = 0.5,
                   lwd = 1.25,
                   position=position_dodge(0.7), outlier.alpha = 0) +
      
      stat_summary(fun = "median", 
                   geom = "point", size = 1.5,
                   position = position_dodge(0.7),
                   color = "floralwhite") +
        stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
        
      #geom_text_repel(data = distinct(select(left_join(metadf, sig), CellType, Type, Wilcox)), 
      #          aes(x = factor(CellType, levels = rev(CT_levels)), label = Wilcox), 
      #          y = 1000, color = "black") +
        
      scale_fill_manual(values = GID_cols) +
      scale_color_manual(values = c("grey30", "grey30")) +
      
      scale_y_log10(limits = lims) +
      
      style + 
      coord_flip() + 
      theme(legend.position = "none", axis.title.y = element_blank())
    
    
  }
  
  NucRatioPlot {
  
  metadf <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/data_resources/metadf.rds")
  
  lims <- c(0.3, 1.5)
  
  #nuc1 <- 
  metadf %>% 
    ggplot(aes(x = SampleType, y = NucleosomeRatio, fill = SampleType)) + 
    
    geom_violin(aes(fill = SampleType),
                scale = "width", 
                trim = F,
                lwd = 1.25, 
                color = "grey30") + 
    
    geom_boxplot(aes(color = SampleType), 
                 fill = "grey30", 
                 width = 0.1, alpha = 1, 
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 3,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    scale_y_log10(limits = lims) +
    
    style + 
    coord_flip() + 
    theme(legend.position = "none", axis.title = element_blank())
  
  
  nuc2 <- 
  metadf %>% 
    ggplot(aes(x = factor(predictedCellType, levels = rev(CT_levels)), y = NucleosomeRatio, fill = SampleType)) + 
    
    geom_boxplot(aes(color = SampleType), 
                 #fill = "grey30", 
                 width = 0.5,
                 lwd = 1.25,
                 position=position_dodge(0.7), outlier.alpha = 0) +
    
    stat_summary(fun = "median", 
                 geom = "point", size = 1.5,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = c("grey30", "grey30")) +
    
    scale_y_log10(limits = lims) +
    
    style + 
    coord_flip() + 
    theme(legend.position = "none", axis.title.y = element_blank())
  
  }
  
  SplicePlot   {
    
    X202108_splicedCountByCell <- read_delim("2021_ReRun/output/202108_splicedCountByCell.tsv", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)
    
    metadf <- X202108_splicedCountByCell %>% mutate(X1 = str_replace(...1, ":", "_")) %>% left_join(rownames_to_column(Sobj@meta.data, "X1"))
    
    lims <- c(0, 50)
    
    
    
    #spl1 <- 
      metadf %>%   
      ggplot(aes(x = Type, y = Percent.spliced*100, fill = Type)) + 
      geom_violin(aes(fill = Type),
                  scale = "width", 
                  trim = F,
                  lwd = 1.25, 
                  color = "grey30") + 
      
      geom_boxplot(aes(color = Type), 
                   fill = "grey30", 
                   width = 0.1, alpha = 1, 
                   lwd = 1.25,
                   position=position_dodge(0.7), outlier.alpha = 0) +
      
      stat_summary(fun = "median", 
                   geom = "point", size = 3,
                   position = position_dodge(0.7),
                   color = "floralwhite") +
      
      scale_fill_manual(values = GID_cols) +
      scale_color_manual(values = c("grey30", "grey30")) +
      
      ylim(limits = lims) +
      
      style + 
      coord_flip() + 
      theme(legend.position = "none", axis.title = element_blank())
    
    
    barcharts { 
      metadf %>% 
      group_by(Type) %>% 
      summarise(mean.spliced = (mean(Percent.spliced)*100)) %>% 
      
      ggplot(aes(x = Type, y = mean.spliced, fill = Type)) + 
      geom_col(color = "grey30", size = 1.25) +
      geom_text(aes(label = round(mean.spliced, digits = 1)), 
                position = position_dodge(width = 1)) +
      
      scale_fill_manual(values = GID_cols) +
      scale_color_manual(values = c("grey30", "grey30")) +
      
      ylim(limits = lims) +
      style + 
      coord_flip() + 
      theme(legend.position = "none", axis.title = element_blank())
    }
    
    spl2 <-
      metadf %>% 
      mutate(Percent.spliced = Percent.spliced*100) %>% 
      
      ggplot(aes(x = factor(CellType, levels = rev(CT_levels)), y = Percent.spliced, fill = Type)) + 
      geom_boxplot(aes(color = Type), 
                   width = 0.5,
                   lwd = 1.25,
                   position=position_dodge(0.7), outlier.alpha = 0) +
      
      stat_summary(fun = "median", 
                   geom = "point", size = 1.5,
                   position = position_dodge(0.7),
                   color = "floralwhite") +
      
      scale_fill_manual(values = GID_cols) +
      scale_color_manual(values = c("grey30", "grey30")) +
      
      ylim(limits = lims) +
      style + 
      coord_flip() + 
      theme(legend.position = "none", axis.title.y = element_blank())
  }
  
  umi1 + nuc1 + spl1 + umi2 + nuc2 + spl2 + patchwork::plot_layout(ncol = 3, heights = c(1,5)) 
  p <- umi1 + nuc1 + spl1 + umi2 + nuc2 + spl2 + patchwork::plot_layout(ncol = 3, heights = c(1,5)) 
  p %>% save_plot(name = "fig_1fgh_global_stats", scale = 5, w = 6.5, h = 3, svg = T)
  
  }

figure_s1b  {
  
  CellType_Markers_Final <- read_excel("2021_ReRun/utilities/CellType_Markers_Final.xlsx", sheet = 2)
  #GeneScoreMatrix <- readRDS("2021_ReRun/data_resources/geneScoreMatrix.rds")
  #test <- subset(GeneScoreMatrix, rowData(GeneScoreMatrix)$name %in% CellType_Markers_Final$Marker)
  
  test <- readRDS("2021_ReRun/data_resources/GeneScores_Subset.rds")
  
  df   <- test@assays@data$GeneScoreMatrix %>% as.data.frame()
  el.meta <- test@elementMetadata %>% as.data.frame()
  CellMeta <- test@colData %>% as.data.frame() %>% rownames_to_column("Cell") %>% select(Cell, predictedCellType)
  
  
  
  AT_mat <- bind_cols(el.meta, df) %>% 
    as_tibble() %>% 
    pivot_longer(colnames(test), names_to = "Cell") %>% 
    left_join(CellMeta) %>% select(Cell, Gene = "name", CellType = "predictedCellType", value) %>% 
    group_by(CellType, Gene) %>% 
    summarise(mean = mean(value)) %>% 
    pivot_wider(names_from = CellType, values_from = mean) %>% 
    column_to_rownames("Gene")
  
  mat <- AverageExpression(Sobj, features = CellType_Markers_Final$Marker, group.by = "CellType", slot = "data", assays = "RNA")$RNA
  
  
  A <- 
    pheatmap(mat, 
             scale = "row",
             show_colnames = T, 
             cluster_cols = F, 
             cluster_rows = F,
             color = colorRampPalette(brewer.pal(11, "BrBG"))(100)
    )
  
  B <- 
    pheatmap(AT_mat[CellType_Markers_Final$Marker, CT_levels], 
             scale = "row",
             show_colnames = T, 
             cluster_cols = F, 
             cluster_rows = F,
             color = colorRampPalette(brewer.pal(11, "BrBG"))(100)
    )
  
  
  p <- as.ggplot(A) + as.ggplot(B) + patchwork::plot_layout(ncol = 2)
  p %>% save_plot(name = "fig_s1b_celltype_calling", scale = 3, w = 2.2, h = 7.5, svg = T)
  
}

figure_s1d  {

Avg <- readRDS("2021_ReRun/output/Avg_RNA_Expression_Sample_CellType.rds")

data <- Avg %>% filter(Gene %in% c("AR", "ESR1", "PGR")) 

p <- 
data %>% 
  ggplot(aes(x = factor(CellType, levels = CT_levels), y = avg_RNA, fill = Type)) + 
  geom_boxplot(outlier.alpha = 0.5) + 
  #geom_point(aes(color = Sample, shape = Type), size = 2, 
  #           position = position_dodge(width = 0.45)) +
  
  scale_color_manual(values = SP_cols) + 
  scale_fill_manual (values = GID_cols) + 
  labs(y = "average expr. in sample", x = NULL) +
  style + theme(legend.position = "none", strip.background = element_blank(), 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  
  facet_wrap("Gene", scales = "free_y", ncol = 1)

p %>% save_plot(name = "fig_s1_HR-expression", scale = 1, w = 9, h = 16, svg = T)
  
}

figure_s1e  {
LibraryMetrics <- read_excel("2021_ReRun/utilities/LibraryMetrics.xlsx")

p <- LibraryMetrics %>% 
  right_join(distinct(select(Sobj@meta.data, Sample, Type))) %>% 
  
  ggplot(aes(x = Type, y = `cDNA size`, fill = Type)) + 
  geom_boxplot(outlier.alpha = 0, size = 1) + 
  geom_point(aes(color = factor(Sample, levels = SP_levels), shape = Type), size = 6,
             position = position_dodge(width = 0.45)) + 
  
  scale_color_manual(values = SP_cols) + 
  scale_fill_manual (values = GID_cols) + 
  
  stat_compare_means(paired = F, method = "wilcox") +
  stat_summary(fun = median, geom = "point", color = "floralwhite", size = 3) +
  
  style + theme(legend.position = "none", 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p %>% save_plot(name = "fig_s1e_cDNA_stats", scale = 3, w = 0.8, h = 2, svg = T)

}

figure_s1f {
  
  DE_Response_by_CellType <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/DE_Response_by_CellType.rds")
  
  
  library(gprofiler2)
  library(viridis)
  
  clust = "UP"
  
  query <- DE_Response_by_CellType %>% filter(avg_log2FC > 0) %>% group_by(Gene) %>% summarise(n = n()) %>% filter(n >= 8) %>% pull(Gene)
  
  Results <- gost(query, organism = "hsapiens", evcodes = F)
  
  p1 <- 
  ggplot(head(filter(Results$result, source == "REAC"), 5), 
         aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
    geom_col() +
    scale_fill_viridis(begin = 0.3, end = 0.8, trans = "reverse") +
    theme(axis.title.y = element_blank(), 
          axis.text  = element_text(size = 18), 
          plot.title = element_text(size = 22, face = "bold")) +
    ylab(label = "Number of Genes in Pathway") +
    ggtitle(paste(unique(head(filter(Results$result, source == "REAC"), 15)$source), "-", clust, "Module")) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
    coord_flip() + style
  
  
  clust = "DOWN"
  
  query <- DE_Response_by_CellType %>% filter(avg_log2FC < 0) %>% group_by(Gene) %>% summarise(n = n()) %>% filter(n >= 7) %>% pull(Gene)
  
  Results <- gost(query, organism = "hsapiens", evcodes = F)
  
  p2 <- 
  ggplot(head(filter(Results$result, source == "REAC"), 5), 
         aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
    geom_col() +
    scale_fill_viridis(begin = 0.3, end = 0.8, trans = "reverse") +
    theme(axis.title.y = element_blank(), 
          axis.text  = element_text(size = 18), 
          plot.title = element_text(size = 22, face = "bold")) +
    ylab(label = "Number of Genes in Pathway") +
    ggtitle(paste(unique(head(filter(Results$result, source == "REAC"), 15)$source), "-", clust, "Module")) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
    coord_flip() + style

  p <- p1 + p2 + plot_layout(ncol = 1)  
  
  p %>% save_plot(name = "fig_s1_global_de_pathways", scale = 1, w = 16, h = 9, svg = T)
  
}

figure_1x {  
  
  #Sobj <- readRDS("2021_ReRun/Seurat_Objects/Sobj_Final_Scaled_MetaUpdt.rds")
  
  ### Adding module scores of DE_Response  
  for (i in 1:length(unique(Sobj$CellType))) {
    
    ### Take up and down genes of each celltype and then add it as a module to the dataset
    do <- DE_Response %>% filter(avg_log2FC < -0.5,  p_val_adj <= 0.05, Cluster == unique(Sobj$CellType)[i]) %>% pull(Gene)
    up <- DE_Response %>% filter(avg_log2FC >  0.5,  p_val_adj <= 0.05, Cluster == unique(Sobj$CellType)[i]) %>% pull(Gene)
    
    
    Sobj <- AddModuleScore  (Sobj, assay = "RNA", 
                             list(do), 
                             nbin = 25, 
                             name = paste0(unique(Sobj$CellType)[i],"_downreg"))
    
    Sobj <- AddModuleScore  (Sobj, assay = "RNA", 
                             list(up), 
                             nbin = 25, 
                             name = paste0(unique(Sobj$CellType)[i],"_upreg"))
  }
  
  
  ### dataframe processing for plotting
  {
    ### processing the score dataframe for plotting by cell
    Sobj@meta.data %>% 
      select(Type, Menopause, Sample, contains("reg1"), CellType) %>% 
      pivot_longer(cols = contains("reg1"), names_to = "Module") %>% 
      mutate(Module = str_remove(Module, "1")) %>%
      mutate(Module = str_replace(Module, "R.", "R-")) #%>% saveRDS("2021_ReRun/output/DE_ResponseScores_For_Menopausal_Impact/DE_ResponseScores_FC_0.5.rds")   
    
    
    ### processing the score dataframe for plotting by sample average
    Sobj@meta.data %>% 
      select(Type, Menopause, Sample, contains("reg1"), CellType) %>% 
      pivot_longer(cols = contains("reg1"), names_to = "Module") %>% 
      mutate(Module = str_remove(Module, "1")) %>% 
      mutate(Module = str_replace(Module, "R.", "R-")) %>% 
      group_by(Sample, CellType, Type, Module, Menopause) %>% 
      summarise(value = mean(value)) #%>% saveRDS("2021_ReRun/output/DE_ResponseScores_For_Menopausal_Impact/DE_ResponseScores_SampleAverage_FC_0.5.rds")   
  }
  
  ### Plotting the scores by cells
  {
    data_cells <- readRDS("2021_ReRun/output/DE_ResponseScores_For_Menopausal_Impact/DE_ResponseScores_FC_0.25.rds")   
    
    ### filtering for modules that are relevant for each celltype i.e.(LUM_HR-pos only keeps scores from LUM_HR-pos DE genes)
    data_cells_clean <- 
      data_cells %>% 
      mutate(x = ifelse(str_detect(Module, as.character(CellType)), "A", "B")) %>% 
      filter(x == "A") %>% 
      select(-x)
    
    ### setting levels of modules for plotting/faceting
    levels <- c(unique(filter(data_cells_clean, str_detect(Module, "upreg"))$Module), unique(filter(data_cells_clean, str_detect(Module, "downreg"))$Module))
    
    ### plot
    data_cells_clean %>% 
      mutate(Module = factor(Module, levels = levels)) %>% 
      
      ggplot(aes(x = factor(Type, levels = c("TM", "CF")), 
                 fill = factor(Menopause, levels = c("NA", "pre", "post")), 
                 y = value)) + 
      geom_violin(lwd = 1.25,
                  color = "grey30") + 
      geom_boxplot(color = "grey30",
                   width = 0.1, alpha = 1, 
                   lwd = 1.25,
                   position=position_dodge(0.9), 
                   outlier.alpha = 0) +
      stat_summary(fun = "median", 
                   geom = "point", size = 2,
                   position = position_dodge(0.9),
                   color = "floralwhite") +
      
      scale_fill_manual(values = Mens_cols) +
      
      labs(title = "Menstrual Status DE_Impact", 
           x     = "Gender", 
           y     = "module score", 
           fill  = "Menopause") +
      
      facet_wrap("CellType", scales = "free_y", ncol = 10) + 
      facet_wrap("Module",   scales = "free_y", nrow = 2) + 
      style + theme(legend.position = "none", 
                    axis.text.y = element_text(size = 10))
    
  }
  
  ### Plotting the scores by sample average
  {
    data_avrg  <- readRDS("2021_ReRun/output/DE_ResponseScores_For_Menopausal_Impact/DE_ResponseScores_SampleAverage_FC_0.25.rds")   
    my_comparisons <- list(c("post", "pre"), c("NA", "pre"), c("post", "NA"))
    ### filtering for modules that are relevant for each celltype i.e.(LUM_HR-pos only keeps scores from LUM_HR-pos DE genes)
    data_avrg_clean <- 
      data_avrg %>% 
      mutate(x = ifelse(str_detect(Module, as.character(CellType)), "A", "B")) %>% 
      filter(x == "A") %>% 
      select(-x)
    
    ### setting levels of modules for plotting/faceting
    levels <- c(unique(filter(data_avrg_clean, str_detect(Module, "upreg"))$Module), unique(filter(data_avrg_clean, str_detect(Module, "downreg"))$Module))
    
    ### plot
    data_avrg_clean %>% 
      mutate(Module = factor(Module, levels = levels)) %>% 
      
      ggplot(aes(x = factor(Type, levels = c("TM", "CF")), 
                 fill = factor(Menopause, levels = c("NA", "pre", "post")), 
                 y = value)) + 
      
      geom_boxplot(color = "grey30",
                   width = 0.5, alpha = 1, 
                   lwd = 1.25, 
                   position = position_dodge(0.9), 
                   outlier.alpha = 0) +
      stat_summary(fun = "median", 
                   geom = "point", size = 1.5,
                   position = position_dodge(0.9),
                   color = "floralwhite") +
      
      stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
      scale_fill_manual(values = Mens_cols) +
      
      labs(title = "Menstrual Status DE_Impact", 
           x     = "Gender", 
           y     = "module score (avg/Sample)", 
           fill  = "Menopause") +
      
      facet_wrap("CellType", scales = "free_y", ncol = 10) + 
      facet_wrap("Module",   scales = "free_y", nrow = 2) + 
      style + theme(legend.position = "bottom", 
                    #axis.text.y = element_text(size = 10)
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank())
    
    
    
    
    data_avrg_clean %>% 
      mutate(Module = factor(Module, levels = levels)) %>% 
      
      ggplot(aes(#x = factor(Type, levels = c("TM", "CF")), 
                 x = Menopause,
                 fill = factor(Menopause, levels = c("NA", "pre", "post")), 
                 y = value)) + 
      
      geom_boxplot(color = "grey30",
                   width = 0.5, alpha = 1, 
                   lwd = 1.25, 
                   position = position_dodge(0.9), 
                   outlier.alpha = 0) +
      stat_summary(fun = "median", 
                   geom = "point", size = 1.5,
                   position = position_dodge(0.9),
                   color = "floralwhite") +
      
      stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
      scale_fill_manual(values = Mens_cols) +
      
      labs(title = "Menstrual Status DE_Impact", 
           x     = "Gender", 
           y     = "module score (avg/Sample)", 
           fill  = "Menopause") +
      
      facet_wrap("CellType", scales = "free_y", ncol = 10) + 
      facet_wrap("Module",   scales = "free_y", nrow = 2) + 
      style + theme(legend.position = "bottom", 
                    #axis.text.y = element_text(size = 10)
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank())
    
    
  }
  
  ### Plotting a heatmap of sample averages
  {  
    
    order <- 
      c(unique(filter(Sobj@meta.data, Menopause == "NA")$Sample),
        unique(filter(Sobj@meta.data, Menopause == "pre")$Sample),
        unique(filter(Sobj@meta.data, Menopause == "post")$Sample)
      )
    
    
    Ann <- Sobj@meta.data %>% select(Sample, Menopause) %>% as_tibble() %>% unique() %>% column_to_rownames("Sample")
    
    
    mat_down <- 
      data_avrg_clean %>% 
      ungroup() %>% 
      filter(str_detect(Module, "downreg")) %>% 
      select(Sample, CellType, value) %>% 
      pivot_wider(names_from = Sample, values_from = value) %>% 
      column_to_rownames("CellType")    
    
    mat_up <- 
      data_avrg_clean %>% 
      ungroup() %>% 
      filter(str_detect(Module, "upreg")) %>% 
      select(Sample, CellType, value) %>% 
      pivot_wider(names_from = Sample, values_from = value) %>% 
      column_to_rownames("CellType")    
    
    
    p1 <- 
      pheatmap(mat_up[, order], 
               cluster_cols = F, 
               cluster_row = F, 
               scale = "row",
               annotation_col = Ann,
               color = colorRampPalette(brewer.pal(11, "BrBG"))(100), 
               main = "Upregulated", fontsize = 15)
    p2 <- 
      pheatmap(mat_down[, order], 
               cluster_cols = F, 
               cluster_row = F, 
               scale = "row",
               annotation_col = Ann,
               color = colorRampPalette(brewer.pal(11, "BrBG"))(100), 
               main = "Downregulated", fontsize = 15)
    
    as.ggplot(p1) + as.ggplot(p2)
    
  }
  
  #my_comparisons <- list(c("post", "pre"), c("NA", "pre"), c("post", "NA"))
  #stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F)
  
}

# Figure 2 ########################################

figure_2a   {
  
  Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/LUM_HR-pos_MetaUpdt.rds")
  Sobj <- AddMetaData(Sobj, as.data.frame(Sobj@reductions$umap@cell.embeddings))
  
  RNAdata  <- Sobj@meta.data
  
  p1 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = Subcluster)) + 
    
    geom_point(size = 3, alpha = 0.5, stroke = 0) +
    scale_color_manual(values = pal_d3("category20")(20)) + 
    style + 
    theme(legend.position = "top", 
          panel.background = element_blank(),
          strip.background = element_blank(), 
          strip.text = element_blank(),
          axis.title = element_blank(), 
          axis.text  = element_blank(), 
          axis.ticks = element_blank())
  
  p2 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = factor(Sample, levels = SP_levels))) + 
    
    geom_point(data = select(RNAdata, UMAP_1, UMAP_2), 
               color = "grey90", size = 4, alpha = 1, stroke = 0) + 
    geom_point(size = 2, alpha = 1, stroke = 0) +
    scale_color_manual(values = SP_cols) + 
    style + 
    facet_wrap("Type", ncol = 1) +

    theme(legend.position = "none", 
          strip.background = element_blank(), 
          panel.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_blank(), 
          axis.text  = element_blank(), 
          axis.ticks = element_blank())
  
  p <- p1 + p2 + patchwork::plot_layout(widths = c(2,1))
  
  p %>% save_plot(name = "fig_2a_type_umaps", scale = 1, w = 16, h = 9, svg = T)

    
}

figure_2b   {
  
  Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/LUM_HR-pos_MetaUpdt.rds")
  
  pt.size = 0.5
  order = T
  
  plots{
    
    p1 <- 
      FeaturePlot(Sobj, features = c("AR"), 
                  cols = viridis(n = 100, option = "A"), 
                  pt.size = pt.size, 
                  order = order) + 
      style + theme(legend.position = "none",
                    panel.background = element_blank(),
                    axis.title = element_blank(), 
                    axis.text  = element_blank(), 
                    axis.ticks = element_blank(), 
                    plot.title = element_blank())
    
    p2 <- FeaturePlot(Sobj, features = c("ESR1"), 
                      cols = viridis(n = 100, option = "A"), 
                      pt.size = pt.size, 
                      order = order) + 
      style + theme(legend.position = "none",
                    panel.background = element_blank(),
                    axis.title = element_blank(), 
                    axis.text  = element_blank(), 
                    axis.ticks = element_blank(),
                    plot.title = element_blank())
    
    p3 <- FeaturePlot(Sobj, features = c("PGR"), 
                      cols = viridis(n = 100, option = "A"), 
                      pt.size = pt.size, 
                      order = order) + 
      style + theme(legend.position = "bottom", 
                    panel.background = element_blank(),
                    axis.title = element_blank(), 
                    axis.text  = element_blank(), 
                    axis.ticks = element_blank(),
                    plot.title = element_blank())
    
    p <- p1 + p2 + p3
   
    p %>% save_plot(name = "fig_2b_HR-expression", scale = 1.5, w = 3, h = 8, svg = F)
     
  }
  
}

figure_2c   {

ar_zoscore_lumhrpos <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/data_resources/ar_zoscore_lumhrpos.rds")

p <- 
ar_zoscore_lumhrpos %>% 
  select(AR_689, Type = "SampleType") %>% 
  ggplot(aes(x = Type, y = AR_689, fill = Type)) + 
  
  geom_violin(color = "grey30", 
              scale = "width", size = 1) + 
  
  geom_boxplot(fill = "grey30", color = "grey30", 
               size = 1,
               width = 0.2, 
               outlier.alpha = 0) + 
  
  stat_summary(fun = median, 
               geom = "point", 
               color = "floralwhite", 
               size = 4) +
  
  scale_fill_manual (values = GID_cols) + 
  ggtitle("AR chromVar z-score") +
  style + 
  theme(legend.position = "none")

p %>% save_plot(name = "fig_2c_chromvar", scale = 1, w = 3, h = 8, svg = T)

}

figure_2b   {
  
  Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/LUM_HR-pos_MetaUpdt.rds")
  
  pt.size = 0.5
  order = T
    
    p <- 
      FeaturePlot(Sobj, features = c("CUX2"), 
                  cols = viridis(n = 100, option = "A"), 
                  pt.size = pt.size, 
                  order = order) + 
      style + theme(legend.position = "none",
                    panel.background = element_blank(),
                    axis.title = element_blank(), 
                    axis.text  = element_blank(), 
                    axis.ticks = element_blank(), 
                    plot.title = element_blank())
    
    p %>% save_plot(name = "fig_2h_cux2-expression", scale = 1.5, w = 3.5, h = 3, svg = F)
    
}

figure_s2b  {
  EnrichmentOf            <- readRDS("2021_ReRun/data_resources/motifEnrichment/EnrichmentOfExpressedTfsOrMotifsInDifferentCellTypes_20210719.RDS") #%>% bind_rows()
  DE_Response_by_CellType <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/DE_Response_by_CellType.rds")
  Motif_Summary           <- EnrichmentOf %>% 
    pull(TF) %>% unique() %>% 
    as_tibble() %>% 
    separate(value, into = c("TF", "B"), sep = "_", remove = F) %>% 
    select(Motif = "value", Gene = "TF") %>% 
    left_join(select(left_join(readRDS("2021_ReRun/output/Percent&Average_GeneExpression_CellType_Type.rds"), DE_Response_by_CellType, by = c("Gene", "CellType" = "Cluster")), -p_val, -pct.1, -pct.2))
  
  Enrichment <- 
    EnrichmentOf %>% select(-Gene) %>% 
    left_join(Motif_Summary, by = c("TF" = "Motif", "CellType")) %>% 
    filter(!is.na("Gene")) %>% 
    separate(Comparison, into = c("A", "Type", "C"), sep = " ", extra = "merge") %>% 
    select(-A, -C)
  
  
  
  celltype = "LUM_HR-pos"
  
  
  motifs <- 
    Enrichment %>% 
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    filter(CellType == celltype, Dataset == "Motif.cisBP", 
           abs(Avg_Div) > 0.5,
           CompareProportion*100 >= 10) %>% 
    group_by(Type) %>% top_n(8, Enrichment) %>% pull(TF) %>% unique()
  
  
  p <- 
    Enrichment %>% filter(CellType == celltype, TF %in% c(motifs, "BATF_529", "BATF_129", "BATF..JUN_485")) %>%
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    
    ggplot(aes(y = reorder(TF, Enrichment_Div), x = Enrichment_Div, fill = Type)) + 
    
    geom_col(width = 0.85, alpha = 1) + 
    geom_text(aes(label = nCompare)) +
    scale_fill_manual(values = GID_cols) +
    
    ggtitle(label = paste0(celltype, " Motif Enrichment")) + 
    scale_y_discrete(position = "right") +
    
    geom_vline(xintercept = 0, lwd = 3, color = "white") +
    
    style +
    theme(legend.position = "none", 
          panel.grid.major.y  = element_line(color = "grey70", 
                                             lineend = "round",
                                             linetype = "dotted", 
                                             size = 1))
  
  
  
  q <- 
    Enrichment %>% filter(CellType == celltype, TF %in% motifs) %>%
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    
    ggplot(aes(y = reorder(TF, Enrichment_Div), x = CompareProportion*100, fill = Type)) + 
    geom_col(width = 0.85, alpha = 1, position = position_dodge()) + 
    scale_fill_manual(values = GID_cols) +
    scale_x_reverse() + 
    style+
    theme(legend.position = "none", 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(), 
          panel.grid.major = element_blank())
  
  
  q + p + patchwork::plot_layout(widths = c(1, 4))
  
}  

figure_s2a_proportion {

  celltype = "Basal"
  
  subset <- filter(Sobj@meta.data, CellType == celltype)
  subset$Subcluster <- as.character(subset$Subcluster)
  
  Prop <- 
    as.data.frame(prop.table(table(subset$Subcluster, 
                                   subset$Sample), 
                             margin = 2)) %>% 
    dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
    mutate(Type = substring(.$Sample, 1, 2)) %>% 
    left_join(distinct(select(subset, Sample, Menopause)))


  ### standart plot without Mens Split but with colored dots
  Prop %>% 
    ggplot(aes(x = Cluster, y = Freq, fill = Type)) + 
    geom_boxplot(size = 1, width = 0.5, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    geom_point(position=position_jitterdodge(jitter.width = 0.1),
               aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
    stat_summary(fun = "median", 
                 geom = "point", size = 2,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle(paste0(celltype)) + 
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  
  ### standart plot with Mens Split
  Prop %>% 
    ggplot(aes(x = Cluster, 
               y = Freq, 
               fill = factor(Menopause, levels = c("NA", "pre", "post"))
               )
           ) + 
    geom_boxplot(size = 1, width = 0.65, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    #geom_point(position=position_jitterdodge(jitter.width = 0.1),
    #           aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
    #           size = 4) + 
    stat_summary(fun = "median", 
                 geom = "point", size = 2,
                 position = position_dodge(0.75),
                 color = "floralwhite") +
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle(paste0(celltype)) + 
    style + 
    theme(legend.position = "bottom", 
          axis.title.x = element_blank())
  
  
  ### standart plot with Mens Split
  Prop %>% 
    ggplot(aes(x = Cluster, 
               y = Freq, 
               fill = factor(Menopause, levels = c("NA", "pre", "post"))
               )
           ) + 
    geom_col(position = position_dodge()) +
    
    
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle(paste0(celltype)) + 
    style + 
    theme(legend.position = "bottom", 
          axis.title.x = element_blank())
  
  
  ### stacked bar plot
  Prop <- 
    as.data.frame(prop.table(table(subset$Subcluster, 
                                   subset$Menopause), 
                             margin = 2)) %>% 
    dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq")
  
  #p1 <- 
  Prop %>% 
    ggplot(aes(x = Cluster, y = Freq, fill = factor(Sample, levels = c("NA", "pre", "post")))) + 
    geom_col(position = position_fill()) +
    scale_fill_manual (values = Mens_cols) +
    ggtitle(paste0(celltype)) + 
    style + coord_flip() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  #p1 %>% save_plot(name = "fig_2_proportions", scale = 1, w = 16, h = 9, svg = T)

}

LUM_hrpos_prepost_test {
  
  making_the_data {
    marks_pre <- FindMarkers(subset(Sobj, cells = pre), 
                             ident.1 = "TM", "CF", test.use = "MAST", 
                             logfc.threshold = 0) %>% 
      rownames_to_column("Gene") %>% 
      mutate(Comp = "pre")
    
    
    marks_post <- FindMarkers(subset(Sobj, cells = post), 
                              ident.1 = "TM", "CF", test.use = "MAST", 
                              logfc.threshold = 0) %>% 
      rownames_to_column("Gene") %>% 
      mutate(Comp = "post")
    
    
    data <- bind_rows(marks_post, marks_pre) %>% 
      select(Gene, avg_log2FC, Comp) %>% 
      pivot_wider(names_from = Comp, values_from = avg_log2FC)
  }
  
  data <- readRDS("2021_ReRun/output/LUM_HR-pos_TMvsPRE+POST_MAST.rds")
  
  label <- c("RYR2", "CUX2", "PGR", "AREG", "EREG", "NR4A1", "ADCY2", "PI15", "ACACA", "FASN")
  
  
  p1 <- 
  data %>% 
    filter(!is.na(pre), !is.na(post)) %>% 
    ggplot(aes(x = pre, y = post)) + 
    stat_cor(method = "pearson") + 
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_point(alpha = 0.5, color = "#A6499B") +
    labs(title = c("DE Correlation Pre & Post")) +
    xlab("TM vs. Pre") +
    ylab("TM vs. Post") +
    geom_smooth(method = "lm", se=FALSE, color="black") +
    geom_text_repel(data = filter(data, Gene %in% label), size = 5,
                    aes(label = Gene), 
                    min.segment.length = 0) + 
    style
  
  p2 <- 
  data %>% 
    mutate(post = replace_na(post, 0), pre = replace_na(pre, 0)) %>% 
    mutate(gene_pool = ifelse(post == 0, "unique-pre", "shared"), gene_pool = ifelse(pre == 0, "unique-post", gene_pool)) %>% 
    group_by(gene_pool) %>% 
    summarise(n_DEgenes = n()) %>% 
    mutate(main = "DE_Genes") %>% 
    ggplot(aes(y = main, x = n_DEgenes, fill = gene_pool)) + 
    geom_col(position = position_stack()) + 
    scale_fill_manual(values = Mens_cols[c(1,3,2)]) +
    style + 
    theme(axis.title = element_blank(), legend.position = "bottom")
  
  p1 + p2 + patchwork::plot_layout(ncol = 1, heights = c(10,1))
  
}

figure_s2_motif_distances {


  AR_coOccurence_distance <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/Submission/Figures/Mehran_PDF/figs2c/AR_coOccurence_distance.RDS")

p <- 
AR_coOccurence_distance %>% 
  filter(Motif != "AR_689") %>% 
  ggplot(aes(x = NumOccurence, y = MeanAbsoluteDistance)) + 
  geom_hline(yintercept=median(AR_coOccurence_distance$MeanAbsoluteDistance), size = 0.5, color = "red") +
  geom_point(size = 3) +
  geom_text_repel(data = filter(AR_coOccurence_distance, 
                                Motif != "AR_689", str_detect(Motif, "PGR|NR3C1|NR3C2|SOX13|NFIC|FOXO3|FOXA1|RUNX1|ZNF148|SP1|PURA|STAT2|BATF|CTCF")),
                                #NumOccurence > 2000 | MeanAbsoluteDistance < 230), 
                  aes(label = Motif), min.segment.length = 0) +
  style

p %>% save_plot(name = "fig_s2_motif_distances", scale = 1, w = 16, h = 9, svg = T)

}

# Figure 3 ########################################

figure_3a {
  
  Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/Fibroblast_MetaUpdt.rds")
  Sobj <- AddMetaData(Sobj, as.data.frame(Sobj@reductions$umap@cell.embeddings))
  
  RNAdata  <- Sobj@meta.data
  
  p1 <- 
    RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = Subcluster)) + 
    
    geom_point(size = 3, alpha = 0.5, stroke = 0) +
    scale_color_manual(values = pal_d3("category20")(20)) + 
    style + 
    theme(legend.position = "top", 
          panel.background = element_blank(),
          strip.background = element_blank(), 
          strip.text = element_blank(),
          axis.title = element_blank(), 
          axis.text  = element_blank(), 
          axis.ticks = element_blank())
  
  p2 <- 
    RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = factor(Sample, levels = SP_levels))) + 
    
    geom_point(data = select(RNAdata, UMAP_1, UMAP_2), 
               color = "grey90", size = 4, alpha = 1, stroke = 0) + 
    geom_point(size = 2, alpha = 1, stroke = 0) +
    scale_color_manual(values = SP_cols) + 
    style + 
    facet_wrap("Type", ncol = 2) +
    
    theme(legend.position = "none", 
          strip.background = element_blank(), 
          panel.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_blank(), 
          axis.text  = element_blank(), 
          axis.ticks = element_blank())
  
  p1 + p2 + patchwork::plot_layout(ncol = 1, height = c(2,1))
  
  #p %>% save_plot(name = "fig_2a_type_umaps", scale = 1, w = 16, h = 9, svg = T)
  
}

figure_3b  {
  EnrichmentOf            <- readRDS("2021_ReRun/data_resources/motifEnrichment/EnrichmentOfExpressedTfsOrMotifsInDifferentCellTypes_20210719.RDS") #%>% bind_rows()
  DE_Response_by_CellType <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/DE_Response_by_CellType.rds")
  Motif_Summary           <- EnrichmentOf %>% 
    pull(TF) %>% unique() %>% 
    as_tibble() %>% 
    separate(value, into = c("TF", "B"), sep = "_", remove = F) %>% 
    select(Motif = "value", Gene = "TF") %>% 
    left_join(select(left_join(readRDS("2021_ReRun/output/Percent&Average_GeneExpression_CellType_Type.rds"), DE_Response_by_CellType, by = c("Gene", "CellType" = "Cluster")), -p_val, -pct.1, -pct.2))
  
  Enrichment <- 
    EnrichmentOf %>% select(-Gene) %>% 
    left_join(Motif_Summary, by = c("TF" = "Motif", "CellType")) %>% 
    filter(!is.na("Gene")) %>% 
    separate(Comparison, into = c("A", "Type", "C"), sep = " ", extra = "merge") %>% 
    select(-A, -C)
  
  
  
  celltype = "Fibroblast"
  
  
  motifs <- 
    Enrichment %>% 
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    filter(CellType == celltype, Dataset == "Motif.cisBP", 
           abs(Avg_Div) > 0.5,
           CompareProportion*100 >= 10) %>% 
    group_by(Type) %>% top_n(8, Enrichment) %>% pull(TF) %>% unique()
  
  
  p <- 
    Enrichment %>% filter(CellType == celltype, TF %in% motifs) %>%
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    
    ggplot(aes(y = reorder(TF, Enrichment_Div), x = Enrichment_Div, fill = Type)) + 
    
    geom_col(width = 0.85, alpha = 1) + 
    geom_text(aes(label = nCompare)) +
    scale_fill_manual(values = GID_cols) +
    
    ggtitle(label = paste0(celltype, " Motif Enrichment")) + 
    scale_y_discrete(position = "right") +
    
    geom_vline(xintercept = 0, lwd = 3, color = "white") +
    
    style +
    theme(legend.position = "none", 
          panel.grid.major.y  = element_line(color = "grey70", 
                                             lineend = "round",
                                             linetype = "dotted", 
                                             size = 1))
  
  
  
  q <- 
    Enrichment %>% filter(CellType == celltype, TF %in% motifs) %>%
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    
    ggplot(aes(y = reorder(TF, Enrichment_Div), x = CompareProportion*100, fill = Type)) + 
    geom_col(width = 0.85, alpha = 1, position = position_dodge()) + 
    scale_fill_manual(values = GID_cols) +
    scale_x_reverse() + 
    style+
    theme(legend.position = "none", 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(), 
          panel.grid.major = element_blank())
  
  
  r <- q + p + patchwork::plot_layout(widths = c(1, 4))
  
  r %>% save_plot(name = "fig_3_motif_enrichment", scale = 1, w = 16, h = 9, svg = T)
  
}  

figure_3d  {
  
  p1 <- 
  VlnPlot(Sobj, idents = c("Matrix_1", "Matrix_2", "Lipo_F"), features = c("LAMA2", "LAMB1"), split.by = "Type", cols = GID_cols, ncol = 1, pt.size = 0)
  
  Searchable_gmt %>% left_join(gene_counts, by = c("Term" = "Module")) %>% filter(Gene %in% c("LAMA2", "LAMB1")) %>% as_tibble() %>% arrange(-nGenes) %>% print(n = 100)
  
  gmt <- read.gmt("~/Downloads/msigdb.v7.4.symbols.gmt")
  genes <- gmt$REACTOME_LAMININ_INTERACTIONS
  
  data <- DE_Response_by_CellType %>% 
    filter(Cluster == "Fibroblast") %>% 
    left_join(filter(Searchable_gmt, Term == "REACTOME_LAMININ_INTERACTIONS")) %>% 
    mutate(Term = replace_na(Term, "none"))
  
  p2 <- 
  data  %>% 
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = Term)) + 
    geom_point(alpha = 0.25, size = 3) + 
    geom_point(data = filter(data, Gene %in% genes), alpha = 0.5, size = 4) + 
    geom_text_repel(data = filter(data, abs(avg_log2FC) > 1, -log10(p_val_adj) > 120), aes(label = Gene)) + 
    scale_color_manual(values = c("black", GID_cols[1])) +
    style + 
    theme(legend.position = "bottom", 
          legend.text = element_text(size = 10))
  
  
  
  
  data <- readRDS("2021_ReRun/output/Fibroblast_TMvsPRE+POST_MAST.rds")
  
  label <- gmt$REACTOME_LAMININ_INTERACTIONS
  
  data <- bind_rows(marks_post, marks_pre) %>% 
    select(Gene, avg_log2FC, Comp) %>% 
    pivot_wider(names_from = Comp, values_from = avg_log2FC)
  
  
  p3 <- 
  data %>% ggplot(aes(x = pre, y = post)) + 
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_smooth(method = "lm", se=FALSE, color="black") +
    geom_point(alpha = 0.25, color = "#A6499B") +
    #stat_regline_equation(label.y = 2, aes(label = ..rr.label..)) +
    geom_text_repel(data = filter(data, Gene %in% label),
                    aes(label = Gene), min.segment.length = 0) + style

  p2 + p1 + p3
    
}  

figure_s3a {

marks_fib  <- readRDS("Output/Previous_Analysis/Fibroblast_Subcluster_Markers.rds")

genes1 <- Subcluster_markers %>% filter(CellType == "Fibroblast", !str_detect(gene, "-AS1|-AS2|\\.|LINC0")) %>% group_by(cluster) %>% top_n(25, pct.1) %>% top_n(3, avg_log2FC) %>% arrange(cluster) %>% pull(gene)


mat <- AverageExpression(Sobj, features = c(genes1[!genes1 %in% c("SVIL", "NEGR1", "SEMA5A", "NFIB", "SLIT2", "FLRT2", "ZBTB7C", "PRKG1")], "PPARG", "LAMA2", "ENTPD1", "EPHA3", "LAMB1", "BMPR1B", "TNXB"), group.by = "Subcluster", slot = "data", assays = "RNA")$RNA

A <- 
  pheatmap(t(mat), 
         scale = "column",
         show_colnames = T, 
         cluster_cols = T, 
         cluster_rows = T,
         color = colorRampPalette(brewer.pal(11, "BrBG"))(100)
         )
A %>% save_plot(name = "Fibroblast_Markers", scale = 1.5, w = 20, h = 10, svg = T)

}

figure_4e {
  
  df_long <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/ChromVar_AllCells_Long.rds")
  df_sum <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/ChromVar_Average_CellType.rds")
  
  VlnPlot(Sobj, idents = c("Matrix_1", "Matrix_2", "Lipo_F", "Vasc_F"), features = c("LAMA2", "LAMB1"), split.by = "Type", cols = GID_cols, ncol = 1, pt.size = 0)
  
  query <- c("BACH2")
  
  p1 <- 
    Avg_RNA_Expression_Sample_CellType %>% 
    filter(Gene %in% query, CellType %in% c("Fibroblast")) %>% 
    #filter(Gene %in% query) %>% 
    
    ggplot(aes(x = CellType, y = avg_RNA, fill = Type)) + 
    
    geom_boxplot(outlier.alpha = 0) +
    
    geom_point(position=position_jitterdodge(jitter.width = 0.5),
               aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
    
    scale_color_manual(values = SP_cols) + 
    scale_fill_manual (values = GID_cols) + 
    scale_shape_manual(values = c(16, 18)) + 
    
    style + theme(legend.position = "none", axis.title.x = element_blank(), 
                  axis.text.x = element_blank()) + 
    facet_wrap("Gene", scales = "free_y")
  
  #p2 <- 
  df_long <- readRDS("2021_ReRun/output/Chromvar_FibroMatrix_Split.rds")
    
   df_long %>% 
    ggplot(aes(x = Subtype, y = value, fill = Type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.25, position = position_dodge(width = 0.9)) +
    scale_fill_manual (values = GID_cols) + 
    labs(y = "avg_ChromVar") +
    style + theme(legend.position = "none", 
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("Motif", scales = "free_y")
  
  p2 <- 
    df_sum %>% 
    filter(CellType == "Fibroblast", str_detect(Motif, "BACH2")) %>% 
    ggplot(aes(x = CellType, y = mean, fill = Type)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position=position_jitterdodge(jitter.width = 0.5),
               aes(shape = Type, color = factor(SampleName, levels = SP_levels)), 
               size = 4) + 
    scale_fill_manual (values = GID_cols) + 
    labs(y = "avg_ChromVar") +
    style + theme(legend.position = "none", 
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("Motif", scales = "free_y")
  
  p1 + p2 + patchwork::plot_layout(ncol = 2) 
  
}



Fibroblast_prepost_test {
  
  making_the_data {
    marks_pre <- FindMarkers(subset(Sobj, cells = pre), 
                             ident.1 = "TM", "CF", test.use = "MAST", 
                             logfc.threshold = 0) %>% 
      rownames_to_column("Gene") %>% 
      mutate(Comp = "pre")
    
    
    marks_post <- FindMarkers(subset(Sobj, cells = post), 
                              ident.1 = "TM", "CF", test.use = "MAST", 
                              logfc.threshold = 0) %>% 
      rownames_to_column("Gene") %>% 
      mutate(Comp = "post")
    
    
    data <- bind_rows(marks_post, marks_pre) %>% 
      select(Gene, avg_log2FC, Comp) %>% 
      pivot_wider(names_from = Comp, values_from = avg_log2FC)
  }
  
  data <- readRDS("2021_ReRun/output/Fibroblast_TMvsPRE+POST_MAST.rds") 
  
  label <- c("LAMB1", "LAMA2", "NR4A1", "ADCY2", "THBS1", "VEGFA", "BACH2", "GREB1L")
  
  p1 <- 
    data %>% 
    filter(!is.na(pre), !is.na(post)) %>% 
    ggplot(aes(x = pre, y = post)) + 
    stat_cor(method = "pearson", size = 5) + 
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_point(alpha = 0.5, color = "#A6499B") +
    labs(title = c("DE Correlation Pre & Post in Fibroblasts")) +
    xlab("TM vs. Pre") +
    ylab("TM vs. Post") +
    geom_smooth(method = "lm", se=FALSE, color="black") +
    geom_text_repel(data = filter(data, Gene %in% label), size = 5,
                    aes(label = Gene), 
                    min.segment.length = 0) + 
    style
  
  p2 <- 
    data %>% 
    mutate(post = replace_na(post, 0), pre = replace_na(pre, 0)) %>% 
    mutate(gene_pool = ifelse(post == 0, "unique-pre", "shared"), gene_pool = ifelse(pre == 0, "unique-post", gene_pool)) %>% 
    group_by(gene_pool) %>% 
    summarise(n_DEgenes = n()) %>% 
    mutate(main = "DE_Genes") %>% 
    ggplot(aes(y = main, x = n_DEgenes, fill = gene_pool)) + 
    geom_col(position = position_stack()) + 
    scale_fill_manual(values = Mens_cols[c(1,3,2)]) +
    style + 
    theme(axis.title = element_blank(), legend.position = "bottom")
  
  p1 + p2 + patchwork::plot_layout(ncol = 1, heights = c(10,1))
  
}


# Figure 4 ########################################

figure_4e {
Avg_RNA_Expression_Sample_CellType %>% 
  filter(Gene %in% c("ENAH", "ITGA2", "ITGB8"), 
         CellType == "LUM_HR-neg") %>% 
  
  ggplot(aes(x = Type, y = avg_RNA, fill = Type)) +
  geom_boxplot(size = 1, width = width, 
               position = position_dodge(width = 0.75),
               color = "grey30",
               outlier.alpha = 0) + 
  
  geom_point(position=position_jitterdodge(jitter.width = 0.1),
             aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
             size = 5) + 
  
  scale_fill_manual (values = GID_cols) + 
  scale_color_manual(values = SP_cols) +
  ggtitle("ENAH") + 
  style + 
  theme(legend.position = "none", 
        axis.title.x = element_blank()) + 
  facet_wrap("Gene", scales = "free_y")
}

figure_4h {
  Avg_RNA_Expression_Sample_CellType %>% 
    filter(Gene %in% c("ACTA2", "OXTR"), 
           CellType == "Basal") %>% 
    
    ggplot(aes(x = Type, y = avg_RNA, fill = Type)) +
    geom_boxplot(size = 1, width = 0.5, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    stat_summary(fun = median, geom = "point", shape = 18, color = "grey30", size = 8) +
    
    geom_point(position=position_jitterdodge(jitter.width = 0.1),
               aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
               size = 5) + 
    
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Basal") + 
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank()) + 
    facet_wrap("Gene", scales = "free_y", ncol = 1)
}

figure_4f {
  
  Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/Basal_MetaUpdt_3.rds")
  Sobj <- AddMetaData(Sobj, as.data.frame(Sobj@reductions$umap@cell.embeddings))
  
  RNAdata  <- Sobj@meta.data
  
  
  p <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = RNA_snn_res.0.1)) + 
    
    geom_point(size = 2, alpha = 0.5, stroke = 0) + 
    
    scale_color_manual(values = pal_d3("category20")(20)) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank()) + coord_fixed()
  
  
  q <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = factor(Sample, levels = SP_levels))) + 
    
    geom_point(data = select(RNAdata, UMAP_1, UMAP_2), 
               color = "grey90", size = 4, alpha = 1, stroke = 0) + 
    
    geom_point(size = 1, alpha = 1, stroke = 0) + 
    
    scale_color_manual(values = SP_cols) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank()) + 
    facet_wrap("Type") + coord_fixed()
  
  
  x <- p + q + plot_layout(ncol = 1, heights = c(2,1))
  
  save_plot(x, "Basal_Umaps", w = 6, h = 9, scale = 1, svg = T)
  
  pt.sze = 0.1
  
  p1 <- FeaturePlot(Sobj, features = c("NRP1"), cols = viridis(n = 100, option = "A"), ncol = 1, pt.size = 0.25, order = T) + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  p1 <- FeaturePlot(Sobj, features = c("KEGG_FOCAL_ADHESION1"), cols = viridis(n = 100, option = "A"), ncol = 1, pt.size = pt.sze, order = T) + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + coord_fixed()
  p2 <- FeaturePlot(Sobj, features = c("REACTOME_CELL_JUNCTION_ORGANIZATION1"), cols = viridis(n = 100, option = "A"), ncol = 1, pt.size = pt.sze, order = T) + theme(legend.position = "none",axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + coord_fixed()
  p3 <- FeaturePlot(Sobj, features = c("REACTOME_SMOOTH_MUSCLE_CONTRACTION1"), cols = viridis(n = 100, option = "A"), ncol = 1, pt.size = pt.sze, order = T) +  theme(legend.position = "none",axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + coord_fixed()
  p2 <- FeaturePlot(Sobj, features = c("ACTA2"), cols = viridis(n = 100, option = "A"), ncol = 1, pt.size = 0.25, order = T) + theme(legend.position = "none",axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  p3 <- FeaturePlot(Sobj, features = c("OXTR"), cols = viridis(n = 100, option = "A"), ncol = 1, pt.size = 0.25, order = T) +  theme(legend.position = "none",axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  
  
  y <- p1 + p2 + p3 + plot_layout(ncol = 1)
  
  save_plot(y, "Basal_Pathways", w = 6, h = 9, scale = 1, svg = T)
  
  
  
}

figure_4g  {
  
  EnrichmentOf            <- readRDS("2021_ReRun/data_resources/motifEnrichment/EnrichmentOfExpressedTfsOrMotifsInDifferentCellTypes_20210719.RDS") #%>% bind_rows()
  DE_Response_by_CellType <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/DE_Response_by_CellType.rds")
  Motif_Summary           <- EnrichmentOf %>% 
    pull(TF) %>% unique() %>% 
    as_tibble() %>% 
    separate(value, into = c("TF", "B"), sep = "_", remove = F) %>% 
    select(Motif = "value", Gene = "TF") %>% 
    left_join(select(left_join(readRDS("2021_ReRun/output/Percent&Average_GeneExpression_CellType_Type.rds"), DE_Response_by_CellType, by = c("Gene", "CellType" = "Cluster")), -p_val, -pct.1, -pct.2))
  
  Enrichment <- 
    EnrichmentOf %>% select(-Gene) %>% 
    left_join(Motif_Summary, by = c("TF" = "Motif", "CellType")) %>% 
    filter(!is.na("Gene")) %>% 
    separate(Comparison, into = c("A", "Type", "C"), sep = " ", extra = "merge") %>% 
    select(-A, -C)
  
  
  
  celltype = "Basal"
  
  
  motifs <- 
    Enrichment %>% 
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    filter(CellType == celltype, 
           Dataset == "Motif.cisBP", 
           abs(Avg_Div) > 3, 
           CompareProportion*100 >= 10) %>% arrange(-Enrichment_Div) %>% 
    group_by(Type) %>% top_n(8, Enrichment) %>% pull(TF) %>% unique()
  
  
  p <- 
    Enrichment %>% filter(CellType == celltype, TF %in% motifs) %>%
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    
    ggplot(aes(y = reorder(TF, Enrichment_Div), x = Enrichment_Div, fill = Type)) + 
    
    geom_col(width = 0.85, alpha = 1) + 
    geom_text(aes(label = nCompare)) +
    scale_fill_manual(values = GID_cols) +
    
    ggtitle(label = paste0(celltype, " Motif Enrichment")) + 
    scale_y_discrete(position = "right") +
    
    geom_vline(xintercept = 0, lwd = 3, color = "white") +
    
    style +
    theme(legend.position = "none", 
          panel.grid.major.y  = element_line(color = "grey70", 
                                             lineend = "round",
                                             linetype = "dotted", 
                                             size = 1))
  
  
  
  q <- 
    Enrichment %>% filter(CellType == celltype, TF %in% motifs) %>%
    mutate(Enrichment_Div = ifelse(Type == "Trans", Enrichment, Enrichment * -1), 
           Avg_Div = ifelse(Type == "Trans", avg_TM, avg_CF * -1)) %>% 
    
    ggplot(aes(y = reorder(TF, Enrichment_Div), x = CompareProportion*100, fill = Type)) + 
    geom_col(width = 0.85, alpha = 1, position = position_dodge()) + 
    scale_fill_manual(values = GID_cols) +
    scale_x_reverse() + 
    style+
    theme(legend.position = "none", 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(), 
          panel.grid.major = element_blank())
  
  
  x <- q + p + patchwork::plot_layout(widths = c(1, 4))
  
  x %>% save_plot(name = "fig_4_enrichment", scale = 1, w = 14, h = 16, svg = T)
  
}  

figure_4_ACTA_OXTR {

query <- c("NRP1", "ACTA2", "OXTR")

p <- Avg_RNA_Expression_Sample_CellType %>% 
    filter(Gene %in% query, CellType %in% c("Basal")) %>% 
    #filter(Gene %in% query) %>% 
    left_join(cluster.translation, by = c(CellType = "Cluster")) %>% 
    
    ggplot(aes(x = factor(Cluster_new, levels = cluster.translation$Cluster_new), y = avg_RNA, fill = Type)) + 
    
    geom_boxplot(outlier.alpha = 0.5) +
    
    geom_point(position=position_jitterdodge(jitter.width = 0.5),
               aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
    
    scale_color_manual(values = SP_cols) + 
    scale_fill_manual (values = GID_cols) + 
    scale_shape_manual(values = c(16, 18)) + 
    labs(x = NULL) +
    
    style + theme(legend.position = "none", strip.background = element_blank(),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("Gene", scales = "free_y", ncol = 3)

p %>% save_plot(name = "fig_4_acta_oxtr", scale = 1, w = 16, h = 9, svg = T)

}



df_sum %>% 
  filter(CellType == "LUM_HR-neg", str_detect(Motif, "ESRR")) %>% 
  ggplot(aes(x = CellType, y = mean, fill = Type)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.5),
             aes(shape = Type, color = factor(SampleName, levels = SP_levels)), 
             size = 4) + 
  scale_fill_manual (values = GID_cols) + 
  labs(y = "avg_ChromVar") +
  style + theme(legend.position = "none", 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  facet_wrap("Motif", scales = "free_y")


figure_4e {

query <- c("NFIC", "BACH2")

p1 <- 
  Avg_RNA_Expression_Sample_CellType %>% 
    filter(Gene %in% query, CellType %in% c("Basal")) %>% 
    #filter(Gene %in% query) %>% 
    
    ggplot(aes(x = CellType, y = avg_RNA, fill = Type)) + 
    
    geom_boxplot(outlier.alpha = 0) +
    
    geom_point(position=position_jitterdodge(jitter.width = 0.5),
               aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
    
    scale_color_manual(values = SP_cols) + 
    scale_fill_manual (values = GID_cols) + 
    scale_shape_manual(values = c(16, 18)) + 
    
    style + theme(legend.position = "none", axis.title.x = element_blank(), 
                  axis.text.x = element_blank()) + 
    facet_wrap("Gene", scales = "free_y")


  #df_long %>% 
  #  filter(CellType == "Basal", str_detect(Motif, "BACH2|NFIC")) %>% 
  #  ggplot(aes(x = Type, y = value)) +
  #  geom_violin() +
  #  geom_boxplot(width = 0.25) +
  #  facet_wrap("Motif", scales = "free_y")
  
  #p2 <- 
  df_sum %>% 
    filter(CellType == "Basal", str_detect(Motif, "BACH2|NFIC")) %>% 
    ggplot(aes(x = CellType, y = mean, fill = Type)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position=position_jitterdodge(jitter.width = 0.5),
               aes(shape = Type, color = factor(SampleName, levels = SP_levels)), 
               size = 4) + 
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) + 
    scale_shape_manual(values = c(16, 18)) +
    labs(y = "avg_ChromVar") +
    style + theme(legend.position = "none", 
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("Motif", scales = "free_y")
  
p <- p1 + p2 + patchwork::plot_layout(ncol = 1) 

p %>% save_plot(name = "fig_4_bach2_nfic", scale = 1, w = 16, h = 9, svg = T)
    
}

figure_s4a {
  
  Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/Epithelial_V2.rds")
  Sobj <- AddMetaData(Sobj, as.data.frame(Sobj@reductions$umap@cell.embeddings))
  
  Modules_List <- readRDS("~/Box/Knott_Lab/Flo/Projects/Organoids/Analysis/2019.05.28_Estrogen_24-48h/utilities/Modules_List.rds")
  Sobj <- AddModuleScore  (Sobj, assay = "RNA", list((top_n(Modules_List$LUMA, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUMA")
  Sobj <- AddModuleScore  (Sobj, assay = "RNA", list((top_n(Modules_List$LUPR, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUPR")
  Sobj <- AddModuleScore  (Sobj, assay = "RNA", list((top_n(Modules_List$MASC, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "MASC")
  Sobj <- AddModuleScore  (Sobj, assay = "RNA", list((top_n(Modules_List$STRM, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "STRM")
  rm(Modules_List)
  
  
  
  RNAdata  <- Sobj@meta.data
  
  p1 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = factor(CellType, levels = CT_levels))) + 
    
    geom_point(size = 2, alpha = 0.5, stroke = 0) + 
    
    scale_color_manual(values = CT_cols) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank())
  
  p2 <- 
  FeaturePlot(Sobj, features = c("LUMA1"), 
              cols = viridis(n = 100, option = "A"), 
              pt.size = 0.5, 
              order = T) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank(), 
                  plot.title = element_blank())
  p3 <- 
  FeaturePlot(Sobj, features = c("LUPR1"), 
              cols = viridis(n = 100, option = "A"), 
              pt.size = 0.5, 
              order = T) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank(), 
                  plot.title = element_blank())
  p4 <- 
  FeaturePlot(Sobj, features = c("MASC1"), 
              cols = viridis(n = 100, option = "A"), 
              pt.size = 0.5, 
              order = T) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank(), 
                  plot.title = element_blank())

  p <- p1 + p2 + p3 + p4 + patchwork::plot_layout(nrow = 1)
  
  p %>% save_plot(name = "fig_4_epithelial_calling", scale = 1, w = 16, h = 4, svg = T)
  
}

figure_4f {
  
  Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/LUM_HR-neg_MetaUpdt.rds")
  Sobj <- AddMetaData(Sobj, as.data.frame(Sobj@reductions$umap@cell.embeddings))
  
  RNAdata  <- Sobj@meta.data
  
  p1 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = Subcluster)) + 
    
    geom_point(size = 2, alpha = 0.5, stroke = 0) + 
    
    scale_color_manual(values = pal_d3("category20")(20)) + 
    style + theme(legend.position = "top",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank())
  
  p2 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = factor(Sample, levels = SP_levels))) + 
    
    geom_point(data = select(RNAdata, UMAP_1, UMAP_2), 
               color = "grey90", size = 4, alpha = 1, stroke = 0) + 
    
    geom_point(size = 1.5, alpha = 1, stroke = 0) + 
    
    scale_color_manual(values = SP_cols) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank()) + 
    facet_wrap("Type", ncol = 1)
  
 p <- p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(2,1))
  
}

figure_s4b_poportion {
  
  Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/Fibroblast_MetaUpdt.rds")
  Idents(Sobj) <- "Subcluster"
  Sobj$Subcluster <- Idents(Sobj)

  Prop <- 
    as.data.frame(prop.table(table(Sobj$Subcluster, 
                                   Sobj$Sample), 
                             margin = 2)) %>% 
    dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
    mutate(Type = substring(.$Sample, 1, 2)) %>% 
    left_join(distinct(select(Sobj@meta.data, Sample, Menopause)))
  
  
  
  
  ### by sample average split by menopause
    Prop %>% 
      ggplot(aes(x = factor(Cluster, levels = c("Matrix_1", "Matrix_2", "Lipo_F", "Vasc_F", "Chondrocyte")), 
                 y = Freq, 
                 fill = factor(Menopause, levels = c("NA", "pre",  "post")))) + 
    
    geom_boxplot(size = 1, width = 0.5, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0.5) +
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Fibroblast Proportions") + 
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
    
    
    ### by sample average split by Type
    Prop %>% 
      ggplot(aes(x = factor(Cluster, levels = c("Matrix_1", "Matrix_2", "Lipo_F", "Vasc_F", "Chondrocyte")), 
                 y = Freq, 
                 fill = Type)) + 
      
      geom_boxplot(size = 1, width = 0.5, 
                   position = position_dodge(width = 0.75),
                   color = "grey30",
                   outlier.alpha = 0.5) +
      
      scale_fill_manual (values = GID_cols) + 
      scale_color_manual(values = SP_cols) +
      ggtitle("Fibroblast Proportions") + 
      style + 
      theme(legend.position = "none", 
            axis.title.x = element_blank())
  
  
  ### by menstrual status summarized
  Prop <- 
    as.data.frame(prop.table(table(Sobj$Subcluster, 
                                   Sobj$Menopause), 
                             margin = 2)) %>% 
    dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq")
  
    p <- Prop %>% 
    ggplot(aes(x = factor(Cluster, levels = rev(c("Matrix_1", "Matrix_2", "Lipo_F", "Vasc_F", "Chondrocyte"))), 
               y = Freq, 
               fill = factor(Sample, levels = c("NA", "pre", "post")))) + 
    geom_col(position = position_fill()) +
    scale_fill_manual (values = Mens_cols) +
    style + coord_flip() +
    theme(legend.position = "bottom", 
          axis.title.y = element_blank())

    p %>% save_plot(name = "fibrblast_prop", scale = 1, w = 16, h = 9, svg = T)
    
    }

RPSA_plot {
  
p1 <- 
Perc_Avg %>% 
  filter(Gene == "RPSA") %>% 
  pivot_longer(cols = c(avg_TM, avg_CF), names_to = "avg_type", values_to = "average") %>% 
  pivot_longer(cols = c(perc_TM, perc_CF), names_to = "perc_type", values_to = "percent") %>% 
  ggplot(aes(x = CellType, y = average, fill = avg_type)) + 
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = GID_cols) +
  coord_flip() +
  style

p2 <- 
Perc_Avg %>% 
  filter(Gene == "RPSA") %>% 
  pivot_longer(cols = c(avg_TM, avg_CF), names_to = "avg_type", values_to = "average") %>% 
  pivot_longer(cols = c(perc_TM, perc_CF), names_to = "perc_type", values_to = "percent") %>% 
  ggplot(aes(x = CellType, y = percent, fill = perc_type)) + 
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = GID_cols) +
  coord_flip() +
  style

p1 + p2 + plot_layout(ncol = 1)

}

# Figure 5 ########################################

figure_5a {
  
  RNAdata  <- readRDS("2021_ReRun/data_resources/Plot_Data/Fig5a_Vascular_Umap.rds")
  
  Sub.levels <- c("Vein", "Capillary", "Artery",  "Lymph_EC",  "Lymph_EC2",  "Pericyte", "vasc_SM1",  "vasc_SM2")
  
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = factor(Subcluster, levels = Sub.levels))) + 
    geom_point(size = 2, alpha = 0.5, stroke = 0) + 
    scale_color_manual(values = pal_d3("category20")(20)[-4]) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank())
}

figure_5b {
  
  Prop  <-  readRDS("2021_ReRun/data_resources/Vascular_Proportions.rds")
  
  #data.subset <- Sobj@meta.data 
  #data.subset$Subcluster <- as.character(data.subset$Subcluster)
  
  #Prop <- 
  #  as.data.frame(prop.table(table(data.subset$Subcluster, 
  #                                 data.subset$Sample), 
  #                           margin = 2)) %>% 
  #  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  #  mutate(Type = substring(.$Sample, 1, 2))# %>% arrange(factor(Cluster, levels = levels), desc(Freq))
  
  Sub.levels <- c("Vein", "Capillary", "Artery")
  
  width = 0.5
  
  #p1 <- 
    Prop %>% 
    filter(Cluster %in% Sub.levels, Sample != "TM-9817") %>% 
    
    #ggplot(aes(x = factor(Cluster, levels = c("Vein", "Capillary", "Artery")), y = Freq, fill = Type)) + 
      ggplot(aes(x = Type, y = Freq, fill = Type)) + 
    
    geom_boxplot(size = 1, width = width, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    geom_point(position=position_jitterdodge(jitter.width = 0.1),
               aes(#shape = Type, 
                   color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
    
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Blood_EC") + 
    #scale_y_sqrt() +
    style + 
      facet_wrap("Cluster") +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  #p2 <- 
    Prop %>% 
    filter(Cluster %in% as.character(unique(filter(Sobj@meta.data, CellType == "Lymph_EC")$Subcluster)), Sample != "TM-9817") %>% 
    
    #ggplot(aes(x = factor(Cluster, levels = c("Lymph_EC", "Lymph_EC2")), y = Freq, fill = Type)) + 
      ggplot(aes(x = Type, y = Freq, fill = Type)) + 
      
    geom_boxplot(size = 1, width = width,
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    geom_point(position=position_jitterdodge(jitter.width = 0.1),
               aes(#shape = Type, 
                   color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
      stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
      
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Lymph_EC") + 
    #scale_y_sqrt() +
    style + 
      facet_wrap("Cluster") +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  p3 <-   
    Prop %>% 
    filter(Cluster %in% as.character(unique(filter(Sobj@meta.data, CellType == "Vasc.Acc.")$Subcluster)), Sample != "TM-9817") %>% 
    
    ggplot(aes(x = factor(Cluster, levels = c("Pericyte", "vasc_SM1", "vasc_SM2")), y = Freq, fill = Type)) + 
    
    geom_boxplot(size = 1, width = width,
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    geom_point(position=position_jitterdodge(jitter.width = 0.1),
               aes(#shape = Type, 
                   color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
    
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Vasc.Acc.") + 
    #scale_y_sqrt() +
    
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  
  p <- p1 + p2 + p3 + patchwork::plot_layout(ncol = 3, widths = c(3, 2, 3))
  
}  

figure_5b_Mens {
  
  Prop  <-  readRDS("2021_ReRun/data_resources/Vascular_Proportions.rds")
  
  #data.subset <- Sobj@meta.data 
  #data.subset$Subcluster <- as.character(data.subset$Subcluster)
  
  #Prop <- 
  #  as.data.frame(prop.table(table(data.subset$Subcluster, 
  #                                 data.subset$Sample), 
  #                           margin = 2)) %>% 
  #  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  #  mutate(Type = substring(.$Sample, 1, 2))# %>% arrange(factor(Cluster, levels = levels), desc(Freq))
  
  Sub.levels <- c("Vein", "Capillary", "Artery")
  
  width = 0.5
  
  p1 <- 
    Prop %>% left_join(select(HRT_data, Sample, Menopause)) %>% 
    filter(Cluster %in% Sub.levels, Sample != "TM-9817") %>% 
    
    ggplot(aes(x = factor(Cluster, levels = c("Vein", "Capillary", "Artery")), 
               y = Freq, 
               fill = factor(Menopause, levels = c("NA", "pre", "post"))
               )
           ) + 
    
    geom_boxplot(size = 1, width = width, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0.5) + 
    
    #geom_point(position=position_jitterdodge(jitter.width = 0.1),
    #           aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
    #           size = 4) + 
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Blood_EC") + 
    #scale_y_sqrt() +
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  p2 <- 
    Prop %>% left_join(select(HRT_data, Sample, Menopause)) %>% 
    filter(Cluster %in% as.character(unique(filter(Sobj@meta.data, CellType == "Lymph_EC")$Subcluster)), Sample != "TM-9817") %>% 
    
    ggplot(aes(x = factor(Cluster, levels = c("Lymph_EC", "Lymph_EC2")),
               y = Freq, 
               fill = factor(Menopause, levels = c("NA", "pre", "post"))
               )
           ) + 
    
    geom_boxplot(size = 1, width = width,
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0.5) + 
    
    #geom_point(position=position_jitterdodge(jitter.width = 0.1),
    #           aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
    #           size = 4) + 
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Lymph_EC") + 
    #scale_y_sqrt() +
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  p3 <-   
    Prop %>% left_join(select(HRT_data, Sample, Menopause)) %>% 
    filter(Cluster %in% as.character(unique(filter(Sobj@meta.data, CellType == "Vasc.Acc.")$Subcluster)), Sample != "TM-9817") %>% 
    
    ggplot(aes(x = factor(Cluster, levels = c("Pericyte", "vasc_SM1", "vasc_SM2")), 
               y = Freq, 
               fill = factor(Menopause, levels = c("NA", "pre", "post"))
               )
           ) + 
    
    geom_boxplot(size = 1, width = width,
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0.5) + 
    
    #geom_point(position=position_jitterdodge(jitter.width = 0.1),
    #           aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
    #           size = 4) + 
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Vasc.Acc.") + 
    #scale_y_sqrt() +
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  
  p1 + p2 + p3 + patchwork::plot_layout(ncol = 3, widths = c(3, 2, 3))
  
}  

figure_5cd {
  
  percentile_all <- 
    read_tsv(paste0("2021_ReRun/output/GRNboost2/mnt/expr_mat.adjacencies.tsv")) %>% 
    group_by(TF) %>% 
    mutate(percentile_rank = ntile(importance, 100))
  
  RNAdata  <- readRDS("2021_ReRun/data_resources/Plot_Data/Fig5c_Blood_EC_Umap.rds")
  
  Sub.levels <- c("Vein", "Capillary", "Artery")
  
  p1 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = factor(Subcluster, levels = Sub.levels))) + 
    geom_point(size = 2, alpha = 0.5, stroke = 0) + 
    scale_color_manual(values = pal_d3("category20")(20)[-4]) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank())
  
  
  query <- percentile_all %>% 
    filter(TF == "PPARG", percentile_rank > 95) %>% pull(target)
  
  Sobj <- AddModuleScore  (Sobj, assay = "RNA", 
                           list(query), 
                           nbin = 25, 
                           name = "PPARG_mod")
  
  
  
  Idents(Sobj) <- "Type"
  pt.size = 1
  
  p2 <- 
  FeaturePlot(subset(Sobj, idents = "CF", downsample = 7250), 
              features = c("PPARG_mod1"), 
              cols = viridis(n = 100, option = "A"), 
              pt.size = pt.size, 
              order = T) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank(), 
                  plot.title = element_blank())
  
  p3 <- 
  FeaturePlot(subset(Sobj, idents = "TM"), 
              features = c("PPARG_mod1"), 
              cols = viridis(n = 100, option = "A"), 
              pt.size = pt.size, 
              order = T) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank(), 
                  plot.title = element_blank())
  
  p4 <- p2 + p3 + plot_layout(nrow = 1)
  p1 + p4 + plot_layout(ncol = 2, widths = c(2,1))
  
  
  p <- percentile_all %>% 
    filter(TF == "PPARG", percentile_rank > 95, !str_detect(target, "-AS1")) %>% 
    top_n(5, importance) %>% 
    ggplot(aes(x = importance, y = reorder(target, importance))) + 
    geom_col(fill = "lightseagreen", 
             color = "grey30") + 
    style + theme(axis.title.y = element_blank())
  
  
}

figure_5e {
  
  Avg_RNA_Expression_Sample_CellType %>% 
    filter(Gene == "PPARG", CellType == "Blood_EC") %>% 
    ggplot(aes(x = Type, y = avg_RNA, fill = Type)) +
    geom_boxplot(size = 1, width = width, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    geom_point(position=position_jitterdodge(jitter.width = 0.1),
               aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
               size = 5) + 
    
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("PPARG") + 
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
}

figure_s5a{
  
  RNAdata  <- readRDS("2021_ReRun/data_resources/Plot_Data/FigS5a_TM9817_Outlier_Umap.rds")
  
  RNAdata[sample(nrow(RNAdata)),] %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2, 
               color = factor(Sample, levels = SP_levels))) + 
    
    geom_point(size = 2, alpha = 1, stroke = 0) + 
    
    scale_color_manual(values = SP_cols) + 
    style + theme(legend.position = "none",
                  axis.title = element_blank(), 
                  axis.text  = element_blank(), 
                  axis.ticks = element_blank())
}

figure_s5b {

library(readxl)
Kalucka_Conserved_Markers <- read_excel("2021_ReRun/utilities/Kalucka_Conserved_Markers.xlsx")

Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/Blood_EC_MetaUpdt.rds")

for (i in 1:3) {
  
  Sobj <- AddModuleScore  (Sobj, assay = "RNA", 
                           list(toupper(filter(Kalucka_Conserved_Markers, 
                                               Vessel == unique(Kalucka_Conserved_Markers$Vessel)[i])$gene)), 
                           nbin = 25, 
                           name = unique(Kalucka_Conserved_Markers$Vessel)[i])
}


pt.size = 0.5

p1 <- 
FeaturePlot(Sobj, features = c("Artery1"), 
            cols = viridis(n = 100, option = "A"), 
            pt.size = pt.size, 
            order = T) + 
  style + theme(legend.position = "none",
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank(), 
                plot.title = element_blank())
p2 <- 
FeaturePlot(Sobj, features = c("Capillary1"), 
            cols = viridis(n = 100, option = "A"), 
            pt.size = pt.size, 
            order = T) + 
  style + theme(legend.position = "none",
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank(), 
                plot.title = element_blank())
p3 <- 
FeaturePlot(Sobj, features = c("Vein1"), 
            cols = viridis(n = 100, option = "A"), 
            pt.size = pt.size, 
            order = T) + 
  style + theme(legend.position = "bottom",
                axis.title = element_blank(), 
                axis.text  = element_blank(), 
                axis.ticks = element_blank(), 
                plot.title = element_blank())

p <- p1 + p2 + p3 + patchwork::plot_layout()

}

figure_s5c {

Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/Vascular_V2.rds")

Idents(Sobj) <- "Subcluster"

query <- c("TLL1", "LNX1", "PCSK5", "RELN", "RADIL", "PDGFRB", "ACTA2", "THSD7B", "CLSTN2", "CACNB2", "CD36", "VWF", "ACTA2")
query <- Marks %>% filter(pct.1 > 0.5, !str_detect(gene, "-AS1|-AS2|\\.|LINC0")) %>% group_by(cluster) %>% top_n(3, avg_log2FC) %>% pull(gene)

mat <- AverageExpression(Sobj, features = query, assays = "RNA", group.by = "Subcluster")$RNA

pheatmap(t(mat)[levels(Sobj), query], 
         fontsize = 20, 
         scale = "column", 
         show_colnames = T,
         cluster_rows = F, 
         cluster_cols = F, 
         color = colorRampPalette(brewer.pal(11, "BrBG"))(100))

}

figure_s5f {

percentile_all <- 
  read_tsv(paste0("2021_ReRun/output/GRNboost2/mnt/expr_mat.adjacencies.tsv")) %>% 
  group_by(TF) %>% 
  mutate(percentile_rank = ntile(importance, 100))


library(gprofiler2)
library(viridis)

clust = "PPARG"

query <- percentile_all %>% filter(TF == clust, percentile_rank > 95) %>% pull(target)

Results <- gost(query, organism = "hsapiens", evcodes = F)

ggplot(head(filter(Results$result, source == "GO:BP"), 6), 
       aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
  geom_col() +
  scale_fill_viridis(begin = 0.3, end = 0.8, trans = "reverse") +
  theme(axis.title.y = element_blank(), 
        axis.text  = element_text(size = 18), 
        plot.title = element_text(size = 22, face = "bold")) +
  ylab(label = "Number of Genes in Pathway") +
  ggtitle(paste(unique(head(filter(Results$result, source == "GO:BP"), 15)$source), "-", clust, "Module")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  coord_flip() + style


}

Vasc_cell_marker_heatmaps {
  
  subs_vas <- Sobj@meta.data %>% filter(CellType %in% c("Blood_EC", "Lymph_EC", "Vasc.Acc.")) %>% pull(Subcluster) %>% unique()
  
  Subcluster_markers %>% filter(cluster == "Lymph_EC2") %>% arrange(-avg_log2FC)
  
  genes <- c("TLL1", "PLCXD3", "ICAM1", "VWF",
             "LNX1", "CD36", "KDR",
             "PCSK5", "NKAIN2", "MECOM",
             #"RELN", 
             #"NRP2", 
             "LYVE1", "PDPN", "FLT4",
             #"RADIL", 
             "SV2C",
             "PDGFRB", "THSD7B", 
             "ACTA2", "CLSTN2", "CACNB2")
  
  
  avg <- AverageExpression(Sobj, features = c(genes), 
                           group.by = "Subcluster")
  
  mat1 <- avg$RNA %>% as.data.frame() 
  
  #breaksList = seq(-1.75, 1.75, by = 0.01)
  
  p <- 
    pheatmap(t(mat1[genes, c("Vein", "Capillary", "Artery", "Lymph_EC", "Lymph_EC2", "Pericyte", "vasc_SM1", "vasc_SM2")]), 
             cluster_rows = F,
             cluster_cols = F,
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation", 
             clustering_method = "ward.D2",
             scale = "column", fontsize = 20, 
             #breaks = breaksList,
             color = colorRampPalette(brewer.pal(11, "BrBG"))(100)
    )
  
  p %>% save_plot(name = "fig_5_vasc_markers", scale = 1, w = 16, h = 9, svg = T)
  
  
  
  
  
  
}

# Figure 6 ########################################

figure_6a {

Sobj <- AddMetaData(Sobj, as.data.frame(Sobj@reductions$umap@cell.embeddings))
RNAdata  <- Sobj@meta.data

#p1 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, 
             color = Subcluster)) + 
  
  geom_point(size = 3, alpha = 0.5, stroke = 0) +
  scale_color_manual(values = pal_d3("category20")(20)) + 
  style + 
  theme(legend.position = "top", 
        panel.background = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.title = element_blank(), 
        axis.text  = element_blank(), 
        axis.ticks = element_blank())

}

figure_6b {
  
  
  Prop  <-  readRDS("2021_ReRun/data_resources/Immune_Proportions.rds")
  
  
  Lymphoid <- c("CD4_T", "CD8_T", "T-Eff", "NK")
  Myeloid  <- c("mono.DC", "Macrophage", "Monocyte", "DC")
  Minor  <- c("B.mzone", "B.plasma", "HSC")
  
  subplot <- Lymphoid
  
  p3 <- 
    Prop %>% filter(Cluster %in% c(Minor)) %>%
    ggplot(aes(x = Type, y = Freq, fill = Type)) + 
    geom_boxplot(size = 1, width = 0.75, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    geom_point(position=position_jitter(width = 0.1),
               aes(#shape = Type, 
                   color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
    stat_summary(fun = "median", 
                 geom = "point", size = 2,
                 position = position_dodge(0.75),
                 color = "floralwhite") +
    
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) +
    #scale_y_sqrt() +
    style + 
    facet_wrap("Cluster", scales = "free_x", nrow = 1) +
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  
  save_plot(p, name = "Immune_Proportions", scale = 1, w = 16, h = 10, svg = T)
  
  
  Prop %>% filter(Cluster %in% Myeloid) %>% 
    left_join(select(HRT_data, Sample, Menopause)) %>% 
    ggplot(aes(x = factor(Cluster, levels = Myeloid), 
               y = Freq, 
               fill = factor(Menopause, levels = c("NA", "pre", "post")))
           ) + 
    geom_boxplot(size = 1, width = 0.5, 
                 position = position_dodge(width = 0.85),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    #geom_point(position=position_jitterdodge(jitter.width = 0.1),
    #           aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
    #           size = 4) + 
    stat_summary(fun = "median", 
                 geom = "point", size = 2,
                 position = position_dodge(0.85),
                 color = "floralwhite") +
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Myeloid Cells") + scale_y_sqrt() +
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  
  Prop %>% filter(Cluster %in% Lymphoid) %>% 
    left_join(select(HRT_data, Sample, Menopause)) %>% 
    ggplot(aes(x = factor(Cluster, levels = Lymphoid), 
               y = Freq, 
               fill = factor(Menopause, levels = c("NA", "pre", "post")))
    ) + 
    geom_boxplot(size = 1, width = 0.65, 
                 position = position_dodge(width = 0.85),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    #geom_point(position=position_jitterdodge(jitter.width = 0.1),
    #           aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
    #           size = 4) + 
    stat_summary(fun = "median", 
                 geom = "point", size = 2,
                 position = position_dodge(0.85),
                 color = "floralwhite") +
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle("Lymphoid Cells") + scale_y_sqrt() +
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
}

figure_6b_proportion {
  
  celltype = c("Myeloid", "Lymphoid")
  
  subset <- filter(Sobj@meta.data, CellType %in% celltype)
  subset$Subcluster <- as.character(subset$Subcluster)
  
  Prop <- 
    as.data.frame(prop.table(table(subset$Subcluster, 
                                   subset$Sample), 
                             margin = 2)) %>% 
    dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
    mutate(Type = substring(.$Sample, 1, 2)) %>% 
    left_join(distinct(select(subset, Sample, Menopause))) %>% 
    mutate(Group = ifelse(Cluster %in% c("CD4_T", "CD8_T", "NK", "T-Eff"), "Lymphoid", "Minor")) %>% 
    mutate(Group = ifelse(Cluster %in% c("Macrophage", "mono.DC", "Monocyte", "DC"), "Myeloid", Group)) %>% 
    mutate(Group = factor(Group, levels = c("Lymphoid", "Myeloid", "Minor")))
  
  
  ### standart plot without Mens Split but with colored dots
  Prop %>% 
    ggplot(aes(x = Cluster, y = Freq, fill = Type)) + 
    geom_boxplot(size = 1, width = 0.5, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0) + 
    
    geom_point(position=position_jitterdodge(jitter.width = 0.1),
               aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
               size = 4) + 
    stat_summary(fun = "median", 
                 geom = "point", size = 2,
                 position = position_dodge(0.7),
                 color = "floralwhite") +
    
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle(paste0(celltype)) + 
    style + 
    facet_wrap("Group", scales = "free") +
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  
  ### standart plot with Mens Split
  Prop %>% 
    ggplot(aes(x = Cluster, 
               y = Freq, 
               fill = factor(Menopause, levels = c("NA", "pre", "post"))
    )
    ) + 
    geom_boxplot(size = 1, width = 0.65, 
                 position = position_dodge(width = 0.75),
                 color = "grey30",
                 outlier.alpha = 0.5) + 
    
    #geom_point(position=position_jitterdodge(jitter.width = 0.1),
    #           aes(shape = Type, color = factor(Sample, levels = SP_levels)), 
    #           size = 4) + 
    stat_summary(fun = "median", 
                 geom = "point", size = 2,
                 position = position_dodge(0.75),
                 color = "floralwhite") +
    
    scale_fill_manual (values = Mens_cols) + 
    scale_color_manual(values = SP_cols) +
    ggtitle(paste0(celltype)) + 
    style + 
    facet_wrap("Group", scales = "free") + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  
  ### stacked bar plot
  Prop <- 
    as.data.frame(prop.table(table(subset$Subcluster, 
                                   subset$Menopause), 
                             margin = 2)) %>% 
    dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq")
  
  #p1 <- 
  Prop %>% 
    ggplot(aes(x = Cluster, y = Freq, fill = Sample)) + 
    geom_col(position = position_fill()) +
    scale_fill_manual (values = Mens_cols) +
    ggtitle(paste0(celltype)) + 
    style + 
    theme(legend.position = "none", 
          axis.title.x = element_blank())
  
  #p1 %>% save_plot(name = "fig_2_proportions", scale = 1, w = 16, h = 9, svg = T)
  
}

Immune_cell_marker_heatmaps {
  
  subs_mye <- Sobj@meta.data %>% filter(CellType %in% c("Myeloid")) %>% pull(Subcluster) %>% unique()
  subs_lym <- Sobj@meta.data %>% filter(CellType %in% c("Lymphoid")) %>% pull(Subcluster) %>% unique()
  
  
  
  Subcluster_markers %>% filter(cluster == "mono.DC") %>% arrange(-avg_log2FC)
  
  
  
  lym_genes <- c("TOX", "DTHD1",
                 "CD4", "IL7R", "CD3E", "LEF1",
                 "CD8A", "LAG3",
                 "NCAM1", "GNLY", "GZMB",
                 "MS4A1", "BANK1",
                 "MZB1", "BMP6")
  
  
  avg <- AverageExpression(Sobj, features = c(lym_genes), 
                           group.by = "Subcluster")
  
  mat1 <- avg$RNA %>% as.data.frame() 
  
  breaksList = seq(-1.75, 1.75, by = 0.01)
  
  p <- 
    pheatmap(t(mat1[,as.character(subs_lym)]), 
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation", 
             clustering_method = "ward.D2",
             scale = "column", fontsize = 20, breaks = breaksList,
             color = colorRampPalette(brewer.pal(11, "BrBG"))(length(breaksList))
    )
  
  p %>% save_plot(name = "fig_6a_lym_markers", scale = 1, w = 16, h = 9, svg = T)
  
  
  
  
  mye_genes <- c("MRC1", "CD163", 
                 "ITGAX", "ITGAM", "CD14",
                 "C3", "PALD1", 
                 "HLA-DRB1", 
                 "CD36", "SNTB1",
                 "CD44", "KIT")
  
  avg <- AverageExpression(Sobj, features = c(mye_genes), 
                           group.by = "Subcluster")
  
  
  
  mat1 <- avg$RNA %>% as.data.frame() 
  
  
  # Sets the minimum (0), the maximum (15), and the increasing steps (+1) for the color scale
  # Note: if some of your genes are outside of this range, they will appear white on the heatmap
  breaksList = seq(-1.75, 1.75, by = 0.01)
  
  p <- 
    pheatmap(t(mat1[,as.character(subs_mye)]), 
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation", 
             clustering_method = "ward.D2", 
             breaks = breaksList,
             scale = "column", fontsize = 20, 
             color = colorRampPalette(brewer.pal(11, "BrBG"))(length(breaksList))
    )
  
  p %>% save_plot(name = "fig_6a_mye_markers", scale = 1, w = 16, h = 9, svg = T)
  
}

# Figure 7 ########################################

adipocyte_UMAP {

Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/Adipocyte_MetaUpdt.rds")
Sobj <- AddMetaData(Sobj, as.data.frame(Sobj@reductions$umap@cell.embeddings))

RNAdata  <- Sobj@meta.data

#p1 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, 
             color = Sample)) + 
  
  geom_point(size = 3, alpha = 0.5, stroke = 0) +
  scale_color_manual(values = pal_d3("category20")(20)) + 
  style + 
  theme(legend.position = "top", 
        panel.background = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.title = element_blank(), 
        axis.text  = element_blank(), 
        axis.ticks = element_blank())

p2 <- 
  RNAdata[sample(nrow(RNAdata)),] %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, 
             color = factor(Sample, levels = SP_levels))) + 
  
  geom_point(data = select(RNAdata, UMAP_1, UMAP_2), 
             color = "grey90", size = 4, alpha = 1, stroke = 0) + 
  geom_point(size = 2, alpha = 1, stroke = 0) +
  scale_color_manual(values = SP_cols) + 
  style + 
  facet_wrap("Type", ncol = 1) +
  
  theme(legend.position = "none", 
        strip.background = element_blank(), 
        panel.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(), 
        axis.text  = element_blank(), 
        axis.ticks = element_blank())

}

pi3k_TFs {

gProfiler_KEGG_PI3K_AKT <- read_csv("utilities/Previous_Analysis/gProfiler_KEGG_PI3K-AKT.csv")

mat <- 
DE_Response %>% 
  filter(Gene %in% gProfiler_KEGG_PI3K_AKT$name, TF == "YES", p_val_adj <= 0.05) %>% 
  group_by(Gene) %>% 
  mutate(n = n()) %>% 
  filter(n > 5) %>% 
  select(Gene, avg_log2FC, Cluster) %>% 
  pivot_wider(names_from = Cluster, 
              values_from = avg_log2FC) %>% 
  column_to_rownames("Gene")

mat[is.na(mat)] <- 0

p <- 
pheatmap(mat, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D2",
         scale = "column", 
         fontsize = 20, 
         color = colorRampPalette(brewer.pal(11, "BrBG"))(100)
         )

p %>% save_plot(name = "fig_7b_pi3k_TFs", scale = 1, w = 16, h = 9, svg = T)

}

Lipolysis_Markers <- read_csv("utilities/Previous_Analysis/Lipolysis_Markers.csv")

Lipolysis_Markers %>% 
  left_join(DE_Response) %>% 
  filter(Cluster == "Adipocyte") %>% 
  select(Gene, LipoFunction, avg_log2FC, p_val_adj) %>% 
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = LipoFunction)) + 
  geom_vline(xintercept = 0, color = "grey80") +
  geom_point() + 
  style + 
  theme(legend.position = "bottom")

pi3k_enrichment {


library(gprofiler2)
library(viridis)

gProfiler_KEGG_PI3K_AKT <- read_csv("utilities/Previous_Analysis/gProfiler_KEGG_PI3K-AKT.csv")
  
query <- 
  DE_Response %>% 
  filter(#Gene %in% gProfiler_KEGG_PI3K_AKT$name, 
         p_val_adj <= 0.05, abs(avg_log2FC) > 0.25,
         Receptor == "YES" | Ligand == "YES") %>% 
  pull(Gene) %>% 
  unique()

clust = "PI3K Signaling"


Results <- gost(query, organism = "hsapiens", evcodes = F)


data <- head(filter(Results$result, source == "KEGG"), 6) %>% 
  mutate(r = as.character(p_value)) %>% separate(r, into = c("A", "B"), sep = "e") %>% 
  mutate(A = round(as.numeric(A), 2)) %>% unite(p_value, "A", "B", sep = "e") %>% 
  mutate(p_value = as.numeric(p_value))


ggplot(data, 
       aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
  geom_col() +
  scale_fill_viridis(begin = 0.2, end = 0.8, trans = "reverse") +
  theme(axis.title.y = element_blank(), 
        axis.text  = element_text(size = 18), 
        plot.title = element_text(size = 22, face = "bold")) +
  ylab(label = "Number of Genes in Pathway") +
  ggtitle(paste(unique(head(filter(Results$result, source == "KEGG"), 15)$source), "-", clust, "Module")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  coord_flip() + 
  style




Results <- gost(query, organism = "hsapiens", evcodes = T)

data <- head(filter(Results$result, source == "GO:BP"), 10) %>% 
  mutate(r = as.character(p_value)) %>% #separate(r, into = c("A", "B"), sep = "e") %>% 
  #mutate(A = round(as.numeric(A), 2)) %>% unite(p_value, "A", "B", sep = "e") %>% 
  mutate(p_value = as.numeric(r))

p <- ggplot(data, 
            aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
  geom_col() +
  scale_fill_viridis(begin = 0.2, end = 0.8, trans = "reverse") +
  theme(axis.title.y = element_blank(), 
        axis.text  = element_text(size = 18), 
        plot.title = element_text(size = 22, face = "bold")) +
  ylab(label = "Number of Genes in Pathway") +
  ggtitle(paste(unique(head(filter(Results$result, source == "GO:BP"), 15)$source), "-", clust, "Module")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  coord_flip() + 
  style


}

figure_s7f_gpam {
  
  percentile_all <- 
    read_tsv(paste0("2021_ReRun/output/GRNboost2/mnt/expr_mat.adjacencies.tsv")) %>% 
    group_by(TF) %>% 
    mutate(percentile_rank = ntile(importance, 100))
  
  
  library(gprofiler2)
  library(viridis)
  
  clust = "GPAM"
  
  query <- percentile_all %>% filter(TF == clust, percentile_rank > 95) %>% pull(target)
  
  Results <- gost(query, organism = "hsapiens", evcodes = T)
  
  ggplot(head(filter(Results$result, source == "WP"), 8), 
         aes(fct_reorder(term_name, p_value, .desc = T), intersection_size, fill = p_value)) + 
    geom_col() +
    scale_fill_viridis(begin = 0.3, end = 0.8, trans = "reverse") +
    theme(axis.title.y = element_blank(), 
          axis.text  = element_text(size = 18), 
          plot.title = element_text(size = 22, face = "bold")) +
    ylab(label = "Number of Genes in Pathway") +
    ggtitle(paste(unique(head(filter(Results$result, source == "WP"), 15)$source), "-", clust, "Module")) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
    coord_flip() + style
  
  
}

figure_7_adipocyte {

p <- VlnPlot(Sobj, idents = "Adipocyte", features = c("PTEN", "PIK3R1", "AKT3"), split.by = "Type", pt.size = 0, cols = GID_cols, ncol = 1)

data <- Sobj@assays$RNA@data[c("TCF7L2", "ADCY2"), ] %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% left_join(rownames_to_column(Sobj@meta.data, "ID")) %>% filter(CellType == "Adipocyte") %>% select(CellType, Type, Sample, TCF7L2, gpam1, Menopause) #ggplot(aes(x = Type, y = TCF7L2, fill = Type)) + geom_violin()

p1 <- 
data %>% 
  ggplot(aes(x = Type, y = TCF7L2, fill = Type)) + 
  geom_violin() + 
  scale_fill_manual (values = GID_cols) +
  style + 
  labs(title = "TCF7L2", x = NULL) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

#p2 <- 
Sobj@meta.data %>% 
  ggplot(aes(x = Type, y = gpam1, fill = Type)) + 
  geom_violin() + 
  scale_fill_manual (values = GID_cols) +
  style + 
  labs(title = "GPAM Mpdule", x = NULL) +
theme(legend.position = "none", 
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 


p2 + p1


}


Sobj <- AddModuleScore(Sobj, assay = "RNA", 
                       list(query), 
                       nbin = 25, 
                       name = "gpam")



Sobj <- AddModuleScore(Sobj, assay = "RNA", 
                       list(filter(Searchable_gmt, Term == terms[15])$Gene), 
                       nbin = 25, 
                       name = "test")


p <- Sobj@meta.data %>% 
  select(Type, CellType, Sample, test1) %>% 
  ggplot(aes(x = CellType, y = test1, fill = Type)) + 
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
  labs(y = "insulin signaling score", x = NULL) +
  style + 
  theme(legend.position = "none", text = element_text(size = 35),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))




figure_GTEX_malefemale {
  
  make_GTEX_Breast_Seurat_Object {
    
    
    # load gtex data and annotations
    dat.gct <- read.delim(file="~/Downloads/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip=2)
    
    GTEx_Analysis_v8_Annotations_SampleAttributesDS  <- read_delim("~/Downloads/SexBias/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", "\t",  escape_double = FALSE, trim_ws = TRUE)
    GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS <- read_delim("~/Downloads/SexBias/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    
    # pull the gtex gene identifiers
    genes <- dat.gct %>% select(1,2) %>% group_by(Description) %>% summarise(n = n()) %>% filter(n == 1) %>% pull(Description)
    
    # create SUBJID from SAMPID
    translate <- 
      GTEx_Analysis_v8_Annotations_SampleAttributesDS %>% 
      select(SAMPID) %>% separate(SAMPID, into = c("A", "B", "C"), sep = "-", extra = "merge", remove = F) %>% 
      unite(SUBJID, A, B, sep = "-") %>% select(-C)
    
    colnames(dat.gct) <- str_replace_all(colnames(dat.gct), "\\.", "-")
    breast_samples <- GTEx_Analysis_v8_Annotations_SampleAttributesDS %>% filter(SMTS == "Breast") %>% pull(SAMPID)
    
    # subset breast samples from GTEX data
    gtex_breast <- dat.gct %>% 
      filter(Description %in% genes) %>% 
      select(Description, intersect(breast_samples, colnames(dat.gct)))
    
    #gtex_breast %>% saveRDS("Output/Previous_Analysis/GTEx_breast_samples.rds")
    
    
    # load GTEX breast data and make Seurat Object
    gtex_breast  <- readRDS("Output/Previous_Analysis/GTEx_breast_samples.rds")
    
    mat <- column_to_rownames(gtex_breast, "Description")
    
    Sobj <- CreateSeuratObject(mat)
    Sobj <- NormalizeData        (Sobj, verbose = T)
    Sobj <- FindVariableFeatures (Sobj, selection.method = "vst", nfeatures = 2000, verbose = T)
    Sobj <- ScaleData     (Sobj, verbose = T)
    Sobj <- RunPCA        (Sobj, npcs = 50, verbose = FALSE)
    Sobj <- RunUMAP       (Sobj, reduction  = "pca", dims = 1:50)
    Sobj <- FindNeighbors (Sobj, reduction  = "pca", dims = 1:50)
    Sobj <- FindClusters  (Sobj, resolution = c(0.2, 0.5, 1))
    
    
    # create metadata
    meta <- 
      translate %>% 
      filter(SAMPID %in% colnames(Sobj)) %>% 
      left_join(GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS) %>% 
      mutate(Gender = ifelse(SEX == "1", "male", "female")) %>% 
      column_to_rownames("SAMPID")
    
    Sobj <- AddMetaData(Sobj, meta)
    
    
    # score data for group and celltype markers from nuclei
    my_markers    <- readRDS("Output/Previous_Analysis/Group_Markers.rds")
    my_ct_markers <- readRDS("Output/Previous_Analysis/CellType_Markers.rds")
    
    
    for (i in 1:length(unique(my_markers$cluster))) {
      
      module <- my_markers    %>% filter(cluster == unique(my_markers$cluster)[i]) %>% top_n(100, avg_log2FC) %>% pull(gene)
      Sobj <- AddModuleScore  (Sobj, assay = "RNA", list(module), nbin = 25, name = paste0(unique(my_markers$cluster)[i], "_"))
      
    }
    
    for (i in 1:length(unique(my_ct_markers$cluster))) {
      
      module <- my_ct_markers    %>% filter(cluster == unique(my_ct_markers$cluster)[i]) %>% top_n(100, avg_log2FC) %>% pull(gene)
      Sobj <- AddModuleScore  (Sobj, assay = "RNA", list(module), nbin = 25, name = paste0(unique(my_ct_markers$cluster)[i], "_"))
      
    }
    
    # select fibroblast enriched samples (two cutoffs, 75th and 50th percentile)
    fibros_meta <- Sobj@meta.data %>% rownames_to_column("ID") %>% select(ID, Fibroblast_1) %>% mutate(FibroRank = percent_rank(Fibroblast_1)) %>% mutate(Fibroblast_top25 = ifelse(FibroRank > 0.75, "Yes", "No"), Fibroblast_top50 = ifelse(FibroRank > 0.5, "Yes", "No")) %>% select(ID, Fibroblast_top25, Fibroblast_top50) %>% column_to_rownames("ID")
    Sobj <- AddMetaData(Sobj, fibros_meta)
    Sobj@meta.data %>% rownames_to_column("SAMPLEID") %>% select(SAMPLEID, SUBJID, AGE, Gender, Group, Fibroblast_top25, Fibroblast_top50, seurat_clusters) %>% saveRDS("Output/Previous_Analysis/GTEX_breast_metadata.rds")
    
    
    Idents(Sobj) <- "seurat_clusters"
    
    Sobj <- RenameIdents(Sobj, 
                         "0" = "Epithelial",
                         "1" = "Adipose",
                         "2" = "Vasculature",
                         "3" = "Epithelial")
    
    Sobj$Group <- Idents(Sobj)
    
    Sobj %>% saveRDS("2021_ReRun/output/GTEx_Breast_Sobj_Processed.rds")
    
    
  }
  
  plot_gtex_breast {
    
    Sobj <- readRDS("2021_ReRun/output/GTEx_Breast_Sobj_Processed.rds")
    totalGeneMappedReads_bySample <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/totalGeneMappedReads_bySample.rds")
    
    
    ## processing and characterization of the data
    p <- DimPlot(Sobj, group.by = "Group",  pt.size = 2, label = F, cols = pal_d3("category20")(20)) + theme(plot.margin=unit(c(0, 0, 0, 0), "lines"), legend.position = "bottom", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    q <- DimPlot(Sobj, group.by = "Gender", pt.size = 2, cols = GID_cols, label = F) +                 theme(plot.margin=unit(c(0, 0, 0, 0), "lines"), legend.position = "bottom", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    
    layout <- theme(plot.margin=unit(c(1, 0, 0, 0), "lines"), legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    
    p1 <- FeaturePlot(Sobj, features = c("Epithelial_1"),  cols = viridis(n = 100, option = "A"), pt.size = 1, order = T) + labs(title = "Epithelial") +  layout
    p2 <- FeaturePlot(Sobj, features = c("Vasculature_1"), cols = viridis(n = 100, option = "A"), pt.size = 1, order = T) + labs(title = "Vasculature") + layout
    p3 <- FeaturePlot(Sobj, features = c("Adipocyte_1"),   cols = viridis(n = 100, option = "A"), pt.size = 1, order = T) + labs(title = "Adipose") +     layout
    p4 <- FeaturePlot(Sobj, features = c("Fibroblast_1"),  cols = viridis(n = 100, option = "A"), pt.size = 1, order = T) + labs(title = "Fibroblast") +  layout
    
    # more vertical version
    x <- p1 + p2 +  patchwork::plot_layout(ncol = 2)
    y <- p3 + p4 +  patchwork::plot_layout(ncol = 2)
    p <- p + q + x + y + patchwork::plot_layout(ncol = 2, nrow = 2, heights = c(2, 1))
    
    # more horizontal version
    x <- p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 4)
    p + q + x + patchwork::plot_layout(ncol = 3, nrow = 1, widths = c(1, 1, 1))
    
  }
  
  plot_gtext_cohort {
    
  cohort <- 
  Sobj@meta.data %>% 
    select(Gender, AGE, DTHHRDY, Group, nCount_RNA, nFeature_RNA) %>% 
    group_by(Gender, AGE) %>% 
    summarise(n = n()) %>% 
    ggplot(aes(x = Gender, y = n, fill = factor(AGE, levels = rev(sort(unique(Sobj@meta.data$AGE)))))) + 
    geom_col(position = position_stack()) + 
    scale_fill_viridis_d(direction = 1, end = 0.9) +
    facet_wrap("Gender", scales = "free") +
    labs(title = "GTEx Breast Sample Cohort", fill = "Age") +
    style + theme(legend.position = "bottom", 
                  axis.title = element_blank(), 
                  axis.text.x = element_blank(), 
                  axis.ticks = element_blank())
  
  
  my_comparisons <- list(c("male", "female"))
  
  nGenes <- 
  Sobj@meta.data %>% 
    rownames_to_column("Sample") %>% 
    left_join(totalGeneMappedReads_bySample) %>% 
    select(Gender, AGE, DTHHRDY, Group, nCount_RNA, nFeature_RNA, TotalGeneMappedReads) %>% 
    ggplot(aes(x = Gender, y = (nFeature_RNA/1000), fill = Gender)) + 
    geom_violin() +
    geom_boxplot(aes(color = Gender), fill = "grey30", 
                 width = 0.2, 
                 size = 0.5) + 
    stat_summary(fun = "median", 
                 geom = "point", size = 2,
                 position = position_dodge(0.9),
                 color = "floralwhite") +
    stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    labs(title = "GTEx Genes Detected", fill = "Age") +
    scale_fill_manual (values = GID_cols) + 
    scale_color_manual(values = c("grey30", "grey30")) +
    style + theme(legend.position = "none")
  
  
  A <- p1 + p2
  B <- p + q
  
  p <- A + B + plot_layout(ncol = 3, widths = c(1,1,4))
  p %>% save_plot(name = "fig_GTEx_Upper", scale = 1, w = 16, h = 9, svg = T)
  x %>% save_plot(name = "fig_GTEx_middle", scale = 1, w = 16, h = 9, svg = T)
  
  }
  
  DE_testing_wilcox {
    Idents(Sobj) <- "Group"
    
    epithelial <- subset(Sobj, idents = c("Epithelial"))
    Idents(epithelial) <- "Gender"
    response_epithelial  <- FindMarkers(epithelial, ident.1 = "male", assay = "RNA", slot = "counts", logfc.threshold = 0.1, test.use = "wilcox") %>% rownames_to_column("Gene")
    
    epithelial <- subset(Sobj, idents = c("Adipose"))
    Idents(epithelial) <- "Gender"
    response_adipose     <- FindMarkers(epithelial, ident.1 = "male", assay = "RNA", slot = "counts", logfc.threshold = 0.1, test.use = "wilcox") %>% rownames_to_column("Gene")
    
    epithelial <- subset(Sobj, idents = c("Vasculature"))
    Idents(epithelial) <- "Gender"
    response_vasculature <- FindMarkers(epithelial, ident.1 = "male", assay = "RNA", slot = "counts", logfc.threshold = 0.1, test.use = "wilcox") %>% rownames_to_column("Gene")
    
    epithelial <- subset(Sobj, cells  =  rownames(filter(Sobj@meta.data, Fibroblast_top25 == "Yes")))
    Idents(epithelial) <- "Gender"
    response_fibros25      <- FindMarkers(epithelial, ident.1 = "male", assay = "RNA", slot = "counts", logfc.threshold = 0.1, test.use = "wilcox") %>% rownames_to_column("Gene")
    
    epithelial <- subset(Sobj, cells  =  rownames(filter(Sobj@meta.data, Fibroblast_top50 == "Yes")))
    Idents(epithelial) <- "Gender"
    response_fibros50      <- FindMarkers(epithelial, ident.1 = "male", assay = "RNA", slot = "counts", logfc.threshold = 0.1, test.use = "wilcox") %>% rownames_to_column("Gene")
    
  }
  DE_testing_Deseq2 {
    response_epithelial <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/Submission/Figures/Mehran_PDF/fig_gtexMaleBreast/data/Epithelial_deseq2_GenderDE.RDS") %>% 
      as_tibble() %>% select(Gene = "GeneSymbol", avg_log2FC = log2FoldChange, p_val_adj = padj)
    
    response_adipose <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/Submission/Figures/Mehran_PDF/fig_gtexMaleBreast/data/Adipose_deseq2_GenderDE.RDS") %>% 
      as_tibble() %>% select(Gene = "GeneSymbol", avg_log2FC = log2FoldChange, p_val_adj = padj)
    
    response_vasculature <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/Submission/Figures/Mehran_PDF/fig_gtexMaleBreast/data/Vasculature_deseq2_GenderDE.RDS") %>% 
      as_tibble() %>% select(Gene = "GeneSymbol", avg_log2FC = log2FoldChange, p_val_adj = padj)
    
    response_fibor25 <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/Submission/Figures/Mehran_PDF/fig_gtexMaleBreast/data/Fibroblast_top25_deseq2_GenderDE.RDS") %>% 
      as_tibble() %>% select(Gene = "GeneSymbol", avg_log2FC = log2FoldChange, p_val_adj = padj)
    
    response_fibor50 <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/Submission/Figures/Mehran_PDF/fig_gtexMaleBreast/data/Fibroblast_top50_deseq2_GenderDE.RDS") %>% 
      as_tibble() %>% select(Gene = "GeneSymbol", avg_log2FC = log2FoldChange, p_val_adj = padj)
    
    response_All <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/Submission/Figures/Mehran_PDF/fig_gtexMaleBreast/data/AllGTExBreasts_deseq2_GenderDE.RDS") %>% 
      as_tibble() %>% select(Gene = "GeneSymbol", avg_log2FC = log2FoldChange, p_val_adj = padj)
  }
  
  
  plot_epithelial {
    
    #labels1 <- c("ACTA2", "AZGP1", "KLF6", "NRP1", "EGFR", "CTNNB1", "TEAD1", "OXTR", "ITGA2", "ITGB8", "ENAH", "NR4A1", "FOXO3", "BACH2")  
    #labels2 <- c("AZGP1", "CUX2", "RYR2", "NR4A1", "FOXO3", "ADCY2", "AREG", "PGR", "EREG", "ACACA")  
    #labels1 <- c("ACTA2", "AZGP1", "KLF6", "NRP1", "EGFR", "CTNNB1", "TEAD1", "OXTR", "ITGA2", "ITGB8", "ENAH", "NR4A1", "FOXO3", "BACH2")  
    lab_lup <- c("AZGP1",  "CUX2", "RYR2", "RXRA", "NR4A1", "FOXO3", "PROS1", "ADCY2", "AREG", "PGR", "ESR1", "EREG", "EGFR", "ADAM17", "ITGB1", "ENAH", "ITGB8", "ITGA2", "ESRRG", "BATF", "BACH2", "TP63", "ACTA2", "OXTR")  
    lab_lup <- c("ESRRG", "CREB5", "EGFR", "ITGA2", "ITGB8", "ENAH", "NR4A1", "FOXO3")  
    lab_lup <- c("ACTA2", "NRP1", "EGFR", "TEAD1", "OXTR", "FOXO3", "BACH2", "PROS1", "ADAM17")  
    
    celltypes <- c("LUM_HR-neg", "Basal", "LUM_HR-pos")
    labels <- c(labels1, labels2)
    
    data <- 
      DE_Response %>% 
      filter(Cluster %in% celltypes, abs(avg_log2FC) > 0.2, p_val_adj <= 0.05) %>%
      select(Gene, Cluster, "logFC_AR" = avg_log2FC) %>% 
      inner_join(filter(response_epithelial), by = "Gene") %>% 
      mutate(color = ifelse(p_val_adj <= 0.05, "p-adj < 0.05", "n.s")) %>% 
      select(Gene, Cluster, color, logFC_AR, "logFC_GTEx" = avg_log2FC)
    
    #p <- 
    data %>% 
      filter(abs(logFC_GTEx) < 4, !is.na(color)) %>% 
      ggplot(aes(x = logFC_AR, y = logFC_GTEx, color = color)) + 
      geom_point(alpha = 0.75) + 
      stat_cor(method = "pearson", size = 5) + 
      geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
      geom_vline(xintercept = 0, color = "grey50", linetype = "dotted") +
      geom_smooth(data = filter(data, color != "n.s"), 
                  method = lm, 
                  se = FALSE, 
                  color = "turquoise4") + 
      geom_label_repel(data = filter(data, Gene %in% lab_lup), 
                       aes(label = Gene), size = 5,
                       min.segment.length = 0) + 
      labs(title = c("Cor. GTEx Male vs. Female and AR-Treatment"), color = "sig. in GTEx") +
      xlab("logFC AR-treatment") + 
      ylab("logFC Male vs. Female in GTEx epithlial breast") +
      #ylim(c(-2.2, 4)) +
      scale_color_manual(values = c("black", "#A6499B")) +
      facet_wrap(facets = "Cluster", ncol =  3, scales = "free") +
      style +
      theme(legend.position = "bottom")
    
  }
  
  plot_adipose {
    
    labels <- c("TCF7L2", "PTEN", "NR4A1", "FOXO3", "ADCY2", "GPAM", "INSR", "THBS1", "VEGFA", "GREB1L")
    celltypes <- c("Adipocyte")
    
    data <- 
      DE_Response %>% 
      filter(Cluster %in% celltypes, abs(avg_log2FC) > 0.3, p_val_adj <= 0.05) %>%
      select(Gene, Cluster, "logFC_AR" = avg_log2FC) %>% 
      inner_join(filter(response_adipose), by = "Gene") %>% 
      mutate(color = ifelse(p_val_adj <= 0.05, "p-adj < 0.05", "n.s")) %>% 
      select(Gene, Cluster, color, logFC_AR, "logFC_GTEx" = avg_log2FC)
    
    #p1 <- 
      data %>% 
      filter(abs(logFC_GTEx) < 4, !is.na(color)) %>% 
      ggplot(aes(x = logFC_AR, y = logFC_GTEx, color = color)) + 
      geom_point(alpha = 0.75) + 
      stat_cor(method = "pearson", size = 5) + 
      geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
      geom_vline(xintercept = 0, color = "grey50", linetype = "dotted") +
      geom_smooth(data = filter(data, color != "n.s"), 
                  method = lm, 
                  se = FALSE, 
                  color = "turquoise4") + 
      geom_label_repel(data = filter(data, Gene %in% labels), 
                       aes(label = Gene), size = 5,
                       min.segment.length = 0) + 
      labs(title = c("Cor. GTEx Male vs. Female and AR-Treatment"), color = "sig. in GTEx") +
      xlab("logFC AR-treatment") + 
      ylab("logFC Male vs. Female in GTEx Adipose breast") +
      #ylim(c(-2.2, 4)) +
      scale_color_manual(values = c("black", "#A6499B")) +
      facet_wrap(facets = "Cluster", ncol =  3, scales = "free") +
      style +
      theme(legend.position = "bottom")
    
  }
  
  plot_vasculature {
    
    labels <- c("PPARG", "CD36", "FOXO3", "FLT4", "KDR", "VEGFA")
    celltypes <- c("Blood_EC")
    
    data <- 
      DE_Response %>% 
      filter(Cluster %in% celltypes, abs(avg_log2FC) > 0.3, p_val_adj <= 0.05) %>%
      select(Gene, Cluster, "logFC_AR" = avg_log2FC) %>% 
      inner_join(filter(response_vasculature), by = "Gene") %>% 
      mutate(color = ifelse(p_val_adj <= 0.05, "p-adj < 0.05", "n.s")) %>% 
      select(Gene, Cluster, color, logFC_AR, "logFC_GTEx" = avg_log2FC)
    
    #p2 <- 
      data %>% 
      filter(abs(logFC_GTEx) < 4, !is.na(color)) %>% 
      ggplot(aes(x = logFC_AR, y = logFC_GTEx, color = color)) + 
      geom_point(alpha = 0.75) + 
      stat_cor(method = "pearson", size = 5) + 
      geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
      geom_vline(xintercept = 0, color = "grey50", linetype = "dotted") +
      geom_smooth(data = filter(data, color != "n.s"), 
                  method = lm, 
                  se = FALSE, 
                  color = "turquoise4") + 
      geom_label_repel(data = filter(data, Gene %in% labels), 
                       aes(label = Gene), size = 5,
                       min.segment.length = 0) + 
      labs(title = c("Cor. GTEx Male vs. Female and AR-Treatment"), color = "sig. in GTEx") +
      xlab("logFC AR-treatment") + 
      ylab("logFC Male vs. Female in GTEx Vascular breast") +
      #ylim(c(-2.2, 4)) +
      scale_color_manual(values = c("black", "#A6499B")) +
      facet_wrap(facets = "Cluster", ncol =  3, scales = "free") +
      style +
      theme(legend.position = "bottom")
    
  }
  
  plot_fibroblast {
    
    labels <- c("LAMA2", "LAMB1", "FOXO3", "INSR", "THBS1", "VEGFA", "NR4A1", "GREB1L")
    celltypes <- c("Fibroblast")
    
    data <- 
      DE_Response %>% 
      filter(Cluster %in% celltypes, abs(avg_log2FC) > 0.2, p_val_adj <= 0.05) %>%
      select(Gene, Cluster, "logFC_AR" = avg_log2FC) %>% 
      inner_join(filter(response_fibor25), by = "Gene") %>% 
      mutate(color = ifelse(p_val_adj <= 0.05, "p-adj < 0.05", "n.s")) %>% 
      select(Gene, Cluster, color, logFC_AR, "logFC_GTEx" = avg_log2FC)
    
    #p3 <- 
      data %>% 
      filter(abs(logFC_GTEx) < 4, !is.na(color)) %>% 
      ggplot(aes(x = logFC_AR, y = logFC_GTEx, color = color)) + 
      geom_point(alpha = 0.75) + 
      stat_cor(method = "pearson", size = 5) + 
      geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
      geom_vline(xintercept = 0, color = "grey50", linetype = "dotted") +
      geom_smooth(data = filter(data, color != "n.s"), 
                  method = lm, 
                  se = FALSE, 
                  color = "turquoise4") + 
      geom_label_repel(data = filter(data, Gene %in% labels), 
                       aes(label = Gene), size = 5,
                       min.segment.length = 0) + 
      labs(title = c("Cor. GTEx Male vs. Female and AR-Treatment"), color = "sig. in GTEx") +
      xlab("logFC AR-treatment") + 
      ylab("logFC Male vs. Female in GTEx fibroblast breast") +
      #ylim(c(-2.2, 4)) +
      scale_color_manual(values = c("black", "#A6499B")) +
      facet_wrap(facets = "Cluster", ncol =  3, scales = "free") +
      style +
      theme(legend.position = "bottom")
    
  }
  
  q <- p1 + p2 + theme(axis.title.y = element_blank()) + p3 + theme(axis.title.y = element_blank()) + plot_layout(nrow = 1)
  x <- p + q + plot_layout(nrow = 2)
  x %>% save_plot(name = "GTEX_scatter", scale = 2, w = 16, h = 10, svg = T)
}

Pre_Post_test {
  
  make_the_data {
    
    ### load main dataset  
    Sobj_Main <- readRDS("2021_ReRun/Seurat_Objects/Sobj_Final_Scaled_MetaUpdt.rds")
    Idents(Sobj_Main) <- "CellType"
    
    ### list for storing the results
    data_list    <- vector("list", length = length(levels(Sobj_Main)))
    
    for (i in 1:length(levels(Sobj_Main))) {
      
      ### subset celltype
      Sobj <- subset(Sobj_Main, idents = levels(Sobj_Main)[i])
      
      Idents(Sobj) <- "Menopause"
      table(Sobj$Menopause)
      
      ### downsample to equal numbers
      pre  <- WhichCells(Sobj, idents = c("NA", "pre"))#,  downsample = min(table(Sobj$Menopause)))
      post <- WhichCells(Sobj, idents = c("NA", "post"))#, downsample = min(table(Sobj$Menopause)))
      
      
      ### perform test
      Idents(Sobj) <- "Type"
      
      marks_pre <- 
        FindMarkers(subset(Sobj, cells = pre), 
                    ident.1 = "TM", 
                    ident.2 = "CF", 
                    test.use = "MAST", 
                    logfc.threshold = 0) %>% 
        rownames_to_column("Gene") %>% 
        mutate(Comp = "pre")
      
      marks_post <- 
        FindMarkers(subset(Sobj, cells = post), 
                    ident.1 = "TM", 
                    ident.2 = "CF", 
                    test.use = "MAST", 
                    logfc.threshold = 0) %>% 
        rownames_to_column("Gene") %>% 
        mutate(Comp = "post")
      
      ### combine two lists and store
      data_list[[i]] <- 
        bind_rows(marks_post, marks_pre) %>% 
        filter(p_val_adj <= 0.05) %>% 
        select(Gene, avg_log2FC, Comp) %>% 
        pivot_wider(names_from = Comp, values_from = avg_log2FC) %>% 
        mutate(CellType = paste0(levels(Sobj_Main)[i])) %>% 
        mutate(post = replace_na(post, 0), pre = replace_na(pre, 0)) %>% 
        mutate(gene_pool = ifelse(post == 0, "unique-pre", "shared"), 
               gene_pool = ifelse(pre == 0, "unique-post", gene_pool))  
      
    }
    
    
  }

  load_data{
  ### data for Plotting the TM vs Pre + Post DE correlation
  load_DE_test_data {
  ### select from subsetted or all cell testing
  #Pre_Post_DEtesting <- do.call(rbind.data.frame, readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/Pre_Post_DEtesting_subsetted.rds"))
  Pre_Post_DEtesting <- do.call(rbind.data.frame, readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/Pre_Post_DEtesting_allcells.rds"))
  }
  ### data for Plotting the MensDE scores by sample average
  load_score_data{
    data_avrg  <- readRDS("2021_ReRun/output/DE_ResponseScores_For_Menopausal_Impact/DE_ResponseScores_SampleAverage_FC_0.25.rds")   
    
    ### filtering for modules that are relevant for each celltype i.e.(LUM_HR-pos only keeps scores from LUM_HR-pos DE genes)
    data_avrg_clean <- 
      data_avrg %>% 
      mutate(x = ifelse(str_detect(Module, as.character(CellType)), "A", "B")) %>% 
      filter(x == "A") %>% 
      select(-x)
    
  }
  }
  
  
  make_label_list {
    label_list        <- vector("list", length = length(unique(Pre_Post_DEtesting$CellType)))
    names(label_list) <- unique(Pre_Post_DEtesting$CellType)
    fill_label_list {
      label_list[["LUM_HR-pos"]] <- c("RYR2", "CUX2", "PGR", "AREG", "EREG", "NR4A1", "AR", "BATF", "JUN", "ESR1", "ACACA", "FASN", "AZGP1", "KCNC2")
      label_list[["LUM_HR-neg"]] <- c("MEF2A", "ITGA1", "ITGB8", "GOLGA8A", "ITGA2", "ENAH", "FOXO3", "CTNNB1", "NR4A1", "AZGP1", "ESRRG")
      label_list[["Basal"]]      <- c("PROS1", "ACTA2", "OXTR", "NRP1", "FOXO3", "AZGP1", "NR4A1", "TP63", "ITGB1")
      label_list[["Fibroblast"]] <- c("LAMA2", "LAMB1", "GREB1L", "THBS1", "VEGFA", "NR4A1", "ADCY2", "FN1")
      label_list[["Adipocyte"]]  <- c("GREB1L", "NR4A1", "ADCY2", "FASN", "INSR", "TCF7L2", "GPAM", "ACACA", "THBS1")
      label_list[["Blood_EC"]]   <- c("PPARG", "CD36", "NR4A1", "FOXO3", "ADAMTSL1")
      label_list[["Lymph_EC"]]   <- c("FLT4", "KDR", "FN1")
      label_list[["Vasc.Acc."]]  <- c("MT2A", "PLCB1")
      label_list[["Myeloid"]]    <- c("AXL", "MERTK", "MRC1", "CD163", "PLAUR")
      label_list[["Lymphoid"]]   <- c("TCF7", "LENG8", "TXK")
      
    }
  }
  
  plotting_functions{
  ### plotting function
  plot_pre_post_single <- function(celltype) {
    
    data <- 
      Pre_Post_DEtesting %>% 
      filter(CellType == celltype) #%>% 
      #mutate(thresh = ifelse(abs(post) > 0.25, "pass", "fail")) %>% 
      #mutate(thresh = ifelse(abs(pre) > 0.25, "pass", thresh))
    
    
    p1 <- 
      data %>% 
      filter(gene_pool == "shared") %>% 
      ggplot(aes(x = pre, 
                 y = post#, color = thresh
      )) + 
      
      geom_hline(yintercept = 0,  color = "grey80", linetype='dashed', size = 1) +
      geom_vline(xintercept = 0,  color = "grey80", linetype='dashed', size = 1) +
      geom_smooth(method = "lm", se = F, color = "lightseagreen") +
      
      geom_point(alpha = 0.5, color = GID_cols[1]) +
      stat_cor(method = "pearson", size = 5) + 
      
      labs(title = paste0("DE Correlation Pre & Post - ", celltype)) +
      xlab("TM vs. Pre") +
      ylab("TM vs. Post") +
      
      geom_text_repel(data = filter(data, Gene %in% label_list[[celltype]]), 
                      size = 5, color = "grey20",
                      aes(label = Gene), 
                      min.segment.length = 0) + 
      style
    
    
    data_bar <- 
      data %>% 
      filter(pre  == 0|abs(pre)  > 0.2 ) %>% 
      filter(post == 0|abs(post) > 0.2 ) %>% 
      group_by(gene_pool) %>% 
      summarise(n_DEgenes = n()) %>% 
      mutate(main = "DE_Genes")
    
    max <- sum(data_bar$n_DEgenes)
    
    p2 <- 
      data_bar %>% 
      ggplot(aes(y = main, x = n_DEgenes, fill = gene_pool)) + 
      geom_col(position = position_stack()) + 
      scale_fill_manual(values = Mens_cols[c(1,3,2)]) +
      scale_x_continuous(breaks = seq(0, max, by = max), limits=c(0, max)) +
      style + 
      theme(axis.title = element_blank(), 
            legend.position = "bottom", 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    
    p1 + p2 + patchwork::plot_layout(ncol = 1, heights = c(10,1))
    
  }
  plot_pre_post_single_vert <- function(celltype) {
    
    data <- 
      Pre_Post_DEtesting %>% 
      filter(CellType == celltype) %>% 
      mutate(thresh = ifelse(abs(post) > 0.25, "pass", "fail")) %>% 
      mutate(thresh = ifelse(abs(pre) > 0.25, "pass", thresh))
    
    
    p1 <- 
      data %>% 
      filter(gene_pool == "shared") %>% 
      ggplot(aes(x = pre, 
                 y = post#, color = thresh
      )) + 
      
      geom_hline(yintercept = 0,  color = "grey80", linetype='dashed', size = 1) +
      geom_vline(xintercept = 0,  color = "grey80", linetype='dashed', size = 1) +
      geom_smooth(method = "lm", se = F, color = "lightseagreen") +
      
      geom_point(alpha = 0.5, color = GID_cols[1]) +
      stat_cor(method = "pearson", size = 5) + 
      
      labs(title = paste0("DE Correlation Pre & Post - ", celltype)) +
      xlab("TM vs. Pre") +
      ylab("TM vs. Post") +
      
      geom_text_repel(data = filter(data, Gene %in% label_list[[celltype]]), 
                      size = 5, color = "grey20",
                      aes(label = Gene), 
                      min.segment.length = 0) + 
      style
    
    
    data_bar <- 
      data %>% 
      filter(pre  == 0|abs(pre)  > 0.2 ) %>% 
      filter(post == 0|abs(post) > 0.2 ) %>% 
      group_by(gene_pool) %>% 
      summarise(n_DEgenes = n()) %>% 
      mutate(main = "DE_Genes")
    
    max <- sum(data_bar$n_DEgenes)
    
    p2 <- 
      data_bar %>% 
      ggplot(aes(x = main, y = n_DEgenes, fill = gene_pool)) + 
      geom_col(position = position_stack()) + 
      scale_fill_manual(values = Mens_cols[c(1,3,2)]) +
      style + 
      scale_y_continuous(breaks = seq(0, max, by = max), limits=c(0, max)) +
      theme(axis.title = element_blank(), 
            legend.position = "none")
    
    p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(10,1))
    
  }
  plot_pre_post_single_dotOnly <- function(celltype) {
    
    data <- 
      Pre_Post_DEtesting %>% 
      filter(CellType == celltype) %>% 
      mutate(thresh = ifelse(abs(post) > 0.25, "pass", "fail")) %>% 
      mutate(thresh = ifelse(abs(pre) > 0.25, "pass", thresh))
    
    
    #p1 <- 
      data %>% 
      filter(gene_pool == "shared") %>% 
      ggplot(aes(x = pre, 
                 y = post#, color = thresh
      )) + 
      
      geom_hline(yintercept = 0,  color = "grey80", linetype='dashed', size = 1) +
      geom_vline(xintercept = 0,  color = "grey80", linetype='dashed', size = 1) +
      geom_smooth(method = "lm", se = F, color = "lightseagreen") +
      
      geom_point(alpha = 0.5, color = GID_cols[1]) +
      stat_cor(method = "pearson", size = 5) + 
      
      labs(title = paste0("Pre & Post - ", celltype)) +
      xlab("TM vs. Pre") +
      ylab("TM vs. Post") +
      
      geom_text_repel(data = filter(data, Gene %in% label_list[[celltype]]), 
                      size = 5, color = "grey20",
                      aes(label = Gene), 
                      min.segment.length = 0) + 
      style + theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
  
    
  }
  plot_pre_post_multi  <- function(celltype) {
    
    data <- 
      Pre_Post_DEtesting %>% 
      filter(CellType == celltype) %>% 
      mutate(thresh = ifelse(abs(post) > 0.25, "pass", "fail")) %>% 
      mutate(thresh = ifelse(abs(pre)  > 0.25, "pass", thresh))
    
    
    p1 <- 
      data %>% 
      filter(gene_pool == "shared") %>% 
      ggplot(aes(x = pre, 
                 y = post#, color = thresh
      )) + 
      
      geom_hline(yintercept = 0,  color = "grey80", linetype='dashed', size = 1) +
      geom_vline(xintercept = 0,  color = "grey80", linetype='dashed', size = 1) +
      geom_smooth(method = "lm", se = F, color = "lightseagreen") +
      
      geom_point(alpha = 0.5, color = GID_cols[1]) +
      stat_cor(method = "pearson", size = 5) + 
      
      labs(title = paste0(celltype)) +
      xlab("TM vs. Pre") +
      ylab("TM vs. Post") +
      
      geom_text_repel(data = filter(data, Gene %in% label_list[[celltype]]), 
                      size = 5, color = "grey20",
                      aes(label = Gene), 
                      min.segment.length = 0) + 
      style 
    
    
    
    data_bar <- 
      data %>% 
      filter(pre  == 0|abs(pre)  > 0.2 ) %>% 
      filter(post == 0|abs(post) > 0.2 ) %>% 
      group_by(gene_pool) %>% 
      summarise(n_DEgenes = n()) %>% 
      mutate(main = "DE_Genes")
    
    max <- sum(data_bar$n_DEgenes)
    
    p2 <- 
      data_bar %>% 
      ggplot(aes(y = main, x = n_DEgenes, fill = gene_pool)) + 
      geom_col(position = position_stack()) + 
      scale_fill_manual(values = Mens_cols[c(1,3,2)]) +
      scale_x_continuous(breaks = seq(0, max, by = max), limits=c(0, max)) +
      style + 
      theme(axis.title = element_blank(), 
            legend.position = "none", 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    
    p1 + p2 + patchwork::plot_layout(ncol = 1, heights = c(10,1))
    
  }
  plot_mens_score <- function(celltype) {
    
    data_avrg_clean %>% 
      filter(CellType == celltype) %>% 
      mutate(Module = str_remove(Module, paste0(celltype, "_"))) %>% 
      mutate(Module = factor(Module, levels = c("upreg", "downreg"))) %>% 
      
      ggplot(aes(x = factor(Type, levels = c("TM", "CF")), 
                 fill = factor(Menopause, levels = c("NA", "pre", "post")), 
                 y = value)) + 
      
      geom_boxplot(color = "grey30",
                   width = 0.75, alpha = 1, 
                   lwd = 1.25, 
                   position = position_dodge(0.9), 
                   outlier.alpha = 0) +
      stat_summary(fun = "median", 
                   geom = "point", size = 1.5,
                   position = position_dodge(0.9),
                   color = "floralwhite") +
      
      geom_point(aes(color = factor(Sample, levels = SP_levels)), 
                 position = position_jitter(width = 0.15), size = 2) +
      
      scale_fill_manual(values = Mens_cols) +
      
      scale_color_manual(values = SP_cols) +
      
      labs(title = paste0(celltype), 
           x     = "Gender", 
           y     = "avg. DE score / Sample", 
           fill  = "Menopause") +
      
      #facet_wrap("CellType", scales = "free_y", ncol = 10) + 
      facet_wrap("Module", scales = "free_y", 
                 nrow = 2) + 
      style + theme(legend.position = "none", axis.text.y = element_text(angle = 90),
                    axis.title.x = element_blank(), 
                    axis.title.y = element_blank())
                    #axis.ticks.y = element_blank())
  }
  }
  
  
  
  
  ### plot around
  celltype <- "Lymphoid"
  
  # plot it
  plot_pre_post_single_dotOnly(celltype)
  plot_mens_score(celltype) + theme(legend.position = "bottom")
  
  ## plot all of them
  plot_list    <- vector("list", length = length(unique(Pre_Post_DEtesting$CellType)))
  
  for (i in 1:length(unique(Pre_Post_DEtesting$CellType))) {
    plot_list[[i]] <- plot_pre_post_multi(unique(Pre_Post_DEtesting$CellType)[i])
  }
  patchwork::wrap_plots(c(plot_list[3:10]), ncol = 4)
  
  
  
  
  ### plot combined stuff for all celltypes
  plot_list    <- vector("list", length = length(unique(Pre_Post_DEtesting$CellType)))
  
  for (i in 1:length(unique(Pre_Post_DEtesting$CellType))) {
    
    p2 <- plot_mens_score(unique(Pre_Post_DEtesting$CellType)[i])
    #p2 <- plot_pre_post_single_dotOnly(unique(Pre_Post_DEtesting$CellType)[i])
    
    #plot_list[[i]] <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,2))
    plot_list[[i]] <- p2 + theme(axis.title = element_blank())
    
  }
  
  p <- patchwork::wrap_plots(c(plot_list[1:10]), ncol = 10)
  
  
  
  as.ggplot(plot_list[[1]]) + plot_spacer() + as.ggplot(plot_list[[2]]) + 
  as.ggplot(plot_list[[3]]) + plot_spacer() + as.ggplot(plot_list[[4]]) + 
    plot_layout(ncol = 3, widths = c(1,0.1,1))
  as.ggplot(plot_list[[5]]) + plot_spacer() + as.ggplot(plot_list[[6]]) + 
  as.ggplot(plot_list[[7]]) + plot_spacer() + as.ggplot(plot_list[[8]]) + 
    plot_layout(ncol = 3, widths = c(1,0.1,1))
  as.ggplot(plot_list[[7]]) + plot_spacer() + as.ggplot(plot_list[[8]]) + 
  as.ggplot(plot_list[[9]]) + plot_spacer() + as.ggplot(plot_list[[10]]) + 
  plot_layout(ncol = 3, widths = c(1,0.1,1))
  
    
    
    
}

blood_pressure_etc {
plot_cohort_details <- function() {

my_comparisons <- list(c("post", "pre"), c("NA", "pre"), c("post", "NA"))

p1 <- 
HRT_data %>% 
  separate(BP, into = c("Cist", "Dia")) %>% 
  mutate(Cist = as.numeric(Cist)) %>% 
  ggplot(aes(x = factor(Menopause, levels = c("NA", "pre", "post")), y = Cist, fill = factor(Menopause, levels = c("NA", "pre", "post")))) + 
  geom_boxplot(size = 1, color = "grey30") + 
  scale_fill_manual(values = Mens_cols) +
  stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) +
  stat_summary(fun = "median", 
               geom = "point", size = 3,
               position = position_dodge(0.75),
               color = "floralwhite") +
  labs(title = "Cistolic BP", fill = "Menopause", x = NULL) + 
  style + theme(legend.position = "none")

p2 <- 
HRT_data %>% 
  separate(BP, into = c("Cist", "Dia")) %>% 
  mutate(Dia = as.numeric(Dia)) %>% 
  ggplot(aes(x = factor(Menopause, levels = c("NA", "pre", "post")), y = Dia, fill = factor(Menopause, levels = c("NA", "pre", "post")))) + 
  geom_boxplot(size = 1, color = "grey30") + 
  scale_fill_manual(values = Mens_cols) +
  stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) +
  stat_summary(fun = "median", 
               geom = "point", size = 3,
               position = position_dodge(0.75),
               color = "floralwhite") +
  labs(title = "Diastolic BP", fill = "Menopause", x = NULL) + 
  style + theme(legend.position = "none")

p3 <- 
HRT_data %>% 
  separate(BP, into = c("Cist", "Dia")) %>% 
  mutate(Cist = as.numeric(Cist)) %>% 
  ggplot(aes(x = factor(Menopause, levels = c("NA", "pre", "post")), y = BMI, fill = factor(Menopause, levels = c("NA", "pre", "post")))) + 
  geom_boxplot(size = 1, color = "grey30") + 
  scale_fill_manual(values = Mens_cols) +
  stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) +
  stat_summary(fun = "median", 
               geom = "point", size = 3,
               position = position_dodge(0.75),
               color = "floralwhite") +
  labs(title = "Body Mass Index", fill = "Menopause", x = NULL) + 
  style + theme(legend.position = "none")

p4 <- 
HRT_data %>% 
  separate(BP, into = c("Cist", "Dia")) %>% 
  mutate(Cist = as.numeric(Cist)) %>% 
  ggplot(aes(x = factor(Menopause, levels = c("NA", "pre", "post")), y = `Gluc (mg/dL)`, fill = factor(Menopause, levels = c("NA", "pre", "post")))) + 
  geom_boxplot(size = 1, color = "grey30") + 
  scale_fill_manual(values = Mens_cols) +
  stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) +
  stat_summary(fun = "median", 
               geom = "point", size = 3,
               position = position_dodge(0.75),
               color = "floralwhite") +
  labs(title = "Body Mass Index", fill = "Menopause", x = NULL) + 
  style + theme(legend.position = "none")

p1 + p2 + p3 + p4

}
plot_cohort_details()
}



Nuclei_Libraries_plus_ATAC_V4 <- 
  read_excel("2021_ReRun/Documents/Nuclei-Libraries_plus_ATAC_V4.xlsx", 
             sheet = "Nuclei - Plotting") %>% 
  separate(Sample, into = c("Type", "Code"), remove = F)

names <- colnames(Nuclei_Libraries_plus_ATAC_V4)
i = 32
Nuclei_Libraries_plus_ATAC_V4 %>% ggplot(aes(x = Type, y = get(names[i]), fill = Type)) + geom_boxplot() + ggtitle(paste0(names[i]))

my_comparisons <- list(c("TM", "CF"))

  
p1 <- Nuclei_Libraries_plus_ATAC_V4 %>% ggplot(aes(x = Type, y = get(names[5]),  fill = Type)) + geom_boxplot() + stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) + ggtitle(paste0(names[5]))  + theme(legend.position = "none") + style
p2 <- Nuclei_Libraries_plus_ATAC_V4 %>% ggplot(aes(x = Type, y = get(names[6]),  fill = Type)) + geom_boxplot() + stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) + ggtitle(paste0(names[6]))  + theme(legend.position = "none") + style
p3 <- Nuclei_Libraries_plus_ATAC_V4 %>% ggplot(aes(x = Type, y = get(names[9]),  fill = Type)) + geom_boxplot() + stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) + ggtitle(paste0(names[9]))  + theme(legend.position = "none") + style + labs(title = c("Genes per Cell"))
p4 <- Nuclei_Libraries_plus_ATAC_V4 %>% ggplot(aes(x = Type, y = get(names[12]), fill = Type)) + geom_boxplot() + stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) + ggtitle(paste0(names[12])) + theme(legend.position = "none") + style + labs(title = c("Seq. Saturation"))
p5 <- Nuclei_Libraries_plus_ATAC_V4 %>% ggplot(aes(x = Type, y = get(names[20]), fill = Type)) + geom_boxplot() + stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) + ggtitle(paste0(names[20])) + theme(legend.position = "none") + style + labs(title = c("Mapped Exon Reads"))
p6 <- Nuclei_Libraries_plus_ATAC_V4 %>% ggplot(aes(x = Type, y = get(names[19]), fill = Type)) + geom_boxplot() + stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) + ggtitle(paste0(names[19])) + theme(legend.position = "none") + style + labs(title = c("Mapped Intron Reads"))

p1 + p2 + p3 + p4 + p5 + p6




Idents(Sobj) <- "Subcluster"
Sobj$Subcluster <- Idents(Sobj)

my_comparisons <- list(c("post", "pre"), c("NA", "pre"), c("post", "NA"))

TM_CF {
my_comparisons <- list(c("CF", "TM"))
as.data.frame(prop.table(table(Sobj$Subcluster, 
                               Sobj$Sample), 
                         margin = 2)) %>% filter(Var2 != "TM-9817") %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2)) %>% 
  left_join( distinct(select(Sobj@meta.data, Sample, Menopause))) %>% 
  #filter(Type == "CF") %>% 
  #ggplot(aes(x = factor(Menopause, levels = c("NA", "pre", "post")), fill = factor(Menopause, levels = c("NA", "pre", "post")), y = Freq)) + 
  ggplot(aes(x = Type, fill = Type, y = Freq)) + 
  geom_boxplot(width = 0.5) + 
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
  #scale_y_sqrt() + 
  #scale_fill_manual(values = Mens_cols) + 
  scale_fill_manual(values = GID_cols) + 
  style + 
  theme(legend.position = "bottom", 
        axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  facet_wrap(facets = "Cluster", scales = "free_y") 

}
Pre_Post {
my_comparisons <- list(c("post", "pre"), c("NA", "pre"), c("post", "NA"))
as.data.frame(prop.table(table(Sobj$Subcluster, 
                               Sobj$Sample), 
                         margin = 2)) %>% 
  dplyr::select(Cluster = "Var1", Sample = "Var2", Freq = "Freq") %>% 
  mutate(Type = substring(.$Sample, 1, 2)) %>% 
  left_join( distinct(select(Sobj@meta.data, Sample, Menopause))) %>% 
  #filter(Type == "CF") %>% 
  ggplot(aes(x = factor(Menopause, levels = c("NA", "pre", "post")), fill = factor(Menopause, levels = c("NA", "pre", "post")), y = Freq)) + 
  #ggplot(aes(x = Type, fill = Type, y = Freq)) + 
  geom_boxplot(width = 0.5) + 
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
  #scale_y_sqrt() + 
  scale_fill_manual(values = Mens_cols) + 
  style + 
  theme(legend.position = "bottom", 
        axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  facet_wrap(facets = "Cluster", scales = "free_y") 
    
}




















Process_ChromVar_Matrix {

Motif.cisBPMatrix <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/Motif.cisBPMatrix.RDS")

df <- Motif.cisBPMatrix@assays@data$z %>% as.data.frame() %>% t()
motifs <- colnames(df) %>% as_tibble() %>% separate(value, into = c("Gene", "add"), sep = "_", remove = F) %>% select(Motif = "value", Gene)
meta.data <- Motif.cisBPMatrix@colData %>% as.data.frame() %>% rownames_to_column("BC") %>% select(BC, SampleName, CellType = "predictedCellType", Type =  "SampleType")


df_sum <- 
  df %>% 
  as.data.frame() %>% 
  rownames_to_column("BC") %>% 
  left_join(meta.data) %>% 
  pivot_longer(cols = motifs$Motif, names_to = "Motif") %>% 
  group_by(Type, SampleName, CellType, Motif) %>% 
  summarise(mean = mean(value)) 


df_long <- 
  df %>% 
  as.data.frame() %>% 
  rownames_to_column("BC") %>% 
  left_join(meta.data) %>% 
  pivot_longer(cols = motifs$Motif, names_to = "Motif") 


#df_long %>% saveRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/ChromVar_AllCells_Long.rds")
#df_sum %>% saveRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/output/ChromVar_Average_CellType.rds")

}





celltype <- "Basal"

query <- c("NR4A1")
query <- 
  DE_Response %>% filter(Cluster == celltype, abs(avg_log2FC) > 0.25, TF == "YES") %>% top_n(30, avg_log2FC) %>% arrange(-avg_log2FC) %>% pull(Gene)
mot <- motifs %>% filter(Gene %in% query)
motis <- mot[match(query, mot$Gene),] %>% filter(!is.na(Motif)) %>% pull(Motif)

data <- df %>% 
  rownames_to_column("BC") %>% 
  left_join(meta.data) %>% 
  select(BC, Type, CellType, motis) %>% 
  pivot_longer(motis, names_to = "Motif", values_to = "z.score") %>% 
  filter(CellType == celltype)

data$Motif <- factor(data$Motif, levels = motis)
  
  data %>% ggplot(aes(x = CellType, y = z.score, fill = Type)) + 
  geom_violin(scale = "width", 
              trim = F,
              lwd = 1.25, 
              color = "grey30") +
  geom_boxplot(aes(color = Type), 
               fill = "grey30", 
               width = 0.25, 
               lwd = 1.25,
               position=position_dodge(0.9), outlier.alpha = 0) + 
  stat_summary(fun = "median", 
               geom = "point", size = 3,
               position = position_dodge(0.9),
               color = "floralwhite") +
  
  scale_fill_manual(values = GID_cols) +
  scale_color_manual(values = c("grey30", "grey30")) +
  style + 
  theme(legend.position = "bottom", 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap("Motif", scales = "free_y", ncol = 7)


  
  
  
  
  
  
  
  test2 <- readRDS("2021_ReRun/output/ChromVar_Average_CellType.rds")
  
  celltype = "Basal"
  
  labels <- c("MEF2A", "KLF6", "NFIC", "TCF7L2", "NR3C1", "NFIB", "NFIA", "NFIX", "EGR1", "KLF12", "TEAD1", "TP63", "FOS", "BACH2")
  
  test2 %>% 
    left_join(filter(DE_Response, Cluster == celltype)) %>% 
    filter(!is.na(avg_log2FC))
    
    
    plot_data <- data %>% 
    separate(Motif, into = c("Gene", "Num"), sep = "_") %>% 
    pivot_wider(names_from = Type, values_from = mean) %>% 
    mutate(delta_TM_CF = `Trans male` - `Cis female`) %>% 
    left_join(filter(DE_Response, Cluster == celltype)) %>% 
    filter(!is.na(avg_log2FC))
    
    p <- plot_data %>% 
    ggplot(aes(x = avg_log2FC, y = delta_TM_CF)) + 
    geom_hline(yintercept = 0, color =   "turquoise") +
    geom_vline(xintercept = 0, color =   "turquoise") +
    geom_point() + 
    geom_text_repel(data = filter(plot_data, Gene %in% labels), 
                    aes(label = Gene)) + 
      style


  
  
  DE_match <- readRDS("2021_ReRun/output/DE_list_Receptor_Ligand_Matched.rds")
  
  DE_match %>% 
    filter(!is.na(PMIDs), 
           Cluster.x == "LUM_HR-neg", 
           Rec_logFC > 0.25, 
           Lig_logFC > 0.25, 
           Cluster.y %in% c("LUM_HR-pos", "Fibroblast", "Adipocyte", "LUM_HR-neg", "Basal")) %>% 
    arrange(-Rec_logFC)
  
  
  DE_match %>% 
    filter(!is.na(PMIDs), 
           Cluster.x == "Basal", 
           Rec_logFC < -0.2, 
           Lig_logFC < -0.25, 
           Cluster.y %in% c("LUM_HR-pos", "Fibroblast", "Adipocyte", "LUM_HR-neg", "Basal")) %>% 
    arrange(Rec_logFC)

  
  

  
  
  
  recps <- 
  DE_Response_by_Subcluster %>% 
    filter(Cluster == "CD4_T", 
           Gene %in% cabello_aguilar$receptor, 
           avg_log2FC > 0.2, p_val_adj <= 0.05) %>% 
    arrange(-avg_log2FC) %>% 
    pull(Gene)
  
  
  ## TFs only up in CD4
  DE_Response_by_Subcluster %>% 
    filter(Cluster %in% c("CD4_T", "CD8_T"), 
           Gene %in% hs_hgnc_tfs$X1, 
           avg_log2FC > 0.2, p_val_adj <= 0.05) %>% 
    arrange(-avg_log2FC) %>% 
    select(Gene, avg_log2FC, Cluster) %>%
    pivot_wider(names_from = Cluster, values_from = avg_log2FC) %>% 
    filter(is.na(CD8_T))
  
  
  ## TFs only up in CD4
  DE_Response_by_Subcluster %>% 
    filter(Cluster %in% c("CD4_T", "CD8_T"), 
           Gene %in% cabello_aguilar$receptor, 
           avg_log2FC > 0.2, p_val_adj <= 0.05) %>% 
    arrange(-avg_log2FC) %>% 
    select(Gene, avg_log2FC, Cluster) %>%
    pivot_wider(names_from = Cluster, values_from = avg_log2FC) %>% 
    filter(is.na(CD8_T))
  
  DE_Response_by_Subcluster %>% 
    filter(Cluster %in% c("Macrophage"), 
           Gene %in% cabello_aguilar$receptor, 
           avg_log2FC < -0.2, p_val_adj <= 0.05) %>% 
    arrange(-avg_log2FC) %>% 
    select(Gene, avg_log2FC, Cluster) %>%
    pivot_wider(names_from = Cluster, values_from = avg_log2FC) %>% 
    filter(is.na(CD8_T))
  
  
  DE_Response_by_Subcluster %>% 
         filter(Cluster == "CD4_T", 
                Gene %in% cabello_aguilar$ligand, 
                avg_log2FC > 0.1, p_val_adj <= 0.05) %>% 
         arrange(-avg_log2FC) 

  
  
  
ligs <- cabello_aguilar %>% filter(!is.na(PMIDs), receptor %in% recps[1:2]) %>% pull(ligand) %>% unique()  

DE_Response %>% filter(Gene %in% ligs, avg_log2FC > 0.2, p_val_adj <= 0.05) %>% arrange(-avg_log2FC)
DE_Response_by_Subcluster %>% filter(Gene %in% ligs, avg_log2FC > 0.2, p_val_adj <= 0.05) %>% arrange(-avg_log2FC)


  






genes1 <- 
  Subcluster_markers %>% 
  filter(CellType == "Lymphoid") %>% 
  filter(!is.na(avg_log2FC), p_val_adj <= 0.05, pct.1 > 0.25) %>% 
  filter(!str_detect(gene, "-AS1|-AS2|\\.|LINC0")) %>% 
  select(-p_val, -p_val_adj) %>% 
  group_by(cluster) %>% 
  top_n(3, avg_log2FC) %>% 
  arrange(cluster) %>% pull(gene)





# Supplementary Tables ##############
writing_supp_DE_tables{
  
  library(xlsx)
  
  ### by celltype
  clusters    <- unique(DE_Response$Cluster)
  sheet.names1 <- c("lum-hr+", "lum-hr-", "basal", "fibroblast", "adipocyte", "blood-ec", "lymph-ec", "vasc.acc.", "myeloid", "lymphoid")
  names(clusters) <- sheet.names1
  
  Response <- 
    DE_Response %>% 
    left_join(select(rownames_to_column(as.data.frame(clusters), "Cluster_new"), 
                     Cluster = clusters, Cluster_new)) %>% 
    left_join(as_tibble(distinct(select(Sobj@meta.data, CellType, Cluster = Subcluster))))
  
  for (i in 1:length(clusters)) {
    
    data <- 
      Response %>% 
      filter(Cluster == clusters[i]) %>% 
      select(Gene, avg_log2FC, p_val_adj, pct.TM = pct.1, pct.CF = pct.2, Cluster = Cluster_new) %>% 
      mutate(Direction = ifelse(avg_log2FC > 0, "up in TM", "up in CF")) %>% 
      arrange(-avg_log2FC)
    
    write.xlsx(as.data.frame(data), 
               file = "~/Documents/Main_DE_Response2.xlsx", 
               sheetName = paste0(sheet.names[i]), 
               row.names=FALSE, 
               append=TRUE
    )
    
  }
  
  
  
  ### by subcluster
  clusters    <- unique(DE_Response_by_Subcluster$Cluster)
  sheet.names2 <- c("lup-1", "lup-2", "lup-3", "lup-4", "lup-ribo", "lup-cycling", 
                    "lun-1", "lun-2", "lun-3", "lun-ribo", "lun-cycling", "lun-4", "lun-5", 
                    "bas-1", "bas-2", "bas-3", "bas-ribo", 
                    "matrix-2", "lipo-f", "matrix-1", "vasc-f", "chondrocyte",
                    "adipocyte", 
                    "vein", "capillary", "artery",
                    "lymph-ec", "lymph-ec2", 
                    "vasc_sm1", "vasc-sm2", "pericyte",
                    "macrophage",  "mono.DC", "monocyte", "DC", "HSC",
                    "CD4", "CD8", "T-eff", "NK", "B-plasma", "B.mzone"
  )
  
  names(clusters) <- sheet.names2
  
  DE_Res_mut <- 
    DE_Response_by_Subcluster %>% 
    left_join(select(rownames_to_column(as.data.frame(clusters), "Cluster_new"), Cluster = clusters, Cluster_new)) %>% 
    left_join(as_tibble(distinct(select(Sobj@meta.data, CellType, Cluster = Subcluster)))) %>% 
    left_join(distinct(select(Response, CellType = "Cluster", CellType_new = "Cluster_new")), by = "CellType")
  
  
  for (i in 1:length(sheet.names1)) {
    
    data <- 
      DE_Res_mut %>% 
      filter(CellType_new == sheet.names1[i], abs(avg_log2FC) > 0.2) %>% 
      select(Gene, avg_log2FC, p_val_adj, pct.TM = pct.1, pct.CF = pct.2, Subcluster = Cluster_new, CellType = CellType_new) %>% 
      arrange(Subcluster, -avg_log2FC) %>% 
      mutate(Direction = ifelse(avg_log2FC > 0, "up in TM", "up in CF"))
    
    write.xlsx(as.data.frame(data), 
               file = "~/Documents/Subcluster_DE_Response.xlsx", 
               sheetName = paste0(sheet.names1[i]), 
               row.names=FALSE, 
               append=TRUE
    )
  }
  
  
}

query <- DE_Response %>% filter(Cluster == "LUM_HR-pos", abs(avg_log2FC) > 0.2) %>% pull(Gene)

GTEx_avg <- 
  dat.gct[ , c("Description", columns[1:5])] %>% 
  filter(Description %in% singles) %>% 
  column_to_rownames("Description") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SAMPID") %>% 
  left_join(pheno) %>% 
  select(SAMPLEID, singles) %>% 
  pivot_longer(query, names_to = "Gene") %>% 
  group_by(SMTS, Gene) %>% 
  mutate(mean = mean(value))





GTEx_avg %>% 
  filter(Gene == "AR") %>%
  select(-SAMPLEID, -value, -SEX) %>% 
  distinct() %>% ungroup() %>% 
  top_n(10, mean) %>% 
  ggplot(aes(x = reorder(SMTS, mean), y = mean)) + 
  geom_col(fill = "lightseagreen") + 
  facet_wrap("Gene") + 
  style + 
  labs(x = "GTEx mean", y = NULL) +
  coord_flip()



Scavenger_Receptors <- read_excel("2021_ReRun/utilities/Scavenger_Receptors.xlsx")

p <- DE_Response_by_Subcluster %>% 
  filter(Subcluster == "Macrophage", Gene %in% c(Scavenger_Receptors$ID...2, "MERTK", "AXL")) %>% 
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
  geom_point(data = filter(DE_Response_by_Subcluster, Subcluster == "Macrophage"), 
             color = "grey70", alpha = 0.5, size = 3) +
  geom_point(color = GID_cols[1], size = 5) + 
  style + xlim(c(-1, 1)) + 
  geom_text_repel(aes(label = Gene), min.segment.length = 0) + 
  geom_vline(xintercept = 0)



p <- Nuclei_Libraries_plus_ATAC_V4 %>% 
  select(Sample, `cDNA size`, `Median Genes per Cell...7`, "Sequencing Saturation...10", "Reads Mapped Confidently to Exonic Regions...18") %>% 
  pivot_longer(2:5) %>% 
  separate(Sample, into = c("Type", "Num"), sep = "-", extra = "merge", remove = F) %>% 
  ggplot(aes(x = Type, y = value, fill = Type)) + 
  geom_boxplot() + 
  geom_point(aes(color = factor(Sample, levels = SP_levels)), 
             position = position_dodge(width = 0.25), 
             size = 3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
  scale_fill_manual(values = GID_cols) +
  scale_color_manual(values = SP_cols) +
  facet_wrap(facets = "name", scales = "free") + 
  style + 
  theme(legend.position = "none")








p <- umap_lumhrpos_originProb %>% 
  filter(!is.na(OriginProb)) %>% 
  arrange(OriginProb) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = OriginProb)) + 
  geom_point(data = filter(umap_lumhrpos_originProb, Subcluster == "lup_4"), 
             color = "grey80", 
             size = 10) +
  geom_point(aes(size = OriginProb)) + 
  scale_color_viridis(option = "D") +
  style




















p <- prop.table(table(Sobj$Subcluster, Sobj$Type), margin = 2) %>% 
  as.data.frame() %>% rename(CellType = "Var1", Type = "Var2") %>% 
  #pivot_wider(names_from = Type, values_from = Freq) %>% 
  #mutate(log2FC = log2(TM/CF)) %>% 
  #left_join(Prop_RNA) %>% 
  #arrange(-row_number()) %>% 
  
  ggplot(aes(x = "", y = Freq, fill = Type)) + 
  geom_bar(stat = "identity", position = position_fill()) +
  coord_polar(theta = "y") +
  facet_wrap("CellType")
  #geom_text(aes(label = round(log2FC, 2)), position = position_fill(vjust = 0.5)) +
  #scale_fill_continuous_divergingx(palette = 'BrBG', mid = 0, limits = c(-1.35, 1.35)) +
  style +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        legend.position = "right")

  
  
  
  
  ATACdata <- ATACdata <- readRDS("2021_ReRun/data_resources/umap_tileMat.rds") %>% 
    rename(UMAP_1 = "UMAP1", UMAP_2 = "UMAP2") %>% 
    mutate(Type = str_replace(SampleType, "Trans male", "TM"), Type = str_replace(Type, "Cis female", "CF")) %>% rownames_to_column("Cell") %>% separate(Cell, into = c("Sample", "ID"), sep = "_", extra = "merge")
  ATACdatafil <- ATACdata %>% filter(CellType %in% c("Lymph_EC", "Blood_EC", "Vasc.Acc."))
  
  p <- prop.table(table(ATACdatafil$CellType, ATACdatafil$Sample), margin = 2) %>% 
    as.data.frame() %>% filter(Var1 == "Lymph_EC") %>% 
    separate(Var2, into = c("Type", "num"), sep = "-", extra = "merge", remove = F) %>% 
    ggplot(aes(x = Type, y = Freq*100, fill = Type)) + 
    geom_boxplot() + 
    geom_point(aes(color = factor(Var2, levels = SP_levels)), 
               position = position_dodge(width = 0.25), 
               size = 3) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = SP_cols) +
    #facet_wrap(facets = "name", scales = "free") + 
    style + 
    ylab("ATAC_proportion_Lymph_EC") +
    theme(legend.position = "none", 
          axis.text.y = element_text(angle = 90, hjust = 0.5))

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
  
  
  
  
  
  
  
  
  
genes <- 
  pgr_ar_annotation %>% 
  as_tibble() %>% 
  filter(Correlation > 0.2, 
         ATAC.DA.Log2FC < -0.1, 
         ATAC.DA.FDR < 0.2, 
         PGR_660 == T,
         AR_689 == F) %>% 
  pull(GeneRna) %>% 
  unique()


DE_Response %>% 
  filter(Gene %in% genes, 
         Cluster == "LUM_HR-pos", 
         avg_log2FC < -0.5) %>% 
  pull(Gene) %>% 
  cat()


















