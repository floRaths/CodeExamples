library(harmony)
library(Seurat)
library(tidyverse)
library(viridis)
library(ggsci)


setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

Meta <- readRDS("2021_ReRun/Seurat_Objects/MetaData_08.25.2021.rds")
#MetaData_list    <- vector("list", length = length(levels(Meta$CellType)))
#names(MetaData_list) <- levels(Meta$CellType)


Idents(Sobj) <- "CellType"

for (i in 1:length(levels(Sobj))) {
        
        Seur <- subset(Sobj, idents = levels(Sobj)[i])
        Seur@reductions <- readRDS(paste0("2021_ReRun/Reductions/", levels(Sobj)[i], ".rds"))
        
        Seur %>% DietSeurat(counts = T, 
                            data = T, 
                            scale.data = F, 
                            assays = "RNA", 
                            dimreducs = c("pca", "harmony", "umap")) %>% 
                
                saveRDS(paste0("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/", levels(Sobj)[i], "_MetaUpdt.rds"))
}



Sobj <- readRDS("2021_ReRun/Seurat_Objects/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/LUM_HR-neg_MetaUpdt.rds")


Sobj <- AddMetaData(Sobj, Meta)
Sobj@reductions <- readRDS(paste0("2021_ReRun/Reductions/", unique(Sobj$CellType), ".rds"))

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Sobj <- CellCycleScoring(Sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
rm(s.genes, g2m.genes)


res <- "RNA_snn_res.0.2"

p2 <-  
DimPlot(Sobj,
        reduction = "umap", 
        #group.by  = "Subcluster",
        group.by  = "Sample",
        split.by  = "Type",
        pt.size = 0.5, 
        label = T,
        repel = T, 
        label.size = 4, cols = pal_d3("category20")(20)
        )

p2 <- DimPlot(subset(Sobj, cells = c(cell_Epi)),
        reduction = "umap", 
        #group.by  = "Subcluster",
        group.by  = res,
        #split.by  = "Type",
        pt.size = 1, 
        label = T,
        repel = T, 
        label.size = 4
        )

p1 + p2

FeaturePlot(Sobj, features = c("AR", "ALDH1A1"), cols = viridis(n = 100, option = "A"), pt.size = 1, order = T)
FeaturePlot(Sobj, features = c("S.Score", "RPLP1", "ESR1"), cols = viridis(n = 100, option = "A"), pt.size = 1, order = T)
FeaturePlot(Sobj, features = c(""), cols = viridis(n = 100, option = "A"), pt.size = 1, order = T)
FeaturePlot(Sobj, features = c("CD44"), cols = viridis(n = 100, option = "A"), pt.size = 0.5, order = T) #+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
VlnPlot(Sobj, group.by = res, features = "Age") #+ ylim(c(0, 5000))


determine.lowqualcells <- function() {

Idents(Sobj) <- res

cell_Epi <- WhichCells(Sobj, idents = c("4", "7", "10"))
cell_Vas <- WhichCells(Sobj, idents = c("8", "5"))
cell_strm <- WhichCells(Sobj, idents = c("4"))
cell_Imu <- WhichCells(Sobj, idents = "5")

cell_lumneg <- WhichCells(Sobj, idents = c("14", "16", "13", "3", "11", "15"))
cell_lumpos <- WhichCells(Sobj, idents = c("9"))
cell_basal  <- WhichCells(Sobj, idents = c("3"))
cell_VasAcc <- WhichCells(Sobj, idents = c("4", "5"))
cell_Blood  <- WhichCells(Sobj, idents = c("3", "4", "5"))
cell_Fibro  <- WhichCells(Sobj, idents = c("5", "2"))
#cell_Adipo  <- WhichCells(Sobj, idents = c("5"))
cell_Lymph  <- WhichCells(Sobj, idents = c("5", "2"))
cell_Myel   <- WhichCells(Sobj, idents = c("3"))


c(cell_lumneg, cell_lumpos, cell_basal, cell_VasAcc, cell_Blood, cell_Fibro, cell_Lymph, cell_Myel)


}



DefaultAssay (Sobj) <- "RNA"

Sobj <- NormalizeData (Sobj, verbose = T)
Sobj <- FindVariableFeatures (Sobj, selection.method = "vst", nfeatures = 5000, verbose = T)
Sobj <- ScaleData     (Sobj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T) 
Sobj <- RunPCA        (Sobj, npcs = 20, verbose = T)
Sobj <- RunHarmony    (Sobj, "batch")
Sobj <- RunUMAP       (Sobj, reduction  = "harmony", dims = 1:20)
Sobj <- FindNeighbors (Sobj, reduction  = "harmony", dims = 1:20)
Sobj <- FindClusters  (Sobj, resolution = c(0.05, 0.1, 0.2))

#saveRDS(Sobj, "Seurat_Objects/Bad_Sample_Exclusion/CellType_Sobjs/2_Final_Harmony/5000_Vargenes/Myeloid_V2.rds")


res <- "RNA_snn_res.0.2"

DimPlot(Sobj,
        reduction = "umap", 
        group.by  = "Type",
        #group.by  = res,
        #split.by  = "Type",
        pt.size = 1, 
        label = T,
        repel = T, 
        label.size = 4
        )


FeaturePlot(Sobj, "test1")

Idents(Sobj) <- res

#Sobj <- RenameIdents(Sobj, "0" = "HSC")
#Sobj <- RenameIdents(Sobj, "0" = "CD4_T", "1" = "CD8_T", "2" = "T-Eff", "3" = "NK", "4" = "B.plasma", "5" = "B.mzone")
#Sobj <- RenameIdents(Sobj, "1" = "Macrophage", "0" = "mono.DC", "2" = "Monocyte", "3" = "DC", "4" = "HSC", "5" = "mono.DC")
#Sobj <- RenameIdents(Sobj, "0" = "vasc_SM1", "1" = "vasc_SM2", "3" = "vasc_SM2", "2" = "Pericyte")
#Sobj <- RenameIdents(Sobj, "0" = "Lymph_EC", "1" = "Lymph_EC2")
#Sobj <- RenameIdents(Sobj, "0" = "Vein", "1" = "Capillary", "2" = "Artery")
#Sobj <- RenameIdents(Sobj, "0" = "Matrix_2", "1" = "Lipo_F", "2" = "Matrix_1", "3" = "Matrix_1", "4" = "Vasc_F", "5" = "Chondrocyte")
#Sobj <- RenameIdents(Sobj, "0" = "Adipocyte", "1" = "Adipocyte")
#Sobj <- RenameIdents(Sobj, "0" = "lun_1", "1" = "lun_2", "2" = "lun_3", "3" = "lun_trans", "5" = "lun_cycling", "4" = "lun_4", "6" = "lun_6")
Sobj <- RenameIdents(Sobj, "0" = "bas_1", "1" = "bas_2", "2" = "bas_3", "3" = "bas_trans")
#Sobj <- RenameIdents(Sobj, "0" = "lup_1", "1" = "lup_2", "2" = "lup_3", "3" = "lup_4", "4" = "lup_trans", "5" = "lup_cycling")


Sobj$Subcluster <- Idents(Sobj)

MetaData_list[[(unique(Sobj$CellType))]] <- select(Sobj@meta.data, Subcluster)

Sobj@reductions %>% saveRDS(paste0("2021_ReRun/Reductions/", unique(Sobj$CellType), ".rds"))

Bind <- do.call(rbind.data.frame, MetaData_list)



terms <- Searchable_gmt %>% filter(str_detect(Term, "BIOCARTA_IL6")) %>% pull(Term) %>% unique()


terms <- grnboost_tumor_responders %>% mutate(perc = percent_rank(importance)) %>% filter(target %notin% remove_genes) %>% top_n(200, perc) %>% pull(target)


for (i in 1:length(terms)) {

Sobj <- AddModuleScore  (Sobj, assay = "RNA", 
                         list(query), 
                         nbin = 25, 
                         name = "ESRRG_Sig")
}

FeaturePlot(Sobj, features = c("SPP1"), cols = viridis(n = 100, option = "A"), pt.size = 0.5, order = T) #+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
FeaturePlot(Sobj, features = c("cxcl13_tumor_cluster_signature1", "CXCL2", "KIT", "ACTA2"), cols = viridis(n = 100, option = "A"), pt.size = 0.5, order = T) #+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
FeaturePlot(Sobj, features = c("BIOCARTA_IL6_PATHWAY1", "BIOCARTA_CDMAC_PATHWAY1", "BIOCARTA_HDAC_PATHWAY1", "REACTOME_RUNX3_REGULATES_NOTCH_SIGNALING1", "REACTOME_SMOOTH_MUSCLE_CONTRACTION1"), cols = viridis(n = 100, option = "A"), ncol = 2, pt.size = 0.5, order = T) #+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())


VlnPlot(Sobj, features = "gpam1", group.by = "CellType", split.by = "Type", pt.size = 0)


"REACTOME_CELL_JUNCTION_ORGANIZATION1"
"HALLMARK_APICAL_JUNCTION1"
"KEGG_FOCAL_ADHESION1"
"REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY1"
"REACTOME_SMOOTH_MUSCLE_CONTRACTION1"
"REACTOME_RUNX3_REGULATES_NOTCH_SIGNALING1"

for (i in 1:length(levels(Marks$cluster))) {
        
        Sobj <- AddModuleScore  (Sobj, assay = "RNA", 
                 list(top_n(filter(Marks, 
                                    cluster == levels(Marks$cluster)[i]), 
                             100, avg_log2FC)$gene), 
                 nbin = 25, 
                 name = levels(Marks$cluster)[i])
        }


FeaturePlot(Sobj, 
            features = str_replace_all(paste0(levels(Marks$cluster), "1"), pattern = "-", replacement = "."),
            cols = viridis(n = 100, option = "A"), pt.size = 1, order = T, ncol = 3)







Plot_list    <- vector("list", length = 5*5)

theta  <- rep(c(0, 1, 2, 3, 4, 5), 5)
lambda <- c(rep(1, 6),rep(2, 6),rep(3, 6),rep(4, 6),rep(5, 6))

for (i in 1:length(theta)) {
                
        Sobj <- RunHarmony    (Sobj, "Batch", theta = theta[i], lambda = lambda[i])
        Sobj <- RunUMAP       (Sobj, reduction  = "harmony", dims = 1:50)
        Sobj <- FindNeighbors (Sobj, reduction  = "harmony", dims = 1:50)
        Sobj <- FindClusters  (Sobj, resolution = c(0.1))

        res <- "RNA_snn_res.0.1"

        Plot_list[[i]] <- 
                DimPlot(Sobj,
                reduction = "umap", 
                group.by  = res,
                pt.size = 1, 
                label = T,
                repel = T, 
                label.size = 4 
                ) + 
                ggtitle(paste0("Theta", theta[i], "_", "Lambda", lambda[i]))
        
        Sobj@reductions %>% saveRDS(paste0("Basal_Reductions/Reductions_Theta", theta[i], "_", "Lambda", lambda[i], ".rds"))
        
        }

Plot_list %>% saveRDS(paste0("Basal_Reductions/Reductions_Plot_List.rds"))


i = 0

patchwork::wrap_plots(Plot_list[[1+i]],  Plot_list[[2+i]],  Plot_list[[3+i]], 
                      Plot_list[[4+i]],  Plot_list[[5+i]],  Plot_list[[6+i]],
                      Plot_list[[7+i]],  Plot_list[[8+i]],  Plot_list[[9+i]],
                      Plot_list[[10+i]], Plot_list[[11+i]], Plot_list[[12+i]],
                      Plot_list[[13+i]], Plot_list[[14+i]], Plot_list[[15+i]],
                      Plot_list[[16+i]], Plot_list[[17+i]], Plot_list[[18+i]]) + 
        patchwork::plot_layout(ncol = 6)







meta %>% select(Sample, Age, Type, Subcluster) %>% group_by(Type, Subcluster) %>% 
        summarise(mean_age = mean(Age))




Idents(Sobj) <- "Subcluster"
Sobj$Subcluster <- Idents(Sobj)


prop.table(table(Sobj$Subcluster, Sobj$Sample), margin = 2) %>% 
        as.data.frame() %>% left_join(select(Sobj@meta.data, Sample, Menopause, Type), by = c("Var2" = "Sample")) %>% 
        distinct() %>% 
        ggplot(aes(x = Var1, y = Freq, fill = Type)) + 
        geom_boxplot() + 
        scale_y_sqrt() + 
        scale_fill_manual(values = GID_cols) +
        style

