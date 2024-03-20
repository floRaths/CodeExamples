library(scDblFinder)
library(DropletUtils)
library(harmony)
library(Seurat)
library(tidyverse)
library(viridis)
library(readxl)


CellRanger_to_Sobj  <- function() {
  
  ### Make Seurat Object List
  rawdata_list <- vector("list", length = length(samples))                 # this will contain our count matrices
  Sobj_list    <- vector("list", length = length(samples))                 # our Seurat objects will be generated in here
  
  
  ### following loop reads in rawdata from cellranger output folders according to the sample name and then creates an Sobj in the list
  
  for (i in seq_along(samples)) {
    rawdata_list[[i]] <- Read10X            (data.dir = file.path("CellRanger_Out/RNA/2021_ReRun_CellR_v6/", samples[i] , "/outs/filtered_feature_bc_matrix/"))
    Sobj_list   [[i]] <- CreateSeuratObject (rawdata_list[[i]], project = samples[i], min.cells = 1, min.features = 1)
  }
  
  names(Sobj_list) <- samples # assign the correct names to the Sobj_List levels
  
  # rename cells using object names as prefix
  for (i in names(Sobj_list)) {
    Sobj_list[[i]] <- RenameCells(Sobj_list[[i]], add.cell.id = i)
  }
  
  # merge into one Seurat Object
  Sobj <<- merge(x = Sobj_list[[1]], y = Sobj_list[-1])
  
}
AssignClusterIdents <- function() {
  
  # make assignment of clusters to celltypes by using the best fit of marker genes
  CellType_Markers_Final <- read_excel("utilities/2021_ReAnalysis/CellType_Markers_Final.xlsx")
  Avg.Expr <- AverageExpression(Sobj, assays = "RNA", slot = "data", features = CellType_Markers_Final$Marker, group.by = res)
  
  
  Assign <- Avg.Expr$RNA %>% as.data.frame() %>% 
    rownames_to_column("Marker") %>% 
    pivot_longer(-1, names_to = res) %>% 
    left_join(CellType_Markers_Final) %>% 
    group_by(Marker) %>% top_n(1, value) %>% 
    group_by(.data[[res]], CellType) %>% 
    summarise(n = n()) %>% group_by(CellType) %>% 
    filter(n > 1) %>% select(-n)
  
  
  # create the new labels
  Idents(Sobj) <- res
  newIdent <- Assign$CellType
  names(newIdent) <- Assign$RNA_snn_res.0.1
  Sobj <- RenameIdents(object = Sobj, newIdent)
  Sobj$CellType <- Idents(Sobj)
  
}
DoubletDetection    <- function() {
  # prepare cluster assignment for merging with sce data
  cluster.assign <- Sobj@meta.data %>% 
    rownames_to_column("Barcode") %>% 
    select(Barcode, Sample = "orig.ident", CellType) %>% 
    mutate(Barcode = substr(Barcode, nchar(Barcode)-17, 40))
  
  #### load 10x data as sce object
  
  # create vector with file destinations
  p = vector()
  for (i in 1:length(samples)) {
    p[i] <- paste0("CellRanger_Out/RNA/2021_ReRun_CellR_v6/", samples[i], "/outs/filtered_feature_bc_matrix/")
  }
  
  sce <- read10xCounts(samples = p, sample.names = samples.short)
  
  # merge sce cellnames (ordering) with cluster assignemnt
  order <- sce@colData %>% as.data.frame() %>% 
    rownames_to_column("Cell") %>% 
    select(Cell, Sample, Barcode) %>% 
    left_join(cluster.assign)
  
  # add metadata column
  sce@colData$CellType <- order$CellType
  
  # run doublet detection
  sce <- scDblFinder(sce, samples = "Sample", clusters = "CellType")
  
  # add doublet detection results to Sobj
  meta <- sce@colData %>% as_tibble() %>% unite(x, Sample, Barcode) %>% column_to_rownames("x")
  Sobj <<- AddMetaData(Sobj, meta)
  
}
SeuratIntegration   <- function() {
  
  # now we integrate ("batch correct") the samples
  ####################################################################################################################################
  
  
  # Some further standard processing
  for (i in 1:length(Sobj_list)) {
    Sobj_list[[i]] <- NormalizeData        (Sobj_list[[i]], verbose = T)
    Sobj_list[[i]] <- FindVariableFeatures (Sobj_list[[i]], selection.method = "vst", nfeatures = 4000, verbose = T)
  }
  
  
  ### Initiating integration procedure
  
  Sobj_anchors    <- FindIntegrationAnchors (object.list = Sobj_list,    dims = 1:50)
  saveRDS(Sobj_anchors, "/home/rathsf/scNuclei-Libraries/Analysis/Sobj_Anchors.rds") # don't have to save each intermediate, but can be nice if things break down
  
  Sobj <- IntegrateData          (anchorset   = Sobj_anchors, dims = 1:50)
  saveRDS(Sobj, "/home/rathsf/scNuclei-Libraries/Analysis/Sobj.rds")
  
  rm(Sobj_anchors, Sobj_list)
  
  DefaultAssay (Sobj) <- "integrated"  # NOTE: there are 2 assays in the resulting object now, RNA and integrated... the integrated assay is ONLY for clustering and visualization purposes. Don't do any DE-testing or other analytical steps on this assay.
  Sobj <- ScaleData  (Sobj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T)
  saveRDS(Sobj, "/home/rathsf/scNuclei-Libraries/Analysis/Sobj_Scaled.rds")
  
  ### continuing with regressed dataset
  Sobj <- RunPCA        (Sobj, npcs = 50, verbose = FALSE)
  Sobj <- RunUMAP       (Sobj, reduction  = "pca", dims = 1:50)
  Sobj <- FindNeighbors (Sobj, reduction  = "pca", dims = 1:50)
  Sobj <- FindClusters  (Sobj, resolution = c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2)) # bit overkill on the resolution
  
  saveRDS(Sobj, "Seurat_Objects/LUMA_Subset_integrated.rds")  # thats you integrated dataset. If the clustering looks good, I eventually remove the integrated assay: Sobj[["integrated"]] <- NULL
  # it reduces memory impact and you don't run into any confusions with the two assays
}
GenerateMetaData    <- function() {
  ### assigning Patient_ID/replicate_ID to the Data
  Idents(Sobj) <- "orig.ident"
  Sobj@meta.data$Library <- Idents(Sobj)
  
  
  Idents(Sobj) <- "orig.ident"
  new.cluster.ids <- str_remove(levels(Sobj), "_2")
  names(new.cluster.ids) <- levels(Sobj)
  Sobj <- RenameIdents(Sobj, new.cluster.ids)
  Sobj@meta.data$Sample <- Idents(Sobj)
  
  
  ### assigning Treatment groups to the Data
  Idents(Sobj) <- "Sample"
  new.cluster.ids <- substring(levels(Sobj), 1, 2) #only keep CF or TM
  names(new.cluster.ids) <- levels(Sobj)
  Sobj <- RenameIdents(Sobj, new.cluster.ids)
  Sobj@meta.data$Type <- Idents(Sobj)
  
  
  ## assigning broader Groups on 0.1 resolution
  Idents(Sobj) <- "CellType"
  new.cluster.ids <- levels(Sobj)
  new.cluster.ids <- recode(new.cluster.ids, 
                            "LUM_HR-pos" = "Epithelial",
                            "LUM_HR-neg" = "Epithelial",
                            "Basal"      = "Epithelial",
                            "Fibroblast" = "Stroma",
                            "Adipocyte"  = "Stroma",
                            "Blood_EC"   = "Vascular",
                            "Lymph_EC"   = "Vascular",
                            "Vasc.Acc."  = "Vascular",
                            "Myeloid"    = "Immune",
                            "Lymphoid"   = "Immune"
                            )
  
  names(new.cluster.ids) <- levels(Sobj)
  Sobj <- RenameIdents(Sobj, new.cluster.ids)
  Sobj@meta.data$Group <- Idents(Sobj)
  
  
  
  LibraryMetrics <- read_excel("utilities/2021_ReAnalysis/LibraryMetrics.xlsx", skip = 0, sheet = 1) %>% filter(Sample %in% Sobj$orig.ident)
  levels <- LibraryMetrics$Sample
  
  # Batch
  Idents(Sobj) <- "orig.ident"
  newIdent <- LibraryMetrics$Batch
  names(newIdent) <- LibraryMetrics$Sample
  Sobj <- RenameIdents(object = Sobj, newIdent)
  Sobj$Batch <- Idents(Sobj)
  
  # PrepMethod
  Idents(Sobj) <- "orig.ident"
  newIdent <- LibraryMetrics$PrepMethod
  names(newIdent) <- LibraryMetrics$Sample
  Sobj <- RenameIdents(object = Sobj, newIdent)
  Sobj$Prep <- Idents(Sobj)
  
  # cDNA Avg
  Idents(Sobj) <- "orig.ident"
  newIdent <- LibraryMetrics$`Lib avg.`
  names(newIdent) <- LibraryMetrics$Sample
  Sobj <- RenameIdents(object = Sobj, newIdent)
  Sobj$cDNA_Avg <- Idents(Sobj)
}


### This workflow automates the processing of multiple 10x runs into one integrated dataset
setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

### Preparation of our file lists that will contain our seurat objects
samples <- list.files("./CellRanger_Out/RNA/2021_ReRun_CellR_v6/")   # folder structure should be:  Cellranger_Out/Sample_Name/filtered_feature_bc_matrix
samples.short <- samples

levels <- c("LUM_HR-pos", "LUM_HR-neg", "Basal", "Fibroblast", "Adipocyte", "Blood_EC", "Lymph_EC", "Vasc.Acc.", "Myeloid", "Lymphoid")
DietSeurat(counts = T, data = T, scale.data = F, assays = "RNA", dimreducs = c("pca", "harmony", "umap"))
### read in 10X data and create a merged Seurat object with unique barcodes
CellRanger_to_Sobj()




Sobj <- readRDS("2021_ReRun/Seurat_Objects/Sobj_Final_Scaled_MetaUpdt.rds")

### standard preprocessing for cluster detection
DefaultAssay (Sobj) <- "RNA"

Sobj <- NormalizeData        (Sobj, verbose = T)
Sobj <- FindVariableFeatures (Sobj, selection.method = "vst", nfeatures = 5000, verbose = T)
Sobj <- ScaleData            (Sobj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T) 
Sobj <- RunPCA               (Sobj, npcs = 50, verbose = T)
Sobj <- RunHarmony           (Sobj, "Batch")
Sobj <- RunUMAP              (Sobj, reduction  = "harmony", dims = 1:50)
Sobj <- FindNeighbors        (Sobj, reduction  = "harmony", dims = 1:50)
Sobj <- FindClusters         (Sobj, resolution = c(0.05, 0.1, 0.2, 0.5))

Sobj %>% DietSeurat(counts = T, data = T, scale.data = F, assays = "RNA", dimreducs = c("pca", "harmony", "umap")) %>% saveRDS("Seurat_Objects/Bad_Sample_Exclusion/2_ManualCell_Removal/Sobj_Manual_Cleaned_Scaled_Harmony.rds")


# identify clusters and then rename
DimPlot(Sobj, 
        reduction = "umap", 
        #group.by  = "CellType",
        #group.by  = res,
        #split.by  = "Type",
        pt.size = 1, 
        label = T, 
        raster = F, shuffle = T,
        repel = T, 
        label.size = 4
        )

# choose the resolution that fits the data best
res <- "RNA_snn_res.0.1"


### Assign cluster identities based on peviously determined markers
AssignClusterIdents()

# confirm that everything looks good
VlnPlot(Sobj, features = c("ANKRD30A", "ELF5", "TP63", "VWF", "FLT4", "PLIN1", "COL6A3", "IL7R", "CD163", "RGS6", "CD69"), 
        pt.size = 0, group.by = "CellType")
### Now we can use this assignment to improve doublet detection in the next step


DoubletDetection()

# Quick check
DimPlot(Sobj, group.by = "scDblFinder.class")
VlnPlot(Sobj, group.by = "scDblFinder.class", features = "nCount_RNA", split.by = "orig.ident")




### Standard preprocessing and filtering 
# mito genes
Sobj[["percent.mt"]] <- PercentageFeatureSet(Sobj, pattern = "^MT-")


# VlnPlot of UMI and mito distribution
VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Sobj, features = c("nCount_RNA"), ncol = 1, group.by = res) + ylim(0, 5000)


# filtering out high mito and high UMI cells
Sobj <- subset(Sobj, subset = percent.mt < 2.5 & nCount_RNA < 30000 & scDblFinder.class == "singlet") # adjust percent.mito and nCount RNA according to VlnPlots

saveRDS(Sobj_list, "Seurat_Objects/Feb15_All_Nuclei_Sobj-List.rds") #this is your pre-procesed list of Seurat objects which will go into the integration process next. You can also go on and do further steps (normalizing, scaling, clustering) on them individually



# Secondary processing
Modules_List <- readRDS("~/Box/Knott_Lab/Flo/Projects/Organoids/Analysis/2019.05.28_Estrogen_24-48h/utilities/Modules_List.rds")
Sobj <- AddModuleScore  (Sobj, assay = "RNA", list((top_n(Modules_List$LUMA, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUMA")
Sobj <- AddModuleScore  (Sobj, assay = "RNA", list((top_n(Modules_List$LUPR, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "LUPR")
Sobj <- AddModuleScore  (Sobj, assay = "RNA", list((top_n(Modules_List$MASC, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "MASC")
Sobj <- AddModuleScore  (Sobj, assay = "RNA", list((top_n(Modules_List$STRM, 50, `Average log fold-change`)$Symbol)), nbin = 25, name = "STRM")
rm(Modules_List)



s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Sobj <- CellCycleScoring(Sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
rm(s.genes, g2m.genes)



GenerateMetaData()



library(harmony)
library(Seurat)
library(tidyverse)

setwd("~/scNuclei-Libraries/Analysis/2021_ReRun/")
Sobj <- readRDS("Seurat_Objects/Bad_Sample_Exclusion/2_ManualCell_Removal/Sobj_Final_Scaled.rds")

#Idents(Sobj) <- "CellType"
Idents(Sobj) <- "Group"

n.var <- 5000

for (i in 1:length(levels(Sobj))) {
  
  DefaultAssay (Sobj) <- "RNA"
  
  subset <- subset(Sobj, idents = levels(Sobj)[i])
  
  subset <- NormalizeData (subset, verbose = T)
  subset <- FindVariableFeatures (subset, selection.method = "vst", nfeatures = n.var, verbose = T)
  subset <- ScaleData     (subset, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = T) 
  subset <- RunPCA        (subset, npcs = 50, verbose = T)
  subset <- RunHarmony    (subset, "SampleID")
  subset <- RunUMAP       (subset, reduction  = "harmony", dims = 1:50)
  subset <- FindNeighbors (subset, reduction  = "harmony", dims = 1:50)
  subset <- FindClusters  (subset, resolution = c(0.05, 0.1, 0.2, 0.5, 1))
  
  subset %>% 
    DietSeurat(counts = T, data = T, scale.data = F, assays = "RNA", dimreducs = c("pca", "harmony", "umap")) %>% 
    saveRDS(paste0("Seurat_Objects/Bad_Sample_Exclusion/CellType_Sobjs/2_Final_Harmony/", n.var , "_Vargenes/", levels(Sobj)[i], "_V2.rds"))

}


DimPlot(Sobj, group.by = "RNA_snn_res.0.1", label = T, split.by = "Type")
DimPlot(Sobj, split.by = "Sample", label = T, ncol = 5)
VlnPlot(Sobj, features = "nCount_RNA", group.by = "RNA_snn_res.0.2") + ylim(c(1, 10000))

FeaturePlot(Sobj, features = c("Cxcl13"), cols = viridis(n = 100, option = "A"), pt.size = 1, order = T)


Idents(Sobj) <- "CellType"

### generate percentage calculations

Perc_list <- vector("list", length = length(levels(Sobj)))
names(Perc_list) <- levels(Sobj)


for (i in 1:length(levels(Sobj))) {
  
  Seur <- subset(Sobj, idents = levels(Sobj)[i], downsample = 1000)
  
  Perc_list[[i]] <- 
    Seur@assays$RNA@counts[1:1000, ] %>% 
    as.data.frame() %>% 
    rownames_to_column("Gene") %>% 
    pivot_longer(colnames(Seur), names_to = "Cell", values_to = "count") %>% 
    left_join(rownames_to_column(select(Seur@meta.data, Type, CellType), "Cell")) %>% 
    mutate(expression = ifelse(count == 0, "absent", "present")) %>% 
    group_by(Gene, Type) %>% 
    mutate(n.Cells = n()) %>% 
    group_by(Gene, expression, Type, CellType, n.Cells) %>% 
    summarise(sum = n()) %>% 
    mutate(perc = (sum/n.Cells)*100)
}

Perc_list %>% saveRDS("Perc_list.rds")




Idents(Sobj) <- "Type"

TM <- Avg_TM$RNA %>% as.data.frame() %>% rownames_to_column("Gene") %>% pivot_longer(levels(Sobj$CellType), names_to = "CellType", values_to = "avg_TM")
CF <- Avg_CF$RNA %>% as.data.frame() %>% rownames_to_column("Gene") %>% pivot_longer(levels(Sobj$CellType), names_to = "CellType", values_to = "avg_CF")

avg_ex <- TM %>% left_join(CF)



percent <- 
Perc_list %>% 
  bind_rows() %>% 
  ungroup() %>% 
  filter(expression == "present") %>% 
  select(Gene, Type, CellType, perc) %>% pivot_wider(names_from = Type, values_from = perc) %>%  
  full_join(avg_ex) %>% 
  mutate(perc_TM = replace_na(TM, 0), perc_CF = replace_na(CF, 0)) %>% 
  select(-CF, -TM)

percent %>% saveRDS("2021_ReRun/output/Percent&Average_GeneExpression_CellType_Type.rds")









prop.table(table(Sobj$CellType, Sobj$Sample)) %>% 
  as.data.frame() %>% as_tibble() %>% mutate(Freq = Freq * 100) %>% 
  ggplot(aes(x = Var1, y = Freq, color = Var2)) + 
  geom_point() + 
  geom_label(aes(label = Var2))




Cohort_demographics %>% left_join(Ann) %>% 
  ggplot(aes(x = Type, y = Age, fill = Mens.stat)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(aes(size = Cells), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + 
  theme(text = element_text(size = 25))