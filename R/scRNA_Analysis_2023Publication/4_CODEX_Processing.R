library(readr)
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(pheatmap)
library(rdist)
library(ggplotify)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(viridis)
library(colorspace)

setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")

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
  
col <- pal_d3("category20")(20)
options(tibble.print_max = 150, tibble.print_min = 25)
`%notin%` <- Negate(`%in%`)


Anno <- 
  read_excel("2021_ReRun/CODEX/TMA_Region_Annotation.xlsx") %>% select("block_region" = 1, Type, Sample)
markers <- readRDS("2021_ReRun/CODEX/markers.rds")
GID_cols  <- c("#A6499B", "#FAA42F")
my_comparisons <- list(c("CF", "TM"))

style <- theme(text = element_text(family = "Lato", size = 25),  
               title = element_text(size = 15, face = "bold"), 
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_rect(fill = "white", color = "grey35", size = 2))


Nuclei_raw <- 
  read_csv("2021_ReRun/CODEX/SetA_B_nuclei.csv") %>% 
  left_join(Anno, by = "block_region") %>% 
  rename(Cell = "...1") %>% 
  separate(spot_name, into = c("block", "punch"), remove = F)

Membrane_raw <- 
  read_csv("2021_ReRun/CODEX/SetA_B_membrane.csv") %>% 
  left_join(Anno, by = "block_region") %>% 
  rename(Cell = "...1") %>% 
  separate(spot_name, into = c("block", "punch"), remove = F)

nucs_long <- 
  Nuclei_raw %>% 
  pivot_longer(all_of(markers), names_to = "marker") %>% 
  select(-X, -Y) %>% 
  mutate(data = "nuc") %>% 
  group_by(marker, block_region) %>% 
  mutate(scaled = scale(value)) %>% 
  ungroup()

mems_long <- 
  Membrane_raw %>% 
  pivot_longer(all_of(markers), names_to = "marker") %>% 
  select(-X, -Y) %>% 
  mutate(data = "mem") %>% 
  group_by(marker, block_region) %>% 
  mutate(scaled = scale(value)) %>% 
  ungroup()



Sobj_All <- readRDS("2021_ReRun/CODEX/CellTyping/S6.2_Sobj_artifacts_cleared_More.rds")


save_plot(p, "PI3K_gprofiler", w = 16, h = 9, scale = 1, svg = T)

### pre process the individual codex runs and prepare for merging
pre_process {
  
  ### make sobjs of each block for filtering                    ###### S1.1_Sobj_List_unfiltered
  cd_blocks <- unique(Nuclei_raw$block)
  sobj_list <- vector("list", length = length(cd_blocks))
  remove_list <- vector("list", length = length(cd_blocks))
  for (i in 1:length(cd_blocks)) {
    
    blocker = cd_blocks[i]
    
    rawA <- 
      nucs_long %>% filter(block == blocker) %>% 
      unite(marker, marker, data, sep = ".") %>% 
      select(Cell, marker, value,) %>% 
      pivot_wider(names_from = marker, values_from = value)
    
    rawB <- 
      mems_long %>% filter(block == blocker) %>% 
      unite(marker, marker, data, sep = ".") %>% 
      select(Cell, marker, value,) %>% 
      pivot_wider(names_from = marker, values_from = value)
    
    raw <- rawA %>% left_join(rawB, by = "Cell") %>% column_to_rownames("Cell") %>% t()
    
    
    
    meta <- Nuclei_raw %>% select(Cell, Size, block_region, block, punch, Type, Sample, DAPI) %>% distinct() %>% column_to_rownames("Cell")
  
    Sobj <- CreateSeuratObject(counts = raw)
    Sobj <- AddMetaData(Sobj, meta)
    
    Sobj <- subset(Sobj, subset = `DAPI.nuc` > 0)
    
    DefaultAssay (Sobj) <- "RNA"
    dims = 1:20
    
    #Sobj <- NormalizeData        (Sobj, verbose = T)
    Sobj <- FindVariableFeatures (Sobj, selection.method = "vst", nfeatures = 74, verbose = T)
    Sobj <- ScaleData            (Sobj, features = rownames(Sobj), vars.to.regress = "DAPI") 
    Sobj <- RunPCA               (Sobj, npcs = 30, verbose = T, features = rownames(Sobj))
    Sobj <- RunHarmony           (Sobj, c("Sample"), theta = 10, lambda = 0.1, dims.use = dims)
    Sobj <- RunUMAP              (Sobj, reduction  = "harmony", dims = dims, n.neighbors = 50, min.dist = 0.01)
    Sobj <- FindNeighbors        (Sobj, reduction  = "harmony", dims = dims)
    Sobj <- FindClusters         (Sobj, resolution = c(0.1, 0.2, 0.5))
    
    
    res <- "RNA_snn_res.0.2"
    Idents(Sobj) <- res
    
    sobj_list[[i]] <- Sobj
  }
  
  ### make list of Dimplot and Dapi combos for reference        ###### S1.2_DimPlot_List_unfiltered
  dimplot_list <- vector("list", length = length(cd_blocks))
  for (i in 1:length(cd_blocks)) {
    
    b <- FeaturePlot(sobj_list[[i]], features = c("DAPI.nuc"), pt.size = 1, order = T, cols = viridis(n = 100, option = "A"), max.cutoff = 5000)
    a <- DimPlot(sobj_list[[i]], 
                 group.by = res,
                 pt.size = 1, label = T)
    
    dimplot_list[[i]] <- a + b
    
  }
  
  ### filter clusters that are low for everything               ###### S2.1_Sobj_List_filter.passed // S2.2_Sobj_List_filter.failed_(removed)
  ### (likely oversegmentation. check thoroughly at some point)
  for (i in 1:length(cd_blocks)) {
    
    #Sobj <- sobj_list[[i]]
    #dimplot_list[[i]]
    
    
    #### Determine clusters that have the highest average expression of any marker
    is_top <- 
      AverageExpression(sobj_list[[i]], slot = "counts")$RNA %>% 
      as.data.frame() %>% 
      rownames_to_column("marker") %>% 
      pivot_longer(levels(sobj_list[[i]])) %>% 
      group_by(marker) %>% 
      top_n(1,value) %>% 
      pull(name) %>% unique()
    
    ### annotation of who is on top for which marker
    top <- 
      AverageExpression(sobj_list[[i]], slot = "counts")$RNA %>% 
      as.data.frame() %>% 
      rownames_to_column("marker") %>% 
      pivot_longer(levels(sobj_list[[i]])) %>% 
      group_by(marker) %>% 
      top_n(1,value) %>% 
      select(-value) %>% 
      separate(marker, into = c("marker", "data")) %>% 
      select(-data) %>% 
      distinct()
    
    ### mark clusters for removal
    ### criteria: cluster is not among "top marker clusters" and ranked lowest for at least 2/3 of all markers
    remove <- 
      AverageExpression(sobj_list[[i]], slot = "counts")$RNA %>% 
      as.data.frame() %>% 
      rownames_to_column("marker") %>% 
      pivot_longer(levels(sobj_list[[i]])) %>% 
      group_by(marker) %>%
      top_n(-5, value) %>% 
      group_by(name) %>% 
      summarise(n = n()) %>% 
      mutate(category = ifelse(name %in% is_top, "include", "remove")) %>% 
      arrange(category) %>% 
      filter(n > 24) %>% 
      left_join(top) %>% 
      filter(category == "remove") %>% 
      pull(name)
    
    remove_list[[i]] <- subset(sobj_list[[i]], idents = remove, invert = F)
    sobj_list[[i]]   <- subset(sobj_list[[i]], idents = remove, invert = T)
    
  }
  
  ### keep a list of cells that passed filter                   ###### S2.3_Pre.Process_Meta.Data
  meta.data_list <- vector("list", length = length(cd_blocks))
  for (i in 1:length(cd_blocks)) {
    
    a <- sobj_list[[i]]@meta.data   %>% select(1:10, -Type, - Sample) %>% rownames_to_column("Cell") %>% mutate(pre_process = "passed") %>% left_join(Anno, by = "block_region")
    b <- remove_list[[i]]@meta.data %>% select(1:10, -Type, - Sample) %>% rownames_to_column("Cell") %>% mutate(pre_process = "failed") %>% left_join(Anno, by = "block_region")
    
    meta.data_list[[i]] <- bind_rows(a, b) 
    
  }
  
  rm(meta.data_list, remove_list, cd_blocks, dimplot_list, nucs_long, mems_long)
  
}
### build main dataset in multiple steps        ##### S3.1_Sobj_Merged_scaled&harmony_across_blocks_Annotated.rds // S3.2_Sobj_Merged_ovseg_cleared.rds                                              
build_main_dataset {
  
  ### merge single sobjs into one
  Sobj <- merge(sobj_list[[1]], sobj_list[2:8])
  ### add some useful metadata for later
  expand_meta.data {
  ### expand meta data which will help with plotting later
  meta_expand <- 
    Sobj@meta.data %>% rownames_to_column("Cell") %>% 
    select(Cell, block_region, block) %>% 
    left_join(Anno) %>% 
    separate(block_region, into = c("trash", "blk_region"), sep = "_TG", remove = F) %>% 
    unite(region_type, blk_region, Type, sep = ".", remove = F) %>% 
    select(-trash) %>% 
    unite(block_type, block, Type, sep = ".", remove = F) %>% column_to_rownames("Cell")
  
  Sobj <- AddMetaData(Sobj, meta_expand)
  rm(meta_expand, sobj_list)
  }
  ### standard preprocessing (scale across block, harmony on block)
  processSobj {
    
    DefaultAssay (Sobj) <- "RNA"
    dims = 1:20
    
    Sobj <- FindVariableFeatures (Sobj, selection.method = "vst", nfeatures = 74, verbose = T)
    Sobj <- ScaleData            (Sobj, features = rownames(Sobj), vars.to.regress = c("DAPI.nuc", "DAPI.mem"), split.by = "block") 
    Sobj <- RunPCA               (Sobj, npcs = 30, verbose = T, features = rownames(Sobj))
    Sobj <- RunHarmony           (Sobj, c("block"), theta = 10, lambda = 0.1, dims.use = dims)
    Sobj <- RunUMAP              (Sobj, reduction  = "harmony", dims = dims, n.neighbors = 50, min.dist = 0.01)
    Sobj <- FindNeighbors        (Sobj, reduction  = "harmony", dims = dims)
    Sobj <- FindClusters         (Sobj, resolution = c(0.1, 0.2))
    
    Sobj$DAPI.gmm <- NULL
    Sobj$RNA_snn_res.0.5 <- NULL
    } 
  ### idendify main clusters and determine further clusters to exclude
  find_main_cluster_identities {
    
    ### set approriate resolution as identity
    res <- "RNA_snn_res.0.2"
    Idents(Sobj) <- res
    
    DimPlot(Sobj, group.by = res,
            cols = col,
            pt.size = 0.25, label = T)
    
    ### average expressions help with celltyping
    Avg <- AverageExpression(Sobj, slot = "scale.data")
    
    Avg$RNA %>% 
      as.data.frame() %>% 
      rownames_to_column("marker") %>% 
      pivot_longer(levels(Sobj), names_to = "cluster", values_to = "avg") %>% 
      group_by(cluster) %>% 
      top_n(2, avg) %>% 
      arrange(cluster)
    
    ### Featureplots to confirm
    FeaturePlot(Sobj, slot = "scale.data", features = c("CD45.mem"), pt.size = 0.25, order = T, cols = viridis(n = 100, option = "A"))
    
    ### rename identities
    Sobj <- RenameIdents(Sobj, 
                         "0" = "Luminal",
                         "1" = "Stroma_Other",
                         "2" = "Basal",
                         "3" = "Endothelial",
                         "4" = "Fibroblast",
                         "5" = "Immune",
                         "6" = "Adipocyte",
                         "7" = "EXCL_ovseg_conn.tissue",
                         "8" = "EXCL_ovseg_bloodclots",
                         "9" = "ENAH",
                         "10" = "Lymph_EC",
                         "11" = "EXCL_ovseg_B2_reg3_specific",
                         "12" = "EXCL_ovseg_B1_reg4_specific",
                         "13" = "EXCL_ovseg_B4_reg15_specific"
    )
    
    
    ### these clusters will be removed because they don't contain actual cells
    EXCL <- levels(Sobj) %>% as_tibble() %>% filter(str_detect(value, "EXCL")) %>% pull(value)
    
    ### following steps help visualize location of EXCL clusters
    Sobj$Group <- "All"
    Sobj$CellType <- Idents(Sobj)
    
    DimPlot(Sobj, group.by = "CellType", 
            cells = sample(colnames(Sobj)),
            cols = col, 
            pt.size = 0.25, label = T)
    
    ### focus on suspicious clusters
    Sobj_All <- Sobj
    Sobj <- subset(Sobj_All, idents = EXCL)
    
    ### copy to celltyping file on cerberus
    cat('   "', levels(Sobj)[1],  '": "',  str_sub(col[1], 1,7), '",\n',
        '   "', levels(Sobj)[2],  '": "',  str_sub(col[2], 1,7), '",\n',
        '   "', levels(Sobj)[3],  '": "',  str_sub(col[3], 1,7), '",\n',
        '   "', levels(Sobj)[4],  '": "',  str_sub(col[4], 1,7), '",\n',
        '   "', levels(Sobj)[5],  '": "',  str_sub(col[5], 1,7), '",\n',
        '   "', levels(Sobj)[6],  '": "',  str_sub(col[6], 1,7), '",\n',
        '   "', levels(Sobj)[7],  '": "',  str_sub(col[7], 1,7), '",\n',
        '   "', levels(Sobj)[8],  '": "',  str_sub(col[8], 1,7), '",\n',
        '   "', levels(Sobj)[9],  '": "',  str_sub(col[9], 1,7), '",\n',
        '   "', levels(Sobj)[10], '": "',  str_sub(col[10], 1,7), '",\n',
        '   "', levels(Sobj)[11], '": "',  str_sub(col[11], 1,7), '"\n',
        sep = "")
    
    
    ### upload annotations to cerberus
    regions <- Sobj@meta.data$block_region %>% unique()
    meta <- rownames_to_column(Sobj@meta.data, "Cell") %>% select(Cell, Group, CellType)
    for (i in 1:length(regions)) {
      Nuclei_raw %>% select(1:3, block_region) %>% 
        inner_join(meta, by = "Cell") %>% filter(block_region == regions[i]) %>%  
        write_csv(paste0("2021_ReRun/CODEX/Cell_Filtering/Annotations/transgenderTMA/", regions[i], "/", regions[i], "_3_annotation.csv"))
    }
    
    "cd Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/CODEX/Cell_Filtering/Annotations"
    "rsync -v -r /Users/rathsf/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/CODEX/Cell_Filtering/Annotations/transgenderTMA rathsf@10.220.237.94:/common/ingn/CODEX_PROCESSED_DATA/"
    
    
    ### find region with high amounts of a certain subcluster
    table(Sobj$CellType, Sobj$block_region) %>% as.data.frame() %>% filter(Var1 == "SLPI_High") %>% arrange(Freq)
  
  }
  ### now we remove the excl cells
  remove_excl_clusters {
    
    Sobj <- Sobj_All
    
    # overview before
    DimPlot(Sobj, group.by = "CellType",
          cols = col,
          pt.size = 0.25, label = T)
  
    ### remove the excluded oversegmentation
    ### we assume that all remaining datapoints are cells and we create a normalized dataset
    Sobj <- subset(Sobj, idents = EXCL, invert = T)
    
    # overview after
    DimPlot(Sobj, group.by = "CellType",
            cols = col,
            pt.size = 0.25, label = T)
  }
  
}  

scale_data_Robust.Scaler {
  
  version = "S6.2"
  
  ### send filtered rawdata tables for scaling
  Nuclei_raw %>% 
    filter(Cell %in% colnames(Sobj)) %>% 
    select(-spot_name) %>% 
    write_tsv(paste0("2021_ReRun/CODEX/CellTyping/Normalization/input_for_scaling/Nuclei_", version, "_raw.tsv.gz"))
  
  Membrane_raw %>% 
    filter(Cell %in% colnames(Sobj)) %>% 
    select(-spot_name) %>% 
    write_tsv(paste0("2021_ReRun/CODEX/CellTyping/Normalization/input_for_scaling/Membrane_", version, "_raw.tsv.gz"))
  
  
  ### shell commands for scaling
  "cd ~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/CODEX/CellTyping/Normalization"
  "conda activate codex"
  paste0("python3 analyze_codex.py input_for_scaling/Nuclei_",   version, "_raw.tsv.gz output_scaled Nuclei_",   version)
  paste0("python3 analyze_codex.py input_for_scaling/Membrane_", version, "_raw.tsv.gz output_scaled Membrane_", version)
  
  
  scale.data_nuc <-
    read_delim(paste0("2021_ReRun/CODEX/CellTyping/Normalization/output_scaled/scaled_data_Nuclei_", version, ".tsv.gz"), 
               delim = "\t", escape_double = FALSE, 
               trim_ws = TRUE)
  
  scale.data_mem <-
    read_delim(paste0("2021_ReRun/CODEX/CellTyping/Normalization/output_scaled/scaled_data_Membrane_", version, ".tsv.gz"), 
               delim = "\t", escape_double = FALSE, 
               trim_ws = TRUE)
  
  
  
  
}

processSobj {
  
  DefaultAssay (Sobj) <- "RNA"
  dims = 1:20
  
  Sobj <- FindVariableFeatures (Sobj, selection.method = "vst", nfeatures = 74, verbose = T)
  Sobj <- ScaleData            (Sobj, features = rownames(Sobj), vars.to.regress = c("DAPI.nuc", "DAPI.mem"), split.by = "block") 
  Sobj <- RunPCA               (Sobj, npcs = 30, verbose = T, features = rownames(Sobj))
  Sobj <- RunHarmony           (Sobj, c("block"), theta = 5, lambda = 1, dims.use = dims)
  Sobj <- RunUMAP              (Sobj, reduction  = "harmony", dims = dims, n.neighbors = 50, min.dist = 0.01)
  Sobj <- FindNeighbors        (Sobj, reduction  = "harmony", dims = dims)
  Sobj <- FindClusters         (Sobj, resolution = c(0.1, 0.2))
  
  Sobj$DAPI.gmm <- NULL
  Sobj$RNA_snn_res.0.5 <- NULL
  
  } 



### the workhorse
operate_cluster_identities {
  
  ### set approriate resolution as identity
  res <- "RNA_snn_res.0.1"
  Idents(Sobj_All) <- res
  
  DimPlot(Sobj, group.by = res,
          cols = col,
          pt.size = 0.25, label = T)
  
  ### average expressions help with celltyping
  Avg <- AverageExpression(Sobj, slot = "scale.data")
  
  Avg$RNA %>% 
    as.data.frame() %>% 
    rownames_to_column("marker") %>% 
    pivot_longer(levels(Sobj), names_to = "cluster", values_to = "avg") %>% 
    group_by(cluster) %>% 
    top_n(2, avg) %>% 
    arrange(cluster)
  
  marks <- FindAllMarkers(Sobj, slot = "scale.data", max.cells.per.ident = 10000)
  marks %>% group_by(cluster) %>% top_n(6, avg_diff)
  
  ### Featureplots to confirm
  FeaturePlot(Sobj, slot = "scale.data", features = c("AR.nuc"), pt.size = 0.25, order = T, cols = viridis(n = 100, option = "A")) + theme(text = element_text(size = 20), axis.text = element_blank(),  axis.ticks = element_blank())
  
  ### rename identities
  Sobj <- RenameIdents(Sobj, 
                       "0" = "Fibro_ARhi",
                       "1" = "Fibro_LAMBhi",
                       "2" = "Fibr_FN1hi",
                       "3" = "Fibr_Immu",
                       "4" = "Fibr_Epi",
                       "5" = "Fibr_Adip",
                       "6" = "Other",
                       "7" = "Other",
                       "8" = "Stroma_Other",
                       "9" = "Stroma_Other",
                       "10" = "Stroma_Other",
                       "11" = "",
                       "12" = "",
                       "13" = ""
                       )
  
  
  ### these clusters will be removed because they don't contain actual cells
  EXCL <- levels(Sobj) %>% as_tibble() %>% filter(str_detect(value, "EXCL")) %>% pull(value)
  
  ### following steps help visualize location of EXCL clusters
  Sobj$Group <- "All"
  Sobj$CellType <- Idents(Sobj)
  
  DimPlot(Sobj, group.by = "CellType", #split.by = "Type",
          cells = sample(colnames(Sobj)),
          cols = col, 
          pt.size = 0.25, label = T, raster = F) + 
    theme(text = element_text(size = 20), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + ggtitle("Basal Cells")
  
  
  ### copy to celltyping file on cerberus
  cat('   "', levels(Sobj)[1],  '": "',  str_sub(col[1], 1,7), '",\n',
      '   "', levels(Sobj)[2],  '": "',  str_sub(col[2], 1,7), '",\n',
      '   "', levels(Sobj)[3],  '": "',  str_sub(col[3], 1,7), '",\n',
      '   "', levels(Sobj)[4],  '": "',  str_sub(col[4], 1,7), '",\n',
      '   "', levels(Sobj)[5],  '": "',  str_sub(col[5], 1,7), '",\n',
      '   "', levels(Sobj)[6],  '": "',  str_sub(col[6], 1,7), '",\n',
      '   "', levels(Sobj)[7],  '": "',  str_sub(col[7], 1,7), '",\n',
      '   "', levels(Sobj)[8],  '": "',  str_sub(col[8], 1,7), '",\n',
      '   "', levels(Sobj)[9],  '": "',  str_sub(col[9], 1,7), '",\n',
      '   "', levels(Sobj)[10], '": "',  str_sub(col[10], 1,7), '",\n',
      '   "', levels(Sobj)[11], '": "',  str_sub(col[11], 1,7), '"\n',
      sep = "")
  
  
  ### upload annotations to cerberus
  regions <- Sobj@meta.data$block_region %>% unique()
  meta <- rownames_to_column(Sobj@meta.data, "Cell") %>% select(Cell, Group, CellType)
  for (i in 1:length(regions)) {
    Nuclei_raw %>% select(1:3, block_region) %>% 
      inner_join(meta, by = "Cell") %>% filter(block_region == regions[i]) %>%  
      write_csv(paste0("2021_ReRun/CODEX/Cell_Filtering/Annotations/transgenderTMA/", regions[i], "/", regions[i], "_3_annotation.csv"))
  }
  
  "cd Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/CODEX/Cell_Filtering/Annotations"
  "rsync -v -r /Users/rathsf/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/CODEX/Cell_Filtering/Annotations/transgenderTMA rathsf@10.220.237.94:/common/ingn/CODEX_PROCESSED_DATA/"
  
  
  ### find region with high amounts of a certain subcluster
  table(Sobj$CellType, Sobj$block_region) %>% as.data.frame() %>% filter(Var1 == "Five") %>% arrange(Freq)
  table(Sobj$CellType, Sobj$block_region) %>% as.data.frame() %>% filter(Var1 == "Five") %>% arrange(Freq) %>% ggplot(aes(x = reorder(Var2, -Freq), y = Freq)) + geom_col()
  
}

### Plotting pie charts
circular_proportion {


Prop_RNA <- 
  prop.table(table(Sobj$CellType)) %>% 
  as.data.frame() %>% rename(CellType = "Var1")


 
  ggplot(data = Prop_RNA, aes(x = "", y = Freq, fill = CellType)) + 
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = round(Freq*100, 1)), 
            size = 5, 
            position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_d3(palette = "category20") +
  style +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(), 
        legend.position = "right")



 
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
  
}


nuclear_size_plot {

Sobj_All@meta.data %>% 
  group_by(blk_region, Type, CellClass) %>% 
  mutate(perc = percent_rank(Size)) %>% 
  filter(perc < 0.99) %>% 
  summarise(n = n(), avg.Size = mean(Size)) %>% 
  group_by(CellClass) %>% filter(n > 10) %>% 
  filter(CellClass %notin% c("Epi_Other", "ENAH")) %>% 
  ggplot(aes(x = CellClass, y = avg.Size, fill = Type)) + 
  scale_fill_manual(values = GID_cols) +
  #stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
  #geom_text_repel(aes(label = blk_region)) +
  geom_boxplot() + 
  #geom_point(position = position_dodge(width = 0.75)) +
  style + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))# +
  facet_wrap("CellClass", nrow = 1)
}


region_scatter {
  
  
  p <- 
    Sobj_All@meta.data %>% 
    filter(blk_region == "B9_reg13", 
           Group %in% c("Epithelial", "Immune")
           ) %>% 
    arrange(Group) %>% 
    mutate(label = ifelse(Group == "Epithelial", as.character(Group), as.character(CellType))) %>% 
    ggplot(aes(x = X, y = Y, color = label)) + 
    geom_point() + 
    scale_color_d3() +
    ylim(c(0, 7577)) +
    xlim(c(0, 7577)) +
    coord_fixed() + 
    style +
    theme(axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()
          )
  
}

vasc_umap {

Sobj <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/CODEX/CellTyping/S5.2_Endothelial_filtered_Annotated.rds")
levels(Sobj) <- c("Endo_Main", "Capillary", "Endo_SMA", "Endo_Immune", "Lymph_EC")
p <- DimPlot(Sobj, cols = col)
save_plot(p, "vascular_codex_umap", w = 15, h = 10, scale = 0.5, svg = T)
}



### Euclidean Distance Plots
distance_boxplots {
bind %>% filter(mins != Inf) %>% 
  group_by(CellType, block_region) %>% summarise(avg = mean(mins)) %>% 
  ggplot(aes(y = CellType, x = avg, fill = CellType)) + 
  geom_boxplot(outlier.size = 0.1) + 
  geom_point(alpha = 0.5, 
             color = "grey20", 
             size = 3, 
             pch = 21,
             position = position_jitterdodge(dodge.width = 0.1, 
                                             jitter.width = 0.35)) +
  scale_x_log10() +
  scale_fill_d3() +
    ggtitle("Distance to Epithelia") +
  guides(fill = guide_legend(nrow = 2)) +
  xlab("avg. distance to Epithelial cells per region") +
  style + theme(legend.position = "bottom", axis.title.y = element_blank())
  
  
  
  Sobj@meta.data %>% filter(Epi.Dist != Inf) %>% 
    group_by(Type, CellType, block_region) %>% 
    summarise(avg = mean(Epi.Dist)) %>% 
    ggplot(aes(y = CellType, x = avg, fill = Type)) + 
    geom_boxplot(outlier.size = 0.1) + 
    geom_point(alpha = 0.5, 
               color = "grey20", 
               size = 3, 
               pch = 21,
               position = position_jitterdodge(dodge.width = 0.1, 
                                               jitter.width = 0.35)) +
    scale_x_log10() +
    scale_fill_d3() +
    guides(fill = guide_legend(nrow = 2)) +
    xlab("avg. distance to Epithelial cells per region") +
    style + theme(legend.position = "bottom", axis.title.y = element_blank())
}

### Plotting intensities from robust scaled
scale.data_plots {
  
  Ann_Add <- Sobj_All@meta.data %>% rownames_to_column("Cell") %>% select(Cell, Group, CellClass, CellType, blk_region)
  data_rob.n <- scale.data_nuc %>% select(-...1, -contains("gmm"), -contains("Quant")) %>% left_join(Ann_Add)
  data_rob.m <- scale.data_mem %>% select(-...1, -contains("gmm"), -contains("Quant")) %>% left_join(Ann_Add)
  rm(scale.data_nuc, scale.data_mem)
  
  #data_rob.n <- scale.data_nuc %>% filter(Cell %in% colnames(Sobj)) %>% left_join(rownames_to_column(Sobj@meta.data, "Cell"))
  #data_rob.m <- scale.data_mem %>% filter(Cell %in% colnames(Sobj)) %>% left_join(rownames_to_column(Sobj@meta.data, "Cell"))
  
  #cellt <- "Luminal_HR+"
  cellt <- levels(Sobj)
  query <- c("AREG")
  
  data_rob.m %>% filter(Gene %in% query) %>% 
    ggplot(aes(x = CellType, y = (1+Robust.Scaled), fill = CellType)) + 
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
    #geom_point(#aes(color = block), alpha = 0.25, position = position_jitterdodge(dodge.width = 0.90, jitter.width = 0.01)) +
    scale_color_d3() +
    scale_fill_d3() +
    #scale_fill_manual(values = GID_cols) +
    scale_y_log10(limits = c(0.1, 15)
    ) +
    facet_wrap("Gene", nrow = 1, scales = "free") +
    style + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  
  
  
  
  
  cellt <- levels(Sobj)[1:3]
  query <- c("AREG")
  
  data_rob.m %>% filter(Gene %in% query, CellType %in% cellt) %>% 
    group_by(Gene, CellType) %>% 
    mutate(perc = percent_rank(Robust.Scaled)) %>% 
    filter(perc < 0.999 & perc > 0.01) %>% 
    ggplot(aes(x = CellType, y = (1+Robust.Scaled), fill = Type)) + 
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
    #geom_point(#aes(color = block), alpha = 0.25, position = position_jitterdodge(dodge.width = 0.90, jitter.width = 0.01)) +
    scale_color_d3() +
    #stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    scale_fill_manual(values = GID_cols) +
    scale_y_log10(#limits = c(0.1, 10)
    ) +
    facet_wrap("Gene", nrow = 1, scales = "free") +
    style + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  
  
  
  

  
  
  
  
  
  
  

  #data.m <- 
  #  data_rob.m %>% 
  #  group_by(CellClass, Gene, block_region) %>% 
  #  mutate(perc = percent_rank(Robust.Scaled)) %>% 
  #  mutate(nCt = n(), data = "membrane")
  
  #data.n <- 
  #  data_rob.n %>% 
  #  group_by(CellClass, Gene, block_region) %>% 
  #  mutate(perc = percent_rank(Robust.Scaled)) %>% 
  #  mutate(nCt = n(), data = "nuclear")
  
  
  data.m <- readRDS("2021_ReRun/CODEX/output/data.m.ready_for_plotting.rds")
  data.n <- readRDS("2021_ReRun/CODEX/output/data.n.ready_for_plotting.rds")
  
  
  celltype = c("Luminal", "Luminal_HR+", "Basal")
  
  data <- 
    data.n %>% 
    filter(CellClass %in% celltype, perc < 0.99) %>% 
    group_by(CellClass, nCt, Gene, block_region, blk_region, block, Type) %>% 
    summarise(Robust.Scaled = mean(Robust.Scaled))
  
  data <- 
    data.m %>% 
    filter(CellClass %in% celltype, perc < 0.99) %>% 
    group_by(CellClass, nCt, Gene, block_region, blk_region, block, Type) %>% 
    summarise(Robust.Scaled = mean(Robust.Scaled))
  
  

  query <- c("AZGP1")
  
  p <- data %>% filter(nCt > 10, Gene %in% query, #Robust.Scaled < 30,
                  #CellType %in% c("Fibro_LAMBhi"),
                  #blk_region %notin% c("B3_reg6", "B3_reg9")
                  ) %>% 
    ggplot(aes(x = Type, y = (1+Robust.Scaled), fill = Type)) + 
    geom_boxplot(width = 0.5, position = position_dodge(width = 0.9)) +
    geom_point(aes(color = block), position = position_jitter(width = 0.1), size = 2) +
    scale_color_d3() +
    #geom_text_repel(aes(label = blk_region)) +
    stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    scale_fill_manual(values = GID_cols) +
    scale_y_log10(#limits = c(1, 8)
    ) +
    #ggtitle(paste0(celltype)) +
    facet_wrap("CellClass", nrow = 1, scales = "fixed") +
    #facet_grid(cols = vars(Gene), rows = vars(CellType), scales = "free_y") +
    style + theme(axis.text.y = element_text(angle = 90))
  
  
  
  
  fibroblast_lamb_staining {
  
  p <- data %>% filter(nCt > 10, Gene %in% query, 
                  CellType %in% c("Fibro_LAMBhi"),
                  #blk_region %notin% c("B3_reg6", "B3_reg9")
  ) %>% 
    ggplot(aes(x = Type, y = (1+Robust.Scaled), fill = Type)) + 
    geom_boxplot(width = 0.5, position = position_dodge(width = 0.9)) +
    geom_point(aes(color = block),
               alpha = 0.75, 
               size = 2,
               position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1)) + 
    scale_color_d3() +
    #geom_text_repel(aes(label = blk_region)) +
    stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    scale_fill_manual(values = GID_cols) +
    scale_y_log10(limits = c(0.25, 5)
    ) +
    ggtitle("Fibro_LAMB1+") +
    facet_wrap("Gene", nrow = 1, scales = "free_x") +
    #facet_grid(rows = vars(Gene), cols = vars(CellType), scales = "free_y") +
    style 
  
  save_plot(p, "fibroblastLAMB_staining", w = 16, h = 10, scale = 1, svg = T)
  
  }
  
  
  
  
  
  query <- c("CD45", "TP63", "CD31", "ACTA2", "LYVE1")
  
  include <- table(Sobj$block_region) %>% as.data.frame() %>% arrange(-Freq) %>% filter(Freq >= 0) %>% pull(Var1)
  
  data_rob %>% filter(Gene %in% query, CellType %in% cellt, block_region %in% include) %>% 
    group_by(block_region, block, Type, CellType, Gene) %>% 
    summarise(Robust.Scaled = mean(Robust.Scaled)) %>% 
    filter(Robust.Scaled > -0.8) %>% 
    ggplot(aes(x = Gene, y = (1+Robust.Scaled), fill = CellType)) + 
    #geom_violin(scale = "width") +
    geom_boxplot(width = 0.5, position = position_dodge(width = 0.9), outlier.size = 0.1) +
    geom_point( #aes(color = block), 
      color = "grey20", 
      pch = 21, 
      alpha = 0.5, 
      position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2)) +
    scale_color_d3() +
    ylab("avg/region") +
    #stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    scale_fill_d3() +
    scale_y_log10(#limits = c(0.1, 20)
    ) + 
    guides(fill = guide_legend(nrow = 2)) +
    facet_wrap("Gene", nrow = 1, scales = "fr") +
    style + theme(legend.position = "bottom") 
  
  }

### the marker heatmap
marker_heatmap {
  
  data_rob <- scale.data_nuc %>% filter(Cell %in% colnames(Sobj)) %>% left_join(rownames_to_column(Sobj@meta.data, "Cell"))

  data <- data_rob %>% select(Cell, X, Y, Size, block, punch, Type, Sample, Gene, Old.Val, Robust.Scaled, blk_region, block_type, Group, CellType)

  cellt_avg <- data_rob.m %>% group_by(Gene, CellType) %>% summarise(avg = mean(Robust.Scaled))
  cellc_avg <- data_rob.m %>% group_by(Gene, CellClass) %>% summarise(avg = mean(Robust.Scaled))
  
  mat <- cellt_avg %>% pivot_wider(names_from = Gene, values_from = avg) %>% column_to_rownames("CellType")
  
  key_markers <- c("ESR1", "KRT8", "KRT23", "SLPI", "TP63", "PDGFRB", "CD45", "CD68",  "CD31", "ACTA2", "LYVE1", "PDPN", "PLIN1", "ENAH")
  cc_levels <- c("Luminal_HR+", "Luminal", "Basal", "Fibroblast",  "Immune", "Endothelial", "Lymph_EC", "Adipocyte", "ENAH", "Stroma_Other")
  
  key_markers <- c("KRT23", "KRT8", "TP63", "PDGFRB", "CD45", "CD68",  "CD31", "ACTA2", "LYVE1", "PDPN", "PLIN1", "ENAH")
  cc_levels <- c("Luminal_HR+", "Luminal", "Basal", "Fibroblast",  "Immune", "Macrophage", "Endothelial", "Lymph_EC", "Adipocyte", "ENAH", "Stroma_Other")
  
  key_markers <- c("CD68", "TP63", "KRT8", "KRT23", "ACTA2",  "CD31")
  cc_levels <- c("Luminal_HR+", "Luminal", "Basal", "Fibroblast",  "Immune", "Endothelial", "Lymph_EC", "Adipocyte", "ENAH", "Stroma_Other")
  
  key_markers <- c("CD68", "TP63", "KRT8", "KRT23", "ACTA2",  "CD31", "LYVE1")
  cc_levels <- c("Immune_Main", "Immune_Epi", "Immune_Endo", "Macrophage")
  
  
  
  p <- pheatmap(mat[cc_levels, key_markers], 
           cluster_rows = F,
           cluster_cols = F,
           scale = "column", 
           fontsize = 25,
           #main = "Marker Classification",
           color = colorRampPalette(brewer.pal(11, "BrBG"))(100))

}

### stacked columns of clusters in all regions
region_stack {
  
  # order via Epithelial content
  order <- 
    prop.table(table(Sobj_All$Group, Sobj$blk_region), margin = 2) %>% 
    as.data.frame() %>% 
    rename(blk_region = "Var2") %>% 
    left_join(Sobj@meta.data, by = "blk_region") %>% 
    as_tibble() %>% 
    filter(!is.na(block_region)) %>% filter(Var1 == "Epithelial") %>% arrange(-Freq) %>% pull(blk_region) %>% unique()
  
  
  prop.table(table(Sobj_All$CellClass, Sobj_All$blk_region), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellClass = "Var1", blk_region = "Var2") %>% 
    left_join(Sobj@meta.data, by = c("blk_region", "CellClass")) %>% 
    as_tibble() %>% 
    #filter(!is.na(blk_region)) %>% 
    ggplot(aes(x = factor(blk_region, levels = order), y = Freq, fill = CellClass)) + 
    geom_col(position = position_fill()) + 
    scale_fill_d3(palette = "category20") +
    style +
    ggtitle("CellType Proportion in all Regions") +
    ylab("Proportion") + xlab("region") +
    scale_y_continuous(n.breaks = 3) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "bottom")
  }

### stacked columns of clusters in all Samples
sample_stack {
  
  p <- prop.table(table(Sobj_All$CellClass, Sobj_All$Sample), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellClass = "Var1", Sample = "Var2") %>% 
    left_join(distinct(select(Sobj@meta.data, Sample, Type)), by = "Sample") %>% 
    as_tibble() %>% 
    ggplot(aes(x = Sample, y = Freq, fill = CellClass)) + 
    geom_col() + 
    #scale_fill_d3(palette = "category20") +
    scale_fill_manual(values = c(col[1:7], col[12], col[8:11])) +
    style +
    ggtitle("CellType Proportion across Samples") +
    ylab("Proportion") + xlab("Sample") +
    scale_y_continuous(n.breaks = 3) +
    theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust = 0.5), 
          #axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(), 
          legend.position = "right") +
    facet_wrap("Type", scales = "free")
}

### Nuclear to membrane plots
nuclear_membrane_ratio {
  
  cells <- WhichCells(Sobj_All, idents = "Fibroblast")
  join <- scale.data_nuc %>% select(Cell, Gene, Type, block, block_region, nuc = "Old.Val") %>% left_join(select(scale.data_mem, Cell, Gene, Type, block, block_region, mem = "Old.Val"))
  data <- join %>% filter(Cell %in% cells, Gene %in% c("AR", "ESR1", "PGR")) 
  
  p <- data %>% 
    mutate(ratio = nuc/mem) %>% 
    group_by(block_region) %>% 
    mutate(nCt = n()) %>% 
    filter(nCt > 1) %>% 
    group_by(Type, block, block_region, Gene) %>% 
    summarise(ratio = mean(ratio)) %>% 
    #filter(ratio < 2 & ratio > 0.9) %>% 
    #filter(Gene %in% c("AR", "ESR1", "PGR")) %>% 
    filter(Gene %in% c("AR")) %>% 
    ggplot(aes(x = Type, y = ratio, fill = Type)) + 
    geom_boxplot(outlier.alpha = 0, width = 0.5) + 
    geom_point(aes(color = block), position = position_jitter(width = 0.1), size = 2) +
    stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    scale_fill_manual(values = GID_cols) + 
    scale_color_d3() +
    facet_wrap("Gene", scales = "free") +
    ggtitle("HR-Rec. Nuclear to Membrane signal") +
    style
  
}

### Plotting proportions
proportion_plot {
  
  Idents(Sobj) <- "CellType"
  Sobj$CellType <- Idents(Sobj)
  
  meta <- Sobj@meta.data %>% unite(block_sample, block, Sample) %>% select(block_sample)
  Sobj <- AddMetaData(Sobj, meta)
  meta <- Sobj@meta.data %>% unite(sample_region, blk_region, Sample, sep = ".") %>% select(sample_region)
  Sobj <- AddMetaData(Sobj, meta)
  
  
  add_ann <- Sobj@meta.data %>% select(region_type, block_type, blk_region, block_region, Type, Sample, block) %>% distinct()
  
  include <- table(Sobj$block_region) %>% as.data.frame() %>% arrange(-Freq) %>% filter(Freq >= 20) %>% pull(Var1)
  #include <- table(Sobj$CellType, Sobj$block_region) %>% as.data.frame() %>% filter(Var1 == "Capillary") %>% filter(Freq > 50) %>% pull(Var2) %>% unique()
  my_comparisons <- list(c("CF", "TM"))
  
  prop.table(table(Sobj$CellType, Sobj$block_region), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellType = "Var1", block_region = "Var2") %>% 
    left_join(add_ann, by = "block_region") %>% 
    filter(block_region %in% include, CellType %in% levels(Sobj)[1:5]) %>% 
    #filter(Sample != "CF-3920", 
    #block_region != "211019_TransgenderTMA_SetA_TGB3_reg4"
    #       ) %>% 
    ggplot(aes(x = Type, y = Freq, fill = Type)) +
    geom_boxplot() +
    geom_point(aes(color = block),
               alpha = 0.75, 
               size = 3,
               position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1)) + 
    scale_fill_manual(values = GID_cols) +
    scale_color_d3() +
    #ggtitle("min 15 Macrophages") +
    stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    facet_wrap("CellType", 
               nrow = if_else(length(levels(Sobj)) > 8, 2, 1), 
               scales = "free") +
    style + theme(legend.position = "none")
  
  

  
  prop.table(table(Sobj$CellType, Sobj$block_sample), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellType = "Var1", block_sample = "Var2") %>% 
    separate(block_sample, into = c("block", "Sample"), sep = "_") %>% 
    left_join(distinct(select(add_ann, Sample, Type)), by = "Sample") %>% 
    filter(Sample != "CF-3920") %>% 
    ggplot(aes(x = Type, y = Freq, fill = Type)) +
    geom_boxplot() +
    geom_point(aes(color = Sample),
               alpha = 0.75, 
               size = 3,
               position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1)) + 
    scale_fill_manual(values = GID_cols) +
    scale_color_d3(palette = "category20") +
    stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) +
    facet_wrap("CellType", nrow = 1, scales = "free") +
    style
  
  
  
  prop.table(table(Sobj$CellType, Sobj$Sample), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellType = "Var1", Sample = "Var2") %>% 
    #separate(block_sample, into = c("block", "Sample"), sep = "_") %>% 
    left_join(distinct(select(add_ann, Sample, Type)), by = "Sample") %>% 
    #filter(Sample != "CF-3920") %>% 
    ggplot(aes(x = Type, y = Freq, fill = Type)) +
    geom_boxplot() +
    geom_point(aes(color = Sample),
               alpha = 0.75, 
               size = 3,
               position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1)) + 
    scale_fill_manual(values = GID_cols) +
    scale_color_d3(palette = "category20") +
    stat_compare_means(comparisons=my_comparisons, method = "t.test", paired = F) +
    facet_wrap("CellType", nrow = 1, scales = "free") +
    style
  
  
  
  
  
  ### stacked columns of clusters in all regions
  prop.table(table(Sobj$CellType, Sobj$block_region), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellType = "Var1", block_region = "Var2") %>% 
    left_join(add_ann, by = "block_region") %>% 
    as_tibble() %>% 
    #filter(block_region %in% include) %>% 
    ggplot(aes(x = reorder(blk_region, Freq), y = Freq, fill = CellType)) + 
    geom_col() + 
    scale_fill_d3(palette = "category20") +
    style +
    theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5))
  
  
  
  prop.table(table(Sobj$CellType, Sobj$Sample), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellType = "Var1", Sample = "Var2") %>% 
    separate(Sample, into = c("Type", "ID"), remove = F) %>% 
    #filter(block_punch %in% include) %>% 
    ggplot(aes(x = Type, y = Freq, fill = Type)) +
    geom_boxplot() +
    geom_point(aes(color = ID),
               alpha = 0.75, 
               size = 3,
               position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1)) + 
    scale_fill_manual(values = GID_cols) +
    #scale_color_d3() +
    stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    facet_wrap("CellType", ncol = 5, scales = "free") +
    style
  
  
  prop.table(table(Sobj$CellType, Sobj$block_type), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellType = "Var1", Sample = "Var2") %>% 
    separate(Sample, into = c("block", "Type"), remove = F) %>% 
    #filter(block_punch %in% include) %>% 
    ggplot(aes(x = Type, y = Freq, fill = Type)) +
    geom_boxplot() +
    geom_point(aes(color = block),
               alpha = 0.75, 
               size = 3,
               position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1)) + 
    scale_fill_manual(values = GID_cols) +
    #scale_color_d3() +
    stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    facet_wrap("CellType", ncol = 5, scales = "free") +
    style
  
}


fibroblast_proportion {
  
  
  Idents(Sobj_All) <- "CellClass"
  Sobj <- subset(Sobj_All, idents = "Fibroblast")
  
  Idents(Sobj) <- "CellType"
  Sobj$CellType <- Idents(Sobj)
  
  
  data <- AverageExpression(Sobj, slot = "scale.data", add.ident = "Type")$RNA %>% t() %>% as.data.frame() %>% select(contains(".mem"))
  
  rows <- c(rownames(data)[1:10])
  cols <- c("AR.mem", "LAMB1.mem", "LAMA2.mem", "FN1.mem", "KRT8.mem", "KRT23.mem", "CD45.mem")
  
  a <- pheatmap(data[rows, cols], 
           cluster_rows = F,
           cluster_cols = F,
           scale = "column", 
           fontsize = 25,
           main = "Fibroblast Staining",
           color = colorRampPalette(brewer.pal(11, "BrBG"))(100))
  



  meta <- Sobj@meta.data %>% unite(block_sample, block, Sample) %>% select(block_sample)
  Sobj <- AddMetaData(Sobj, meta)
  meta <- Sobj@meta.data %>% unite(sample_region, blk_region, Sample, sep = ".") %>% select(sample_region)
  Sobj <- AddMetaData(Sobj, meta)
  
  
  add_ann <- Sobj@meta.data %>% select(region_type, block_type, blk_region, block_region, Type, Sample, block) %>% distinct()
  
  include <- table(Sobj$block_region) %>% as.data.frame() %>% arrange(-Freq) %>% filter(Freq >= 10) %>% pull(Var1)
  #include <- table(Sobj$CellType, Sobj$block_region) %>% as.data.frame() %>% filter(Var1 == "Capillary") %>% filter(Freq > 50) %>% pull(Var2) %>% unique()
  my_comparisons <- list(c("CF", "TM"))
  
  
  prop.table(table(Sobj$CellType, Sobj$block_region), margin = 2) %>% 
    as.data.frame() %>% 
    rename(CellType = "Var1", block_region = "Var2") %>% 
    left_join(add_ann, by = "block_region") %>% 
    filter(block_region %in% include, CellType %in% levels(Sobj)[1:5]
           #CellType == "Macrophage"
    ) %>% 
    #filter(Sample != "CF-3920", 
    #block_region != "211019_TransgenderTMA_SetA_TGB3_reg4"
    #       ) %>% 
    ggplot(aes(x = Type, y = Freq, fill = Type)) +
    geom_boxplot() +
    
    scale_fill_manual(values = GID_cols) +
    scale_color_d3() +
    #ggtitle("min 15 Macrophages") +
    stat_compare_means(comparisons=my_comparisons, method = "wilcox", paired = F) +
    #facet_wrap("CellType", 
    #           nrow = if_else(length(levels(Sobj)) > 8, 2, 1), 
    #           scales = "free") +
    style + theme(legend.position = "none") +
    facet_wrap("CellType")
  
    #coord_flip()
  
  
  p <- as.ggplot(a) + b
  
  save_plot(p, "fibroblast proportions", w = 16, h = 10, scale = 1, svg = T)



  }

immune_distance_matrix {
  
  
  Idents(Sobj_All) <- "Group"
  
  region_stats <- 
    table(Sobj_All$Group, Sobj_All$block_region) %>% 
    as.data.frame() %>% 
    filter(Var1 == "Immune") %>% 
    pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
    rename("block_region" = 1, "n.CellType" = Immune) %>% 
    filter(n.CellType > 10)
  
  regions <- region_stats %>% pull(block_region) %>% unique()  

  
  minima_list   <- vector("list", length = length(regions))
  for (i in 1:length(regions)) {
    
    print(paste("working on region ", i, " of ", length(regions)))
  
    ct_mat <- Sobj_All@meta.data %>% filter(Group == "Immune", block_region == regions[i]) %>% select(X, Y)
    tg_mat <- Sobj_All@meta.data %>% filter(block_region == regions[i]) %>% select(X, Y)
    
    ### calculate distance matrix and grad minima
    dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
    mins <- matrixStats::rowMins(dist)
    
    rownames(dist) <- rownames(ct_mat)
    colnames(dist) <- rownames(tg_mat)
    
    ct_ann <- select(rownames_to_column(Sobj_All@meta.data, "Cell"), Cell, CellType)
    tg_ann <- select(rownames_to_column(Sobj_All@meta.data, "target"), target, Group)
    
    minima_list[[i]] <- 
      dist %>% as.data.frame() %>% 
      rownames_to_column("Cell") %>% 
      pivot_longer(cols = rownames(tg_mat), names_to = "target") %>% 
      filter(value != 0) %>% 
      left_join(ct_ann, by = "Cell") %>% 
      left_join(tg_ann, by = "target") %>% 
      mutate(region = regions[i]) %>% 
      group_by(CellType, Group, region) %>%
      top_n(-10, value) %>% 
      #summarise(min = min(value))
      summarise(min = mean(value))
    
  }
  
  data <- bind_rows(minima_list) %>% left_join(ann, by = c("region" = "block_region")) %>% distinct()
  
  data %>% group_by(CellType, Group, Type, region) %>% 
    summarise(mean = mean(min)) %>% 
    ggplot(aes(y = mean, x = CellType, fill = Type)) + 
    geom_boxplot() + 
    facet_wrap(facets = "Group", nrow = 1) + 
    scale_y_log10() + 
    style



}

epithelial_perimeter_plot {
  
  minima_list  <-  readRDS("2021_ReRun/CODEX/output/Epithelial_Perimeter.rds")

bind_rows(minima_list) %>% left_join(Anno) %>% 
  filter(comparison == "CellClass", n.Ct.OI > 10, 
         #Target.ID == "Macrophage"
  ) %>% 
  filter(Target.ID %notin% c("Fibr_Adip", "Fibr_Immu", "Endo_Immu", "Immune_Endo", "ENAH")) %>% 
  
  ggplot(aes(x = Target.ID, y = peri.prop, fill = Type)) + 
  geom_boxplot() + 
  geom_point(pch = 21, alpha = 0.5, size = 2, 
             position = position_dodge(width = 0.75)) +
  style + 
  ylab(paste0("Absolute Neigborhood - 100px")) +
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
  #ylim(c(0, 10)) +
  scale_fill_manual(values = GID_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  #facet_wrap("Target.ID", nrow = 1) +
  facet_wrap("Ct.OI")
  

blocks <- Sobj_All@meta.data %>% select(block_region, block) %>% distinct()
  
bind_rows(minima_list) %>% left_join(Anno) %>% left_join(blocks) %>% 
  filter(comparison == "CellType", n.Ct.OI > 10, 
         str_detect(Target.ID, "Macro")
         #Target.ID %in% c("Capillary", "Endo_Main", "Endo_Immu", "Endo_SMA", "Lymph_EC")
         #Target.ID %in% c("Macrophage", "")
  ) %>% 
  #filter(Target.ID %notin% c("Fibr_Adip", "Fibr_Immu", "Endo_Immu", "Immune_Endo", "ENAH")) %>% 
  ggplot(aes(x = Type, y = peri.prop, fill = Type)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
  geom_point(aes(color = block), #pch = 21, 
             alpha = 1, size = 3, 
             position = position_dodge(width = 0.25)) +
  style + 
  ylab(paste0("Neigborhood Proportion")) +
  scale_y_log10() +
  scale_color_d3() +
  #ylim(c(0, 10)) +
  scale_fill_manual(values = GID_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  #facet_wrap("Target.ID", nrow = 1)
  facet_wrap("Target.ID")



}

capillary_fibro_distance_matrix {
  
  
  Idents(Sobj_All) <- "CellType"
  
  region_stats <- 
    table(Sobj_All$CellType, Sobj_All$block_region) %>% 
    as.data.frame() %>% 
    filter(Var1 == "Capillary") %>% 
    pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
    rename("block_region" = 1, "n.CellType" = Capillary) %>% 
    filter(n.CellType > 10)
  
  regions <- region_stats %>% pull(block_region) %>% unique()  
  
  
  minima_list   <- vector("list", length = length(regions))
  for (i in 1:length(regions)) {
    
    print(paste("working on region ", i, " of ", length(regions)))
    
    ct_mat <- Sobj_All@meta.data %>% filter(CellType == "Capillary", block_region == regions[i]) %>% select(X, Y)
    tg_mat <- Sobj_All@meta.data %>% filter(CellClass == "Fibroblast", block_region == regions[i]) %>% select(X, Y)
    
    ### calculate distance matrix and grad minima
    dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
    mins <- matrixStats::rowMins(dist)
    
    rownames(dist) <- rownames(ct_mat)
    colnames(dist) <- rownames(tg_mat)
    
    ct_ann <- select(rownames_to_column(Sobj_All@meta.data, "Cell"), Cell, CellType)
    tg_ann <- select(rownames_to_column(Sobj_All@meta.data, "target"), target, CellType)
    
    #minima_list[[i]] <- 
      dist %>% as.data.frame() %>% 
      rownames_to_column("Cell") %>% 
      pivot_longer(cols = rownames(tg_mat), names_to = "target") %>% 
      filter(value != 0) %>% 
      #left_join(ct_ann, by = "Cell") %>% 
      left_join(tg_ann, by = "target") %>% 
      mutate(region = regions[i]) %>% 
        left_join(ann, by = "region")
      #group_by(CellType, Group, region) %>%
      #top_n(-10, value) %>% 
      #summarise(min = min(value))
      #summarise(min = mean(value))
    
  }
  
  
  ann <- Sobj_All@meta.data %>% select(region = block_region, Type)
  data <- bind_rows(minima_list) %>% left_join(ann, by = c("region")) %>% distinct()
  
  data %>% group_by(CellType, Group, Type, region) %>% 
    summarise(mean = mean(min)) %>% 
    ggplot(aes(y = mean, x = CellType, fill = Type)) + 
    geom_boxplot() + 
    facet_wrap(facets = "Group", nrow = 1) + 
    scale_y_log10() + 
    style
  
  
  
}

endothelial_matrix {
  
  
  data <- data_rob.m %>% filter(CellClass %in% c("Lymph_EC", "Endothelial"), Gene %in% c("ACTA2", "LNX1", "CD45", "CD36", "LYVE1", "PDPN"))
  mat <- data %>% group_by(CellType, Gene) %>% summarise(mean = mean(Robust.Scaled)) %>% pivot_wider(names_from = Gene, values_from = mean) %>% column_to_rownames("CellType")
  p <- pheatmap(mat, 
           cluster_rows = T,
           cluster_cols = T,
           scale = "column", 
           fontsize = 25,
           main = "Marker Classification",
           color = colorRampPalette(brewer.pal(11, "BrBG"))(100))

 }

capillary_epi_distance_dotplots {
  
  Sobj_All@meta.data %>% as_tibble() %>% 
    filter(blk_region == "B9_reg13", Group %in% c("Vasculature", "Epithelial")) %>% 
    mutate(label = ifelse(Group == "Vasculature", as.character(CellType), as.character(Group))) %>% 
    arrange(Group) %>% ggplot(aes(x = X, y = Y, color = label)) + 
    geom_point() + style + 
    scale_color_d3()
  
}




### Euclidean distances

calc_euclideans {

  ### prepare the Sobjs and Frequency resource
  Idents(Sobj_All) <- "CellType"
  region_num <- table(Sobj_All$block_region, Sobj_All$CellClass) %>% as.data.frame() ### we will filter based on region frequencies later
  
  celltype <- "Fibroblast" ### this is the celltype we are interested in
  exclude <- c("Fibr_Adip", "Endo_Immu", "Immune_Endo", "Epi_Endo", "ENAH")
  targets  <- levels(Sobj_All)[!levels(Sobj_All) %in% c(celltype, exclude)]    ### these are the celltypes we are querying for distance
  
  compare_to_all_others {
    
    Sobj      <- subset(Sobj_All, idents = celltype)
    
    ### this loop iterates over all target cell types and finds the minimum distances in each region
    target_dist.mats <- vector("list", length = length(targets))
    target_mini_list <- vector("list", length = length(targets))
    for (j in 1:length(targets)) {
        
      query_celltypes <- targets[j]
      
      ### here we are making sure that we only look at region that have the query and target celltypes in them
      region_stats <- 
        region_num %>% filter(Var2 %in% c(celltype, c("Basal", "Luminal", "Luminal_HR+"))) %>% 
        #pivot_wider(names_from = "Var2", values_from = "Freq") %>% 
        #rename("block_region" = 1, "n.CellType" = 2, "n.Target" = 3) %>% 
        #filter(n.CellType > 0, n.Target > 0)
        filter(Freq > 10)
      
      regions <- region_stats %>% pull(Var1) %>% unique()
      
      
      ### cellnames for getting the X Y coordinates (could be stored in Sobj for effciency)
      immune_cells <- colnames(Sobj)
      epithe_cells <- WhichCells(Sobj_All, idents = query_celltypes)
      
      ### grabbing X Y coordinates and block-ID from raw data
      Immune <- Nuclei_raw %>% filter(Cell %in% immune_cells) %>% select(Cell, X, Y, block_region)
      Epith  <- Nuclei_raw %>% filter(Cell %in% epithe_cells) %>% select(Cell, X, Y, block_region)
      
      
      minima_list   <- vector("list", length = length(regions))
      dist.mat_list <- vector("list", length = length(regions))
      ### iterate over all included regions and make a table with minum distances
      for (i in 1:length(regions)) {
        
        ### turn the XY data into a matrix for the current region
        mat_I <- Immune %>% filter(block_region == regions[i]) %>% select(-block_region) %>% column_to_rownames("Cell")
        mat_E <- Epith  %>% filter(block_region == regions[i]) %>% select(-block_region) %>% column_to_rownames("Cell")
        
        ### calculate distance matrix and grad minima
        dist <- cdist(mat_I, mat_E, metric = "euclidean", p = 2)
        mins <- matrixStats::rowMins(dist)
        
        rownames(dist) <- rownames(mat_I)
        colnames(dist) <- rownames(mat_E)
        
        dist.mat_list[[i]]      <- dist
        names(dist.mat_list)[i] <- paste0(regions[i])
        
        ### store the minima alongside their cell names
        minima_list[[i]] <- 
          Immune %>% 
          filter(block_region == regions[i]) %>% 
          bind_cols(as.data.frame(mins)) %>% 
          left_join(region_stats, by = c("block_region" = "Var1"))
        
        names(minima_list)[i] <- paste0(regions[i])
        
      }
      
    target_dist.mats[[j]] <- dist.mat_list
    target_mini_list[[j]] <- bind_rows(minima_list) %>% mutate(CellType = celltype, "TargetCells" = targets[j])
    
    names(target_dist.mats)[j] <- paste0(targets[j])
    names(target_mini_list)[j] <- paste0(targets[j])
    
    }
    
    bind <- bind_rows(target_mini_list)
  
  }
  
  compare_to_rnd.sample {
  
    rnd1 <- subset(Sobj_All, idents = celltype, invert = T, downsample = 1000)
    Idents(rnd1) <- "sample.1k"
    rnd1$CellType <- Idents(rnd1)
    
    targets <- levels(rnd1)[!levels(rnd1) %in% levels(Sobj)]
    
    region_num1 <- table(rnd1$block_region, rnd1$CellType) %>% as.data.frame()
    region_num2 <- table(Sobj$block_region) %>% as.data.frame() %>% mutate(Var2 = levels(Sobj))
    region_num  <- region_num1 %>% bind_rows(region_num2)
    
    bind_list <- vector("list", length = length(levels(rnd1)))
    
    for (j in 1:length(targets)) {
      
      query_celltypes <- targets[j]
      
      region_stats <- region_num %>% filter(Var2 %in% c(query_celltypes, levels(Sobj))) %>% 
        pivot_wider(names_from = "Var2", values_from = "Freq") %>% 
        select("block_region" = 1, "n.CellType" = 3, "n.Target" = 2) %>% 
        filter(n.CellType > 0, n.Target > 0)
      
      regions <- region_stats %>% pull(block_region) %>% unique()
      
      
      immune_cells <- colnames(Sobj)
      epithe_cells <- WhichCells(rnd1, idents = query_celltypes)
      
      Immune <- Nuclei_raw %>% filter(Cell %in% immune_cells) %>% select(Cell, X, Y, block_region) #%>% left_join(imu_meta)
      Epith  <- Nuclei_raw %>% filter(Cell %in% epithe_cells) %>% select(Cell, X, Y, block_region)
      
      
      df_list <- vector("list", length = length(regions))
      
      ### iterate over all included regions and make a table with minum distances
      for (i in 1:length(regions)) {
        
        mat_I <- Immune %>% filter(block_region == regions[i]) %>% select(-block_region) %>% column_to_rownames("Cell")
        mat_E <- Epith  %>% filter(block_region == regions[i]) %>% select(-block_region) %>% column_to_rownames("Cell")
        
        dist <- cdist(mat_I, mat_E, metric = "euclidean", p = 2)
        mins <- matrixStats::rowMins(dist)
        
        df_list[[i]] <- 
          Immune %>% 
          filter(block_region == regions[i]) %>% 
          bind_cols(as.data.frame(mins)) %>% left_join(region_stats, by = "block_region")
        
      }
      
      bind_list[[j]] <- bind_rows(df_list) %>% mutate(CellType = levels(Sobj), "TargetCells" = targets[j])
      
    }
    
    bind_rnd <- bind_rows(bind_list)
    
  }
  
  
  saveRDS(target_dist.mats, paste0("2021_ReRun/CODEX/output/Distance_Matrices/Distance_Matrices_Per_CellType/CellType_", celltype, "_Matrices.rds"))
  saveRDS(target_mini_list, paste0("2021_ReRun/CODEX/output/Distance_Matrices/Minimal_Distances_Per_CellType//CellType_", celltype, "_MinDist.rds"))
  
  
  bind %>% bind_rows(bind_rnd) %>% 
    left_join(select(rownames_to_column(Sobj_All@meta.data, "Cell"), Cell, Type)) %>% 
    filter(mins != Inf, n.CellType > 10, n.Target > 10) %>% 
    group_by(CellType, TargetCells, block_region, Type) %>% 
    summarise(mean = mean(mins)) %>% 
    group_by( Type, CellType, TargetCells) %>% 
    mutate(rank = mean(mean)) %>% 
    ggplot(aes(y = reorder(TargetCells, rank), x = mean+1, fill = Type)) + 
    geom_boxplot() + 
    scale_x_sqrt() + 
    scale_fill_manual(values = GID_cols) +
    #geom_point(pch = 21, size = 2, alpha = 0.25, fill = "grey50") +
    ggtitle(paste0(celltype, " - avg. distance to all others")) +
    xlab("average distance") +
    style + 
    theme(legend.position = "left", axis.title.y = element_blank())
  
  
  
}


more_distance_crap {
  
  
  
  Idents(Sobj_All) <- "CellType"
  levels <- levels(Sobj_All)
  
  final_list <- vector("list", length = length(levels))
  names(final_list) <- levels
  for (k in 1:length(levels)) {
    
    celltype <- levels[k]
    target   <- levels[!levels %in% levels[k]]
    
    bind_list <- vector("list", length = length(target))
    for (j in 1:length(target)) {
      
      region_stats <- 
        table(Sobj_All$CellType, Sobj_All$block_region) %>% 
        as.data.frame() %>% 
        filter(Var1 %in% c(celltype, target[j])) %>% 
        pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
        rename("block_region" = 1, "n.CellType" = all_of(celltype), "n.Target" = all_of(target[j])) %>% 
        filter(n.CellType > 0, n.Target > 0)
      
      regions <- region_stats %>% pull(block_region) %>% unique()  
      
      
      minima_list   <- vector("list", length = length(regions))
      for (i in 1:length(regions)) {
        
        ct_mat <- Sobj_All@meta.data %>% filter(CellType == celltype,  block_region == regions[i]) %>% select(X, Y)
        tg_mat <- Sobj_All@meta.data %>% filter(CellType == target[j], block_region == regions[i]) %>% select(X, Y)
        
        ### calculate distance matrix and grad minima
        dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
        mins <- matrixStats::rowMins(dist)
        
        rownames(dist) <- rownames(ct_mat)
        colnames(dist) <- rownames(tg_mat)
        
        minima_list[[i]] <- 
          ct_mat %>% mutate(mins = mins, Target = target[j]) %>% 
          rownames_to_column("Cell") %>% 
          left_join(select(rownames_to_column(Sobj_All@meta.data, "Cell"), -X, -Y), by = "Cell") %>% 
          left_join(region_stats, by = "block_region")
        
      } # end loop i
      
      bind_list[[j]] <- bind_rows(minima_list) %>% mutate(Cell.OI = levels[k])
    } # end loop j
    
    final_list[[k]] <- bind_rows(bind_list)
  } # end loop k
  
  
  
  
  
  Idents(Sobj_All) <- "CellType"
  levels <- levels(Sobj_All)
  
  ### the whole shabang with top_N(-10)
  final_list <- vector("list", length = length(levels))
  names(final_list) <- levels
  for (k in 1:length(levels)) {
    
    celltype <- levels[k]
    target   <- levels[!levels %in% levels[k]]
    #target   <- levels
    
    print(paste0("Working on ", celltype, " #"  , k, " of ", length(levels)))
    
    bind_list <- vector("list", length = length(target))
    for (j in 1:length(target)) {
      
      print(paste0("Working on ", target[j], " #"  , k, " of ", length(target)))
      
      region_stats <- 
        table(Sobj_All$CellType, Sobj_All$block_region) %>% 
        as.data.frame() %>% 
        filter(Var1 %in% c(celltype, target[j])) %>% 
        pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
        rename("block_region" = 1, "n.CellType" = all_of(celltype), "n.Target" = all_of(target[j])) %>% 
        filter(n.CellType > 0, n.Target > 0)
      
      regions <- region_stats %>% pull(block_region) %>% unique()  
      
      
      minima_list   <- vector("list", length = length(regions))
      for (i in 1:length(regions)) {
        
        print(paste("working on region ", i, " of ", length(regions)))
        
        ct_mat <- Sobj_All@meta.data %>% filter(CellType == celltype,  block_region == regions[i]) %>% select(X, Y)
        tg_mat <- Sobj_All@meta.data %>% filter(CellType == target[j], block_region == regions[i]) %>% select(X, Y)
        
        ### calculate distance matrix and grad minima
        dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
        
        rownames(dist) <- rownames(ct_mat)
        colnames(dist) <- rownames(tg_mat)
        
        minima_list[[i]] <- 
          dist %>% as.data.frame() %>% 
          rownames_to_column("Cell") %>% 
          pivot_longer(cols = rownames(tg_mat)) %>% 
          group_by(Cell) %>% 
          top_n(-10, value) %>% 
          summarise(mins = mean(value)) %>% 
          mutate(Target = target[j]) %>% 
          left_join(select(rownames_to_column(Sobj_All@meta.data, "Cell"), -X, -Y), by = "Cell") %>% 
          left_join(region_stats, by = "block_region")
      } # end loop i
      
      bind_list[[j]] <- bind_rows(minima_list) %>% mutate(Cell.OI = levels[k])
    } # end loop j
    
    final_list[[k]] <- bind_rows(bind_list)
  } # end loop k
  final_list %>% saveRDS("2021_ReRun/CODEX/output/AverageDistances_List.rds")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  multiplicative_neighbor_ratio {
    
    ###################################
    ##### make the mutiplicative distance ratio set
    Idents(Sobj_All) <- "Group"
    levels <- levels(Sobj_All)
    
    nn = 10
    
    ### the whole shabang with top_N(-10)
    final_list <- vector("list", length = length(levels))
    names(final_list) <- levels
    for (k in 1:length(levels)) {
      
      celltype <- levels[k]
      print(paste0("Working on ", celltype, " #"  , k, " of ", length(levels)))
      
      region_stats <- 
        table(Sobj_All$Group, Sobj_All$block_region) %>% 
        as.data.frame() %>% 
        filter(Var1 %in% c(celltype)) %>% 
        pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
        rename("block_region" = 1, "n.Ct.OI" = all_of(celltype)) %>% 
        filter(n.Ct.OI > 0)
      
      regions <- region_stats %>% pull(block_region) %>% unique()  
      
      minima_list   <- vector("list", length = length(regions))
      for (i in 1:length(regions)) {
        
        print(paste("working on region ", i, " of ", length(regions)))
        
        ct_mat <- Sobj_All@meta.data %>% filter(Group     == celltype, block_region == regions[i]) %>% select(X, Y)
        tg_mat <- Sobj_All@meta.data %>% filter(block_region == regions[i]) %>% select(X, Y)
        
        ### calculate distance matrix and grad minima
        dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
        
        rownames(dist) <- rownames(ct_mat)
        colnames(dist) <- rownames(tg_mat)
        
        dist_long <- 
          dist %>% 
          as.data.frame() %>% 
          rownames_to_column("Cell") %>% 
          pivot_longer(cols = rownames(tg_mat), names_to = "Targets") %>% 
          left_join(select(rownames_to_column(Sobj_All@meta.data, "Targets"), 
                           "Targets", "CellType", "CellClass", "Group"), 
                    by = "Targets")
        
        
        nn.stats <- 
          dist_long %>% 
          group_by(Cell) %>% 
          filter(value != 0) %>% 
          
          top_n(-nn, value) %>%
          
          group_by(CellType) %>% 
          mutate(n.N.CellType = n()) %>% 
          group_by(CellClass) %>% 
          mutate(n.N.CellClass = n()) %>% 
          group_by(Group) %>% 
          mutate(n.N.Group = n()) %>% 
          select(CellType, CellClass, Group , n.N.CellType, n.N.CellClass, n.N.Group) %>% 
          distinct() %>% ungroup()
        
        
        df_celltype <- 
          dist_long %>% 
          select(Targets, CellType) %>% distinct() %>% 
          group_by(CellType) %>% 
          summarise(abundance = n()) %>% 
          left_join(select(nn.stats, CellType, n.N.CellType), 
                    by = "CellType") %>% 
          mutate(freq.NN = replace_na(n.N.CellType, 0)) %>% 
          select(-n.N.CellType) %>% 
          
          mutate(ratio = (freq.NN) / (length(rownames(ct_mat)) * abundance)) %>% 
          
          arrange(-ratio) %>% 
          mutate(Ct.OI = celltype, block_region = regions[i], comparison = "CellType") %>% 
          left_join(region_stats, by = "block_region") %>% 
          rename("Target.ID" = "CellType")
        
        df_cellclass <- 
          dist_long %>% 
          select(Targets, CellClass) %>% distinct() %>% 
          group_by(CellClass) %>% 
          summarise(abundance = n()) %>% 
          left_join(distinct(select(nn.stats, CellClass, n.N.CellClass)), 
                    by = "CellClass") %>% 
          mutate(freq.NN = replace_na(n.N.CellClass, 0)) %>% 
          select(-n.N.CellClass) %>% 
          
          mutate(ratio = (freq.NN) / (length(rownames(ct_mat)) * abundance)) %>% 
          
          arrange(-ratio) %>% 
          mutate(Ct.OI = celltype, block_region = regions[i], comparison = "CellClass") %>% 
          left_join(region_stats, by = "block_region")%>% 
          rename("Target.ID" = "CellClass")
        
        df_group <- 
          dist_long %>% 
          select(Targets, Group) %>% distinct() %>% 
          group_by(Group) %>% 
          summarise(abundance = n()) %>% 
          left_join(distinct(select(nn.stats, Group, n.N.Group)), 
                    by = "Group") %>% 
          mutate(freq.NN = replace_na(n.N.Group, 0)) %>% 
          select(-n.N.Group) %>% 
          
          mutate(ratio = (freq.NN) / (length(rownames(ct_mat)) * abundance)) %>% 
          
          arrange(-ratio) %>% 
          mutate(Ct.OI = celltype, block_region = regions[i], comparison = "Group") %>% 
          left_join(region_stats, by = "block_region") %>% 
          rename("Target.ID" = "Group")
        
        minima_list[[i]] <- bind_rows(df_cellclass, df_celltype, df_group)
        
      } # end loop i
      
      final_list[[k]] <- bind_rows(minima_list)
    } # end loop k
    #final_list %>% saveRDS("2021_ReRun/CODEX/output/Neighborhood_Ratios_multiplied_List.rds")
    
    
    
    ###################################
    ##### make the random sampled set
    Idents(Sobj_All) <- "CellType"
    levels <- levels(Sobj_All)
    nn = 10
    
    final_list <- vector("list", length = length(levels))
    names(final_list) <- levels
    for (k in 1:length(levels)) {
      
      celltype <- levels[k]
      
      region_stats <- 
        table(Sobj_All$CellType, Sobj_All$block_region) %>% 
        as.data.frame() %>% 
        filter(Var1 %in% c(celltype)) %>% 
        pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
        rename("block_region" = 1, "n.Ct.OI" = all_of(celltype)) %>% 
        filter(n.Ct.OI > 0)
      
      regions <- region_stats %>% pull(block_region) %>% unique()  
      #regions <- Sobj_All@meta.data %>% pull(block_region) %>% unique()
      
      minima_list   <- vector("list", length = length(regions))
      for (i in 1:length(regions)) {
        
        ct_mat <- Sobj_All@meta.data %>% filter(CellType     == celltype, block_region == regions[i]) %>% select(X, Y)
        tg_mat <- Sobj_All@meta.data %>% filter(block_region == regions[i]) %>% select(X, Y)
        
        
        iter = 10
        
        sam_list   <- vector("list", length = iter)
        for (p in 1:iter) {
          
          sam <- rownames(tg_mat) %>% sample(0.2*length(rownames(tg_mat)))
          
          dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
          rownames(dist) <- rownames(ct_mat)
          colnames(dist) <- rownames(tg_mat)
          
          dist_long <- 
            dist %>% 
            as.data.frame() %>% 
            rownames_to_column("Cell") %>% 
            pivot_longer(cols = rownames(tg_mat), names_to = "Targets") %>% 
            left_join(select(rownames_to_column(Sobj_All@meta.data, "Targets"), 
                             "Targets", "CellType"), 
                      by = "Targets") %>% 
            mutate(CellType = ifelse(Targets %in% sam, "rnd", "other"))
          
          
          
          nn.stats <- 
            dist_long %>% 
            group_by(Cell) %>% 
            filter(value != 0) %>% 
            
            top_n(-nn, value) %>%
            
            group_by(CellType) %>% 
            summarise(n.N.CellType = n())
          
          
          sam_list[[p]] <- 
            dist_long %>% 
            select(Targets, CellType) %>% distinct() %>% 
            group_by(CellType) %>% 
            summarise(abundance = n()) %>% 
            left_join(nn.stats, by = "CellType") %>% 
            mutate(freq.NN = replace_na(n.N.CellType, 0)) %>% 
            select(-n.N.CellType) %>% 
            
            mutate(ratio = (freq.NN) / (length(rownames(ct_mat)) * abundance)) %>% 
            
            mutate(Ct.OI = celltype, block_region = regions[i], comparison = paste0("iteration.", p)) %>% 
            left_join(region_stats, by = "block_region") %>% 
            rename("Target.ID" = "CellType")
          
        }
        
        minima_list[[i]] <- bind_rows(sam_list)
        
      } # end loop i
      
      final_list[[k]] <- bind_rows(minima_list)
    } # end loop k
    #final_list %>% saveRDS("2021_ReRun/CODEX/output/Neighborhood_Ratios_multiplied_Random_List.rds")
    
    
    
    
    
    
    ratio_data <- 
      bind_rows(final_list) %>% 
      left_join(Anno) 
    
    inclusion_of_random_Set {
      
      #rnd.set <- minima_list_rnd %>% bind_rows() %>% filter(Target.ID == "rnd") %>% group_by(block_region) %>% mutate(ratio = mean(ratio), freq.NN = mean(freq.NN), comparison = "random.set") %>% distinct()
      
      ratio_data <- 
        bind_rows(final_list) %>% 
        left_join(Anno) %>% 
        ratio_data %>% bind_rows(rnd_data)
      
      rnd_data <- filter(bind_rows(final_list_rnd), Target.ID ==  "rnd") %>% 
        group_by(block_region, Target.ID, abundance, Ct.OI, n.Ct.OI) %>% 
        summarise(ratio = mean(ratio), freq.NN = mean(freq.NN)) %>% left_join(Anno) %>% 
        mutate(comparison = "random.sample")
      
    }
    
    
    exclude <- c("Fibr_Adip", "Endo_Immu", "Immune_Endo", "Epi_Endo", "ENAH", "Fibr_Immu", "Epi_Other")
    include <- c("Luminal_HR+", "Luminal_K8hi", "Luminal_K8lo", "Basal_Lo", "Basal_Hi", "rnd")
    
    ratio_data %>% 
      filter(Ct.OI == "Adipose", 
             comparison %in% c("CellClass"),
             n.Ct.OI > 10, 
             abundance > 10,
             #Target.ID %notin% "Endothelial",
      ) %>% 
      ggplot(aes(x = Target.ID, y = (ratio + 1), fill = Type)) + 
      geom_boxplot() + 
      geom_point(pch = 21, alpha = 0.5, size = 2, 
                 position = position_dodge(width = 0.75)) +
      #ylim(c(0, 0.05)) +
      scale_y_log10(limits = c(1, 1.05)) +
      ggtitle("Neighborhood Ratio (#N Y.of.X) / (#X * #Y)"  ) +
      style + 
      ylab(paste0("Neighborhood Ratio - nn=", nn)) +
      scale_fill_manual(values = GID_cols) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap("Ct.OI")
    
  }
  
  
  
  
  neighborhood_proportions {
    
    ##############################################
    ##### Neighborhood Proportions and Aggregation
    Idents(Sobj_All) <- "Group"
    levels <- levels(Sobj_All)
    nn = 25 #n.Neighbors
    
    # storage for final results 
    # "2021_ReRun/CODEX/output/Neighborhood_Proportions_and_Aggregation_List.rds"
    final_list <- vector("list", length = length(levels)) 
    names(final_list) <- levels
    for (k in 1:length(levels)) {
      
      celltype <- levels[k]
      print(paste0("Working on ", celltype, " #"  , k, " of ", length(levels)))
      
      # proportion of Celltype of interest in all regions and filter for at least 1 cells to avoid errors
      region_stats <- 
        table(Sobj_All$Group, Sobj_All$block_region) %>% 
        as.data.frame() %>% 
        filter(Var1 %in% c(celltype)) %>% 
        pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
        rename("block_region" = 1, "n.Ct.OI" = all_of(celltype)) %>% 
        filter(n.Ct.OI > 0)
      
      # list of regions to iterate on
      regions <- region_stats %>% pull(block_region) %>% unique()  
      minima_list   <- vector("list", length = length(regions))
      for (i in 1:length(regions)) {
        
        print(paste("working on region ", i, " of ", length(regions)))
        
        ct_mat <- Sobj_All@meta.data %>% filter(Group     == celltype, block_region == regions[i]) %>% select(X, Y)
        tg_mat <- Sobj_All@meta.data %>% filter(block_region == regions[i]) %>% select(X, Y)
        
        
        dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
        rownames(dist) <- rownames(ct_mat)
        colnames(dist) <- rownames(tg_mat)
        
        
        dist_long <- 
          dist %>% 
          as.data.frame() %>% 
          rownames_to_column("Cell") %>% 
          pivot_longer(cols = rownames(tg_mat), names_to = "Targets") %>% 
          left_join(select(rownames_to_column(Sobj_All@meta.data, "Targets"), 
                           "Targets", "CellType", "CellClass", "Group"), 
                    by = "Targets")
        
        
        cc_stats <- 
          dist_long %>%
          group_by(Cell, CellClass) %>%
          mutate(n.Target = n()) %>% 
          group_by(Cell) %>% 
          filter(value != 0) %>% 
          top_n(-nn, value) %>%
          group_by(CellClass, n.Target) %>% 
          summarise(in.hood = n(), .groups =  "drop") %>% 
          mutate(hood.size = (nn*length(rownames(ct_mat)))) %>% 
          #mutate(adjacency = (in.hood/n.Target)*100, 
          #mutate(aggregation = in.hood / (n.Target * length(rownames(ct_mat))),        
          mutate(hood.perc   = (in.hood/hood.size)*100,
                 aggregation = hood.perc / ((n.Target / length(rownames(tg_mat))*100))
          ) %>% 
          mutate(Ct.OI = celltype,
                 block_region = regions[i],
                 comparison = "CellClass") %>% 
          rename(Target.ID  = "CellClass") %>% 
          left_join(region_stats, by = "block_region")
        
        ct_stats <- 
          dist_long %>% 
          group_by(Cell, CellType) %>%
          mutate(n.Target = n()) %>% 
          group_by(Cell) %>% 
          filter(value != 0) %>% 
          top_n(-nn, value) %>%
          group_by(CellType, n.Target) %>% 
          summarise(in.hood = n(), .groups =  "drop") %>% 
          mutate(hood.size = (nn*length(rownames(ct_mat)))) %>% 
          #mutate(adjacency = (in.hood/n.Target)*100, 
          #mutate(aggregation = in.hood / (n.Target * length(rownames(ct_mat))),        
          mutate(hood.perc   = (in.hood/hood.size)*100,
                 aggregation = hood.perc / ((n.Target / length(rownames(tg_mat))*100))
          ) %>%
          mutate(Ct.OI = celltype,
                 block_region = regions[i],
                 comparison = "CellType") %>% 
          rename(Target.ID  = "CellType") %>% 
          left_join(region_stats, by = "block_region")
        
        gr_stats <- 
          dist_long %>% 
          group_by(Cell, Group) %>%
          mutate(n.Target = n()) %>% 
          group_by(Cell) %>% 
          filter(value != 0) %>% 
          top_n(-nn, value) %>%
          group_by(Group, n.Target) %>% 
          summarise(in.hood = n(), .groups =  "drop") %>% 
          mutate(hood.size = (nn*length(rownames(ct_mat)))) %>% 
          #mutate(adjacency = (in.hood/n.Target)*100, 
          #mutate(aggregation = in.hood / (n.Target * length(rownames(ct_mat))),        
          mutate(hood.perc   = (in.hood/hood.size)*100,
                 aggregation = hood.perc / ((n.Target / length(rownames(tg_mat))*100))
          ) %>%
          mutate(Ct.OI = celltype,
                 block_region = regions[i],
                 comparison = "Group") %>% 
          rename(Target.ID  = "Group") %>% 
          left_join(region_stats, by = "block_region")
        
        
        minima_list[[i]] <- bind_rows(ct_stats, cc_stats, gr_stats)
      } # end loop i
      
      final_list[[k]] <- bind_rows(minima_list)
    } # end loop k
    
    # random sampling as control
    # "2021_ReRun/CODEX/output/Neighborhood_Proportions_and_Aggregation_List.rds"
    final_list <- vector("list", length = length(levels)) 
    names(final_list) <- levels
    for (k in 1:length(levels)) {
      
      celltype <- levels[k]
      print(paste0("Working on ", celltype, " #"  , k, " of ", length(levels)))
      
      # proportion of Celltype of interest in all regions and filter for at least 1 cells to avoid errors
      region_stats <- 
        table(Sobj_All$CellType, Sobj_All$block_region) %>% 
        as.data.frame() %>% 
        filter(Var1 %in% c(celltype)) %>% 
        pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
        rename("block_region" = 1, "n.Ct.OI" = all_of(celltype)) %>% 
        filter(n.Ct.OI > 0)
      
      # list of regions to iterate on
      regions <- region_stats %>% pull(block_region) %>% unique()  
      minima_list   <- vector("list", length = length(regions))
      for (i in 1:length(regions)) {
        
        print(paste("working on region ", i, " of ", length(regions)))
        
        ct_mat <- Sobj_All@meta.data %>% filter(CellType     == celltype, block_region == regions[i]) %>% select(X, Y)
        tg_mat <- Sobj_All@meta.data %>% filter(block_region == regions[i]) %>% select(X, Y)
        
        
        dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
        rownames(dist) <- rownames(ct_mat)
        colnames(dist) <- rownames(tg_mat)
        
        dist_long <- 
          dist %>% 
          as.data.frame() %>% 
          rownames_to_column("Cell") %>% 
          pivot_longer(cols = rownames(tg_mat), names_to = "Targets")
        
        iter = 10
        
        sam_list   <- vector("list", length = iter)
        for (p in 1:iter) {
          
          print(paste("sampling ", p, " of ", iter))
          
          sam <- rownames(tg_mat) %>% sample(0.2*length(rownames(tg_mat)))
          
          dist_rnd <- dist_long %>% mutate(CellType = ifelse(Targets %in% sam, "rnd", "other"))
          
          sam_list[[p]] <- 
            dist_rnd %>%
            group_by(Cell, CellType) %>%
            mutate(n.Target = n()) %>% 
            group_by(Cell) %>% 
            filter(value != 0) %>% 
            top_n(-nn, value) %>%
            group_by(CellType, n.Target) %>% 
            summarise(in.hood = n(), .groups =  "drop") %>% 
            mutate(hood.size = (nn*length(rownames(ct_mat)))) %>% 
            #mutate(adjacency = (in.hood/n.Target)*100, 
            #mutate(aggregation = in.hood / (n.Target * length(rownames(ct_mat))),        
            mutate(hood.perc   = (in.hood/hood.size)*100,
                   aggregation = hood.perc / ((n.Target / length(rownames(tg_mat))*100))
            ) %>% 
            mutate(Ct.OI = celltype,
                   block_region = regions[i], comparison = paste0("iteration.", p)) %>% 
            rename(Target.ID  = "CellType") %>% 
            left_join(region_stats, by = "block_region") %>% 
            filter(Target.ID == "rnd")  
        }
        
        minima_list[[i]] <- bind_rows(sam_list)
        
      } # end loop i
      
      final_list[[k]] <- bind_rows(minima_list)
    } # end loop k
    #final_list %>% saveRDS("2021_ReRun/CODEX/output/Neighborhood_Proportions_and_Aggregation_Random_List.rds")
    
    
    
    
    
    
    data <- final_list %>% #[[k]] %>% 
      bind_rows() %>% 
      bind_rows(minima_list) %>% 
      filter(str_detect(comparison, "CellClass|iteration"),
             Ct.OI == "Luminal_HR+",
             n.Ct.OI  > 10,
             n.Target > 10) %>% 
      left_join(Anno)
    
    
    data <- readRDS("2021_ReRun/CODEX/output/Neighborhood_Proportions_and_Aggregation_Group.Data.rds")
    
    data_prop <- 
      data %>% 
      filter(comparison == "CellType",
             Ct.OI == "Epithelial",
             n.Ct.OI  > 10,
             n.Target > 10) %>% 
      left_join(Anno)
    
    
    data_prop %>% #filter(Target.ID == "Macrophage") %>% 
      ggplot(aes(x = Target.ID, y = aggregation, fill = Type)) + 
      geom_boxplot() + 
      geom_point(pch = 21, alpha = 0.5, size = 2, 
                 position = position_dodge(width = 0.75)) +
      style + 
      ylab(paste0("Aggregation - nn.", nn)) +
      ylim(c(0, 10)) +
      scale_fill_manual(values = GID_cols) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap("Ct.OI")
    
    
    data_prop %>% filter(Target.ID == "Stroma_Other") %>% 
      ggplot(aes(x = Target.ID, y = hood.perc, fill = Type)) + 
      geom_boxplot() + 
      geom_point(pch = 21, alpha = 0.5, size = 2, 
                 position = position_dodge(width = 0.75)) +
      style + 
      ylab(paste0("Neighborhood Content - nn.", nn)) +
      ylim(c(0,6)) +
      scale_fill_manual(values = GID_cols) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap("Ct.OI")
    
    
  }
  
  
  
  Sobj_All@meta.data[inperimeter, ] %>% ggplot(aes(x = X, y = Y, color = Group)) + geom_point() + ylim(c(1800, 3500)) + xlim(c(1200, 4000))
  Sobj_All@meta.data[rownames(tg_mat), ] %>% ggplot(aes(x = X, y = Y, color = Group)) + geom_point() + ylim(c(1800, 3500)) + xlim(c(1200, 4000))
  
  
  
  
  neighborhood_absolute_distance {
    
    ##############################################
    ##### Neighborhood Absolute Distance
    ap_group <- readRDS("2021_ReRun/CODEX/output/AffinityPropagation_Group_cellClusterIDs.RDS")
    
    Idents(Sobj_All) <- "Group"
    levels <- levels(Sobj_All)
    distance = 100 #n.Neighbors
    
    # storage for final results 
    # "2021_ReRun/CODEX/output/Neighborhood_Proportions_and_Aggregation_List.rds"
    final_list <- vector("list", length = length(levels)) 
    names(final_list) <- levels
    for (k in 1:length(levels)) {
      
      celltype <- levels[k]
      print(paste0("Working on ", celltype, " #"  , k, " of ", length(levels)))
      
      # proportion of Celltype of interest in all regions and filter for at least 1 cells to avoid errors
      region_stats <- 
        table(Sobj_All$Group, Sobj_All$block_region) %>% 
        as.data.frame() %>% 
        filter(Var1 %in% c(celltype)) %>% 
        pivot_wider(names_from = "Var1", values_from = "Freq") %>% 
        rename("block_region" = 1, "n.Ct.OI" = all_of(celltype)) %>% 
        filter(n.Ct.OI > 0)
      
      ap_cells <- ap_group %>% filter(Group == celltype, AP.CT.Count > 35) %>% rownames()
      
      # list of regions to iterate on
      regions <- region_stats %>% pull(block_region) %>% unique()  
      minima_list   <- vector("list", length = length(regions))
      for (i in 1:length(regions)) {
        
        print(paste("working on region ", i, " of ", length(regions)))
        
        ct_mat <- Sobj_All@meta.data %>% filter(Group == celltype, block_region == regions[i]) %>% select(X, Y)
        tg_mat <- Sobj_All@meta.data %>% filter(block_region == regions[i]) %>% select(X, Y)
        
        
        dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
        rownames(dist) <- rownames(ct_mat)
        colnames(dist) <- rownames(tg_mat)
        
        
        dist_long <- 
          dist %>% 
          as.data.frame() %>% 
          rownames_to_column("Cell") %>% 
          pivot_longer(cols = rownames(tg_mat), names_to = "Targets") %>% 
          left_join(select(rownames_to_column(Sobj_All@meta.data, "Targets"), 
                           "Targets", "CellType", "CellClass", "Group"), 
                    by = "Targets")
        
        
        inperimeter <- 
          dist_long %>% 
          filter(Cell %in% ap_cells, 
                 Targets %notin% rownames(ct_mat)) %>% 
          #group_by(Cell, CellClass) %>%
          #mutate(n.Target = n()) %>% 
          group_by(Cell) %>% 
          #filter(value != 0) %>% 
          filter(value < distance) %>% 
          pull(Targets) %>% 
          unique()
        
        
        cc_stats <- 
          dist_long %>% select(-CellType, -Cell, -value) %>% 
          distinct() %>% 
          filter(Targets %in% inperimeter) %>% 
          group_by(Group, CellClass) %>% 
          summarise(inperi = n(), .groups = "drop") %>% 
          mutate(population = sum(inperi), 
                 n.Ct.OI = length(intersect(rownames(ct_mat), ap_cells))) %>% 
          mutate(peri.prop = (inperi/population)*100, 
                 ct.peri.ratio = n.Ct.OI/population) %>% 
          mutate(Ct.OI = celltype,
                 block_region = regions[i],
                 comparison = "CellClass") %>% 
          rename(Target.ID  = "CellClass")
        
        
        ct_stats <- 
          dist_long %>% select(-CellClass, -Cell, -value) %>% 
          distinct() %>% 
          filter(Targets %in% inperimeter) %>% 
          group_by(Group, CellType) %>% 
          summarise(inperi = n(), .groups = "drop") %>% 
          mutate(population = sum(inperi), 
                 n.Ct.OI = length(intersect(rownames(ct_mat), ap_cells))) %>% 
          mutate(peri.prop = (inperi/population)*100, 
                 ct.peri.ratio = n.Ct.OI/population) %>% 
          mutate(Ct.OI = celltype,
                 block_region = regions[i],
                 comparison = "CellType") %>% 
          rename(Target.ID  = "CellType")
        
        
        minima_list[[i]] <- bind_rows(ct_stats, cc_stats)
      } # end loop i
      
      final_list[[k]] <- bind_rows(minima_list)
    } # end loop k
    
    
    
    
    
    
    
    bind_rows(minima_list) %>% left_join(Anno) %>% 
      filter(comparison == "CellClass", n.Ct.OI > 10, 
             #Target.ID == "Macrophage"
      ) %>% 
      filter(Target.ID %notin% c("Fibr_Adip", "Fibr_Immu", "Endo_Immu", "Immune_Endo", "ENAH")) %>% 
      ggplot(aes(x = Target.ID, y = peri.prop, fill = Type)) + 
      geom_boxplot() + 
      geom_point(pch = 21, alpha = 0.5, size = 2, 
                 position = position_dodge(width = 0.75)) +
      style + 
      ylab(paste0("Absolute Neigborhood - 100px")) +
      #ylim(c(0, 10)) +
      scale_fill_manual(values = GID_cols) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap("Ct.OI")
    
    
    
    # plot the perimeter ratios
    bind_rows(minima_list) %>% 
      mutate(peri.ct.ratio = population/n.Ct.OI) %>% 
      left_join(Anno) %>% 
      filter(n.Ct.OI > 10) %>% 
      select(Ct.OI, Type, block_region, ct.peri.ratio) %>% 
      distinct() %>%
      ggplot(aes(x = Type, y = ct.peri.ratio, fill = Type)) + 
      geom_boxplot() + 
      geom_point(pch = 21, alpha = 0.5, size = 2, 
                 position = position_dodge(width = 0.75)) +
      style + 
      ylab(paste0("Ratio of Epith. to Non-Epith. cells in perimeter")) +
      stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
      #ylim(c(0, 10)) +
      scale_fill_manual(values = GID_cols) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap("Ct.OI")
    
    p <- bind_rows(minima_list) %>% 
      mutate(peri.ct.ratio = population/n.Ct.OI) %>% 
      left_join(Anno) %>% 
      left_join(banno) %>% 
      filter(n.Ct.OI > 10) %>% 
      select(Ct.OI, Type, block_region, peri.ct.ratio, block) %>% 
      distinct() %>%
      ggplot(aes(x = Type, y = peri.ct.ratio, fill = Type)) + 
      geom_boxplot() + 
      geom_point(aes(color = block), 
                 #pch = 21, 
                 #alpha = 0.5, 
                 size = 2, 
                 position = position_dodge(width = 0.25)) +
      style + 
      ylab(paste0("Ratio of Non-Epith. to Epith. cells in perimeter")) +
      stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
      #ylim(c(0, 10)) +
      scale_color_d3() +
      scale_y_log10() +
      scale_fill_manual(values = GID_cols) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap("Ct.OI")
    
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
  
  data <- ap_group %>% filter(Group == celltype, AP.CT.Count >30 & AP.CT.Count < 200) %>% group_by(blk_region, AP.CT) %>% summarise(cent.X = mean(X), cent.Y = mean(Y)) %>% mutate(nAP = n())
  
  
  regions <- data %>% filter(nAP > 10) %>% pull(blk_region) %>% unique()
  final_list <- vector("list", length = length(regions)) 
  for (i in 1:length(regions)) {
    
    ct_mat <- data %>% ungroup() %>% filter(blk_region == regions[i]) %>% select(-nAP, -blk_region) %>% column_to_rownames("AP.CT")
    tg_mat <- data %>% ungroup() %>% filter(blk_region == regions[i]) %>% select(-nAP, -blk_region) %>% column_to_rownames("AP.CT")
    
    dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
    rownames(dist) <- rownames(ct_mat)
    colnames(dist) <- rownames(tg_mat)
    
    final_list[[i]] <- 
      dist %>% as.data.frame() %>% 
      rownames_to_column("AP.CT") %>% 
      pivot_longer(cols = rownames(tg_mat), names_to = "Targets") %>% 
      filter(value != 0) %>% 
      group_by(AP.CT) %>% 
      top_n(-1, value) %>% 
      mutate(blk_region = regions[i])
    
  }
  
  final_list %>% bind_rows() %>% 
    left_join(ann) %>%
    group_by(Type, blk_region) %>% 
    summarise(mean = mean(value)) %>% 
    ggplot(aes(x = Type, y = mean, fill = Type)) + 
    geom_boxplot() + 
    style
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### Plotting Average Distances
  ###############################################################
  cf_tm_block_Split {
    
    query   <- c("Capillary")
    include <- c("Basal_Hi", "Basal_Lo", "Luminal_HR+", "Luminal_K8hi")
    
    sampling <- prop.table(table(Sobj_All$Group, Sobj_All$block_region), margin = 2) %>% as.data.frame() %>% filter(Var1 == "Epithelial") %>% filter(Freq > 0.3) %>% pull(Var2)
    
    
    cf_blocks <- Sobj_All@meta.data %>% filter(Type == "CF") %>% pull(block_region) %>% unique()
    tm_blocks <- Sobj_All@meta.data %>% filter(Type == "TM") %>% pull(block_region) %>% unique()
    
    
    data_cf <- 
      bind_rows(final_list_avg[query]) %>%
      filter(block_region %in% cf_blocks) %>% 
      mutate(Type = ifelse(block_region %in% sample(cf_blocks, length(cf_blocks)/2), "mock.A", "mock.B"))
    #mutate(Type = ifelse(block_region %in% sample(tm_blocks, 31), "mock.A", "mock.B")) %>% 
    
    data_tm <- 
      bind_rows(final_list_avg[query]) %>%
      filter(block_region %in% tm_blocks) %>% 
      mutate(Type = ifelse(block_region %in% sample(tm_blocks, length(tm_blocks)/2), "mock.C", "mock.D"))
    
    data <- bind_rows(data_tm, data_cf) %>% 
      #group_by(block_region, Target) %>% 
      filter(mins != Inf, n.CellType > 10, n.Target > 10) %>% 
      group_by(CellType, Target, block_region, Type) %>% 
      summarise(mean = mean(mins)) %>% 
      group_by(Type, CellType, Target) %>% 
      mutate(rank = mean(mean))
    
    
    data %>% 
      filter(Target %in% include) %>% 
      ggplot(aes(x = Target, 
                 y = mean+1, 
                 fill = Type)) + 
      geom_boxplot() + 
      scale_y_sqrt() + 
      #scale_fill_manual(values = GID_cols) +
      geom_point(pch = 21, size = 2, alpha = 0.5, #fill = "grey50", 
                 position = position_dodge(width = 0.75)) +
      ggtitle(paste0(query, " - avg. distances")) +
      xlab("average distance") +
      style + 
      #scale_fill_manual(values = GID_cols) +
      #coord_flip() +
      theme(legend.position = "top", 
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)) #+
    #facet_wrap("Cell.OI")
    
  }
  
  
  query   <- c("Capillary", "Endo_Main")
  query <- c("Fibro_ARhi", "Fibro_LAMBhi", "Fibr_Epi")
  exclude <- c("Fibr_Adip", "Endo_Immu", "Immune_Endo", "Epi_Endo", "ENAH", "Fibr_Immu")
  exclude <- c("Epi_Other", "Stroma_Other", "Macrophage", "ENAH")
  include <- c("Luminal_HR+", "Basal_Lo", "Basal_Hi", "Fibro_ARhi", "Lymph_EC", "Capillary", "Macrophage", "Immune_Main", "Luminal_K8hi", "Luminal_K8lo", "Stroma_Other", "Endo_Main")
  include <- c("Luminal_HR+")
  
  
  data <- 
    bind_rows(final_list_avg[query]) %>%
    filter(mins != Inf, 
           n.CellType > 10, 
           n.Target > 10) %>% 
    group_by(Cell.OI, CellType, Target, block_region, Type) %>% 
    summarise(mean = mean(mins)) %>% 
    group_by(Type, CellType, Target) %>% 
    mutate(rank = mean(mean))
  
  data %>% 
    filter(Target %in% include) %>% 
    ggplot(aes(x = Target, 
               y = mean+1, 
               fill = Type)) + 
    geom_boxplot() + 
    scale_y_sqrt() + 
    #scale_fill_manual(values = GID_cols) +
    geom_point(pch = 21, size = 2, alpha = 0.5, #fill = "grey50", 
               position = position_dodge(width = 0.75)) +
    ggtitle(paste0(query, " - avg. distances")) +
    xlab("average distance") +
    style + 
    scale_fill_manual(values = GID_cols) +
    #coord_flip() +
    theme(legend.position = "top", 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap("Cell.OI")
  
  
  
}


calc_fibro_epi_distance {
  
  ### prepare the Sobjs and Frequency resource
  Idents(Sobj_All) <- "CellClass"
  ann <- Sobj_All@meta.data %>% rownames_to_column("Cell") %>% select(Cell, CellClass, CellType, Type)
  
  celltype <- "Fibroblast"
    
  celltype_regions <- 
    table(Sobj_All$CellClass, Sobj_All$blk_region) %>% 
    as.data.frame() %>% 
    filter(Var1 == celltype, Freq > 10) %>% 
    pull(Var2) %>% unique()
  
  regions <- 
    table(Sobj_All$Group, Sobj_All$blk_region) %>% 
    as.data.frame() %>% 
    filter(Var2 %in% celltype_regions, 
           Var1 == "Epithelial", Freq > 10) %>% 
    pull(Var2) %>% unique()
      
      
  minima_list   <- vector("list", length = length(regions))
  ### iterate over all included regions and make a table with minum distances
  for (i in 1:length(regions)) {
    
    ### turn the XY data into a matrix for the current region
    mat_I <- Sobj_All@meta.data %>% filter(CellClass == celltype, blk_region == regions[i]) %>% select(X, Y) #%>% column_to_rownames("Cell")
    mat_E <- Sobj_All@meta.data %>% filter(Group == "Epithelial", blk_region == regions[i]) %>% select(X, Y)
    
    ### calculate distance matrix and grad minima
    dist <- cdist(mat_I, mat_E, metric = "euclidean", p = 2)
    mins <- matrixStats::rowMins(dist)
    
    rownames(dist) <- rownames(mat_I)
    colnames(dist) <- rownames(mat_E)
    
    dist.mat_list[[i]]      <- dist
    names(dist.mat_list)[i] <- paste0(regions[i])
    
    
    minima_list[[i]] <- 
      dist %>% as.data.frame() %>% 
      rownames_to_column("Cell") %>% 
      pivot_longer(cols = rownames(mat_E), names_to = "Target") %>% 
      group_by(Cell) %>% 
      top_n(-1, value) %>% 
      left_join(ann) %>% mutate(blk_region = regions[i])
    
  
  }
      

  data <- minima_list %>% bind_rows()
  
  
  p <- data %>% filter(CellType != "Fibr_Adip") %>% 
    group_by(CellType, blk_region) %>% mutate(n = n()) %>% arrange(n) %>% 
    filter(n > 5) %>% group_by(blk_region, CellType) %>% 
    summarise(mean = mean(value*0.325)) %>% 
    ggplot(aes(x = mean, y = CellType, fill = CellType)) + 
    geom_boxplot() + 
    geom_point(pch = 21, size = 3, alpha = 0.5,
               position = position_dodge(width = 0.75)) +
    scale_x_log10() +
    style + 
    scale_fill_d3()
  
  
  
}

calc_Immu_epi_distance {
  
  Sobj_Immu <- readRDS("2021_ReRun/CODEX/CellTyping/S5.3_Immune_filtered_Annotated.rds")
  
  ### prepare the Sobjs and Frequency resource
  Idents(Sobj_All) <- "CellClass"
  ann <- Sobj_All@meta.data %>% rownames_to_column("Cell") %>% select(Cell, CellClass, CellType, Type)
  
  celltype1 <- "Macrophage"
  celltype2 <- "Immune"
  
  celltype_regions1 <- 
    table(Sobj_All$CellClass, Sobj_All$blk_region) %>% 
    as.data.frame() %>% 
    filter(Var1 == celltype1, Freq > 5) %>% 
    pull(Var2) %>% unique()
  
  celltype_regions2 <- 
    table(Sobj_All$CellClass, Sobj_All$blk_region) %>% 
    as.data.frame() %>% 
    filter(Var1 == celltype2, Freq > 10) %>% 
    pull(Var2) %>% unique()
  
  celltype_regions <- c(celltype_regions1, celltype_regions2) %>% unique()
  
  regions <- 
    table(Sobj_All$Group, Sobj_All$blk_region) %>% 
    as.data.frame() %>% 
    filter(Var2 %in% celltype_regions, 
           Var1 == "Epithelial", Freq > 10) %>% 
    pull(Var2) %>% unique()
  
  
  minima_list   <- vector("list", length = length(regions))
  ### iterate over all included regions and make a table with minum distances
  for (i in 1:length(regions)) {
    
    ### turn the XY data into a matrix for the current region
    #mat_I <- Sobj_All@meta.data %>% filter(CellClass %in% c(celltype1, celltype2), blk_region == regions[i]) %>% select(X, Y) #%>% column_to_rownames("Cell")
    mat_I <- Sobj_All@meta.data[colnames(Sobj_Immu), ] %>% filter(blk_region == regions[i]) %>% select(X, Y) #%>% column_to_rownames("Cell")
    mat_E <- Sobj_All@meta.data %>% filter(Group == "Epithelial", blk_region == regions[i]) %>% select(X, Y)
    
    ### calculate distance matrix and grad minima
    dist <- cdist(mat_I, mat_E, metric = "euclidean", p = 2)
    mins <- matrixStats::rowMins(dist)
    
    rownames(dist) <- rownames(mat_I)
    colnames(dist) <- rownames(mat_E)
    
    minima_list[[i]] <- 
      dist %>% as.data.frame() %>% 
      rownames_to_column("Cell") %>% 
      pivot_longer(cols = rownames(mat_E), names_to = "Target") %>% 
      group_by(Cell) %>% 
      top_n(-1, value) %>% 
      left_join(ann) %>% mutate(blk_region = regions[i])
    
    
  }
  
  data <- minima_list %>% bind_rows()
  
  p1 <- data %>% 
    group_by(CellType, blk_region) %>% mutate(n = n()) %>% arrange(n) %>% 
    filter(n > 5) %>% 
    group_by(blk_region, CellType) %>% 
    summarise(mean = mean(value*0.325)) %>% 
    ggplot(aes(x = mean, y = CellType, fill = CellType)) + 
    geom_boxplot() + 
    geom_point(pch = 21, size = 3, alpha = 0.5,
               position = position_jitterdodge(jitter.width = 2.5)) +
    scale_x_log10() +
    style + 
    scale_fill_d3() + coord_flip() + theme(legend.position = "bottom")
  
  
  
  meta <- data %>% group_by(Cell) %>% summarise(value = mean(value)) %>% column_to_rownames("Cell")
  Sobj_Immu <- AddMetaData(Sobj_Immu, meta)
  Sobj_Immu <- AddMetaData(Sobj_Immu, as.data.frame(Sobj_Immu@reductions$umap@cell.embeddings))
  

  
  p2 <- Sobj_Immu@meta.data %>% 
    ggplot(aes(x = UMAP_1, 
               y = UMAP_2, 
               size = (value*0.325), 
               color = CellType)) + 
    geom_point(alpha = 0.5, pch = 21, 
               fill = "white", 
               stroke = 1) + 
    scale_size(range = c(0,20)) +
    #scale_fill_d3() +
    scale_color_d3() +
    style + 
    theme(axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          legend.position = "none")
}
 



Macrophage_Proportion_Plot {

nums <- table(Sobj_All$CellType, Sobj_All$block_region) %>% as.data.frame() %>% rename(n = "Freq")

prop.table(table(Sobj_All$CellType, Sobj_All$block_region), margin = 2) %>% 
  as.data.frame() %>% 
  left_join(Anno, by = c("Var2" = "block_region")) %>% 
  left_join(nums) %>% 
  left_join(distinct(select(Sobj_All@meta.data, Var2 = block_region, block))) %>% 

  filter(Var1 == "Macrophage" & n > 15) %>% 
  
  ggplot(aes(x = Type, y = Freq, fill = Type)) + 
  geom_boxplot(outlier.size = 0) + 
  ggtitle("Macrophage Proportion") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
  geom_point(aes(color = block), 
             position = position_dodge(width = 0.25), 
             size = 4) + 
  scale_y_log10() +
  scale_color_d3() +
  facet_wrap("Var1", nrow = 1) +
  scale_fill_manual(values = GID_cols) +
  style +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
  
}



Epithelial_Proportion_Plot {
  
  nums <- table(Sobj_All$CellClass, Sobj_All$block_region) %>% as.data.frame() %>% rename(n = "Freq")
  
  p <- prop.table(table(Sobj_All$CellClass, Sobj_All$block_region), margin = 2) %>% 
    as.data.frame() %>% 
    left_join(Anno, by = c("Var2" = "block_region")) %>% 
    left_join(nums) %>% 
    left_join(distinct(select(Sobj_All@meta.data, Var2 = block_region, block))) %>% 
    
    #filter(Var1 %in% c("Basal", "Luminal", "Luminal_HR+") & n > 10) %>% 
    filter(#Var1 %in% c("Basal", "Luminal", "Luminal_HR+") & 
             n > 10) %>% 
    
    ggplot(aes(x = Type, y = Freq, fill = Type)) + 
    geom_boxplot(outlier.size = 0) + 
    ggtitle("Epithelial Proportion") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
    geom_point(aes(color = block), 
               position = position_dodge(width = 0.25), 
               size = 4) + 
    scale_y_log10() +
    scale_color_d3() +
    facet_wrap("Var1", nrow = 1) +
    scale_fill_manual(values = GID_cols) +
    style +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
  
}


vasculature_Proportion_Plot {
  
  nums <- table(Sobj_All$CellType, Sobj_All$block_region) %>% as.data.frame() %>% rename(n = "Freq")
  
  Idents(Sobj_All) <- "Group"
  Sobj <- subset(Sobj_All, idents = "Vasculature")
  Idents(Sobj) <- "CellType"
  Sobj$CellType <- Idents(Sobj)
  
  include <- table(Sobj$CellClass, Sobj$block_region) %>% as.data.frame() %>% filter(Var1 %in% c("Endothelial", "Lymph_EC")) %>% 
    filter(Freq > 50) %>% pull(Var2) %>% unique()
  #include <- table(Sobj$block_region) %>% as.data.frame() %>% filter(Freq > 10) %>% pull(Var1) %>% unique()
  
  p <- prop.table(table(Sobj$CellType, Sobj$block_region), margin = 2) %>% 
    as.data.frame() %>% 
    left_join(Anno, by = c("Var2" = "block_region")) %>% 
    #left_join(nums) %>% 
    left_join(distinct(select(Sobj_All@meta.data, Var2 = block_region, block))) %>% 
    filter(Var2 %in% include) %>% 
    #filter(Var1 %in% c("Basal", "Luminal", "Luminal_HR+") & n > 10) %>% 
    #filter(#Var1 %in% c("Basal", "Luminal", "Luminal_HR+") & n > 10) %>% 
    
    ggplot(aes(x = Type, y = Freq, fill = Type)) + 
    geom_boxplot(outlier.size = 0) + 
    ggtitle("Vasculature Proportion") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
    geom_point(aes(color = block), 
               position = position_dodge(width = 0.25), 
               size = 4) + 
    #scale_y_log10() +
    scale_color_d3() +
    facet_wrap("Var1", nrow = 1) +
    scale_fill_manual(values = GID_cols) +
    style +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
  
}
     
      

nums <- table(Sobj_Immu$CellType, Sobj_Immu$block_region) %>% as.data.frame() %>% rename(n = "Freq")

p <- prop.table(table(Sobj_Immu$CellType, Sobj_Immu$block_region), margin = 2) %>% 
  as.data.frame() %>% 
  left_join(Anno, by = c("Var2" = "block_region")) %>% 
  left_join(nums) %>% 
  left_join(distinct(select(Sobj_Immu@meta.data, Var2 = block_region, block))) %>% 
  
  filter(n > 5) %>% 
  
  ggplot(aes(x = Type, y = Freq, fill = Type)) + 
  geom_boxplot(outlier.size = 0) + 
  ggtitle("Immune Proportion") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
  geom_point(aes(color = block), 
             position = position_dodge(width = 0.25), 
             size = 4) + 
  #scale_y_log10() +
  scale_color_d3() +
  facet_wrap("Var1", nrow = 1) +
  scale_fill_manual(values = GID_cols) +
  style +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))







plot_data <- Sobj_All@meta.data %>% filter(blk_region == "B8_reg15")

plot_data %>% 
  filter(Group == "Epithelial") %>% 
  ggplot(aes(x = X, y = Y)) + 
  geom_point(color = "grey60", size = 1) + 
  geom_point(data = filter(plot_data, Group == "Immune"), 
             aes(color = CellType)) +
  scale_color_d3() +
  style











#### working with affinity propagation

ap_group <- readRDS("2021_ReRun/CODEX/output/AffinityPropagation_Group_cellClusterIDs.RDS")

data <- 
  ap_group %>% filter(str_detect(block_region, "TGB3_reg10")#, Group == "Epithelial", AP.CT.Count > 35
                      ) #%>% 
    #group_by(block_region, AP.CT) %>% mutate(n = n()) %>% filter(n > 10)


data[ , ]                                %>% ggplot(aes(x = X, y = Y, color = Group)) + geom_point() + facet_wrap("block_region") + scale_color_d3() + style + theme(legend.position = "bottom", axis.title = element_blank(), axis.text.y = element_text(angle = 90)) + coord_fixed()
p <- data[c(rownames(ct_mat), inperimeter), ] %>% ggplot(aes(x = X, y = Y, color = Group)) + geom_point() + facet_wrap("block_region") + scale_color_d3() + style + theme(legend.position = "bottom", axis.title = element_blank(), axis.text.y = element_text(angle = 90)) + coord_fixed()



data[rownames(ct_mat) , ] %>% ggplot(aes(x = X, y = Y, color = AP.CT)) + geom_point() + geom_point() + facet_wrap("block_region") + style


data[rownames(ct_mat) , ] %>% ggplot(aes(x = X, y = Y, color = AP.CT)) + geom_point() + geom_point() + geom_point(data = filter(data, AP.CT.Examplar == "TRUE"), color = "black") + facet_wrap("block_region") + style+ theme(legend.position = "bottom")
data[c(inperimeter, rownames(ct_mat)), ] %>% ggplot(aes(x = X, y = Y, color = Group)) + geom_point() + geom_point() + facet_wrap("block_region") + scale_color_d3() + style + theme(legend.position = "bottom")











Anno <- Anno %>% filter(Type != "orientation")
banno <- Sobj_All@meta.data %>% select(region = "block_region", block) %>% distinct()



acini_morphology <- read_csv("2021_ReRun/IHC/acini_modfhology_v4.csv")
acini_morphology <- read_csv("2021_ReRun/IHC/acini_morphology_v5.csv")
acini_morphology <- read_csv("2021_ReRun/IHC/acini_morphology_v6.csv")
acini_morphology <- read_csv("2021_ReRun/CODEX/output/acini_morphology_v8.csv")

a <- acini_morphology %>% 
  left_join(Anno, by = c("region" = "block_region")) %>% 
  group_by(region, Type) %>% 
  filter(!is.na(Type)) %>% 
  left_join(banno) %>% 
  group_by(block, region, Type, Sample) %>% 
  mutate(nAcini = n()) %>% 
  filter(#acta2_pct > 0.05, 
    nAcini > 10) %>% 
  summarise(mean = mean(acta2_pct)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Type, y = mean, fill = Type)) +
  geom_boxplot() +
  geom_point(aes(color = block),
             #pch = 21, 
             size = 2, 
             #alpha = 0.5, #fill = "grey50", 
             position = position_dodge(width = 0.25)) +
  scale_color_d3() +
  ylab(paste0("avg. ACTA2 coverage")) +
  stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
  scale_fill_manual(values = GID_cols) +
  ggtitle("ACTA2 coverage") +
  style + 
  theme(axis.text.y = element_text(angle = 90), 
        legend.position = "none")


b <- acini_morphology %>% 
  left_join(Anno, by = c("region" = "block_region")) %>% 
  group_by(region, Type) %>% 
  filter(!is.na(Type)) %>% 
  left_join(banno) %>% 
  group_by(block, region, Type, Sample) %>% 
  mutate(nAcini = n()) %>% 
  filter(#acta2_pct > 0.05, 
    nAcini > 10) %>% 
  summarise(mean = mean(area*0.325)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Type, y = mean, fill = Type)) +
  geom_boxplot() +
  geom_point(aes(color = block),
             #pch = 21, 
             size = 2, 
             #alpha = 0.5, #fill = "grey50", 
             position = position_dodge(width = 0.25)) +
  scale_color_d3() +
  ylab(paste0("avg. Acini KRT Area")) +
  stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
  scale_fill_manual(values = GID_cols) +
  ggtitle("Acini Area") +
  style + 
  theme(axis.text.y = element_text(angle = 90), 
        legend.position = "none")



acini_morphology %>% 
  left_join(Anno, by = c("region" = "block_region")) %>% 
  group_by(region, Type) %>% 
  filter(!is.na(Type)) %>% 
  group_by(region, Type, Sample) %>% 
  mutate(nAcini = n()) %>% 
  filter(#acta2_pct > 0.05, 
    nAcini > 10) %>% 
  summarise(mean = mean(Immune)) %>% 
  
  ggplot(aes(x = Type, y = mean, fill = Type)) +
  geom_boxplot() +
  geom_point(pch = 21, size = 3, alpha = 0.5, #fill = "grey50", 
             position = position_dodge(width = 0.75)) +
  ylab(paste0("percent Immune Cells")) +
  stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
  scale_fill_manual(values = GID_cols) +
  ggtitle("Immune Cells in Acini") +
  style








p <- Sobj_All@meta.data %>% 
  filter(Group == "Epithelial", CellClass != "Epi_Other") %>% 
  group_by(Type, block, block_region, CellClass) %>% 
  summarise(mean = mean(Size*0.325), n = n()) %>% 
  filter(n > 10, mean < 283) %>% 
  ggplot(aes(x = Type, y = mean, fill = Type)) + 
  geom_boxplot(outlier.size = 0) + 
  geom_point(aes(color = block),
    #pch = 21, 
    size = 2, 
    #alpha = 0.5, #fill = "grey50", 
    position = position_dodge(width = 0.25)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
  scale_color_d3() +
  ylab("avg. Nucleus Area") +
  scale_fill_manual(values = GID_cols) +
  facet_wrap("CellClass", nrow = 1) + 
  style




