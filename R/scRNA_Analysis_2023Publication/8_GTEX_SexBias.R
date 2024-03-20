# align the Sample attributes with the SMTS identifier
pheno <- 
  GTEx_Analysis_v8_Annotations_SampleAttributesDS %>% 
  select(1, SMTS, SMTSD) %>% 
  separate(SAMPID, into = c("A", "B", "C"), extra = "merge", remove = F) %>% 
  unite(SAMPLEID, A, B, sep = "-") %>% 
  select(-C) %>% 
  left_join(GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS, by = c("SAMPLEID" = "SUBJID")) %>% 
  mutate(SEX = ifelse(SEX == 1, "M", "F"), SAMPID = str_replace_all(SAMPID, "-", "\\."))


# take mean of effect sizes to accomodate SMTS identifier
effect_size_collapsed <- 
  dat.gct %>% select(1,2) %>% 
  right_join(effect_size, by = c("Name" = "gene")) %>% 
  select(-Name) %>% pivot_longer(2:45) %>% 
  left_join(collapse_effsize_columns , by = c("name" = "value")) %>% mutate(value = replace_na(value, 0)) %>% 
  group_by(SMTS, Description) %>% 
  summarise(eff_size = mean(value))


AR_median <- 
  dat.gct %>% 
  select(Description, intersect(colnames(dat.gct), pheno$SAMPID)) %>% 
  filter(Description == "AR") %>% 
  pivot_longer(intersect(colnames(dat.gct), pheno$SAMPID), names_to = "SAMPID") %>% 
  left_join(pheno) %>% 
  group_by(SMTS) %>% 
  summarise(median = median(value)) %>% 
  arrange(-median) %>% 
  filter(SMTS %notin% c("Bladder", "Cervix Uteri", 
                        "Fallopian Tube", "Ovary", 
                        "Prostate", "Testis", "Uterus", 
                        "Vagina", "Cells - Leukemia cell line (CML)", 
                        "Kidney - Medulla"))

p1 <- 
  ggplot(AR_median, aes(x = factor(SMTS, levels = SMTS), y = median)) + geom_col() + 
  style + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



fc = 0.25

regulation_grouping <- 
  DE_Response_by_CellType %>% 
  filter(Cluster == "LUM_HR-pos", abs(avg_log2FC) > fc) %>% 
  mutate(Dir = ifelse(avg_log2FC > 0, "UP", "DO")) %>% 
  select(Gene, Dir)


mat <- 
  effect_size_collapsed %>% 
  left_join(regulation_grouping, by = c("Description" = "Gene")) %>% 
  filter(!is.na(Dir)) %>% 
  group_by(SMTS, Dir) %>% 
  summarise(mean = mean(eff_size)) %>% 
  pivot_wider(names_from = Dir, values_from = mean) %>% 
  column_to_rownames("SMTS")


breaksList = seq(-max(mat), max(mat), by = 2*max(mat)/100)

p5 <- 
pheatmap(t(mat)[ ,AR_median$SMTS], 
         cluster_cols = F,
         breaks = breaksList,
         color = colorRampPalette(brewer.pal(11, "BrBG"))(101)
         )

p <- as.ggplot(p5) + p1 + plot_layout(ncol = 1)



p %>% save_plot(name = "sexbias_heatmap", scale = 1, w = 16, h = 9, svg = T)









Nat_comm_sobj_Epi_Only <- readRDS("~/Box/Knott_Lab/Flo/Projects/Breast_Epithelial/Nat_comm_sobj_Epi_Only.rds")
Sentinel_Cell_Markers_from_Nuclei <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/Sentinel_Cells/Sentinel_Cell_Markers_from_Nuclei.rds")

Idents(Nat_comm_sobj_Epi_Only) <- "RNA_snn_res.0.05"
DimPlot(Nat_comm_sobj_Epi_Only, group.by = "RNA_snn_res.0.05")

Sobj <- subset(Nat_comm_sobj_Epi_Only, idents = 2)


{

cells     <- WhichCells(Sobj, expression = Cxcl13 > 0)
cells_neg <- WhichCells(Sobj, expression = Cxcl13 == 0, downsample = 22)


A <- 
Sobj@meta.data %>% 
  rownames_to_column("ID") %>% 
  filter(ID %in% cells) %>% 
  mutate(Cxcl13_Expr = "YES") %>% 
  select(ID, Cxcl13_Expr) %>% 
  column_to_rownames("ID")

B <- 
Sobj@meta.data %>% 
  rownames_to_column("ID") %>% 
  filter(ID %notin% c(cells, cells_neg)) %>% 
  mutate(Cxcl13_Expr = "NO") %>% 
  select(ID, Cxcl13_Expr) %>% 
  column_to_rownames("ID")

C <- 
  Sobj@meta.data %>% 
  rownames_to_column("ID") %>% 
  filter(ID %in% cells_neg) %>% 
  mutate(Cxcl13_Expr = "CNTR") %>% 
  select(ID, Cxcl13_Expr) %>% 
  column_to_rownames("ID")

meta <- bind_rows(A, B, C)

Sobj <- AddMetaData(Sobj, meta)



Idents(Sobj) <- "Cxcl13_Expr"

marks_cx  <- FindMarkers(Sobj, ident.1 = "YES",  ident.2 = "NO", test.use = "MAST") %>% rownames_to_column("gene")
marks_ctr <- FindMarkers(Sobj, ident.1 = "CNTR", ident.2 = "NO", test.use = "MAST") %>% rownames_to_column("gene")

nuclei_sentinel_lower <- Sentinel_Cell_Markers_from_Nuclei %>% select(gene = "Gene", avg_log2FC_nuclei = avg_log2FC) %>% mutate(gene = str_to_sentence(gene))

marks_cx %>% left_join(nuclei_sentinel_lower, by = "gene") %>% as_tibble() %>% filter(!is.na(avg_log2FC_nuclei)) %>% arrange(-avg_log2FC)
marks_ctr %>% left_join(nuclei_sentinel_lower, by = "gene") %>% as_tibble() %>% filter(!is.na(avg_log2FC_nuclei)) %>% arrange(-avg_log2FC)

}