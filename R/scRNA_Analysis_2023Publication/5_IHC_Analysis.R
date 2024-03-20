library(readxl)
library(tidyverse)
library(rdist)
library(readr)
library(ggpubr)
library(viridis)


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
  

  
  setwd("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/")
  
  Sample_Code <- read_excel("2021_ReRun/IHC/Sample_Code.xlsx") %>% separate(Sample, into = c("Type", "Num"), sep = "-", extra = "merge", remove = F) %>% select(-Num)
  
  `%notin%` <- Negate(`%in%`)
  
  SP_levels <- c("CF-3920", "CF-1380", "CF-7780", "CF-2797", "CF-318-813", "CF-0404", "CF-428-112", "CF-2099", "CF-4014", "TM-9469", "TM-1956", "TM-6544", "TM-6477", "TM-8249", "TM-7567", "TM-2768", "TM-3937", "TM-9817")
  SP_cols   <- c("#F26C66", "#FFBB78", "#F57F46", "#7FC210", "#2CA02C",    "#BCBD22", "#EDBD1F",    "#FF9896", "#C49C94", "#9467BD", "#C5B0D5", "#17BECF", "#C5D4D2", "#438ABB", "#FA7AA3", "#F7B6D2", "#A8486A", "#3FB8AF")
  GID_cols  <- c("#A6499B", "#FAA42F")
  
  my_comparisons <- list(c("CF", "TM"))
  
  style <- theme(text = element_text(family = "Lato", size = 25),  
                 title = element_text(size = 15, face = "bold"), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = "white", color = "grey35", size = 2))

}

save_plot(p, "Adipocyte_Vacuole_Size", w = 5, h = 9, scale = 1, svg = T)

Scene_List <- read_excel("2021_ReRun/IHC/Scene_List.xlsx", col_names = FALSE)

Scene_Clean <- 
  Scene_List %>% 
  filter(str_detect(...1, ".czi."), !str_detect(...1, "-CT"), !str_detect(...1, "-Ct")) %>% 
  separate(...1, into = c("Sample", "Scene"), sep = ".czi.") %>%
  mutate(Scene = str_remove(Scene, ".tif")) %>% 
  mutate(set = substr(Sample, nchar(Sample)-0 , nchar(Sample)), 
         Sample = substr(Sample, 1 , nchar(Sample)-2)) %>% 
  mutate(Sample = str_replace(Sample, "CF318-813", "CF-318-813"), 
         Sample = str_replace(Sample, "CF428-112", "CF-428-112"), 
         Sample = str_replace(Sample, "S1-934120", "S19-34120")) %>% 
  mutate(x = "set") %>% unite(set, x, set, sep = "") %>% 
  rename(Scene_Czi = "Scene")

table(Scene_Clean$Sample, Scene_Clean$set)





sets <- list.files("2021_ReRun/IHC/tables/")

set_list    <- vector("list", length = length(sets))
for (i in 1:length(sets)) {
  
  print(paste("working on:", sets[i]))
  
  files <- list.files(paste0("2021_ReRun/IHC/tables/", sets[i])) %>% 
    as_tibble() %>% 
    mutate(set = sets[i]) %>% 
    mutate(sample = str_sub(value, 1, nchar(value)-10)) %>% 
    separate(set, into = c("set", "markers"), extra = "merge") %>% 
    distinct() %>% 
    rename(file = "value")
  
  file_list    <- vector("list", length = length(files$file))
  for (j in 1:length(files$file)) {
    
    print(paste("working on:", sets[i], "and", files$file[j]))

    table <- read_csv(paste0("2021_ReRun/IHC/tables/", sets[i], "/", files$file[j]), show_col_types = FALSE)

    ann <- files[j ,] %>% mutate(file_name = unique(table$file_name))
    file_list[[j]] <- table %>% select(-sample) %>% left_join(ann, by = "file_name") %>% select(-file_name) %>% rename(Cell = "...1")
  
  }
  set_list[[i]] <- file_list %>% bind_rows()
}
#set_list %>% saveRDS("2021_ReRun/IHC/set_list.rds")  
set_list <- readRDS("2021_ReRun/IHC/set_list.rds")  


filter_small_scenes_and_prep_for_scaling {
table <- bind_rows(set_list) %>% filter(!is.na(DAPI)) %>% group_by(set, sample, scene) %>% mutate(nNucs = n()) %>% ungroup() %>% filter(nNucs > 500)

table %>% filter(set == "set1") %>% select(-nNucs) %>% pivot_longer(cols = c("DAPI", "AF488", "Cy3", "Cy5"), names_to = "channel") %>% write_tsv("2021_ReRun/IHC/ihc_set1_table_long.tsv.gz")
table %>% filter(set == "set2") %>% select(-nNucs) %>% pivot_longer(cols = c("DAPI", "AF488", "Cy3", "Cy5"), names_to = "channel") %>% write_tsv("2021_ReRun/IHC/ihc_set2_table_long.tsv.gz")
table %>% filter(set == "set3") %>% select(-nNucs) %>% pivot_longer(cols = c("DAPI", "AF488", "Cy3"),        names_to = "channel") %>% write_tsv("2021_ReRun/IHC/ihc_set3_table_long.tsv.gz")
table %>% filter(set == "set4") %>% select(-nNucs) %>% pivot_longer(cols = c("DAPI", "AF488", "Cy3"),        names_to = "channel") %>% write_tsv("2021_ReRun/IHC/ihc_set4_table_long.tsv.gz")
table %>% filter(set == "set5") %>% select(-nNucs) %>% pivot_longer(cols = c("DAPI", "AF488", "Cy3"),        names_to = "channel") %>% write_tsv("2021_ReRun/IHC/ihc_set5_table_long.tsv.gz")
table %>% filter(set == "set6") %>% select(-nNucs) %>% pivot_longer(cols = c("DAPI", "AF488", "Cy3", "Cy5"), names_to = "channel") %>% write_tsv("2021_ReRun/IHC/ihc_set6_table_long.tsv.gz")
table %>% filter(set == "set7") %>% select(-nNucs) %>% pivot_longer(cols = c("DAPI", "AF488", "Cy3", "Cy5"), names_to = "channel") %>% write_tsv("2021_ReRun/IHC/ihc_set7_table_long.tsv.gz")

}


"cd ~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/IHC"
"conda activate codex"
paste0("python3 analyze_ihc.py ihc_set3_table_long.tsv.gz scaled_output set3")

scaled_data_set3 <- 
  read_delim("2021_ReRun/IHC/scaled_output/scaled_data_set5.tsv.gz", 
             delim = "\t", escape_double = FALSE, 
             col_types = cols(...1 = col_skip()), 
             trim_ws = TRUE) %>% 
  select(-value, -Robust.Scaled.DAPI, -QuantileTransformed) %>% 
  left_join(Sample_Code, by = c("sample" = "Sample_IHC"))


set_TCF4  <- scaled_data_set3 %>% pivot_wider(names_from = "channel", values_from = "Robust.Scaled") %>% mutate(Class = ifelse(AF488 > 1, "Adipose", "Other"))
set_INSR  <- scaled_data_set3 %>% pivot_wider(names_from = "channel", values_from = "Robust.Scaled") %>% mutate(Class = ifelse(AF488 > 1, "Adipose", "Other"))
set_NR4A1 <- scaled_data_set3 %>% pivot_wider(names_from = "channel", values_from = "Robust.Scaled") %>% mutate(Class = ifelse(Cy3 > 1, "Adipose", "Other"))


adipose_boxplots {

  
  sum_TCF4 <- 
    set_TCF4 %>% mutate(ratio = Cy3 / DAPI) %>% 
    group_by(set, Sample, Type, scene, Class) %>% 
    summarise(n = n(), median_TCF4 = median(Cy3))
  
  sum_TCF4 %>% ungroup() %>% 
    filter(n > 100, 
           Sample %notin% c("CF-318-813", "CF-428-112")
    ) %>% 
    ggplot(aes(x = Type, y = 1+median_TCF4, fill = Type)) + 
    geom_boxplot(outlier.alpha = 0) + 
    geom_point(aes(color = factor(Sample, levels = SP_levels[c(1:4,6,8:18)])), 
               size = 2, 
               position = position_jitter(width = 0.1)) + 
    
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = SP_cols[c(1:4,6,8:18)]) +
    
    stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
    
    facet_wrap("Class") + 
    scale_y_log10() +
    style + 
    theme(legend.position = "none")
  
  
  
  sum_INSR <- 
    set_INSR %>% group_by(set, Sample, Type, scene, Class) %>% 
    summarise(n = n(), median_INSR = median(Cy3))
  
  b <- sum_INSR %>% ungroup() %>% 
    filter(n > 100, 
           Sample %notin% c("CF-318-813", "CF-428-112")
    ) %>% 
    ggplot(aes(x = Type, y = 1+median_INSR, fill = Type)) + 
    geom_boxplot(outlier.alpha = 0) + 
    geom_point(aes(color = factor(Sample, levels = SP_levels[c(1:4,6,8:18)])), 
               size = 2, 
               position = position_jitter(width = 0.1)) + 
    
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = SP_cols[c(1:4,6,8:18)]) +
    
    stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
    
    facet_wrap("Class") + 
    scale_y_log10() +
    style + 
    theme(legend.position = "none")
  
  
  sum_NR4A <- 
    set_NR4A1 %>% group_by(set, Sample, Type, scene, Class) %>% 
      summarise(n = n(), median_N4A = median(AF488))
  
  c <- sum_NR4A %>% ungroup() %>% 
    filter(n > 100, 
           Sample %notin% c("CF-318-813", "CF-428-112")
           ) %>% 
    ggplot(aes(x = Type, y = 1+median_N4A, fill = Type)) + 
    geom_boxplot(outlier.alpha = 0) + 
    geom_point(aes(color = factor(Sample, levels = SP_levels[c(1:4,6,8:18)])), 
               size = 2, 
               position = position_jitter(width = 0.1)) + 
    
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = SP_cols[c(1:4,6,8:18)]) +
    
    stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F) +
    
    facet_wrap("Class") + 
    scale_y_log10() +
    style + 
    theme(legend.position = "none")


  p <- a + b + c
  
}




AffinityPropagation_IHC %>% 
  as_tibble() %>% 
  filter(Sample == "TM-9817", 
         scene == "ScanRegion18", 
         Class == "Epithelial") %>% 
  ggplot(aes(x = X, y = Y, color = AP.CT)) + 
  geom_point() +
  facet_wrap("scene")
  
  
  
  





neighborhood_absolute_distance {
  
  ##############################################
  ##### Neighborhood Absolute Distance
  library(tidyverse)
  library(rdist)
  set3_wide <- readRDS("2021_ReRun/IHC/set7_wide.rds")
  AffinityPropagation_IHC <- readRDS("2021_ReRun/IHC/AffinityPropagation_Class_cellClusterIDs_IHC.RDS")
  
  AffinityPropagation_IHC <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/IHC/AffinityPropagation_Class_cellClusterIDs_IHC.RDS")
  set3_wide <- scaled_data_set3 %>% pivot_wider(names_from = "channel", values_from = "Robust.Scaled") %>% unite(sample_scene, Sample, scene, remove = F)
  
  distance = 50
  
  Anno <- dist_data %>% select(Targets = "Cell", Class)
    
  ap_num = 20
  
  ap_cells <- 
    AffinityPropagation_IHC %>% 
    as_tibble() %>%
    filter(Class == "Epithelial") %>% 
    group_by(sample_scene, Class, AP.CT) %>% 
    mutate(n = n()) %>% 
    arrange(n) %>% 
    filter(n > ap_num) %>% 
    pull(Cell)
  
  regions <- 
    AffinityPropagation_IHC %>% 
    as_tibble() %>%
    filter(Class == "Epithelial") %>% 
    group_by(sample_scene, Class, AP.CT) %>% 
    mutate(n = n()) %>% 
    arrange(-n) %>% 
    filter(n > ap_num) %>% 
    pull(sample_scene) %>% 
    unique()
  
  
  celltype <- "Epithelial"
  
  minima_list   <- vector("list", length = length(regions))
  for (i in 1:length(regions)) {
    
    print(paste("working on region ", i, " of ", length(regions)))
    
    ct_mat <- dist_data %>% filter(Class == celltype, sample_scene == regions[i]) %>% select(Cell, X, Y) %>% column_to_rownames("Cell")
    tg_mat <- dist_data %>% filter(Class %in% c("Immune", "T-Cell"), sample_scene == regions[i]) %>% select(Cell, X, Y) %>% column_to_rownames("Cell")
    
    
    dist <- cdist(ct_mat, tg_mat, metric = "euclidean", p = 2)
    rownames(dist) <- rownames(ct_mat)
    colnames(dist) <- rownames(tg_mat)
    
    
    dist_long <- 
      dist %>% 
      as.data.frame() %>% 
      rownames_to_column("Cell") %>% 
      pivot_longer(cols = rownames(tg_mat), names_to = "Targets") %>% 
      left_join(Anno, by = "Targets")
    
    
    inperimeter <- 
      dist_long %>% 
      filter(Cell %in% ap_cells) %>% 
      #group_by(Cell, CellClass) %>%
      #mutate(n.Target = n()) %>% 
      group_by(Cell) %>% 
      #filter(value != 0) %>% 
      filter(value < distance) %>% 
      pull(Targets) %>% 
      unique()
    
    
    stats <- 
      dist_long %>% select(-Cell, -value) %>% 
      distinct() %>% 
      filter(Targets %in% inperimeter) %>% 
      group_by(Class) %>% 
      summarise(inperi = n(), .groups = "drop") %>% 
      mutate(population = sum(inperi), 
             n.Ct.OI = length(intersect(rownames(ct_mat), ap_cells))) %>% 
      mutate(peri.prop = (inperi/population)*100, 
             ct.peri.ratio = n.Ct.OI/population) %>% 
      mutate(Ct.OI = celltype,
             sample_scene = regions[i],
             comparison = "Class") %>% 
      rename(Target.ID  = "Class")
    
    
    minima_list[[i]] <- stats
  } # end loop i
    
  
  minima_list_dist50 <- readRDS("~/Box/Knott_Lab/Flo/Projects/Transgender_Nuclei/2021_ReRun/IHC/minima_list_dist50.rds")
  
  anno <- dist_data %>% select(sample_scene, Type) %>% distinct()
  
  # plot the perimeter ratios
  p <- bind_rows(minima_list_dist50) %>% 
    left_join(anno) %>% 
    select(sample_scene, Type, population, Target.ID, inperi) %>% 
    pivot_wider(names_from = Target.ID, values_from = inperi) %>% 
    rename(TCell = `T-Cell`) %>% 
    filter(!is.na(TCell), !is.na(Immune)) %>% 
    filter(Immune > 20, TCell > 10) %>% 
    mutate(ratio = TCell/Immune) %>% 
    separate(sample_scene, into = c("Sample", "scene"), sep = "_", remove = F) %>% 
    ggplot(aes(x = Type, y = ratio, fill = Type)) + 
    geom_boxplot(outlier.alpha = 0) + 
    geom_point(aes(color = factor(Sample, levels = SP_levels[c(1:4,6,8:18)])), size = 2,
               position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.25)) +
    ylab(paste0("Ratio of CD3+ / CD3- immune cells in perimeter")) +
    stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
    ggtitle("CD3 enrichment in Epithelium") +
    style +
    scale_fill_manual(values = GID_cols) +
    scale_color_manual(values = SP_cols[c(1:4,6,8:18)]) + 
    theme(legend.position = "none")
  
  
  
  
  
  

  
}






Lymphoid_Cell_Ratio { 

set3_wide <- scaled_data_set3 %>% pivot_wider(names_from = "channel", values_from = "Robust.Scaled") %>% unite(sample_scene, Sample, scene, remove = F)

dist_data <- set3_wide %>% 
  filter(percent_rank(AF488) < 0.999) %>% 
  #filter(percent_rank(Cy3)   < 0.999) %>% 
  filter(percent_rank(Cy5)   < 0.99) %>% 
  mutate(Class = ifelse(Cy5 > 1, "Epithelial", "Other")) %>% 
  #mutate(Class = ifelse(Class == "Other" & AF488 > 1.25, "Immune", Class)) %>% 
  mutate(Class = ifelse(AF488 > 1, "Immune", Class)) %>% 
  mutate(Class = ifelse(Class == "Immune" & Cy3 > 3, "T-Cell", Class))


p <- dist_data %>% 
  group_by(Sample, Type, scene, Class) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = Class, values_from = n) %>% 
  mutate(ratio = `T-Cell`/Immune) %>% 
  filter(Immune > 10) %>% 
  ggplot(aes(x = Type, y = ratio, fill = Type)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(aes(color = factor(Sample, levels = SP_levels[c(1:4,6,8:18)])), 
             position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.25)) +
  stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
  style +
  scale_fill_manual(values = GID_cols) +
  scale_color_manual(values = SP_cols[c(1:4,6,8:18)]) + 
  theme(legend.position = "none")
  
  
  
}


Macrophage_Ratio { 
  
  macro_wide <- scaled_data_set3 %>% pivot_wider(names_from = "channel", values_from = "Robust.Scaled") %>% unite(sample_scene, Sample, scene, remove = F)
  
  dist_data_mac <- macro_wide %>% 
    filter(percent_rank(AF488) < 0.999) %>% 
    #filter(percent_rank(Cy3)   < 0.999) %>% 
    filter(percent_rank(Cy5)   < 0.99) %>% 
    mutate(Class = ifelse(Cy5 > 1, "Epithelial", "Other")) %>% 
    mutate(Class = ifelse(Class == "Other" & AF488 > 1, "Immune", Class)) %>% 
    #mutate(Class = ifelse(AF488 > 1, "Immune", Class)) %>% 
    #mutate(Class = ifelse(Class == "Immune" & Cy3 > 1, "Macrophage", Class))
    mutate(Class = ifelse(Cy3 > 2, "Macrophage", Class))
  
  
  dist_data_mac %>% 
    group_by(Sample, Type, scene, Class) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = Class, values_from = n) %>% 
    mutate(ratio = Macrophage/Immune) %>% 
    filter(Macrophage > 10, ratio < 15) %>% 
    #filter(Immune > 100, ratio < 15) %>% 
    ggplot(aes(x = Type, y = ratio, fill = Type)) + 
    geom_boxplot() + 
    geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.25)) +
    stat_compare_means(comparisons = my_comparisons, paired = F, method = "wilcox") +
    style
  
}
  
  
  
  
  
  
  
  
  
p <- adipocyte_features_v2_plt_src %>% 
  filter(!is.na(type)) %>% 
  ggplot(aes(x = type, y= sqrt_area, fill = type)) + 
  geom_boxplot() + 
  geom_point(aes(color = factor(sample_name, levels = SP_levels[c(1:4,6,8:18)])), 
             position = position_jitter(width = 0.15)
             ) +
  scale_fill_manual(values = GID_cols) +
  scale_color_manual(values = SP_cols[c(1:4,6,8:18)]) +  
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", paired = F)+
  style + 
  theme(legend.position = "none")
  
  