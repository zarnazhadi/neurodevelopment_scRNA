library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS("hpf20_tSNE.cluster.rds")
sample.name <- "20hpf"

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
sample.markers <- FindAllMarkers(object = sample, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
markers <- sample.markers %>% group_by(cluster) %>% top_n(21915, avg_log2FC)
write.csv(markers, file = paste0(sample.name, "_all_markers.csv"))

# top two genes for each cluster
markers <- sample.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
write.csv(markers, file = paste0(sample.name, "_top2_markers.csv"))

# find min p_val_adj
#non_zero <- markers[apply(markers!=0, 1, all),]
#min(non_zero) # 6.939244e-293

# top ten genes for each cluster
markers_10 <- sample.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(markers_10, file = paste0(sample.name, "_top10_markers.csv"))

# find min p_val_adj
#non_zero <- markers[apply(markers!=0, 1, all),]
#min(non_zero) # 6.939244e-293

for(c in seq(from= 0, to= 40)){
  top2 <- markers[markers$cluster == c,]$gene
  print(top2)
  
  pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/cluster", c,"_top2_markers_VLN_plot.pdf"))
  p <- VlnPlot(object = sample, features = top2[1]) 
  p <- p + theme(axis.text=element_text(size=9),axis.title=element_text(size=14))
  print(p)
  p <- VlnPlot(object = sample, features = top2[2]) 
  p <- p + theme(axis.text=element_text(size=9),axis.title=element_text(size=14))
  print(p)
  dev.off()
  
  pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/cluster", c, "_top2_markers_feature_plot.pdf"))
  p <- FeaturePlot(object = sample, features = top2)
  print (p)
  dev.off()
}