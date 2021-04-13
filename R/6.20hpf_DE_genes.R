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
sample.markers <- FindAllMarkers(object = sample, only.pos = TRUE, min.pct = 0.18,  min.diff.pct = 0.15)

# top two genes for each cluster
markers <- sample.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
write.csv(markers, file = paste0(sample.name, "_top2_markers.csv"))

for(c in seq(from= 0, to= 28)){
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