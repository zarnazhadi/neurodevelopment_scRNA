library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS("hpf18_tSNE.cluster.rds")
sample.name <- "18hpf"

# find markers for the first cluster 
cluster1.markers <- FindMarkers(object = sample, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = sample, ident.1 = 2, ident.2 = c(0, 3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
sample.markers <- FindAllMarkers(object = sample, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

# top two genes for each cluster
markers <- sample.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
write.csv(markers, file = paste0(sample.name, "_top2_markers.csv"))

for(c in seq(from= 0, to= 33)){
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

# find classification power
cluster1.markers <- FindMarkers(object = sample, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)

# top10 markers per cluster 
top10 <- sample.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(markers, file = paste0(sample.name, "_top10_markers.csv"))

pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/top10_markers_feature_plot.pdf"))
p <- DoHeatmap(object = sample, features = top10$gene, label = TRUE) + 
  theme(text = element_text(size = 5))
print(p)
dev.off()