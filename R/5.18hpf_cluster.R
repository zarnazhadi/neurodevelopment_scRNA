library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)
library(patchwork)

#------------------------------------———CLUSTER CELLS—————————————————————---------------------------- 
setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS("hpf18_jackstraw.cluster.rds")
sample.name <- "18hpf"

# find k-nearest neighbours
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:80)

# cluster cells
sample <- FindClusters(sample, resolution = 0.6, algorithm = 1)

# loop to find the best resolution
#res <- 0
#for(r in seq(1:15)){
  #res <- res + 0.2
  #sample <- FindClusters(sample, resolution = res, algorithm = 1)  
  #sample <- RunTSNE(object = sample, dims.use = 1:80, do.fast = TRUE) 
  #pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/tSNE_", r, ".pdf"))
  #p <- DimPlot(object = sample, reduction = "tsne") + plot_annotation(title = paste0("T-SNE Clusters for 80 PCs (18 hpf)", " res = ", res))
  #print(p)
  #dev.off()
#}

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_jack_straw_plot.pdf")
JackStrawPlot(object = sample, dims = 1:80, reduction = "pca", xmax = 0.0025) + NoLegend()
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_elbow_plot.pdf")
ElbowPlot(object = sample, ndims= 80)
dev.off()

sample <- RunTSNE(object = sample, dims.use = 1:80, do.fast = TRUE)

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_tSNE_plot.pdf")
DimPlot(object = sample, reduction = "tsne") + plot_annotation(title = 'T-SNE Clusters for 80 PCs (18 hpf), resolution= 0.6')
dev.off()

# save as .rds file
saveRDS(sample, file = "hpf18_tSNE.cluster.rds")