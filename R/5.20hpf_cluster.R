library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)
library(patchwork)

#------------------------------------———CLUSTER CELLS—————————————————————---------------------------- 
setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS("hpf20_jackstraw.cluster.rds")
sample.name <- "20hpf"

# find k-nearest neighbours
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:80)

# cluster cells
sample <- FindClusters(sample, resolution = 0.8, algorithm = 1)

#res <- 0
#for(r in seq(1:15)){
  #res <- res + 0.2
  #sample <- FindClusters(sample, resolution = res, algorithm = 1)
  #sample <- RunTSNE(object = sample, dims.use = 1:80, do.fast = TRUE) # Louvain algorithm
  #pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/tSNE_", r, ".pdf"))
  #p <- DimPlot(object = sample, reduction = "tsne") + plot_annotation(title = paste0("T-SNE Clusters for 80 PCs (20 hpf)", " res = ", res))
  #print(p)
  #dev.off()
#}

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-clustering/20hpf_jack_straw_plot.pdf")
JackStrawPlot(object = sample, dims = 1:80, reduction = "pca", xmax = 0.0025) + NoLegend()
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-clustering/20hpf_elbow_plot.pdf")
ElbowPlot(object = sample, ndims= 80)
dev.off()

sample <- RunTSNE(object = sample, dims.use = 1:80, do.fast = TRUE)

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-clustering/20hpf_tSNE_plot.pdf")
DimPlot(object = sample, reduction = "tsne") + plot_annotation(title = 'T-SNE Clusters for 80 PCs (20 hpf), resolution= 0.8')
dev.off()

# save as .rds file
saveRDS(sample, file = "hpf20_tSNE.cluster.rds")