library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
s <- readRDS("hpf18.cluster.rds")

s <- RunTSNE(object = s, dims.use = 1:10, do.fast = TRUE)

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_tSNE_plot.pdf")
DimPlot(object = s, reduction = "tsne")
dev.off()

s <- RunUMAP(s, reduction = "pca", dims = 1:20)

# this doesnt work :(
#pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_UMAP.pdf")
#DimPlot(s, reduction = "umap", split.by = "seurat_clusters")
#dev.off()

# this works :)
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_UMAP.pdf")
DimPlot(s, reduction = "umap", split.by = "orig.ident")
dev.off()