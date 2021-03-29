library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
sample <- readRDS("hpf18.cluster.rds")

sample <- RunTSNE(object = sample, dims.use = 1:10, do.fast = TRUE)

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_tSNE_plot.pdf")
DimPlot(object = sample, reduction = "tsne")
dev.off()

sample <- RunUMAP(sample, reduction = "pca", dims = 1:20)

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_UMAP.pdf")
DimPlot(sample, reduction = "umap", split.by = "seurat_clusters")
dev.off()