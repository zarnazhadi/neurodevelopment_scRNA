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

# this doesnt work :(
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-clustering/18hpf_UMAP.pdf")
DimPlot(sample, reduction = "umap", split.by = "seurat_clusters") + RotatedAxis()
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/other/other.pdf")
#DotPlot(sample, features = "nCount_RNA", split.by = "seurat_clusters") + RotatedAxis()

#DoHeatmap(sample, features = VariableFeatures(sample)[1:100], cells = 1:500, size = 1, 
          #angle = 90) + NoLegend()

#DimPlot(sample)

#sample.no.tsne <- sample
#sample.no.tsne[["tsne"]] <- NULL
#DimPlot(sample.no.tsne) + RotatedAxis()

DimPlot(sample, reduction = "umap")
dev.off()