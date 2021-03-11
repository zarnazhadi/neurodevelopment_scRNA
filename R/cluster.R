library("Matrix")
library(tidyverse)
library("Seurat")
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/appDir/data")

hpf18 <- readRDS("hpf18_combined_normalised.rds")
hpf20 <- readRDS("hpf20_combined_normalised.rds")

# find variable genes
hpf18 <- FindVariableFeatures(object = hpf18, mean.function = ExpMean, dispersion.function = LogVMR)
hpf20 <- FindVariableFeatures(object = hpf20, mean.function = ExpMean, dispersion.function = LogVMR)

head(x = HVFInfo(object = hpf18))
head(x = HVFInfo(object = hpf20))

# scaling data - worked after making all samples have the same feature name "percent.mito"
hpf18 <- ScaleData(object = hpf18, vars.to.regress = c("nCount_RNA", "percent.mito")) 
hpf20 <- ScaleData(object = hpf20, vars.to.regress = c("nCount_RNA", "percent.mito"))

# dimensional reduction
hpf18 <- RunPCA(object = hpf18,  npcs = 30, verbose = FALSE)
hpf20 <- RunPCA(object = hpf20,  npcs = 30, verbose = FALSE)

# 18 hpf
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_hmgn2_dimplot.pdf")
FeaturePlot(object = hpf18, features = "hmgn2")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_variable_feature_plot.pdf")
VariableFeaturePlot(object = hpf18)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_heat_map.pdf")
DimHeatmap(object = hpf18, reduction = "pca", cells = 200, balanced = TRUE)
dev.off()

# 20 hpf
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-normalisation/20hpf_hmgn2_dimplot.pdf")
FeaturePlot(object = hpf20, features = "hmgn2")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-normalisation/20hpf_variable_feature_plot.pdf")
VariableFeaturePlot(object = hpf20)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-normalisation/20hpf_heat_map.pdf")
DimHeatmap(object = hpf20, reduction = "pca", cells = 200, balanced = TRUE)
dev.off()

# determine statistically significant genes
hpf18 <- JackStraw(object = hpf18, reduction = "pca", dims = 20, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)
hpf18 <- ScoreJackStraw(object = hpf18, dims = 1:20, reduction = "pca")

hpf20 <- JackStraw(object = hpf20, reduction = "pca", dims = 20, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)
hpf20 <- ScoreJackStraw(object = hpf20, dims = 1:20, reduction = "pca")

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_jack_straw_plot.pdf")
JackStrawPlot(object = hpf18, dims = 1:20, reduction = "pca", xmax = 0.0025)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-normalisation/20hpf_jack_straw_plot.pdf")
JackStrawPlot(object = hpf20, dims = 1:20, reduction = "pca", xmax = 0.01)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/post-normalisation/18hpf_elbow_plot.pdf")
ElbowPlot(object = hpf18)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/post-normalisation/20hpf_elbow_plot.pdf")
ElbowPlot(object = hpf20)
dev.off()

# find k-nearest neighbours
hpf18 <- FindNeighbors(hpf18, reduction = "pca", dims = 1:20)
hpf20 <- FindNeighbors(hpf20, reduction = "pca", dims = 1:20)

# cluster cells
hpf18 <- FindClusters(hpf18, resolution = 0.5, algorithm = 1)
hpf20 <- FindClusters(hpf20, resolution = 0.5, algorithm = 1)

# save as .rds file
saveRDS(hpf18, file = "hpf18.cluster.rds")
saveRDS(hpf18, file = "hpf20.cluster.rds")

hpf18 <- readRDS("hpf18.cluster.rds")
hpf20 <- readRDS("hpf20.cluster.rds")
