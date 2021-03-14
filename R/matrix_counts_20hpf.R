# script which identifies all 20hpf raw samples, cleans, normalises 
# and merges them to prepare for analysis

library("Matrix")
library(tidyverse)
library("Seurat")
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/appDir/data")

#------------------------------------———RAW DATA PREP—————————————————————---------------------------- 
# make list of all 20hpf files
S1a.list <- list.files(pattern=glob2rx("*20hpf1_S1*")) # [1] barcodes [2] genes [3] matrix
S1b.list <- list.files(pattern=glob2rx("*20hpf3_S1*")) # [1] barcodes [2] genes [3] matrix
S2.list <- list.files(pattern=glob2rx("*20hpf2_S2*")) # [1] barcodes [2] genes [3] matrix

# read and prepare all S1a files
S1a.cell.ids <- read_tsv(S1a.list[1], col_names = FALSE)$X1
S1a.gene.ids <- read_tsv(S1a.list[2], col_names = FALSE)$X2 # was x1 before and was getting the ensembl ids instead of gene names
S1a.counts <- readMM(S1a.list[3]) 

rownames(S1a.counts) <- S1a.gene.ids
colnames(S1a.counts) <- S1a.cell.ids

# read and prepare S1b files
S1b.cell.ids <- read_tsv(S1b.list[1], col_names = FALSE)$X1
S1b.gene.ids <- read_tsv(S1b.list[2], col_names = FALSE)$X2
S1b.counts <- readMM(S1b.list[3]) 

rownames(S1b.counts) <- S1b.gene.ids
colnames(S1b.counts) <- S1b.cell.ids

# read and prepare S2 files
S2.cell.ids <- read_tsv(S2.list[1], col_names = FALSE)$X1
S2.gene.ids <- read_tsv(S2.list[2], col_names = FALSE)$X2
S2.counts <- readMM(S2.list[3]) 

rownames(S2.counts) <- S2.gene.ids
colnames(S2.counts) <- S2.cell.ids

#------------------------------------———QC—————————————————————---------------------------- 
# examine the memory savings between regular and sparse matrices
S1a.dense.size <- object.size(x = as.matrix(x = S1a.counts)) # 710762024 bytes
S1b.dense.size <- object.size(x = as.matrix(x = S1b.counts)) # 2045470520 bytes
S2.dense.size <- object.size(x = as.matrix(x = S2.counts)) # 669801080 bytes

S1a.sparse.size <- object.size(x = S1a.counts) # 66976824 bytes
S1b.sparse.size <- object.size(x = S1b.counts) # 196081160 bytes
S2.sparse.size <- object.size(x = S2.counts) # 64607120 bytes

S1a.dense.size/S1a.sparse.size # 10.6 bytes
S1b.dense.size/S1b.sparse.size # 10.4 bytes
S2.dense.size/S2.sparse.size # 10.4 bytes

# initialize the Seurat object with the raw (non-normalized data)
# exclusion criteria: genes expressed in less than three cells and cells with less than 500 genes
S1a <- CreateSeuratObject(counts = S1a.counts, min.cells = 3, min.features  = 500, project = "S1a", assay = "RNA")
S1b <- CreateSeuratObject(counts = S1b.counts, min.cells = 3, min.features  = 500, project = "S1b", assay = "RNA")
S2 <- CreateSeuratObject(counts = S2.counts, min.cells = 3, min.features  = 500, project = "S2", assay = "RNA")

# selection of mitochondiral genes
S1a.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(S1a@assays[["RNA"]]), value = TRUE)
S1b.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(S1b@assays[["RNA"]]), value = TRUE)
S2.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(S2@assays[["RNA"]]), value = TRUE)

# finding percentage of mitochondrial genes 
S1a.percent.mito <- Matrix::colSums(S1a@assays[["RNA"]][S1a.mito.genes, ])/Matrix::colSums(S1a@assays[["RNA"]])
S1b.percent.mito <- Matrix::colSums(S1b@assays[["RNA"]][S1b.mito.genes, ])/Matrix::colSums(S1b@assays[["RNA"]])
S2.percent.mito <- Matrix::colSums(S2@assays[["RNA"]][S2.mito.genes, ])/Matrix::colSums(S2@assays[["RNA"]])

# adding percent.mito to the Seurat object
S1a <- AddMetaData(object = S1a, metadata = S1a.percent.mito, col.name = "percent.mito")
S1b <- AddMetaData(object = S1b, metadata = S1b.percent.mito, col.name = "percent.mito")
S2 <- AddMetaData(object = S2, metadata = S2.percent.mito, col.name = "percent.mito")

# selection of ribosomal genes
S1a.ribo.genes <- grep(pattern = "^rps|^rpl", x = rownames(S1a@assays[["RNA"]]), value = TRUE)
S1b.ribo.genes <- grep(pattern = "^rps|^rpl", x = rownames(S1b@assays[["RNA"]]), value = TRUE)
S2.ribo.genes <- grep(pattern = "^rps|^rpl", x = rownames(S2@assays[["RNA"]]), value = TRUE)

# finding percentage of mitochondrial genes 
S1a.percent.ribo <- Matrix::colSums(S1a@assays[["RNA"]][S1a.ribo.genes, ])/Matrix::colSums(S1a@assays[["RNA"]])
S1b.percent.ribo <- Matrix::colSums(S1b@assays[["RNA"]][S1b.ribo.genes, ])/Matrix::colSums(S1b@assays[["RNA"]])
S2.percent.ribo <- Matrix::colSums(S2@assays[["RNA"]][S2.ribo.genes, ])/Matrix::colSums(S2@assays[["RNA"]])

# adding percent.mito to the Seurat object
S1a <- AddMetaData(object = S1a, metadata = S1a.percent.ribo, col.name = "percent.ribo")
S1b <- AddMetaData(object = S1b, metadata = S1b.percent.ribo, col.name = "percent.ribo")
S2 <- AddMetaData(object = S2, metadata = S2.percent.ribo, col.name = "percent.ribo")

# selection of hsp genes
S1a.hsp.genes <- grep(pattern = "^hsp", x = rownames(S1a@assays[["RNA"]]), value = TRUE)
S1b.hsp.genes <- grep(pattern = "^hsp", x = rownames(S1b@assays[["RNA"]]), value = TRUE)
S2.hsp.genes <- grep(pattern = "^^hsp", x = rownames(S2@assays[["RNA"]]), value = TRUE)

# finding percentage of mitochondrial genes 
S1a.percent.hsp <- Matrix::colSums(S1a@assays[["RNA"]][S1a.hsp.genes, ])/Matrix::colSums(S1a@assays[["RNA"]])
S1b.percent.hsp <- Matrix::colSums(S1b@assays[["RNA"]][S1b.hsp.genes, ])/Matrix::colSums(S1b@assays[["RNA"]])
S2.percent.hsp <- Matrix::colSums(S2@assays[["RNA"]][S2.hsp.genes, ])/Matrix::colSums(S2@assays[["RNA"]])

# adding percent.mito to the Seurat object
S1a <- AddMetaData(object = S1a, metadata = S1a.percent.hsp, col.name = "percent.hsp")
S1b <- AddMetaData(object = S1b, metadata = S1b.percent.hsp, col.name = "percent.hsp")
S2 <- AddMetaData(object = S2, metadata = S2.percent.hsp, col.name = "percent.hsp")

#------------------------------------———PLOTS—————————————————————---------------------------- 
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1a_mito_vln_plot.pdf")
VlnPlot(object = S1a, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1a_ribo_vln_plot.pdf")
VlnPlot(object = S1a, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1a_hsp_vln_plot.pdf")
VlnPlot(object = S1a, features = c("nFeature_RNA", "nCount_RNA", "percent.hsp"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1b_mito_vln_plot.pdf")
VlnPlot(object = S1b, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1b_ribo_vln_plot.pdf")
VlnPlot(object = S1b, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1b_hsp_vln_plot.pdf")
VlnPlot(object = S1b, features = c("nFeature_RNA", "nCount_RNA", "percent.hsp"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S2_mito_vln_plot.pdf")
VlnPlot(object = S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S2_ribo_vln_plot.pdf")
VlnPlot(object = S2, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S2_hsp_vln_plot.pdf")
VlnPlot(object = S2, features = c("nFeature_RNA", "nCount_RNA", "percent.hsp"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1a_mito_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1a, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1a_ribo_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1a, feature1 = "nCount_RNA", feature2 = "percent.ribo")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1a_hsp_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1a, feature1 = "nCount_RNA", feature2 = "percent.hsp")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1a_RNA_scatter_plot.pdf")
FeatureScatter(object = S1a, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1b_mito_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1b, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1b_ribo_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1b, feature1 = "nCount_RNA", feature2 = "percent.ribo")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1b_hsp_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1b, feature1 = "nCount_RNA", feature2 = "percent.hsp")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1b_RNA_scatter_plot.pdf")
FeatureScatter(object = S1b, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S1a_mito_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1a, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S2_ribo_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S2, feature1 = "nCount_RNA", feature2 = "percent.ribo")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S2_hsp_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S2, feature1 = "nCount_RNA", feature2 = "percent.hsp")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/20hpf/pre-normalisation/S2_RNA_scatter_plot.pdf")
FeatureScatter(object = S2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# filter
S1a <- subset(x = S1a, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & S1a.percent.mito >  -Inf & S1a.percent.mito < 0.05 )
S1b <- subset(x = S1b, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & S1b.percent.mito >  -Inf & S1b.percent.mito < 0.05 )
S2 <- subset(x = S2, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & S2.percent.mito >  -Inf & S2.percent.mito < 0.05 )

# normalise
# S1a <- NormalizeData(object = S1a, normalization.method = "LogNormalize", scale.factor = 10000)
# S1b <- NormalizeData(object = S1b, normalization.method = "LogNormalize", scale.factor = 10000)
# S2 <- NormalizeData(object = S2, normalization.method = "LogNormalize", scale.factor = 10000)

# merge data
hpf20.combined <- merge(S1a, y = c(S1b, S2), add.cell.ids = c("S1a", "S1b", "S2"), project = "hpf20")

# cell names have an added identifier 
head(colnames(hpf20.combined))

# visualise
table(hpf20.combined$orig.ident)

# normalise
hpf20.combined.normalised <- NormalizeData(object = hpf20.combined, normalization.method = "LogNormalize", scale.factor = 10000)
GetAssayData(hpf20.combined.normalised[1:10, 1:15])

# output file as .rds
saveRDS(hpf20.combined.normalised, file = "hpf20_combined_normalised.rds")