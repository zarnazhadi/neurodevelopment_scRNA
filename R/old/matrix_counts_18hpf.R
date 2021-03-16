# script which identifies all 18hpf raw samples, cleans, normalises 
# and merges them to prepare for analysis

library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/appDir/data")

#------------------------------------———RAW DATA PREP—————————————————————---------------------------- 
# make list of all 18hpf files
S3.list <- list.files(pattern=glob2rx("*18hpf1_18S_S3*")) # [1] barcodes [2] genes [3] matrix
S4.list <- list.files(pattern=glob2rx("*18hpf2_18S_S4*")) # [1] barcodes [2] genes [3] matrix

# read and prepare S3 files
S3.cell.ids <- read_tsv(S3.list[1], col_names = FALSE)$X1
S3.gene.ids <- read_tsv(S3.list[2], col_names = FALSE)$X2
S3.counts <- readMM(S3.list[3]) 

rownames(S3.counts) <- S3.gene.ids
colnames(S3.counts) <- S3.cell.ids

# read and prepare S4 files
S4.cell.ids <- read_tsv(S4.list[1], col_names = FALSE)$X1
S4.gene.ids <- read_tsv(S4.list[2], col_names = FALSE)$X2
S4.counts <- readMM(S4.list[3]) 

rownames(S4.counts) <- S4.gene.ids
colnames(S4.counts) <- S4.cell.ids

#------------------------------------———QC—————————————————————---------------------------- 
# examine the memory savings between regular and sparse matrices
S3.dense.size <- object.size(x = as.matrix(x = S3.counts)) # 1338056984 bytes
S4.dense.size <- object.size(x = as.matrix(x = S4.counts)) # 1476654392 bytes

S3.sparse.size <- object.size(x = S3.counts) # 147442200 bytes
S4.sparse.size <- object.size(x = S4.counts) # 153873288 bytes

S3.dense.size/S3.sparse.size # 9.1 bytes
S4.dense.size/S4.sparse.size # 9.6 bytes

# initialize the Seurat object with the raw (non-normalized data)
# exclusion criteria: genes expressed in less than three cells and cells with less than 500 genes
S3 <- CreateSeuratObject(counts = S3.counts, min.cells = 3, min.features  = 500, project = "S3", assay = "RNA")
S4 <- CreateSeuratObject(counts = S4.counts, min.cells = 3, min.features  = 500, project = "S4", assay = "RNA")

# selection of mitochondiral genes
S3.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(S3@assays[["RNA"]]), value = TRUE)
S4.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(S4@assays[["RNA"]]), value = TRUE)

# finding percentage of mitochondrial genes 
S3.percent.mito <- Matrix::colSums(S3@assays[["RNA"]][S3.mito.genes, ])/Matrix::colSums(S3@assays[["RNA"]])
S4.percent.mito <- Matrix::colSums(S4@assays[["RNA"]][S4.mito.genes, ])/Matrix::colSums(S4@assays[["RNA"]])

# adding percent.mito to the Seurat object
S3 <- AddMetaData(object = S3, metadata = S3.percent.mito, col.name = "percent.mito")
S4 <- AddMetaData(object = S4, metadata = S4.percent.mito, col.name = "percent.mito")

# selection of ribosomal genes
S3.ribo.genes <- grep(pattern = "^rps|^rpl", x = rownames(S3@assays[["RNA"]]), value = TRUE)
S4.ribo.genes <- grep(pattern = "^rps|^rpl", x = rownames(S4@assays[["RNA"]]), value = TRUE)

# finding percentage of ribsomal genes 
S3.percent.ribo <- Matrix::colSums(S3@assays[["RNA"]][S3.ribo.genes, ])/Matrix::colSums(S3@assays[["RNA"]])
S4.percent.ribo <- Matrix::colSums(S4@assays[["RNA"]][S4.ribo.genes, ])/Matrix::colSums(S4@assays[["RNA"]])

# adding percent.ribo to the Seurat object
S3 <- AddMetaData(object = S3, metadata = S3.percent.ribo, col.name = "percent.ribo")
S4 <- AddMetaData(object = S4, metadata = S4.percent.ribo, col.name = "percent.ribo")

# selection of hsp genes 
S3.hsp.genes <- grep(pattern = "^hsp", x = rownames(S3@assays[["RNA"]]), value = TRUE)
S4.hsp.genes <- grep(pattern = "^hsp", x = rownames(S4@assays[["RNA"]]), value = TRUE)

# finding percentage of hsp genes 
S3.percent.hsp <- Matrix::colSums(S3@assays[["RNA"]][S3.hsp.genes, ])/Matrix::colSums(S3@assays[["RNA"]])
S4.percent.hsp <- Matrix::colSums(S4@assays[["RNA"]][S4.hsp.genes, ])/Matrix::colSums(S4@assays[["RNA"]])

# adding percent.hsp to the Seurat object
S3 <- AddMetaData(object = S3, metadata = S3.percent.hsp, col.name = "percent.hsp")
S4 <- AddMetaData(object = S4, metadata = S4.percent.hsp, col.name = "percent.hsp")

#------------------------------------———PLOTS—————————————————————---------------------------- 
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S3_mito_vln_plot.pdf")
VlnPlot(object = S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S3_ribo_vln_plot.pdf")
VlnPlot(object = S3, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S3_hsp_vln_plot.pdf")
VlnPlot(object = S3, features = c("nFeature_RNA", "nCount_RNA", "percent.hsp"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S3_mito_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S3, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S3_ribo_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S3, feature1 = "nCount_RNA", feature2 = "percent.ribo")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S3_hsp_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S3, feature1 = "nCount_RNA", feature2 = "percent.hsp")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S3_RNA_scatter_plot.pdf")
FeatureScatter(object = S3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S4_mito_vln_plot.pdf")
VlnPlot(object = S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S4_ribo_vln_plot.pdf")
VlnPlot(object = S4, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S4_mito_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S4, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S4_ribo_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S4, feature1 = "nCount_RNA", feature2 = "percent.ribo")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S4_hsp_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S4, feature1 = "nCount_RNA", feature2 = "percent.hsp")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/18hpf/pre-normalisation/S4_RNA_scatter_plot.pdf")
FeatureScatter(object = S4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

#------------------------------------———NORMALISATION—————————————————————---------------------------- 
# filter low quality cells and duplets
# also remove duplicate genes which only occur in fish
S3 <- subset(x = S3, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & 
               S3.percent.mito >  -Inf & S3.percent.mito < 0.09)

S4 <- subset(x = S4, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & S4.percent.mito >  -Inf & S4.percent.mito < 0.09 )

# merge
hpf18.combined <- merge(S3, y = S4, add.cell.ids = c("S3", "S4"), project = "hpf18")

# cell names have an added identifier 
head(colnames(hpf18.combined))

# visualise
table(hpf18.combined$orig.ident)

# normalise merged data
hpf18.combined.normalised <- NormalizeData(object = hpf18.combined, normalization.method = "LogNormalize", scale.factor = 10000)
GetAssayData(hpf18.combined.normalised[1:10, 1:15])

# output file as .rds
saveRDS(hpf18.combined.normalised, file = "hpf18_combined_normalised.rds")
