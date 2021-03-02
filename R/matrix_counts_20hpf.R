library("Matrix")
library(tidyverse)
library("Seurat")
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/appDir/data")

S1a.list <- list.files(pattern=glob2rx("*20hpf1_S1*")) # [1] barcodes [2] genes [3] matrix
S1b.list <- list.files(pattern=glob2rx("*20hpf3_S1*")) # [1] barcodes [2] genes [3] matrix
S2.list <- list.files(pattern=glob2rx("*20hpf2_S2*")) # [1] barcodes [2] genes [3] matrix

# read S1a files
S1a.cell.ids <- read_tsv(S1a.list[1], col_names = FALSE)$X1
S1a.gene.ids <- read_tsv(S1a.list[2], col_names = FALSE)$X1
S1a.counts <- readMM(S1a.list[3]) 

rownames(S1a.counts) <- S1a.gene.ids
colnames(S1a.counts) <- S1a.cell.ids
#S1a.counts

# read S1b files
S1b.cell.ids <- read_tsv(S1b.list[1], col_names = FALSE)$X1
S1b.gene.ids <- read_tsv(S1b.list[2], col_names = FALSE)$X1
S1b.counts <- readMM(S1b.list[3]) 

rownames(S1b.counts) <- S1b.gene.ids
colnames(S1b.counts) <- S1b.cell.ids
#S1b.counts

# read S2 files
S2.cell.ids <- read_tsv(S2.list[1], col_names = FALSE)$X1
S2.gene.ids <- read_tsv(S2.list[2], col_names = FALSE)$X1
S2.counts <- readMM(S2.list[3]) 

rownames(S2.counts) <- S2.gene.ids
colnames(S2.counts) <- S2.cell.ids
#S2.counts

# combine counts by creating seurat objects
S1a <- CreateSeuratObject(counts = S1a.counts, project = "S1a")
S1b <- CreateSeuratObject(counts = S1b.counts, project = "S1b")
S2 <- CreateSeuratObject(counts = S2.counts, project = "S2")

# merge multiple counts
hpf20.combined <- merge(S1a, y = c(S1b, S2), add.cell.ids = c("S1a", "S1b", "S2"), project = "hpf20")
hpf20.combined

# cell names have an added identifier 
head(colnames(hpf20.combined))

# visualise
table(hpf20.combined$orig.ident)

#------------------------------------———NORMALISATION—————————————————————---------------------------- 

# standard pre processing workflow
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

# initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 500 detected genes

S1a <- CreateSeuratObject(counts = S1a.counts, min.cells = 3, min.features  = 500, project = "S1a", assay = "RNA")
S1b <- CreateSeuratObject(counts = S1b.counts, min.cells = 3, min.features  = 500, project = "S1b", assay = "RNA")
S2 <- CreateSeuratObject(counts = S2.counts, min.cells = 3, min.features  = 500, project = "S2", assay = "RNA")

# QC and selecting cells for further analysis
# % mitochondrial genes

S1a.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(S1a@assays[["RNA"]]), value = TRUE)
S1b.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(S1b@assays[["RNA"]]), value = TRUE)
S2.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(S2@assays[["RNA"]]), value = TRUE)

S1a.percent.mito <- Matrix::colSums(S1a@assays[["RNA"]][S1a.mito.genes, ])/Matrix::colSums(S1a@assays[["RNA"]])
S1b.percent.mito <- Matrix::colSums(S1b@assays[["RNA"]][S1b.mito.genes, ])/Matrix::colSums(S1b@assays[["RNA"]])
S2.percent.mito <- Matrix::colSums(S2@assays[["RNA"]][S2.mito.genes, ])/Matrix::colSums(S2@assays[["RNA"]])

S1a <- AddMetaData(object = S1a, metadata = S1a.percent.mito, col.name = "S1a.percent.mito")
S1b <- AddMetaData(object = S1b, metadata = S1b.percent.mito, col.name = "S1b.percent.mito")
S2 <- AddMetaData(object = S2, metadata = S2.percent.mito, col.name = "S2.percent.mito")

# vln plots
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/UPDATE_20hpf_S1a_vln_plot.pdf")
VlnPlot(object = S1a, features = c("nFeature_RNA", "nCount_RNA", "S1a.percent.mito"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/UPDATE_20hpf_S1b_vln_plot.pdf")
VlnPlot(object = S1b, features = c("nFeature_RNA", "nCount_RNA", "S1b.percent.mito"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/UPDATE_20hpf_S2_vln_plot.pdf")
VlnPlot(object = S2, features = c("nFeature_RNA", "nCount_RNA", "S2.percent.mito"), ncol = 3)
dev.off()

# scatter plots
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/20hpf_S1a_mito_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1a, feature1 = "nCount_RNA", feature2 = "S1a.percent.mito")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/20hpf_S1a_RNA_scatter_plot.pdf")
FeatureScatter(object = S1a, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/18hpf_S1b_mito_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S1b, feature1 = "nCount_RNA", feature2 = "S1b.percent.mito")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/18hpf_S1b_RNA_scatter_plot.pdf")
FeatureScatter(object = S1b, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/18hpf_S2_mito_scatter_plot.pdf")
par(mfrow = c(1, 2))
FeatureScatter(object = S2, feature1 = "nCount_RNA", feature2 = "S2.percent.mito")
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/18hpf_S2_RNA_scatter_plot.pdf")
FeatureScatter(object = S2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# filter
S1a <- subset(x = S1a, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & S1a.percent.mito >  -Inf & S1a.percent.mito < 0.05 )
S1b <- subset(x = S1b, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & S1b.percent.mito >  -Inf & S1b.percent.mito < 0.05 )
S2 <- subset(x = S2, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & S2.percent.mito >  -Inf & S2.percent.mito < 0.05 )

# normalise
S1a <- NormalizeData(object = S1a, normalization.method = "LogNormalize", scale.factor = 10000)
S1b <- NormalizeData(object = S1b, normalization.method = "LogNormalize", scale.factor = 10000)
S2 <- NormalizeData(object = S2, normalization.method = "LogNormalize", scale.factor = 10000)

# merge based on normalised data
hpf20.normalised <- merge(S1a, y = c(S1b, S2), add.cell.ids = c("S1a", "S1b", "S2"), project = "hpf20", merge.data = TRUE)
GetAssayData(hpf20.combined[1:10, 1:15])
GetAssayData(hpf20.normalised[1:10, 1:15])

# normalise merged data
hpf20.combined.normalised <- NormalizeData(object = hpf20.combined, normalization.method = "LogNormalize", scale.factor = 10000)
GetAssayData(hpf20.combined.normalised[1:10, 1:15])

# output as matrix - need to finish this
#as.matrix(object@hpf20.combined.normalised)
#write.csv2(as.matrix(hpf20.combined.normalised), file= "hpf20_combined_normalised.csv")
