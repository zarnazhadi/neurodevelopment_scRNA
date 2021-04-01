library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)

setwd("/rds/projects/v/vianaj-development-rna/appDir/data")

#------------------------------------———RAW DATA PREP—————————————————————---------------------------- 
# make list of file names
sample.list <- list.files(pattern=glob2rx("*20hpf3_S1*")) # [1] barcodes [2] genes [3] matrix

# read and prepare files
sample.cell.ids <- read_tsv(sample.list[1], col_names = FALSE)$X1
sample.gene.ids <- read_tsv(sample.list[2], col_names = FALSE)$X2
sample.counts <- readMM(sample.list[3]) 

rownames(sample.counts) <- sample.gene.ids
colnames(sample.counts) <- sample.cell.ids

#------------------------------------———QC—————————————————————---------------------------- 
# initialize the Seurat object with the raw (non-normalized data)
# exclusion criteria: genes expressed in less than three cells and cells with less than 500 genes
sample <- CreateSeuratObject(counts = sample.counts, min.cells = 3, min.features  = 500, project = "sample", assay = "RNA")

# selection of mitochondiral genes
sample.mito.genes <- grep(pattern = "mt-|^AC0", x = rownames(sample@assays[["RNA"]]), value = TRUE)

# finding percentage of mitochondrial genes 
sample.percent.mito <- Matrix::colSums(sample@assays[["RNA"]][sample.mito.genes, ])/Matrix::colSums(sample@assays[["RNA"]])

# adding percent.mito to the Seurat object
sample <- AddMetaData(object = sample, metadata = sample.percent.mito, col.name = "percent.mito")

# filter out mitochondrial genes
sample <- subset(x = sample, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & 
                   sample.percent.mito >  -Inf & sample.percent.mito < 0.09 & 
                   nCount_RNA >  -Inf & nCount_RNA < 10000)

# selection of ribosomal genes
sample.ribo.genes <- grep(pattern = "^rps|^rpl", x = rownames(sample@assays[["RNA"]]), value = TRUE)

# filter out ribosomal genes
ribo.index <- which(rownames(sample.counts) %in% sample.ribo.genes)
sample.counts <- sample.counts[-ribo.index, ]
sample <- subset(sample, features = rownames(sample.counts))

# selection of hsp genes
sample.hsp.genes <- grep(pattern = "^hsp", x = rownames(sample@assays[["RNA"]]), value = TRUE)

# filter out hsp genes
hsp.index <- which(rownames(sample.counts) %in% sample.hsp.genes)
sample.counts <- sample.counts[-hsp.index, ]
sample <- subset(sample, features = rownames(sample.counts))

# selection of duplicate genes
sample.dup.genes <- grep(pattern = "1 of many", x = rownames(sample@assays[["RNA"]]), value = TRUE)

# filter out duplicate genes
dup.index <- which(rownames(sample.counts) %in% sample.dup.genes)
sample.counts <- sample.counts[-dup.index, ]
sample <- subset(sample, features = rownames(sample.counts))

# make sure nFeatures has subsetted properly
length(sample$nFeature_RNA)
sample <- subset(x = sample, subset = nFeature_RNA > 500)
length(sample$nFeature_RNA)

# save clean data 
setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data")
saveRDS(sample, file = "S1b_hpf20_clean.rds")