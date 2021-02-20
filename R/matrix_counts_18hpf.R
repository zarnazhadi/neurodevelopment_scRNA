library("Matrix")
library(tidyverse)
library("Seurat")

setwd("/rds/projects/v/vianaj-development-rna/appDir/data")

S3.list <- list.files(pattern=glob2rx("*18hpf1_18S_S3*")) # [1] barcodes [2] genes [3] matrix
S4.list <- list.files(pattern=glob2rx("*18hpf2_18S_S4*")) # [1] barcodes [2] genes [3] matrix

# read S3 files
S3.cell.ids <- read_tsv(S3.list[1], col_names = FALSE)$X1
S3.gene.ids <- read_tsv(S3.list[2], col_names = FALSE)$X1
S3.counts <- readMM(S3.list[3]) 

rownames(S3.counts) <- S3.gene.ids
colnames(S3.counts) <- S3.cell.ids
S3.counts

# read S4 files
S4.cell.ids <- read_tsv(S4.list[1], col_names = FALSE)$X1
S4.gene.ids <- read_tsv(S4.list[2], col_names = FALSE)$X1
S4.counts <- readMM(S4.list[3]) 

rownames(S4.counts) <- S4.gene.ids
colnames(S4.counts) <- S4.cell.ids
S4.counts

# combine counts by creating two seurat objects
S3 <- CreateSeuratObject(counts = S3.counts, project = "S3")
S4 <- CreateSeuratObject(counts = S4.counts, project = "S4")

# merge counts
hpf18.combined <- merge(S3, y = S4, add.cell.ids = c("S3", "S4"), project = "hpf18")
hpf18.combined

# cell names have an added identifier 
head(colnames(hpf18.combined))

# visualise
table(hpf18.combined$orig.ident)
