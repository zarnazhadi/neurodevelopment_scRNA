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
#S3.counts

# read S4 files
S4.cell.ids <- read_tsv(S4.list[1], col_names = FALSE)$X1
S4.gene.ids <- read_tsv(S4.list[2], col_names = FALSE)$X1
S4.counts <- readMM(S4.list[3]) 

rownames(S4.counts) <- S4.gene.ids
colnames(S4.counts) <- S4.cell.ids
#S4.counts

# combine counts by creating two seurat objects
S3 <- CreateSeuratObject(counts = S3.counts, project = "S3")
S4 <- CreateSeuratObject(counts = S4.counts, project = "S4")

# merge counts
hpf18.combined <- merge(S3, y = S4, add.cell.ids = c("S3", "S4"), project = "hpf18")
#hpf18.combined

# cell names have an added identifier 
head(colnames(hpf18.combined))

# visualise
table(hpf18.combined$orig.ident)

#------------------------------------———NORMALISATION—————————————————————---------------------------- 

# standard pre processing workflow
# examine the memory savings between regular and sparse matrices
S3.dense.size <- object.size(x = as.matrix(x = S3.counts)) # 1338056984 bytes
S4.dense.size <- object.size(x = as.matrix(x = S4.counts)) # 1476654392 bytes

S3.sparse.size <- object.size(x = S3.counts) # 147442200 bytes
S4.sparse.size <- object.size(x = S4.counts) # 153873288 bytes

S3.dense.size/S3.sparse.size # 9.1 bytes
S4.dense.size/S4.sparse.size # 9.6 bytes

# initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes

S3 <- CreateSeuratObject(counts = S3.counts, min.cells = 3, min.features  = 200, project = "S3", assay = "RNA")
S4 <- CreateSeuratObject(counts = S4.counts, min.cells = 3, min.features  = 200, project = "S4", assay = "RNA")

# QC and selecting cells for further analysis
# % mitochondrial genes

S3.mito.genes <- grep(pattern = "^MT-", x = rownames(S3@assays[["RNA"]]), value = TRUE)
S4.mito.genes <- grep(pattern = "^MT-", x = rownames(S4@assays[["RNA"]]), value = TRUE)

S3.percent.mito <- Matrix::colSums(S3@assays[["RNA"]][S3.mito.genes, ])/Matrix::colSums(S3@assays[["RNA"]])
S4.percent.mito <- Matrix::colSums(S4@assays[["RNA"]][S4.mito.genes, ])/Matrix::colSums(S4@assays[["RNA"]])

S3 <- AddMetaData(object = S3, metadata = S3.percent.mito, col.name = "S3.percent.mito")
S4 <- AddMetaData(object = S4, metadata = S4.percent.mito, col.name = "S4.percent.mito")

# plots
pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/18hpf_S3_vln_plot.pdf")
VlnPlot(object = S3, features = c("nFeature_RNA", "nCount_RNA", "S3.percent.mito"), ncol = 3)
dev.off()

pdf("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/18hpf_S4_vln_plot.pdf")
VlnPlot(object = S4, features = c("nFeature_RNA", "nCount_RNA", "S4.percent.mito"), ncol = 3)
dev.off()
