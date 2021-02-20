library("Matrix")
library(tidyverse)

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
S1a.counts

# read S1b files
S1b.cell.ids <- read_tsv(S1b.list[1], col_names = FALSE)$X1
S1b.gene.ids <- read_tsv(S1b.list[2], col_names = FALSE)$X1
S1b.counts <- readMM(S1b.list[3]) 

rownames(S1b.counts) <- S1b.gene.ids
colnames(S1b.counts) <- S1b.cell.ids
S1b.counts

# read S2 files
S2.cell.ids <- read_tsv(S2.list[1], col_names = FALSE)$X1
S2.gene.ids <- read_tsv(S2.list[2], col_names = FALSE)$X1
S2.counts <- readMM(S2.list[3]) 

rownames(S2.counts) <- S2.gene.ids
colnames(S2.counts) <- S2.cell.ids
S2.counts
