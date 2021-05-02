library("Matrix")
library("Seurat")
library(tidyverse)
library(dplyr)
library(cowplot)
library(data.table)

# import Seurat object and markers
setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data/markers")
sample <- readRDS("hpf18_tSNE.cluster.rds")
sample.name <- "18hpf"
top2 <- read.csv(file = paste0(sample.name, "_top2_markers.csv"), header = TRUE)
top10 <- read.csv(file = paste0(sample.name, "_top10_markers.csv"), header = TRUE)

# import gene annotation resource
setwd("/rds/projects/v/vianaj-development-rna/appDir/data")
atlas <- as.data.frame(fread("ZFBrainAtlasMaster.tsv", header = TRUE))
atlas <- atlas %>% filter(STAGE == "18 hpf")

# assign cell identities to markers
identities.list <- list()
for(row in seq(from= 1, to = dim(top2)[1])){
  identity <- toString(atlas[atlas$`ENRICHED MARKERS` %like% top2[row,8], ]$`ASSIGNED CELL TYPE/STATE`)
  if(identity == ""){
    identity <- toString(top2[row, 7])
    print(identity)
  }
  top2$cell.identity[[row]] <- identity
  identities.list[[row]] <- identity
}

top2 <- apply(top2,2,as.character)
write.csv(as.data.frame(top2), file = paste0(sample.name, "_top2_marker_names.csv"))

identities.list <- list()
for(row in seq(from= 1, to = dim(top10)[1])){
  identity <- toString(atlas[atlas$`ENRICHED MARKERS` %like% top10[row,8], ]$`ASSIGNED CELL TYPE/STATE`)
  if(identity == ""){
    identity <- toString(top10[row, 7])
    print(identity)
  }
  top10$cell.identity[[row]] <- identity
  identities.list[[row]] <- identity
}

top10 <- top10 %>%
  filter(!(str_detect(cell.identity, "\\d")))
top10 <- apply(top10,2,as.character)
write.csv(top10, file = paste0(sample.name, "_top10_marker_names.csv"))

# manual insertion of cluster identities 
current.cluster.ids <- seq(from=0, to=25)
new.cluster.ids <- c("Hindbrain (Floor plate/Boundary- R3),\nVentral Diencephalon, Dienc-Midbrain Boundary,\nDiencephalon", "Hindbrain (Ant. Dorsal, Rhombic Lip)", "2", 
                     "Telencephalon (I, II)", "Optic Vesicle (Dorsal)", "Hindbrain- R7, Pharyngeal Arch (Posterior),\nVentral Midbrain",
                     "Floor Plate", "Ventral Diencephalon", "Epidermis", 
                     "Placode (Epibranchial/Lens),\nEpidermis", "Midbrain, Midbrain Neural Rod", "Neural Crest", 
                     "Hindbrain- R5/6", "Diencephalon", "Placode (Epibranchial)",
                     "Commited Progenitors,\nGanglion/Glutamatergic Neurons", "Neural Crest*", "Heart Primordium", 
                     "Placode", "Telencephalon (Pallium, Neuron)", "Lens (Differentiating)", 
                     "Rostral Blood Island (Myeloid)", "Placode (Otic)", "Periderm", 
                     "Prechordal Plate/Polster", "Blood Precursors")
sample@active.ident <- plyr::mapvalues(x = sample@active.ident, from = current.cluster.ids, to = new.cluster.ids)

# output t-SNE plot
pdf(paste0("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/plots/", sample.name, "/post-clustering/NEW_tSNE_plot.pdf"))
pp <- DimPlot(object = sample, reduction = "tsne", label = FALSE, pt.size = 0.5, label.size = 2) + theme_classic(base_size = 10) + 
  guides(color = guide_legend(override.aes = list(size=1), ncol=1))
print(pp)
dev.off()

