library(tidyverse)
library(dplyr)

# Import marker datasets
#setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data/markers")
setwd("/Volumes/rdsprojects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data/markers")
markers.18hpf <- read.csv(file = "18hpf_all_markers.csv") # n = 10431 
markers.20hpf <- read.csv(file = "20hpf_all_markers.csv") # n = 21915 

# Remove clusters belonging to non-brain tissue
non.brain.18hpf <- c(2, 6, 8, 17, 21, 23, 24, 25)
brain.18hpf <- filter(markers.18hpf, !(cluster %in% non.brain.18hpf)) # n = 6847

non.brain.20hpf <- c(7, 16, 14, 20, 21, 23, 24, 25, 26, 29, 30, 32, 39)
brain.20hpf <- filter(markers.20hpf, !(cluster %in% non.brain.20hpf)) # n = 16113 

# Filter genes to include significant genes (p < 0.01)
genes.18hpf <- filter(brain.18hpf, p_val_adj < 0.01)$gene # n = 5238
genes.20hpf <- filter(brain.20hpf, p_val_adj < 0.01)$gene # n = 12507 

# Merge gene names together and remove duplicates
g <- union(genes.18hpf, genes.20hpf) # n = 5235
genes <- Filter(function(x) length(grep(x, g))==1, g) # n = 4996

# Add all cluster, avg2logFC and p_val_adj information for each gene 
hpf18 <- brain.18hpf %>% 
  select(gene, avg_log2FC, p_val_adj, cluster) %>%
  filter(gene %in% genes) %>% 
  rename(avg_log2FC_18hpf = avg_log2FC, p_val_adj_18hpf = p_val_adj, cluster_18hpf = cluster) # n = 6284 
  
hpf20 <- brain.20hpf %>% 
  select(gene, avg_log2FC, p_val_adj, cluster) %>%
  filter(gene %in% genes) %>% 
  rename(avg_log2FC_20hpf = avg_log2FC, p_val_adj_20hpf = p_val_adj, cluster_20hpf = cluster) # n = 14902  

# Merge into one dataset
merged <- merge(hpf18, hpf20) # n = 26738 
write.csv(merged, file = "merged_markers.csv")

# Identify genes significant in one time point but not the other
insignificant.18hpf <- filter(merged, p_val_adj_18hpf > 0.05)
insignificant.18hpf <- filter(insignificant.18hpf, p_val_adj_20hpf < 0.05) # n = 3213
write.csv(insignificant.18hpf, file = "insignificant_markers_18hpf.csv")
write.csv(unique(insignificant.18hpf$gene), file = "insignificant_markers_18hpf.txt")

insignificant.20hpf <- filter(merged, p_val_adj_20hpf > 0.05)
insignificant.20hpf <- filter(insignificant.20hpf, p_val_adj_18hpf < 0.05) # n = 3771 
write.csv(insignificant.20hpf, file = "insignificant_markers_20hpf.csv")
write.csv(unique(insignificant.20hpf$gene), file = "insignificant_markers_20hpf.txt")

# and same clusters but different genes
# diencephalon
diencephalon.20hpf <- insignificant.18hpf[insignificant.18hpf$cluster_20hpf == insignificant.18hpf$cluster_18hpf, ]
diencephalon.20hpf <- diencephalon.20hpf[diencephalon.20hpf$cluster_20hpf == 13,] # 8
diencephalon.18hpf <- insignificant.20hpf[insignificant.20hpf$cluster_20hpf == insignificant.20hpf$cluster_18hpf, ]
diencephalon.18hpf <- diencephalon.18hpf[diencephalon.18hpf$cluster_18hpf == 13,] # n = 4
diencephalon <- rbind(diencephalon.18hpf, diencephalon.20hpf) # n = 12
write.csv(diencephalon, file = "diencephalon_markers.csv")

# ventral diencephalon
vd.20hpf <- insignificant.18hpf[insignificant.18hpf$cluster_20hpf == 13, ]
vd.20hpf <- vd.20hpf[vd.20hpf$cluster_18hpf == 7,] # n = 14
vd.18hpf <- insignificant.20hpf[insignificant.20hpf$cluster_18hpf == 7, ]
vd.18hpf <- vd.18hpf[vd.18hpf$cluster_20hpf == 13,] # n = 4
vd <- rbind(vd.18hpf, vd.20hpf) # n = 18
write.csv(vd, file = "vd_markers.csv")

# telencephalon (pallium, neuron)
tpn_20hpf <- insignificant.18hpf[insignificant.18hpf$cluster_20hpf == 27,]
tpn_20hpf <- tpn_20hpf[tpn_20hpf$cluster_18hpf == 19,] # n = 8
tpn_18hpf <- insignificant.20hpf[insignificant.20hpf$cluster_18hpf == 19,]
tpn_18hpf <- tpn_18hpf[tpn_18hpf$cluster_20hpf == 27,] # n = 72
tpn <- rbind(tpn_18hpf, tpn_20hpf) # n = 80
write.csv(tpn, file = "tpn_markers.csv")

# neural crest
neural.crest.20hpf <- insignificant.18hpf[insignificant.18hpf$cluster_20hpf == 4, ] 
neural.crest.20hpf <- filter(neural.crest.20hpf, cluster_18hpf %in% c(11,16)) # n = 32
neural.crest.18hpf <- filter(insignificant.20hpf, cluster_18hpf %in% c(11,16)) 
neural.crest.18hpf <- neural.crest.18hpf[neural.crest.18hpf$cluster_20hpf == 4,] # n = 0 
neural.crest <- rbind(neural.crest.18hpf, neural.crest.20hpf) # n = 32
write.csv(neural.crest, file = "neural_crest_markers.csv")

# hindbrain - R5/6
r56.20hpf <- insignificant.18hpf[insignificant.18hpf$cluster_20hpf == 6, ]
r56.20hpf <- filter(r56.20hpf, cluster_18hpf %in% 12) # n = 7
r56.18hpf <- filter(insignificant.20hpf, cluster_18hpf %in% 12) 
r56.18hpf <- r56.18hpf[r56.18hpf$cluster_20hpf == 6, ] # n = 0
r56 <- rbind(r56.18hpf, r56.20hpf) # n = 7
write.csv(r56, file = "hindbrain_r56_markers.csv")

# lens 
lens.20hpf <- insignificant.18hpf[insignificant.18hpf$cluster_20hpf == 22, ]
lens.20hpf <- filter(lens.20hpf, cluster_18hpf %in% 20) # 27
lens.18hpf <- filter(insignificant.20hpf, cluster_18hpf %in% 20) 
lens.18hpf <- lens.18hpf[lens.18hpf$cluster_20hpf == 22, ] # 9
lens <- rbind(lens.18hpf, lens.20hpf) # 36
write.csv(lens, file = "lens_markers.csv")

# neurons
neurons.20hpf <- filter(insignificant.18hpf, cluster_20hpf %in% c(2,3))
neurons.20hpf <- filter(neurons.20hpf, cluster_18hpf %in% 15) # 47
neurons.18hpf <- filter(insignificant.20hpf, cluster_18hpf %in% 15) 
neurons.18hpf <- filter(neurons.18hpf, cluster_20hpf %in% c(2,3)) # 4
neurons <- rbind(neurons.18hpf, neurons.20hpf) # 51
write.csv(neurons, file = "neurons_markers.csv")

# midbrain
midbrain.20hpf <- filter(insignificant.18hpf, cluster_20hpf %in% c(11,12))
midbrain.20hpf <- filter(midbrain.20hpf, cluster_18hpf %in% 10) # 4
midbrain.18hpf <- filter(insignificant.20hpf, cluster_18hpf %in% 10) 
midbrain.18hpf <- filter(midbrain.18hpf, cluster_20hpf %in% c(11,12)) # 32
midbrain <- rbind(midbrain.18hpf, midbrain.20hpf) # 36
write.csv(midbrain, file = "midbrain_markers.csv")
