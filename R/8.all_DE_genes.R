library(tidyverse)
library(dplyr)

# Import marker datasets
setwd("/rds/projects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data/markers")
setwd("/Volumes/rdsprojects/v/vianaj-development-rna/Zarnaz/neurodevelopment_scRNA/data/markers")
markers.18hpf <- read.csv(file = "18hpf_all_markers.csv", stringsAsFactors=FALSE) # n = 10431
markers.20hpf <- read.csv(file = "20hpf_all_markers.csv", stringsAsFactors=FALSE) # n = 21915

# midbrain
#cluster 10 in 18hpf and 11&12 in 20hpf

midbrain18 <-markers.18hpf[markers.18hpf$cluster == 10 & markers.18hpf$p_val_adj, ] # n = 486 

midbrain20 <-markers.20hpf[markers.20hpf$cluster == 11 | markers.20hpf$cluster == 12 , ] # n = 1226 

cols_final <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") # names of the columns we want to keep


midbrain1 <- merge(midbrain18[, cols_final], midbrain20[, cols_final], by="gene")

midbrain2 <- midbrain1[midbrain1$p_val_adj.x < 0.01 | midbrain1$p_val_adj.y < 0.01 , ] # n = 347 

dim(midbrain2[midbrain2$p_val_adj.x < 0.01 | midbrain2$p_val_adj.y < 0.01 , ])
dim(midbrain2[midbrain2$p_val_adj.x > 0.01 & midbrain2$p_val_adj.y < 0.01 , ])
dim(midbrain2[midbrain2$p_val_adj.x < 0.01 & midbrain2$p_val_adj.y > 0.01 , ])
dim(midbrain2[midbrain2$p_val_adj.x < 0.01 & midbrain2$p_val_adj.y < 0.01 , ])

colnames(midbrain2) <- c("Gene", "p-value 18hpf", "Avg Log2FC 18hpf", "pct 1 18hpf", "pct 2 18hpf", "Bonferroni p-value 18hpf", "cluster 18hpf", "p-value 20hpf", "Avg Log2FC 20hpf", "pct 1 20hpf", "pct 2 20hpf", "Bonferroni p-value 20hpf", "cluster 20hpf")
write.csv(midbrain2, file = "midbrain_markers.csv")


# telencephalon
#cluster 10 in 18hpf and 11&12 in 20hpf
telencephalon18 <-markers.18hpf[markers.18hpf$cluster == 3 & markers.18hpf$p_val_adj, ] # n = 191 

telencephalon20 <-markers.20hpf[markers.20hpf$cluster == 5, ] # n = 468 

cols_final <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") # names of the columns we want to keep


telencephalon1 <- merge(telencephalon18[, cols_final], telencephalon20[, cols_final], by="gene")

telencephalon2 <- telencephalon1[telencephalon1$p_val_adj.x < 0.01 | telencephalon1$p_val_adj.y < 0.01 , ] # n = 131

dim(telencephalon2[telencephalon2$p_val_adj.x < 0.01 | telencephalon2$p_val_adj.y < 0.01 , ])
dim(telencephalon2[telencephalon2$p_val_adj.x > 0.01 & telencephalon2$p_val_adj.y < 0.01 , ])
dim(telencephalon2[telencephalon2$p_val_adj.x < 0.01 & telencephalon2$p_val_adj.y > 0.01 , ])
dim(telencephalon2[telencephalon2$p_val_adj.x < 0.01 & telencephalon2$p_val_adj.y < 0.01 , ])

colnames(telencephalon2) <- c("Gene", "p-value 18hpf", "Avg Log2FC 18hpf", "pct 1 18hpf", "pct 2 18hpf", "Bonferroni p-value 18hpf", "cluster 18hpf", "p-value 20hpf", "Avg Log2FC 20hpf", "pct 1 20hpf", "pct 2 20hpf", "Bonferroni p-value 20hpf", "cluster 20hpf")
write.csv(telencephalon2, file = "telencephalon_markers.csv")



# neural crest
nc18 <-markers.18hpf[markers.18hpf$cluster == 11 | markers.18hpf$cluster == 16 & markers.18hpf$p_val_adj, ] # n= 760

nc20 <-markers.20hpf[markers.20hpf$cluster == 4 | markers.20hpf$cluster == 33, ] # n = 1030

cols_final <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") # names of the columns we want to keep


nc1 <- merge(nc18[, cols_final], nc20[, cols_final], by="gene")

nc2 <- nc1[nc1$p_val_adj.x < 0.01 | nc1$p_val_adj.y < 0.01 , ] # n = 322

dim(nc2[nc2$p_val_adj.x < 0.01 | nc2$p_val_adj.y < 0.01 , ]) 
dim(nc2[nc2$p_val_adj.x > 0.01 & nc2$p_val_adj.y < 0.01 , ])
dim(nc2[nc2$p_val_adj.x < 0.01 & nc2$p_val_adj.y > 0.01 , ])
dim(nc2[nc2$p_val_adj.x < 0.01 & nc2$p_val_adj.y < 0.01 , ])

colnames(nc2) <- c("Gene", "p-value 18hpf", "Avg Log2FC 18hpf", "pct 1 18hpf", "pct 2 18hpf", "Bonferroni p-value 18hpf", "cluster 18hpf", "p-value 20hpf", "Avg Log2FC 20hpf", "pct 1 20hpf", "pct 2 20hpf", "Bonferroni p-value 20hpf", "cluster 20hpf")
write.csv(nc2, file = "neural_crest_markers.csv")


# ventral diencephalon
vd18 <-markers.18hpf[markers.18hpf$cluster == 7 & markers.18hpf$p_val_adj, ] #all genes  in cluster 10 

vd20 <-markers.20hpf[markers.20hpf$cluster == 13, ] # all genes in either cluster 11 OR 12 

cols_final <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") # names of the columns we want to keep


vd1 <- merge(vd18[, cols_final], vd20[, cols_final], by="gene")

vd2 <- vd1[vd1$p_val_adj.x < 0.01 | vd1$p_val_adj.y < 0.01 , ] # n = 73

dim(vd2[vd2$p_val_adj.x < 0.01 | vd2$p_val_adj.y < 0.01 , ]) 
dim(vd2[vd2$p_val_adj.x > 0.01 & vd2$p_val_adj.y < 0.01 , ])
dim(vd2[vd2$p_val_adj.x < 0.01 & vd2$p_val_adj.y > 0.01 , ])
dim(vd2[vd2$p_val_adj.x < 0.01 & vd2$p_val_adj.y < 0.01 , ])

colnames(vd2) <- c("Gene", "p-value 18hpf", "Avg Log2FC 18hpf", "pct 1 18hpf", "pct 2 18hpf", "Bonferroni p-value 18hpf", "cluster 18hpf", "p-value 20hpf", "Avg Log2FC 20hpf", "pct 1 20hpf", "pct 2 20hpf", "Bonferroni p-value 20hpf", "cluster 20hpf")
write.csv(vd2, file = "ventral_diencephalon_markers.csv")

# tpn
tpn18 <-markers.18hpf[markers.18hpf$cluster == 19 & markers.18hpf$p_val_adj, ] #all genes  in cluster 10 

tpn20 <-markers.20hpf[markers.20hpf$cluster == 27, ] # all genes in either cluster 11 OR 12 

cols_final <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") # names of the columns we want to keep

tpn1 <- merge(tpn18[, cols_final], tpn20[, cols_final], by="gene")

tpn2 <- tpn1[tpn1$p_val_adj.x < 0.01 | tpn1$p_val_adj.y < 0.01 , ] # n = 161

dim(tpn2[tpn2$p_val_adj.x < 0.01 | tpn2$p_val_adj.y < 0.01 , ]) 
dim(tpn2[tpn2$p_val_adj.x > 0.01 & tpn2$p_val_adj.y < 0.01 , ])
dim(tpn2[tpn2$p_val_adj.x < 0.01 & tpn2$p_val_adj.y > 0.01 , ])
dim(tpn2[tpn2$p_val_adj.x < 0.01 & tpn2$p_val_adj.y < 0.01 , ])


colnames(tpn2) <- c("Gene", "p-value 18hpf", "Avg Log2FC 18hpf", "pct 1 18hpf", "pct 2 18hpf", "Bonferroni p-value 18hpf", "cluster 18hpf", "p-value 20hpf", "Avg Log2FC 20hpf", "pct 1 20hpf", "pct 2 20hpf", "Bonferroni p-value 20hpf", "cluster 20hpf")
write.csv(tpn2, file = "tpn_markers.csv")


# hindbrain - R5/6
hb18 <-markers.18hpf[markers.18hpf$cluster == 12 & markers.18hpf$p_val_adj, ] #all genes  in cluster 10 

hb20 <-markers.20hpf[markers.20hpf$cluster == 6, ] # all genes in either cluster 11 OR 12 

cols_final <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") # names of the columns we want to keep


hb1 <- merge(hb18[, cols_final], hb20[, cols_final], by="gene")

hb2 <- hb1[hb1$p_val_adj.x < 0.01 | hb1$p_val_adj.y < 0.01 , ] #now filter by significance

dim(hb2[hb2$p_val_adj.x < 0.01 | hb2$p_val_adj.y < 0.01 , ]) 
dim(hb2[hb2$p_val_adj.x > 0.01 & hb2$p_val_adj.y < 0.01 , ])
dim(hb2[hb2$p_val_adj.x < 0.01 & hb2$p_val_adj.y > 0.01 , ])
dim(hb2[hb2$p_val_adj.x < 0.01 & hb2$p_val_adj.y < 0.01 , ])

colnames(hb2) <- c("Gene", "p-value 18hpf", "Avg Log2FC 18hpf", "pct 1 18hpf", "pct 2 18hpf", "Bonferroni p-value 18hpf", "cluster 18hpf", "p-value 20hpf", "Avg Log2FC 20hpf", "pct 1 20hpf", "pct 2 20hpf", "Bonferroni p-value 20hpf", "cluster 20hpf")
write.csv(hb2, file = "hindbrain_markers.csv")


# neurons 
neuron18 <-markers.18hpf[markers.18hpf$cluster == 15 & markers.18hpf$p_val_adj, ] #all genes  in cluster 10 

neuron20 <-markers.20hpf[markers.20hpf$cluster == 2 | markers.20hpf$cluster == 3, ] # all genes in either cluster 11 OR 12 

cols_final <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") # names of the columns we want to keep


neuron1 <- merge(neuron18[, cols_final], neuron20[, cols_final], by="gene")

neuron2 <- neuron1[neuron1$p_val_adj.x < 0.01 | neuron1$p_val_adj.y < 0.01 , ] #now filter by significance

dim(neuron2[neuron2$p_val_adj.x < 0.01 | neuron2$p_val_adj.y < 0.01 , ]) 
dim(neuron2[neuron2$p_val_adj.x > 0.01 & neuron2$p_val_adj.y < 0.01 , ])
dim(neuron2[neuron2$p_val_adj.x < 0.01 & neuron2$p_val_adj.y > 0.01 , ])
dim(neuron2[neuron2$p_val_adj.x < 0.01 & neuron2$p_val_adj.y < 0.01 , ])

colnames(neuron2) <- c("Gene", "p-value 18hpf", "Avg Log2FC 18hpf", "pct 1 18hpf", "pct 2 18hpf", "Bonferroni p-value 18hpf", "cluster 18hpf", "p-value 20hpf", "Avg Log2FC 20hpf", "pct 1 20hpf", "pct 2 20hpf", "Bonferroni p-value 20hpf", "cluster 20hpf")
write.csv(neuron2, file = "neuron_markers.csv")


