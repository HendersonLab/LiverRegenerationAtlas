library(Seurat)
library(SeuratPipe)
library(dplyr)
library(ggplot2)
source("feature_cor_plot.r")
###########################
# Age Analysis
############################

####### Hepatoyctes

# Load Data 
obj <- readRDS("/Human/Data/Hepatocytes/seu_harmony_seed12_npcs40_r14.rds")
seu <- obj$seu
seu <- add_umap_embedding(seu, "/Human/Data/Hepatocytes/umap_embedding.csv")
anno <- read.csv("/Human/Data/Hepatocytes/seu_meta_res0.4_anno.csv", row.names = 1)
seu <- AddMetaData(seu, anno)
seu$seurat_clusters <- as.factor(seu$seurat_clusters)
Idents(seu) <- "seurat_clusters"

# PCs
png("/Age_Analysis/Hepatocytes_Age_PCs.png", width = 15, height = 6, res = 500, units = "in")
print(feature_cor_plot(seu, features = "age", reduction = "pca"))
dev.off()

# Harmony
png("/Age_Analysis/Hepatocytes_Age_HCs.png", width = 15, height = 6, res = 500, units = "in")
print(feature_cor_plot(seu, features = "age", reduction = "harmony"))
dev.off()


