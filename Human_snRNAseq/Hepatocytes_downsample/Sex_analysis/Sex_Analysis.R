library(Seurat)
library(SeuratPipe)
library(dplyr)
library(ggplot2)
source("feature_cor_plot.r")
###########################
# Sex Analysis
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


png("/Sex_Analysis/Hepatocytes_Sex_Split_Plot.png", width=10,height=5,res=500,units="in")
print(dim_plot(seu, group.by = "Sex", split.by = "Sex"))
dev.off()


# Contribution plots

cp <- scales::hue_pal()(2)
cl_sample <- seu@meta.data |> group_by(seurat_clusters, Sex) |> summarise(n = n()) |> mutate(freq = n / sum(n))
png("/Sex_Analysis/Hepatocytes_Sex_Contribution_plot.png", width = 14, height = 7, res = 500, units = "in")
print(ggplot(data = cl_sample, aes(x = seurat_clusters, y = freq, fill = Sex)) +
        geom_bar(stat = "identity", color="black", size = 0.05) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
        xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
        scale_fill_manual(values = cp))
dev.off()



# PCs
png("/Sex_Analysis/Hepatocytes_Sex_PCs.png", width = 15, height = 6, res = 500, units = "in")
print(feature_cor_plot(seu, features = "Sex_num", reduction = "pca"))
dev.off()

# Harmony
png("/Sex_Analysis/Hepatocytes_Sex_HCs.png", width = 15, height = 6, res = 500, units = "in")
print(feature_cor_plot(seu, features = "Sex_num", reduction = "harmony"))
dev.off()



