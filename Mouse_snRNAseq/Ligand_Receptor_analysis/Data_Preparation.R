library(Seurat)
library(SeuratPipe)
library(dplyr)

###-------------------------------------------------
# Data Preparation for Lingand Receptor Analysis
###-------------------------------------------------

# Load Data

## All Lineage

obj <- readRDS("/Mouse/Data/All_Lineages/seu_harmony_npcs80.rds")
anno <- read.csv("/Mouse/Data/All_Lineages/seu_meta_res0.4_anno.csv", row.names = 1)
all <- obj$seu
all <- add_umap_embedding(all, "/Mouse/Data/All_Lineages/umap_embedding.csv")
all <- AddMetaData(all, anno)
dim_plot(all, group.by = "anno", label = T)

## Lineage Annotation Files


Endothelia <- read.csv("/Mouse/Data/Endothelia/seu_meta_res0.3_anno.csv", row.names = 1)
Hepatocytes <- read.csv("/Mouse/Data/Hepatocytes/seu_meta_res0.2.csv", row.names = 1)
Mesenchyme <- read.csv("/Mouse/Data/Mesenchyme/seu_meta_res0.8.csv", row.names = 1)
MPs <- read.csv("/Mouse/Data/MPs/seu_meta_res0.2.csv", row.names = 1)

### Add Cell Type to annotation Names
Endothelia$seurat_clusters <- paste0("End_",Endothelia$seurat_clusters)
Hepatocytes$seurat_clusters <- paste0("Hep_",Hepatocytes$seurat_clusters)
Mesenchyme$seurat_clusters <- paste0("Mes_",Mesenchyme$seurat_clusters)
MPs$seurat_clusters <- paste0("MP_",MPs$seurat_clusters)

# Removing Cells

## Remove Cells that were removed during each lineage analysis

all$keep <- "TRUE" # Set gobal variable that can be changed to false for cell removal
all$keep[all$anno == "Endothelia"] <- colnames(all)[all$anno == "Endothelia"] %in% rownames(Endothelia) # Change variable to False for cells not in the annotation file
all$keep[all$anno == "Hepatocytes"] <- colnames(all)[all$anno == "Hepatocytes"] %in% rownames(Hepatocytes) # Change variable to False for cells not in the annotation file
all$keep[all$anno == "Mesenchyme"] <- colnames(all)[all$anno == "Mesenchyme"] %in% rownames(Mesenchyme) # Change variable to False for cells not in the annotation file
all$keep[all$anno == "MPs"] <- colnames(all)[all$anno == "MPs"] %in% rownames(MPs) # Change variable to False for cells not in the annotation file
# Cells from the ILCs, Plasma, B cell and T cell, Cholangioytes, pDCs,Mesothelia lineages will be kept as is.
dim_plot(all, group.by = "keep")
all <- subset(all, subset = keep == "TRUE")
dim_plot(all, group.by = "keep", label = T)

# Add Cluster Annotations

all <- AddMetaData(all, metadata = Endothelia[,"seurat_clusters",drop=FALSE], col.name = "Endothelia") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = Hepatocytes[,"seurat_clusters",drop=FALSE], col.name = "Hepatocytes") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = Mesenchyme[,"seurat_clusters",drop=FALSE], col.name = "Mesenchyme") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = MPs[,"seurat_clusters",drop=FALSE], col.name = "MPs") # Add cluster specific metadata to seurat object
df <- all@meta.data %>% mutate(cluster_anno = coalesce(Endothelia,Hepatocytes,Mesenchyme,MPs)) %>% # Merge all these columns into one and add back in as merge cluster annoations
  select(cluster_anno)
all <- AddMetaData(all, df)
all$cluster_anno[is.na(all$cluster_anno)] <- all$anno[is.na(all$cluster_anno)] # add in names of clusters for those lineages that werent subset and so dont have cluster annoations
dim_plot(all, group.by = "cluster_anno")
table(all$cluster_anno, all$anno)


table(all$condition)

# Subset apap - For CellChat 

apap <- subset(all, subset = condition == "apap")
table(apap$condition)
saveRDS(apap, "~/Mouse/apap.RDS")


