library(Seurat)
library(SeuratPipe)
library(dplyr)

###-------------------------------------------------
# Data Preparation for Lingand Receptor Analysis
###-------------------------------------------------

# Load Data

## All Lineage

obj <- readRDS("/Human/Data/All_Lineages/seu_harmony_npcs80.rds")
anno <- read.csv("/Human/Data/All_Lineages/seu_meta_res0.3_anno.csv", row.names = 1)
all <- obj$seu
all <- add_umap_embedding(all, "/Human/Data/All_Lineages/umap_embedding.csv")
all <- AddMetaData(all, anno)
dim_plot(all, group.by = "anno_broad", label = T)

## Lineage Annotation Files


Cholangiocytes <- read.csv("/Human/Data/Cholangiocytes/seu_meta_res0.3_anno.csv", row.names = 1)
Endothelia <- read.csv("/Human/Data/Endothelia/seu_meta_res0.7_anno.csv", row.names = 1)
Hepatocytes <- read.csv("/Human/Data/Hepatocytes/seu_meta_res0.4_anno.csv", row.names = 1)
Mesenchyme <- read.csv("/Human/Data/Mesenchyme/seu_meta_res0.1_anno.csv", row.names = 1)
MPs <- read.csv("/Human/Data/MPs/seu_meta_res0.2_anno.csv", row.names = 1)

### Add Cell Type to annotation Names
Cholangiocytes$anno_name <- paste0("Cho_",Cholangiocytes$anno)
Endothelia$anno_name <- paste0("End_",Endothelia$anno)
Hepatocytes$anno_name <- paste0("Hep_",Hepatocytes$anno)
Mesenchyme$anno_name <- paste0("Mes_",Mesenchyme$anno)
MPs$anno_name <- paste0("MP_",MPs$anno)

### Add Cell Type to annotation Numbers
Cholangiocytes$anno_number <- paste0("Cho_",Cholangiocytes$seurat_clusters)
Endothelia$anno_number <- paste0("End_",Endothelia$seurat_clusters)
Hepatocytes$anno_number <- paste0("Hep_",Hepatocytes$seurat_clusters)
Mesenchyme$anno_number <- paste0("Mes_",Mesenchyme$seurat_clusters)
MPs$anno_number <- paste0("MP_",MPs$seurat_clusters)

# Removing Cells

##  Remove Unidentified Cells

Idents(all) <- "anno"
all <- subset(all, idents = "Unidentified", invert = TRUE)
dim_plot(all, group.by = "anno", label = T)

## Remove Cells that were removed during each lineage analysis

all$keep <- "TRUE" # Set gobal variable that can be changed to false for cell removal
all$keep[all$anno == "Cholangiocytes"] <- colnames(all)[all$anno == "Cholangiocytes"] %in% rownames(Cholangiocytes) # Change variable to False for cells not in the annotation file
all$keep[all$anno == "Endothelia"] <- colnames(all)[all$anno == "Endothelia"] %in% rownames(Endothelia) # Change variable to False for cells not in the annotation file
all$keep[all$anno == "Hepatocytes"] <- colnames(all)[all$anno == "Hepatocytes"] %in% rownames(Hepatocytes) # Change variable to False for cells not in the annotation file
all$keep[all$anno == "Mesenchyme"] <- colnames(all)[all$anno == "Mesenchyme"] %in% rownames(Mesenchyme) # Change variable to False for cells not in the annotation file
all$keep[all$anno == "MPs"] <- colnames(all)[all$anno == "MPs"] %in% rownames(MPs) # Change variable to False for cells not in the annotation file
# Cells from the ILCs, Plasma, B cell and T cell lineages will be kept as is.
dim_plot(all, group.by = "keep")
all <- subset(all, subset = keep == "TRUE")
dim_plot(all, group.by = "keep", label = T)

# Add Cluster Annotations by name 

all <- AddMetaData(all, metadata = Cholangiocytes[,"anno_name",drop=FALSE], col.name = "Cholangiocytes_name") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = Endothelia[,"anno_name",drop=FALSE], col.name = "Endothelia_name") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = Hepatocytes[,"anno_name",drop=FALSE], col.name = "Hepatocytes_name") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = Mesenchyme[,"anno_name",drop=FALSE], col.name = "Mesenchyme_name") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = MPs[,"anno_name",drop=FALSE], col.name = "MPs_name") # Add cluster specific metadata to seurat object
df <- all@meta.data %>% mutate(cluster_anno_name = coalesce(Cholangiocytes_name,Endothelia_name,Hepatocytes_name,Mesenchyme_name,MPs_name)) %>% # Merge all these columns into one and add back in as merge cluster annoations
  select(cluster_anno_name)
all <- AddMetaData(all, df)
all$cluster_anno_name[is.na(all$cluster_anno_name)] <- all$anno[is.na(all$cluster_anno_name)] # add in names of clusters for those lineages that werent subset and so dont have cluster annoations
dim_plot(all, group.by = "cluster_anno_name")
table(all$cluster_anno_name, all$anno)



# Add Cluster Annotations by Number 

all <- AddMetaData(all, metadata = Cholangiocytes[,"anno_number",drop=FALSE], col.name = "Cholangiocytes_number") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = Endothelia[,"anno_number",drop=FALSE], col.name = "Endothelia_number") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = Hepatocytes[,"anno_number",drop=FALSE], col.name = "Hepatocytes_number") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = Mesenchyme[,"anno_number",drop=FALSE], col.name = "Mesenchyme_number") # Add cluster specific metadata to seurat object
all <- AddMetaData(all, metadata = MPs[,"anno_number",drop=FALSE], col.name = "MPs_number") # Add cluster specific metadata to seurat object
df <- all@meta.data %>% mutate(cluster_anno_number = coalesce(Cholangiocytes_number,Endothelia_number,Hepatocytes_number,Mesenchyme_number,MPs_number)) %>% # Merge all these columns into one and add back in as merge cluster annoations
  select(cluster_anno_number)
all <- AddMetaData(all, df)
all$cluster_anno_number[is.na(all$cluster_anno_number)] <- all$anno[is.na(all$cluster_anno_number)] # add in names of clusters for those lineages that werent subset and so dont have cluster annoations
dim_plot(all, group.by = "cluster_anno_number")
table(all$cluster_anno_number, all$anno)



# Subset apap - For CellChat 

apap <- subset(all, subset = condition == "apap")
table(apap$condition)
saveRDS(apap, "~/Human/apap.RDS")


