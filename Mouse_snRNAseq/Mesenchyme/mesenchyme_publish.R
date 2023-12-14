library(Seurat)
library(SeuratPipe)
library(dplyr)
library(tibble)

destination <- "/Analysis/Mouse_Nuclei_apap"

iter_npcs <- 20
iter_dims <- 10
iter_res <- 0.8
obj <- readRDS(paste0(destination, "/mesenchyme/ml_apap_mes_harmony_npcs", iter_npcs, ".rds"))
ml_apap_mes <- obj$seu
opts <- obj$opts
rm(obj)

# Read cluster annotation
ml_apap_mes_anno <- read.csv(
  file = paste0(destination, "/mesenchyme/ml_apap_mes_harmony_npcs", iter_npcs, "/dim", iter_dims,
                "/clusters/seu_meta_res", iter_res, ".csv")) %>%
  tibble::column_to_rownames(var = "X") %>%
  select(seurat_clusters) %>%
  mutate_at(c("seurat_clusters"), factor)

# Add metadata information to current Seurat object
ml_apap_mes <- AddMetaData(ml_apap_mes, metadata = ml_apap_mes_anno)
Idents(object = ml_apap_mes) <- "seurat_clusters"

# Perform data integration using Harmony
ml_apap_mes <- run_harmony_pipeline(
  seu_obj = ml_apap_mes,
  out_dir = "/Analysis/Mouse_Nuclei_apap/mesenchyme_publish/",
  batch_id = "sample",
  npcs = iter_npcs,
  ndims = iter_dims,
  res = iter_res,
  modules_group = list(),
  metadata_to_plot = c("sample", "condition", "doublet_prediction", "timepoint"),
  qc_to_plot = c("nFeature_RNA", "percent.mito", "doublet_score"),
  logfc.threshold = 0.7,
  min.pct = 0.25,
  only.pos = TRUE,
  topn_genes = 20,
  diff_cluster_pct = 0.1,
  pval_adj = 0.05,
  pcs_to_remove = NULL,
  obj_filename = "ml_apap_mes_harmony",
  force_reanalysis = FALSE,
  plot_cluster_markers = FALSE,
  max.cutoff = "q98",
  n_hvgs = 3000,
  seed = 1,
  discrete_col_pal = SeuratPipe:::discrete_col_pal,
  label = TRUE,
  label.size = 6,
  pt.size = 1.4,
  fig.res = 200
)


### TIFFs

library(ggplot2)
library(RColorBrewer)
library(magick)

ml_apap_mes.markers <- find_all_markers(seu=ml_apap_mes, random.seed=1, only.pos=T, min.pct=0.25, logfc.threshold=0.7)
write.csv(ml_apap_mes.markers, file=paste0(destination,"/mesenchyme_publish/csv_mes_markers.csv"))

colours_clus <- SeuratPipe:::discrete_col_pal[1:nlevels(ml_apap_mes$seurat_clusters)]
names(colours_clus) <- levels(ml_apap_mes$seurat_clusters)
colours_timepoint <- viridis::viridis(nlevels(ml_apap_mes$timepoint))
names(colours_timepoint) <- levels(ml_apap_mes$timepoint)

xlim_umap <- c(min(ml_apap_mes@reductions$umap@cell.embeddings[,1])-1, max(ml_apap_mes@reductions$umap@cell.embeddings[,1])+1)
ylim_umap <- c(min(ml_apap_mes@reductions$umap@cell.embeddings[,2])-1, max(ml_apap_mes@reductions$umap@cell.embeddings[,2])+1)


t <- UMAPPlot(ml_apap_mes, pt.size=0.1)
tiff(paste0(destination,"/mesenchyme_publish/umap_mes_clusters_clean.tiff"), width=7,height=7,res=300,units="in")
t <- t + theme_classic() + 
  theme(legend.position="right", legend.key.size=unit(1.5,"line"), legend.text=element_text(size=12), 
        axis.line=element_line(colour="black", size=1.25, linetype="solid"), 
        axis.text.x=element_blank(), axis.text.y=element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.ticks.length=unit(0.2,"cm")) +
  geom_point(aes(t$data$UMAP_1, t$data$UMAP_2, fill=ml_apap_mes$seurat_clusters), shape=21, size=1, stroke=0.1) +
  guides(size=FALSE, colour=FALSE, fill=FALSE) +
  scale_fill_manual(name="", values=colours_clus) +
  scale_x_continuous(name="", minor_breaks=NULL, limits=xlim_umap) + 
  scale_y_continuous(name="", minor_breaks=NULL, limits=ylim_umap)
print(t)
dev.off()
png(paste0(destination,"/mesenchyme_publish/umap_mes_clusters_clean.png"), width=7,height=7,res=300,units="in")
print(t)
dev.off()


dir.create(paste0(destination,"/mesenchyme_publish/timepoints"))
for (time in levels(ml_apap_mes$timepoint)) {
  t <- UMAPPlot(ml_apap_mes, cells=Cells(ml_apap_mes)[ml_apap_mes$timepoint==time], pt.size=0.1, group.by="timepoint")
  tiff(paste0(destination,"/mesenchyme_publish/timepoints/umap_mes_timepoint_",time,".tiff"), width=5,height=5,res=300,units="in")
  t <- t + theme_classic() + ggtitle(time) +
    theme(legend.position="right", legend.key.size=unit(1.5,"line"), legend.text=element_text(size=12),
          axis.line=element_line(colour="black", size=1.25, linetype="solid"), 
          axis.text.x=element_blank(), axis.text.y=element_blank(), 
          axis.title.x=element_blank(), axis.title.y=element_blank(), 
          axis.ticks.length=unit(0.2,"cm")) +
    geom_point(aes(t$data$UMAP_1, t$data$UMAP_2, fill=ml_apap_mes$timepoint[ml_apap_mes$timepoint==time]), shape=21, size=1, stroke=0.1) +
    guides(size=FALSE, colour=FALSE, fill=FALSE) +
    scale_fill_manual(name="", values=colours_timepoint) +
    scale_x_continuous(name="", minor_breaks=NULL, limits=xlim_umap) + 
    scale_y_continuous(name="", minor_breaks=NULL, limits=ylim_umap)
  print(t)
  dev.off()
}
imgs <- list.files(paste0(destination,"/mesenchyme_publish/timepoints/"))
idx <- gtools::mixedorder(imgs)
imgs <- imgs[idx]
imgs <- paste0(destination,"/mesenchyme_publish/timepoints/",imgs)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps=1, delay=70)
image_write(image=img_animated, path=paste0(destination,"/mesenchyme_publish/umap_mes_timepoints.gif"))


markers <- ml_apap_mes.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
markers <- markers[!duplicated(markers$gene),]
markers <- markers[markers$gene %in% rownames(ml_apap_mes@assays$RNA@scale.data),]
cell_order <- unlist(lapply(split(Cells(ml_apap_mes), ml_apap_mes$seurat_clusters), FUN=sample))
annotation_col <- data.frame(Cluster=ml_apap_mes$seurat_clusters[cell_order], 
                             Timepoint=ml_apap_mes$timepoint[cell_order])
rownames(annotation_col) <- cell_order
annotation_row <- data.frame(Cluster=markers$cluster)
rownames(annotation_row) <- markers$gene
annotation_colours <- list(Cluster=colours_clus,
                           Timepoint=colours_timepoint)
tmp <- ml_apap_mes@assays$RNA@scale.data[as.character(markers$gene[order(annotation_row[,1])]),
                                        rownames(annotation_col)[order(annotation_col[,1])]]
tmp[tmp>2.5] <- 2.5
tmp[tmp<-2.5] <- -2.5
gaps_col <- cumsum(summary(ml_apap_mes$seurat_clusters))
gaps_row <- cumsum(summary(markers$cluster))
pheatmap::pheatmap(tmp, color=colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(1000), 
                   breaks=seq(-2.5,2.5,0.005), gaps_row=gaps_row, gaps_col=gaps_col,
                   tree_height_row=NA, tree_height_col=NA, cluster_cols=FALSE, cluster_rows=FALSE,
                   annotation_col=annotation_col, annotation_row=annotation_row, annotation_colors=annotation_colours, 
                   annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=TRUE,
                   show_rownames=TRUE, show_colnames=FALSE, legend=TRUE)
dev.print(tiff, paste0(destination,"/mesenchyme_publish/heatmap_mes_labelled.tiff"), width=20,height=nrow(markers)*3/13,res=300,units="in")
dev.print(png, paste0(destination,"/mesenchyme_publish/heatmap_mes_labelled.png"), width=20,height=nrow(markers)*3/13,res=300,units="in")
