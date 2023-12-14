
###-------------------------------------------------------------
# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(ggplot2))

# Source utils scripts
source("zzz_utils.R")


###---------------------------------------------------------------------------
# Global settings
###---------------------------------------------------------------------------

# Global colour pallete
cp <- SeuratPipe:::discrete_col_pal
condition_col_pal <- c("lightskyblue3", "indianred3")

# Output directory for all figures
fig_dir <- paste0(results_dir, "/zzz_manuscript_figs/")
if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive = TRUE)}


###---------------------------------------------------------------------------
# Whole lineage analysis
###---------------------------------------------------------------------------

npcs <- 80
data_dir <- paste0(results_dir, "/01_integrate/00_iter1/")
obj <- readRDS(paste0(data_dir, "/seu_harmony_npcs", npcs, ".rds"))
seu <- obj$seu
opts <- obj$opts
rm(obj)
gc(verbose = FALSE)
opts$npcs <- npcs
opts$ndims <- 70 # Top PC dimensions to perform UMAP and clustering
opts$res <- 0.4 # Clustering resolution
opts$diff_cluster_pct <- 0.1
opts$topn_genes <- 10
out_dir <- paste0(data_dir, "/seu_harmony_npcs", npcs, "/dim", opts$ndims, "/")

## Colour palette for clusters
cluster_col_pal <- c(
  cp[1], "indianred2", "lightsalmon3", cp[5], cp[7], cp[27], cp[4], "orangered3", cp[29],
  "indianred1", cp[10], cp[31], cp[22], "brown2", "burlywood2", "rosybrown2", "darksalmon",
  cp[6], "salmon2", cp[2], "plum", "lightpink2", cp[3], cp[9], "maroon", cp[26])

anno_col_pal <- cp[c(1, 2, 5, 7, 4, 6, 3, 9, 10)]

# Run UMAP
seu <- add_umap_embedding(seu = seu, embedding = paste0(out_dir, "umap_embedding.csv"))
# Marker genes per cluster
mark <- read.csv(file = paste0(out_dir, "/clusters/seu_markers_res", opts$res, ".csv"))
# Make cluster a factor from character
mark$cluster <- factor(mark$cluster)

# Add cluster info
seu_anno <- read.csv(
  file = paste0(out_dir, "/clusters/seu_meta_res", opts$res, ".csv")) |>
  tibble::column_to_rownames(var = "X") |>
  dplyr::select(seurat_clusters) |>
  dplyr::mutate_at("seurat_clusters", factor)
# Add metadata information to current Seurat object and filter
seu <- AddMetaData(seu, metadata = seu_anno)

# Add anno info
seu_anno <- read.csv(
  file = paste0(out_dir, "/clusters/seu_meta_res", opts$res, "_anno.csv")) |>
  tibble::column_to_rownames(var = "X") |>
  dplyr::select(anno) |>
  dplyr::mutate_at("anno", factor)
# Add metadata information to current Seurat object and filter
seu <- AddMetaData(seu, metadata = seu_anno)
seu$anno <- factor(
  seu$anno,
  levels = c("Hepatocytes", "Cholangiocytes", "Endothelia", "Mesenchyme",
             "MPs", "T cells/ILCs", "Mesothelia", "B cells", "pDCs"),
  labels = c("Hepatocytes", "Cholangiocytes", "Endothelia", "Mesenchyme",
             "MPs", "T cells/ILCs", "Mesothelia", "B cells", "pDCs"),
  ordered = TRUE)
# Idents(object = seu) <- "anno"


## UMAP anno
png(paste0(fig_dir, "SFig4_mouse_all_anno.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = anno_col_pal,
              label = FALSE, label.size = 4, legend.position = "none",
              pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_all_anno_label.png"), width = 12, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = anno_col_pal,
                       label = FALSE, label.size = 4, legend.position = "right",
                       pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + NoAxes())
dev.off()


## UMAP spliut by condition
png(paste0(fig_dir, "SFig4_mouse_all_condition_split.png"), width = 13, height = 6, units = "in", res = 1000)
print(subset_dim_plot(
  seu = seu, subset.by = "condition", reduction = "umap",
  ncol = 2, col_pal = condition_col_pal, pt.size = 1.4, pt.stroke = 0.04,
  back.pt.size = 0.3, back.pt.alpha = 0.1,
  back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes() & ggplot2::labs(title = NULL))
dev.off()

png(paste0(fig_dir, "SFig4_mouse_all_condition_split_label.png"), width = 13, height = 6, units = "in", res = 1000)
print(subset_dim_plot(
  seu = seu, subset.by = "condition", reduction = "umap",
  ncol = 2, col_pal = condition_col_pal, pt.size = 1.4, pt.stroke = 0.04,
  back.pt.size = 0.3, back.pt.alpha = 0.1,
  back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes())
dev.off()


# UMAP cluster
png(paste0(fig_dir, "SFig4_mouse_all_cluster.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cluster_col_pal,
                       label = FALSE, label.size = 4, legend.position = "none",
                       pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_all_cluster_label.png"), width = 11, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cluster_col_pal,
                       label = TRUE, label.size = 7, legend.position = "right",
                       pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + NoAxes())
dev.off()


# ANXA2 expression
png(paste0(fig_dir, "SFig6_mouse_all_Anxa2.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "Anxa2", col_pal = "RdYlBu",
                            legend.position = "none", pt.size = 0.8, pt.alpha = 1,
                            pt.shape = 21, pt.stroke = 0.03, order_points_by_value = TRUE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig6_mouse_all_Anxa2_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "Anxa2", col_pal = "RdYlBu",
                            legend.position = "none", pt.size = 0.8, pt.alpha = 1,
                            pt.shape = 21, pt.stroke = 0.03, order_points_by_value = TRUE) + NoAxes())
dev.off()

# QC
png(paste0(fig_dir, "SFig4_mouse_all_qc_nfeature.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 0.8, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.03, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_all_qc_nfeature_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 0.8, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.03, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_all_qc_mito.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 0.8, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.03, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_all_qc_mito_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                           legend.position = "top", pt.size = 0.8, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.03, order_points_by_value = FALSE) + NoAxes())
dev.off()


# Integration mixing of timepoints
pdf(paste0(fig_dir, "SFig4_mouse_all_timepoint_mixing.pdf"), width = 10, height = 4.5)
cl_sample <- seu@meta.data |> group_by(anno, timepoint) |> summarise(n = n()) |> mutate(freq = n / sum(n))
print(ggplot(data = cl_sample, aes(x = anno, y = freq, fill = timepoint)) +
  geom_bar(stat = "identity", color="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
  scale_fill_manual(values = cp))
dev.off()



# Clustering heatmap
Idents(object = seu) <- "seurat_clusters"
sub <- subset(seu, downsample = 1500)
local_heatmap_plot(seu = sub, markers = mark, topn_genes = 5,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig4_mouse_all_heatmap.png"),
                   fig.res = 600,
                   cluster_col_pal = cluster_col_pal,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = TRUE)

local_heatmap_plot(seu = sub, markers = mark, topn_genes = 5,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig4_mouse_all_heatmap_label.png"),
                   fig.res = 600,
                   cluster_col_pal = cluster_col_pal,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = FALSE)

# Dot plot module scores
for (m in names(mouse_liver_lineage_markers)) {
  seu <- compute_module_score(
    seu = seu, features = mouse_liver_lineage_markers[[m]], name = m)
}
feats <- c("lin_cycling", "lin_mesothelia", "lin_pdc", "lin_bcell", "lin_ilc", "lin_tcell", "lin_mp",
           "lin_mesenchyme", "lin_endo", "lin_cholangiocyte", "lin_hepatocyte")
feats_labels <- c("Cycling cells", "Mesothelia", "pDCs", "B cells", "ILCs", "T cells", "MPs",
                  "Mesenchyme", "Endothelia", "Cholangiocytes", "Hepatocytes")
seu$clusters <- seu$seurat_clusters
seu$clusters <- factor(
  seu$clusters,
  levels = c(0, 1, 2, 5, 7, 8, 9, 11, 12, 13, 16, 18, 24, 19, 3, 15, 21, 4, 6, 10, 20, 17, 23, 25, 22, 14),
  labels = c(0, 1, 2, 5, 7, 8, 9, 11, 12, 13, 16, 18, 24, 19, 3, 15, 21, 4, 6, 10, 20, 17, 23, 25, 22, 14))
Idents(object = seu) <- "clusters"

pdf(paste0(fig_dir, "SFig4_mouse_all_dotplot.pdf"), width = 11, height = 4.5)
print(dot_plot(seu = seu, features = feats, group.by = "clusters",
         labels = feats_labels, xlab = "Signature", ylab = "Cluster",
         legend.position = "right", col_pal = "OrRd",
         col.min = 0, col.max = 2, direction = 1))
dev.off()




###---------------------------------------------------------------------------
# Hepatocyte all conditions analysis
###---------------------------------------------------------------------------
lineage <- "Hepatocytes"
npcs <- 80
seed <- 1
pcs_to_remove <- NULL
npcs_dir <- npcs
if(!is.null(pcs_to_remove)) {
  npcs_dir <- paste0(npcs_dir, "_r", paste(pcs_to_remove, collapse = ""))
}
data_dir <- paste0(results_dir, lineage, "/01_integrate/00_iter1/")
obj <- readRDS(paste0(data_dir, "/seu_harmony_seed", seed, "_npcs", npcs_dir, ".rds"))
seu <- obj$seu
opts <- obj$opts
rm(obj)
gc(verbose = FALSE)
opts$npcs <- npcs
opts$ndims <- 80 # Top PC dimensions to perform UMAP and clustering
opts$res <- 0.2 # Clustering resolution
opts$diff_cluster_pct <- 0.1
opts$topn_genes <- 35
out_dir <- paste0(data_dir, "/seu_harmony_seed",
                  seed, "_npcs", npcs_dir, "/dim", opts$ndims, "/")

## Colour palette for clusters
cluster_col_pal <- cp[c(1, 3, 9, 5, 2, 7, 8, 6, 4)]
anno_col_pal <- cluster_col_pal

# Run UMAP
seu <- add_umap_embedding(seu = seu, embedding = paste0(out_dir, "umap_embedding.csv"))
# Marker genes per cluster
mark <- read.csv(file = paste0(out_dir, "/clusters/seu_markers_res", opts$res, ".csv"))
# Make cluster a factor from character
mark$cluster <- factor(mark$cluster)
mark$anno <- plyr::mapvalues(
  mark$cluster,
  from = seq(0, 8),
  to = c("Portal", "MitoHigh", "Central2", "EarlyResponse", "Central1",
         "Migratory", "Stressed", "Cycling", "STAT1+"))
migratory <- mark |> filter(anno == "Migratory") |>
  dplyr::arrange(-avg_log2FC) |>
  top_n(50, avg_log2FC) |> dplyr::pull(gene)
saveRDS(object = migratory, file = paste0(anno_dir, "/hep_mouse_markers/mouse_hep_migratory_markers.rds"))


# Add cluster info
seu_anno <- read.csv(
  file = paste0(out_dir, "/clusters/seu_meta_res", opts$res, ".csv")) |>
  tibble::column_to_rownames(var = "X") |>
  dplyr::select(seurat_clusters) |>
  dplyr::mutate_at("seurat_clusters", factor)
# Add metadata information to current Seurat object and filter
seu <- AddMetaData(seu, metadata = seu_anno)


# Add anno info
seu$anno <- plyr::mapvalues(
  seu$seurat_clusters,
  from = seq(0, 8),
  to = c("Portal", "MitoHigh", "Central2", "EarlyResponse", "Central1",
         "Migratory", "Stressed", "Cycling", "STAT1+"))
seu$anno <- factor(
  seu$anno,
  levels = c("Portal", "MitoHigh", "Central2", "EarlyResponse", "Central1",
             "Migratory", "Stressed", "Cycling", "STAT1+"),
  labels = c("Portal", "MitoHigh", "Central2", "EarlyResponse", "Central1",
             "Migratory", "Stressed", "Cycling", "STAT1+"),
  ordered = TRUE)
seu$anno_ordered <- factor(
  seu$anno,
  levels = c("Portal", "Central1", "Central2", "MitoHigh", "EarlyResponse",
             "Migratory", "Cycling", "Stressed", "STAT1+"),
  labels = c("Portal", "Central1", "Central2", "MitoHigh", "EarlyResponse",
             "Migratory", "Cycling", "Stressed", "STAT1+"),
  ordered = TRUE)
# Idents(object = seu) <- "anno"


## UMAP anno
png(paste0(fig_dir, "Fig2_mouse_Hep_anno.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = anno_col_pal,
                        label = FALSE, label.size = 4, legend.position = "none",
                        pt.size = 2, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "Fig2_mouse_Hep_anno_label.png"), width = 12, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = anno_col_pal,
                        label = FALSE, label.size = 4, legend.position = "right",
                        pt.size = 2, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
dev.off()


# UMAP cluster
png(paste0(fig_dir, "SFig4_mouse_Hep_cluster.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = anno_col_pal,
                        label = FALSE, label.size = 4, legend.position = "none",
                        pt.size = 2, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_Hep_cluster_label.png"), width = 11, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = anno_col_pal,
                        label = TRUE, label.size = 7, legend.position = "right",
                        pt.size = 2, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + NoAxes())
dev.off()


# QC
png(paste0(fig_dir, "SFig4_mouse_Hep_qc_nfeature.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                            legend.position = "none", pt.size = 1.5, pt.alpha = 1,
                            pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_Hep_qc_nfeature_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                            legend.position = "top", pt.size = 1.5, pt.alpha = 1,
                            pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_Hep_qc_mito.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                            legend.position = "none", pt.size = 1.5, pt.alpha = 1,
                            pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig4_mouse_Hep_qc_mito_legend.png"),
    width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                            legend.position = "top", pt.size = 1.5, pt.alpha = 1,
                            pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
dev.off()



# Integration mixing of donors
cl_sample <- seu@meta.data |> group_by(anno_ordered, donor) |>
  summarise(n = n()) |> mutate(freq = n / sum(n))
pdf(paste0(fig_dir, "SFig4_mouse_Hep_donor_mixing_legend.pdf"), width = 10, height = 5)
print(ggplot(data = cl_sample, aes(x = anno_ordered, y = freq, fill = donor)) +
        geom_bar(stat = "identity", color="black", size = 0.1) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
        xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
        scale_fill_manual(values = cp))
dev.off()

# Integration mixing of donors
pdf(paste0(fig_dir, "SFig4_mouse_Hep_donor_mixing.pdf"), width = 10, height = 5)
print(ggplot(data = cl_sample, aes(x = anno_ordered, y = freq, fill = donor)) +
        geom_bar(stat = "identity", color="black", size = 0.1) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
        xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
        scale_fill_manual(values = cp))
dev.off()


# Clustering heatmap
Idents(object = seu) <- "seurat_clusters"
sub <- subset(seu, downsample = 1500)
local_heatmap_plot(seu = sub, markers = mark, topn_genes = 35,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig4_mouse_Hep_heatmap.png"),
                   fig.res = 300,
                   cluster_col_pal = cluster_col_pal,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = TRUE)

local_heatmap_plot(seu = sub, markers = mark, topn_genes = 35,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig4_mouse_Hep_heatmap_label.png"),
                   fig.res = 300,
                   cluster_col_pal = cluster_col_pal,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = FALSE)


# GO enrichment analysis
go_dir <- paste0(results_dir, lineage, "/geneset_analysis/", "npcs", npcs_dir,
                 "_seed", seed, "_dim", opts$ndims, "_res", opts$res, "/go_enrichment/outs/")
cls <- seu[[c("seurat_clusters", "anno")]] %>% dplyr::distinct() %>% dplyr::arrange(seurat_clusters)
for (cl in 1:NROW(cls)) {
  df <- read.csv(file = paste0(go_dir, "/Cluster", cls$seurat_clusters[cl], "_GO_BP_enrichment_level6.csv"))
  pdf(paste0(fig_dir, "Fig2_mouse_Hep_GO_", cls$anno[cl], ".pdf"), width = 5, height = 10)
  print(plot_enrichment(df = df, terms = NULL, pval_thresh = 0.05, max_terms = 50,
                        text_size = 0.7, dot_size = 2.0, label_format = 45))
  dev.off()
}

cl <- 6
df <- read.csv(file = paste0(go_dir, "/Cluster", cls$seurat_clusters[cl], "_GO_BP_enrichment_level6.csv"))
terms <- c("actin filament organization", "regulation of cell morphogenesis", "ameboidal-type cell migration",
           "focal adhesion assembly", "regulation of actin cytoskeleton organization",
           "epithelial cell migration", "cell-cell junction adhesion",
           "protein localization to plasma membrane")
df <- df |> dplyr::filter(Description %in% terms)
pdf(paste0(fig_dir, "Fig2_mouse_Hep_GO_", cls$anno[cl], "_selected_terms.pdf"), width = 5, height = 2.5)
print(plot_enrichment(df = df, terms = NULL, pval_thresh = 0.05, max_terms = 50,
                      text_size = 1, dot_size = 3.0, label_format = 35))
dev.off()

# Feature plot module scores and marker genes
seu <- compute_module_score(seu = seu, features = hep_migratory_livreg_top25_markers$hep_migratory_top25,
                            name = "Migratory_top25")
seu <- compute_module_score(seu = seu, features = hep_mouse_zonation_markers$hep_zonation_central,
                            name = "Zonation_central")
seu <- compute_module_score(seu = seu, features = hep_mouse_zonation_markers$hep_zonation_portal,
                            name = "Zonation_portal")
seu <- compute_module_score(seu = seu, features = hep_migratory_livreg_top25_markers$hep_migratory_top25,
                            name = "Migratory_Signature")

features <- c("Migratory_Signature", "Migratory_top25", "Zonation_central", "Zonation_portal", "Anxa2")
for (f in features) {
  png(paste0(fig_dir, "Fig2_mouse_Hep_score_order_", f, "_legend.png"),
      width = 7, height = 6, units = "in", res = 1000)
  print(feature_plot_tailored(seu = seu, feature = f, col_pal = "RdYlBu",
                              max.cutoff = "q98", min.cutoff = 0,
                              legend.position = "top", pt.size = 1.5, pt.alpha = 1,
                              pt.shape = 21, pt.stroke = 0.04) + NoAxes())
  dev.off()
}


for (f in features) {
  png(paste0(fig_dir, "Fig2_mouse_Hep_score_", f, "_legend.png"),
      width = 7, height = 6, units = "in", res = 1000)
  print(feature_plot_tailored(seu = seu, feature = f, col_pal = "RdYlBu",
                              max.cutoff = "q98", min.cutoff = 0,
                              legend.position = "top", pt.size = 1.5, pt.alpha = 1,
                              pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
  dev.off()
}


gg <- list()
for (f in features) {
  gg[[f]] <- subset_feature_plot(
    seu = seu, subset.by = "timepoint", feature = f, min.cutoff = 0,
    max.cutoff = "q98", reduction = "umap", slot = "data",
    ncol = NULL, col_pal = NULL, pt.size = 1.8, stroke = 0.05,
    legend.position = "top", back.pt.size = 0.3, back.alpha = 0.1,
    back.color = "grey", combine = FALSE)
}

library(magick)
for (f in features) {
  gif_out_dir <- paste0(fig_dir, "/zzz_gif_", f, "/")
  if (!dir.exists(gif_out_dir)) {dir.create(gif_out_dir, recursive = TRUE)}
  for (s in levels(seu$timepoint)) {
      png(filename = paste0(gif_out_dir, f, "_", s, ".png"),
          width = 7, height = 6, res = 200, units = "in")
      print(gg[[f]][[s]] + ggplot2::ggtitle(s) + Seurat::NoAxes())
      dev.off()
  }
  ## list file names and read in
  imgs <- list.files(gif_out_dir, pattern = "^.*png$", full.names = TRUE)
  idx <- gtools::mixedorder(imgs)
  imgs <- imgs[idx]
  img_list <- lapply(imgs, magick::image_read)
  ## join the images together
  img_joined <- magick::image_join(img_list)
  ## animate at 2 frames per second
  img_animated <- magick::image_animate(img_joined, fps = 1, delay = 70)
  ## save to disk
  magick::image_write(image = img_animated,
                      path = paste0(gif_out_dir, f, ".gif"))
}



## Joint Central and Portal Hep scores
gif_out_dir <- paste0(fig_dir, "/zzz_gif_Zonation/")
if (!dir.exists(gif_out_dir)) {dir.create(gif_out_dir, recursive = TRUE)}
for (s in levels(seu$timepoint)) {
  png(filename = paste0(gif_out_dir,  "Zonation_", s, ".png"), width = 14,
      height = 6, res = 200, units = "in")
  g1 <- gg[["Zonation_portal"]][[s]] + ggplot2::ggtitle(s) + Seurat::NoAxes()
  g2 <- gg[["Zonation_central"]][[s]] + ggplot2::ggtitle(s) + Seurat::NoAxes()
  to_plot <- list(g2, g1)
  print(patchwork::wrap_plots(to_plot, ncol = 2))
  dev.off()
}
## list file names and read in
imgs <- list.files(gif_out_dir, pattern = "^.*png$", full.names = TRUE)
idx <- gtools::mixedorder(imgs)
imgs <- imgs[idx]
img_list <- lapply(imgs, magick::image_read)
## join the images together
img_joined <- magick::image_join(img_list)
## animate at 2 frames per second
img_animated <- magick::image_animate(img_joined, fps = 1, delay = 70)
## save to disk
magick::image_write(image = img_animated,
                    path = paste0(gif_out_dir, "Zonation.gif"))




###---------------------------------------------------------------------------
# NPC analysis
###---------------------------------------------------------------------------

all_lineages <- c("Mesenchyme", "Endothelia", "MPs")
all_npcs <- c(20, 40, 30)
all_ndims <- c(10, 30, 20)
all_res <- c(0.8, 0.3, 0.2)

for (i in 1:length(all_lineages)) {
  lineage <- all_lineages[i]

  npcs <- all_npcs[i]
  data_dir <- paste0(results_dir, lineage, "/01_integrate/00_iter2/")
  obj <- readRDS(paste0(data_dir, "/seu_harmony_npcs", npcs, ".rds"))
  seu <- obj$seu
  opts <- obj$opts
  rm(obj)
  gc(verbose = FALSE)
  opts$npcs <- npcs
  opts$ndims <- all_ndims[i] # Top PC dimensions to perform UMAP and clustering
  opts$res <- all_res[i] # Clustering resolution
  opts$diff_cluster_pct <- 0.1
  opts$topn_genes <- 10
  out_dir <- paste0(data_dir, "/seu_harmony_npcs", npcs, "/dim", opts$ndims, "/")

  # Add UMAP
  seu <- SeuratPipe::add_umap_embedding(seu = seu, embedding = paste0(out_dir, "/umap_embedding.csv"))
  # Marker genes per cluster
  mark <- read.csv(file = paste0(out_dir, "/clusters/seu_markers_res", opts$res, ".csv"))
  # Make cluster a factor from character
  mark$cluster <- factor(mark$cluster)

  # Add cluster info
  seu_anno <- read.csv(
    file = paste0(out_dir, "/clusters/seu_meta_res", opts$res, ".csv")) |>
    tibble::column_to_rownames(var = "X") |>
    dplyr::select(seurat_clusters) |>
    dplyr::mutate_at("seurat_clusters", factor)
  # Add metadata information to current Seurat object and filter
  seu <- AddMetaData(seu, metadata = seu_anno)


  ## UMAP split by condition
  png(paste0(fig_dir, "SFig5_", lineage, "_condition_split.png"),
      width = 6, height = 10, units = "in", res = 1000)
  print(subset_dim_plot(
    seu = seu, subset.by = "condition", reduction = "umap",
    ncol = 1, col_pal = condition_col_pal, pt.size = 3, pt.stroke = 0.04,
    back.pt.size = 0.5, back.pt.alpha = 0.2,
    back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes() & ggplot2::labs(title = NULL))
  dev.off()

  png(paste0(fig_dir, "SFig5_", lineage, "_condition_split_label.png"),
      width = 6, height = 11, units = "in", res = 1000)
  print(subset_dim_plot(
    seu = seu, subset.by = "condition", reduction = "umap",
    ncol = 1, col_pal = condition_col_pal, pt.size = 3, pt.stroke = 0.04,
    back.pt.size = 0.5, back.pt.alpha = 0.2,
    back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes())
  dev.off()


  # UMAP cluster
  png(paste0(fig_dir, "SFig5_", lineage, "_cluster.png"), width = 10, height = 8, units = "in", res = 1000)
  print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cp,
                          label = FALSE, label.size = 4, legend.position = "none",
                          pt.size = 3, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
  dev.off()


  png(paste0(fig_dir, "SFig5_", lineage, "_cluster_label.png"), width = 11, height = 8, units = "in", res = 1000)
  print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cp,
                          label = FALSE, label.size = 7, legend.position = "right",
                          pt.size = 3, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
  dev.off()


  # Clustering heatmap
  Idents(object = seu) <- "seurat_clusters"
  #sub <- subset(seu, downsample = 1500)
  local_heatmap_plot(seu = seu, markers = mark, topn_genes = 10,
                     diff_cluster_pct = opts$diff_cluster_pct,
                     pval_adj = 0.05, filename = paste0(fig_dir, "SFig5_", lineage, "_heatmap.png"),
                     fig.res = 600,
                     cluster_col_pal = cp,
                     condition_col_pal = condition_col_pal,
                     donor_col_pal = cp, is_clean_plot = TRUE)

  local_heatmap_plot(seu = seu, markers = mark, topn_genes = 10,
                     diff_cluster_pct = opts$diff_cluster_pct,
                     pval_adj = 0.05, filename = paste0(fig_dir, "SFig5_", lineage, "_heatmap_label.png"),
                     fig.res = 600,
                     cluster_col_pal = cp,
                     condition_col_pal = condition_col_pal,
                     donor_col_pal = cp, is_clean_plot = FALSE)
}

