
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
condition_col_pal <- c("lightskyblue3", "indianred3", "navajowhite2")

# Output directory for all figures
fig_dir <- paste0(results_dir, "/zzz_manuscript_figs_new/")
if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive = TRUE)}


# Read all sample metadata, including sample anonymised IDs
sample_meta <- read.csv(file = paste0(raw_data_dir, "meta_hln_livreg.csv")) |>
  dplyr::filter(pass_qc == TRUE)
idx <- gtools::mixedorder(sample_meta$sample)
sample_meta <- sample_meta[idx, ]
rm(idx)

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
opts$ndims <- 80 # Top PC dimensions to perform UMAP and clustering
opts$res <- 0.3 # Clustering resolution
opts$diff_cluster_pct <- 0.1
opts$topn_genes <- 10
out_dir <- paste0(data_dir, "/seu_harmony_npcs", npcs, "/dim", opts$ndims, "/")



## Colour palette for clusters
cluster_col_pal <- c(
  cp[1], cp[5], cp[4], "indianred2", cp[2], cp[6], cp[7], "lightsalmon3",
  "burlywood2",  cp[8], "orangered3", cp[3], "rosybrown2", cp[10], "lightpink3",
  cp[9], "yellowgreen", "dodgerblue", "indianred1", cp[11])

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
  levels = c("Hepatocytes", "Cholangiocytes", "Endothelia", "Mesenchyme", "MPs",
             "T cells", "ILCs", "B cells", "Plasma", "Unidentified"),
  labels = c("Hepatocytes", "Cholangiocytes", "Endothelia", "Mesenchyme", "MPs",
             "T cells", "ILCs", "B cells", "Plasma", "Unidentified"),
  ordered = TRUE)
# Idents(object = seu) <- "anno"


## UMAP anno
anno_col_pal <- cp[c(1, 2, 5, 7, 4, 6, 3, 9, 10, 8)]
png(paste0(fig_dir, "Fig1_all_anno.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = anno_col_pal,
              label = FALSE, label.size = 4, legend.position = "none",
              pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.03) + NoAxes())
dev.off()

png(paste0(fig_dir, "Fig1_all_anno_label.png"), width = 12, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = anno_col_pal,
                       label = FALSE, label.size = 4, legend.position = "right",
                       pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.03) + NoAxes())
dev.off()


## UMAP split by condition
png(paste0(fig_dir, "Fig1_all_condition_split.png"), width = 19, height = 6, units = "in", res = 1000)
print(subset_dim_plot(
  seu = seu, subset.by = "condition", reduction = "umap",
  ncol = 3, col_pal = condition_col_pal, pt.size = 1.4, pt.stroke = 0.03,
  back.pt.size = 0.3, back.pt.alpha = 0.1,
  back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes() & ggplot2::labs(title = NULL))
dev.off()

png(paste0(fig_dir, "Fig1_all_condition_split_label.png"), width = 19, height = 6, units = "in", res = 1000)
print(subset_dim_plot(
  seu = seu, subset.by = "condition", reduction = "umap",
  ncol = 3, col_pal = condition_col_pal, pt.size = 1.4, pt.stroke = 0.03,
  back.pt.size = 0.3, back.pt.alpha = 0.1,
  back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes())
dev.off()


# UMAP cluster
png(paste0(fig_dir, "SFig1_all_cluster.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cluster_col_pal,
                       label = FALSE, label.size = 4, legend.position = "none",
                       pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.03) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig1_all_cluster_diverging_cols.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cp,
                       label = FALSE, label.size = 4, legend.position = "none",
                       pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.03) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig1_all_cluster_label.png"), width = 11, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cluster_col_pal,
                       label = TRUE, label.size = 7, legend.position = "right",
                       pt.size = 1.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.03) + NoAxes())
dev.off()


# QC
png(paste0(fig_dir, "SFig1_all_qc_nfeature.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 0.6, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.03, order_points_by_value = FALSE) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "SFig1_all_qc_nfeature_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 0.6, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.03, order_points_by_value = FALSE) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "SFig1_all_qc_mito.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 0.6, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.03, order_points_by_value = FALSE) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "SFig1_all_qc_mito_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                           legend.position = "top", pt.size = 0.6, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.03, order_points_by_value = FALSE) + Seurat::NoAxes())
dev.off()


# Integration mixing of donors
cl_sample <- seu@meta.data |> group_by(anno, donor) |> summarise(n = n()) |> mutate(freq = n / sum(n))
pdf(paste0(fig_dir, "SFig1_all_donor_mixing_legend.pdf"), width = 10, height = 4.5)
print(ggplot(data = cl_sample, aes(x = anno, y = freq, fill = donor)) +
  geom_bar(stat = "identity", color="black", size = 0.05) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
  scale_fill_manual(values = cp))
dev.off()

pdf(paste0(fig_dir, "SFig1_all_donor_mixing.pdf"), width = 8, height = 4.5)
print(ggplot(data = cl_sample, aes(x = anno, y = freq, fill = donor)) +
        geom_bar(stat = "identity", color="black", size = 0.05) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
        xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
        scale_fill_manual(values = cp))
dev.off()


# Integration mixing of donors
cl_sample <- seu@meta.data |> group_by(seurat_clusters, donor) |> summarise(n = n()) |> mutate(freq = n / sum(n))
pdf(paste0(fig_dir, "SFig1_all_donor_mixing_cluster_legend.pdf"), width = 14, height = 4.5)
print(ggplot(data = cl_sample, aes(x = seurat_clusters, y = freq, fill = donor)) +
        geom_bar(stat = "identity", color="black", size = 0.05) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
        xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
        scale_fill_manual(values = cp))
dev.off()

pdf(paste0(fig_dir, "SFig1_all_donor_mixing_cluster.pdf"), width = 12, height = 4.5)
print(ggplot(data = cl_sample, aes(x = seurat_clusters, y = freq, fill = donor)) +
        geom_bar(stat = "identity", color="black", size = 0.05) +
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
local_heatmap_plot(seu = sub, markers = mark, topn_genes = 5,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig1_all_heatmap.png"),
                   fig.res = 600,
                   cluster_col_pal = cluster_col_pal,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = TRUE)

local_heatmap_plot(seu = sub, markers = mark, topn_genes = 5,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig1_all_heatmap_label.png"),
                   fig.res = 600,
                   cluster_col_pal = cluster_col_pal,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = FALSE)


# Dot plot module scores
for (m in names(human_liver_lineage_markers)) {
  seu <- compute_module_score(
    seu = seu, features = human_liver_lineage_markers[[m]], name = m)
}
feats <- c("lin_cycling", "lin_plasma", "lin_bcell", "lin_ilc", "lin_tcell", "lin_mp",
           "lin_mesenchyme", "lin_endo", "lin_cholangiocyte", "lin_hepatocyte")
feats_labels <- c("Cycling cells", "Plasma cells", "B cells", "ILCs", "T cells",
                  "MPs", "Mesenchyme", "Endothelia", "Cholangiocytes", "Hepatocytes")
seu$clusters <- seu$seurat_clusters
seu$clusters <- factor(
  seu$clusters,
  levels = c(0, 3, 7, 10, 18, 4, 17, 1, 12, 14, 6, 16, 2, 5, 11, 15, 13, 19, 8, 9),
  labels = c(0, 3, 7, 10, 18, 4, 17, 1, 12, 14, 6, 16, 2, 5, 11, 15, 13, 19, 8, 9))
Idents(object = seu) <- "clusters"

pdf(paste0(fig_dir, "Fig1_all_dotplot.pdf"), width = 9, height = 3.5)
print(dot_plot(seu = seu, features = feats, group.by = "clusters",
         labels = feats_labels, xlab = "Signature", ylab = "Cluster",
         legend.position = "right", col_pal = "OrRd",
         col.min = 0, col.max = 2, limits = c(0, 2), direction = 1))
dev.off()






###---------------------------------------------------------------------------
# NPC analysis
###---------------------------------------------------------------------------

all_lineages <- c("Mesenchyme", "Cholangiocytes", "Endothelia", "MPs")
all_npcs <- c(10, 60, 70, 50)
all_ndims <- c(10, 50, 60, 40)
all_res <- c(0.1, 0.3, 0.6, 0.2)

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
  png(paste0(fig_dir, "SFig2_", lineage, "_condition_split.png"),
      width = 6, height = 18, units = "in", res = 1000)
  print(subset_dim_plot(
    seu = seu, subset.by = "condition", reduction = "umap",
    ncol = 1, col_pal = condition_col_pal, pt.size = 3, pt.stroke = 0.04,
    back.pt.size = 0.5, back.pt.alpha = 0.2,
    back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes() & ggplot2::labs(title = NULL))
  dev.off()

  png(paste0(fig_dir, "SFig2_", lineage, "_condition_split_label.png"),
      width = 6, height = 19, units = "in", res = 1000)
  print(subset_dim_plot(
    seu = seu, subset.by = "condition", reduction = "umap",
    ncol = 1, col_pal = condition_col_pal, pt.size = 3, pt.stroke = 0.04,
    back.pt.size = 0.5, back.pt.alpha = 0.2,
    back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes())
  dev.off()


  # UMAP cluster
  png(paste0(fig_dir, "SFig2_", lineage, "_cluster.png"), width = 10, height = 8, units = "in", res = 1000)
  print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cp,
                         label = FALSE, label.size = 4, legend.position = "none",
                         pt.size = 3, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
  dev.off()


  png(paste0(fig_dir, "SFig2_", lineage, "_cluster_label.png"), width = 11, height = 8, units = "in", res = 1000)
  print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cp,
                         label = FALSE, label.size = 7, legend.position = "right",
                         pt.size = 3, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
  dev.off()


  # Clustering heatmap
  Idents(object = seu) <- "seurat_clusters"
  #sub <- subset(seu, downsample = 1500)
  local_heatmap_plot(seu = seu, markers = mark, topn_genes = 10,
                     diff_cluster_pct = opts$diff_cluster_pct,
                     pval_adj = 0.05, filename = paste0(fig_dir, "SFig2_", lineage, "_heatmap.png"),
                     fig.res = 600,
                     cluster_col_pal = cp,
                     condition_col_pal = condition_col_pal,
                     donor_col_pal = cp, is_clean_plot = TRUE)

  local_heatmap_plot(seu = seu, markers = mark, topn_genes = 10,
                     diff_cluster_pct = opts$diff_cluster_pct,
                     pval_adj = 0.05, filename = paste0(fig_dir, "SFig2_", lineage, "_heatmap_label.png"),
                     fig.res = 600,
                     cluster_col_pal = cp,
                     condition_col_pal = condition_col_pal,
                     donor_col_pal = cp, is_clean_plot = FALSE)
}




###---------------------------------------------------------------------------
# Healthy Hepatocytes
###---------------------------------------------------------------------------
lineage <- "Hepatocytes_Healthy"
npcs <- 40
pcs_to_remove <- c(1, 2)
npcs_dir <- npcs
if(!is.null(pcs_to_remove)) {
  npcs_dir <- paste0(npcs_dir, "_r", paste(pcs_to_remove, collapse = ""))
}
data_dir <- paste0(results_dir, lineage, "/01_integrate/00_iter2/")
obj <- readRDS(paste0(data_dir, "/seu_harmony_npcs", npcs_dir, ".rds"))
seu <- obj$seu
opts <- obj$opts
rm(obj)
gc(verbose = FALSE)
opts$npcs <- npcs
opts$ndims <- 30 # Top PC dimensions to perform UMAP and clustering
opts$res <- 0.2 # Clustering resolution
opts$diff_cluster_pct <- 0.1
opts$topn_genes <- 10
out_dir <- paste0(data_dir, "/seu_harmony_npcs", npcs_dir, "/dim", opts$ndims, "/")



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
  levels = c("PortalHep", "CentralHep", "MitoHighHep", "StressedHep"),
  labels = c("PortalHep", "CentralHep", "MitoHighHep", "StressedHep"),
  ordered = TRUE)
# Idents(object = seu) <- "anno"


## UMAP anno
png(paste0(fig_dir, "SFig3_Healthy_Hep_anno.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = cp,
                       label = FALSE, label.size = 4, legend.position = "none",
                       pt.size = 2.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "SFig3_Healthy_Hep_anno_label.png"), width = 12, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = cp,
                       label = FALSE, label.size = 4, legend.position = "right",
                       pt.size = 2.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
dev.off()

# QC
png(paste0(fig_dir, "SFig3_Healthy_Hep_qc_nfeature.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 1.6, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "SFig3_Healthy_Hep_qc_nfeature_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 1.6, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "SFig3_Healthy_Hep_qc_mito.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 1.6, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "SFig3_Healthy_Hep_qc_mito_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                           legend.position = "top", pt.size = 1.6, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + Seurat::NoAxes())
dev.off()


# Integration mixing of donors
pdf(paste0(fig_dir, "SFig3_Healthy_Hep_donor_mixing.pdf"), width = 10, height = 4.5)
cl_sample <- seu@meta.data |> group_by(anno, donor) |> summarise(n = n()) |> mutate(freq = n / sum(n))
print(ggplot(data = cl_sample, aes(x = anno, y = freq, fill = donor)) +
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
local_heatmap_plot(seu = sub, markers = mark, topn_genes = 10,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig3_Healthy_Hep_heatmap.png"),
                   fig.res = 600,
                   cluster_col_pal = cp,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = TRUE)

local_heatmap_plot(seu = sub, markers = mark, topn_genes = 10,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig3_Healthy_Hep_heatmap_label.png"),
                   fig.res = 600,
                   cluster_col_pal = cp,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = FALSE)


# Portal/central signature scores
feats <- list("hep_zonation_central" = "Central_Hep",
              "hep_zonation_portal" = "Portal_Hep")
for (m in names(feats)) {
  seu <- compute_module_score(
    seu = seu, features = hep_zonation_markers[[m]], name = feats[[m]])
}

png(paste0(fig_dir, "SFig3_Healthy_Hep_spata_zon_score.png"),
    width = 14, height = 6, units = "in", res = 1000)
p1 <- feature_plot_tailored(seu = seu, feature = "Portal_Hep", max.cutoff = "q98", reduction = "umap",
                            pt.size = 1.7, legend.position = "none", order_points_by_value = FALSE) & Seurat::NoAxes()
p2 <- feature_plot_tailored(seu = seu, feature = "Central_Hep", max.cutoff = "q98", reduction = "umap",
                            pt.size = 1.7, legend.position = "none", order_points_by_value = FALSE) & Seurat::NoAxes()
print(p1 | p2)
dev.off()

png(paste0(fig_dir, "SFig3_Healthy_Hep_spata_zon_score_legend.png"),
    width = 13, height = 6, units = "in", res = 1000)
p1 <- feature_plot_tailored(seu = seu, feature = "Portal_Hep", max.cutoff = "q98", reduction = "umap",
                            pt.size = 1.7, legend.position = "top", order_points_by_value = FALSE) & Seurat::NoAxes()
p2 <- feature_plot_tailored(seu = seu, feature = "Central_Hep", max.cutoff = "q98", reduction = "umap",
                            pt.size = 1.7, legend.position = "top", order_points_by_value = FALSE) & Seurat::NoAxes()
print(p1 | p2)
dev.off()



png(paste0(fig_dir, "SFig3_Healthy_Hep_spata_zon_score_order_and_cutoff.png"),
    width = 14, height = 6, units = "in", res = 1000)
p1 <- feature_plot_tailored(seu = seu, feature = "Portal_Hep", max.cutoff = "q98", min.cutoff = 0, reduction = "umap",
                            pt.size = 1.7, legend.position = "none", order_points_by_value = TRUE) & Seurat::NoAxes()
p2 <- feature_plot_tailored(seu = seu, feature = "Central_Hep", max.cutoff = "q98", min.cutoff = 0, reduction = "umap",
                            pt.size = 1.7, legend.position = "none", order_points_by_value = TRUE) & Seurat::NoAxes()
print(p1 | p2)
dev.off()

png(paste0(fig_dir, "SFig3_Healthy_Hep_spata_zon_score_legend_order_and_cutoff.png"),
    width = 13, height = 6, units = "in", res = 1000)
p1 <- feature_plot_tailored(seu = seu, feature = "Portal_Hep", max.cutoff = "q98", min.cutoff = 0, reduction = "umap",
                            pt.size = 1.7, legend.position = "top", order_points_by_value = TRUE) & Seurat::NoAxes()
p2 <- feature_plot_tailored(seu = seu, feature = "Central_Hep", max.cutoff = "q98", min.cutoff = 0, reduction = "umap",
                            pt.size = 1.7, legend.position = "top", order_points_by_value = TRUE) & Seurat::NoAxes()
print(p1 | p2)
dev.off()




###---------------------------------------------------------------------------
# Hepatocyte all conditions analysis
###---------------------------------------------------------------------------
lineage <- "Hepatocytes_downsample"
npcs <- 40
seed <- 12
pcs_to_remove <- c(1, 4)
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
opts$ndims <- 30 # Top PC dimensions to perform UMAP and clustering
opts$res <- 0.4 # Clustering resolution
opts$diff_cluster_pct <- 0.1
opts$topn_genes <- 50
out_dir <- paste0(data_dir, "/seu_harmony_seed",
                  seed, "_npcs", npcs_dir, "/dim", opts$ndims, "/")

## Colour palette for clusters
cluster_col_pal <- cp[c(8, 1, 7, 3, 2, 5, 6, 4, 32)]
anno_col_pal <- cluster_col_pal



# Add UMAP
seu <- SeuratPipe::add_umap_embedding(seu = seu, embedding = paste0(out_dir, "/umap_embedding.csv"))
# Marker genes per cluster
mark <- read.csv(file = paste0(out_dir, "/clusters/seu_markers_res", opts$res, ".csv"))
# Make cluster a factor from character
mark$cluster <- factor(mark$cluster)
mark$anno <- plyr::mapvalues(
  mark$cluster,
  from = seq(0, 8),
  to = c("DAH1", "Portal", "Migratory", "MitoHigh", "Central", "DAH2",
         "Cycling", "STAT1+", "FTH1+"))
migratory <- mark |> filter(anno == "Migratory") |>
  dplyr::arrange(-avg_log2FC) |>
  top_n(50, avg_log2FC) |> dplyr::pull(gene)
# saveRDS(object = migratory, file = paste0(anno_dir, "/hep_human_markers/human_hep_migratory_markers.rds"))


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
  to = c("DAH1", "Portal", "Migratory", "MitoHigh", "Central", "DAH2",
         "Cycling", "STAT1+", "FTH1+"))
seu$anno <- factor(
  seu$anno,
  levels = c("DAH1", "Portal", "Migratory", "MitoHigh", "Central", "DAH2",
             "Cycling", "STAT1+", "FTH1+"),
  labels = c("DAH1", "Portal", "Migratory", "MitoHigh", "Central", "DAH2",
             "Cycling", "STAT1+", "FTH1+"),
  ordered = TRUE)
seu$anno_ordered <- factor(
  seu$anno,
  levels = c("Central", "Portal", "MitoHigh", "DAH1", "DAH2", "Migratory",
             "Cycling", "STAT1+", "FTH1+"),
  labels = c("Central", "Portal", "MitoHigh", "DAH1", "DAH2", "Migratory",
             "Cycling", "STAT1+", "FTH1+"),
  ordered = TRUE)
# Idents(object = seu) <- "anno"


## UMAP anno
png(paste0(fig_dir, "Fig2_Hep_anno.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = anno_col_pal,
                       label = FALSE, label.size = 4, legend.position = "none",
                       pt.size = 2.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
dev.off()

png(paste0(fig_dir, "Fig2_Hep_anno_label.png"), width = 12, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "anno", col_pal = anno_col_pal,
                       label = FALSE, label.size = 4, legend.position = "right",
                       pt.size = 2.5, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + Seurat::NoAxes())
dev.off()


## UMAP spliut by condition
png(paste0(fig_dir, "SFig3_Hep_condition_split.png"), width = 19, height = 6, units = "in", res = 1000)
print(subset_dim_plot(
  seu = seu, subset.by = "condition", reduction = "umap",
  ncol = 3, col_pal = condition_col_pal, pt.size = 2.4, pt.stroke = 0.04,
  back.pt.size = 0.5, back.pt.alpha = 0.1,
  back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes() & ggplot2::labs(title = NULL))
dev.off()

png(paste0(fig_dir, "SFig3_Hep_condition_split_label.png"), width = 19, height = 6, units = "in", res = 1000)
print(subset_dim_plot(
  seu = seu, subset.by = "condition", reduction = "umap",
  ncol = 3, col_pal = condition_col_pal, pt.size = 2.4, pt.stroke = 0.04,
  back.pt.size = 0.5, back.pt.alpha = 0.1,
  back.pt.color = "grey", combine = TRUE) & Seurat::NoAxes())
dev.off()


# UMAP cluster
png(paste0(fig_dir, "SFig3_Hep_cluster.png"), width = 10, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cluster_col_pal,
                       label = FALSE, label.size = 4, legend.position = "none",
                       pt.size = 2.4, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig3_Hep_cluster_label.png"), width = 11, height = 8, units = "in", res = 1000)
print(dim_plot_tailored(seu = seu, group.by = "seurat_clusters", col_pal = cluster_col_pal,
                       label = TRUE, label.size = 7, legend.position = "right",
                       pt.size = 2.4, pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04) + NoAxes())
dev.off()


# QC
png(paste0(fig_dir, "SFig3_Hep_qc_nfeature.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 1.7, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig3_Hep_qc_nfeature_legend.png"), width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "nFeature_RNA", col_pal = "RdYlBu",
                           legend.position = "top", pt.size = 1.7, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig3_Hep_qc_mito.png"), width = 8, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                           legend.position = "none", pt.size = 1.7, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
dev.off()

png(paste0(fig_dir, "SFig3_Hep_qc_mito_legend.png"),
    width = 7, height = 6, units = "in", res = 1000)
print(feature_plot_tailored(seu = seu, feature = "percent.mito", col_pal = "RdYlBu",
                           legend.position = "top", pt.size = 1.7, pt.alpha = 1,
                           pt.shape = 21, pt.stroke = 0.04, order_points_by_value = FALSE) + NoAxes())
dev.off()



# Integration mixing of donors
cl_sample <- seu@meta.data |> group_by(anno_ordered, donor) |>
  summarise(n = n()) |> mutate(freq = n / sum(n))
pdf(paste0(fig_dir, "SFig3_Hep_donor_mixing_legend.pdf"), width = 10, height = 4.5)
print(ggplot(data = cl_sample, aes(x = anno_ordered, y = freq, fill = donor)) +
  geom_bar(stat = "identity", color="black", size = 0.05) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
  scale_fill_manual(values = cp))
dev.off()

pdf(paste0(fig_dir, "SFig3_Hep_donor_mixing.pdf"), width = 8, height = 4.5)
print(ggplot(data = cl_sample, aes(x = anno_ordered, y = freq, fill = donor)) +
        geom_bar(stat = "identity", color="black", size = 0.05) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
        xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
        scale_fill_manual(values = cp))
dev.off()


# Condition contribution per cluster
pdf(paste0(fig_dir, "Fig2_Hep_condition_rel_contribution.pdf"), width = 8, height = 3)
cl_condition <- seu@meta.data |> group_by(anno_ordered, condition) |>
  summarise(n = n()) |> mutate(freq = n / sum(n))
print(ggplot(data = cl_condition, aes(x = anno_ordered, y = freq, fill = condition)) +
  geom_bar(stat = "identity", color="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  xlab(NULL) + ylab("Relative contribution") + ggtitle(NULL) +
  scale_fill_manual(values = condition_col_pal))
dev.off()

pdf(paste0(fig_dir, "Fig2_Hep_condition_contribution.pdf"), width = 8, height = 3)
cl_condition <- seu@meta.data |> group_by(condition, anno_ordered) |>
  summarise(n = n()) |> mutate(freq = n / sum(n))
print(ggplot(data = cl_condition, aes(x = anno_ordered, y = freq, fill = condition)) +
  geom_bar(stat = "identity", color="black", position = position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  xlab(NULL) + ylab("Contribution") + ggtitle(NULL) +
  scale_fill_manual(values = condition_col_pal) +
  scale_y_continuous(limits=c(0, 0.40), breaks=seq(0, 0.40, 0.05)))
dev.off()




# Clustering heatmap
Idents(object = seu) <- "seurat_clusters"
local_heatmap_plot(seu = seu, markers = mark, topn_genes = 50,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig3_Hep_heatmap.png"),
                   fig.res = 300,
                   cluster_col_pal = cluster_col_pal,
                   condition_col_pal = condition_col_pal,
                   donor_col_pal = cp, is_clean_plot = TRUE)

local_heatmap_plot(seu = seu, markers = mark, topn_genes = 50,
                   diff_cluster_pct = opts$diff_cluster_pct,
                   pval_adj = 0.05, filename = paste0(fig_dir, "SFig3_Hep_heatmap_label.png"),
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
  pdf(paste0(fig_dir, "Fig2_Hep_GO_", cls$anno[cl], ".pdf"), width = 5, height = 10)
  print(plot_enrichment(df = df, terms = NULL, pval_thresh = 0.05, max_terms = 50,
                       text_size = 0.7, dot_size = 2.0, label_format = 45))
  dev.off()
}

cl <- 3
df <- read.csv(file = paste0(go_dir, "/Cluster", cls$seurat_clusters[cl], "_GO_BP_enrichment_level6.csv"))
terms <- c("actin filament organization", "regulation of cell morphogenesis", "ameboidal-type cell migration",
           "focal adhesion assembly", "regulation of cell-matrix adhesion",
           "regulation of actin cytoskeleton organization", "epithelial cell migration", "regulation of cell shape")
df <- df |> dplyr::filter(Description %in% terms)
pdf(paste0(fig_dir, "Fig2_Hep_GO_", cls$anno[cl], "_selected_terms.pdf"), width = 5, height = 2.6)
print(plot_enrichment(df = df, terms = NULL, pval_thresh = 0.05, max_terms = 50,
                      text_size = 1.1, dot_size = 3, label_format = 35))
dev.off()

# Feature plot module scores and marker genes
seu <- compute_module_score(seu = seu, features = hep_migratory_livreg_markers$hep_migratory,
                            name = "Migratory_hep")
seu <- compute_module_score(seu = seu, features = hep_zonation_markers$hep_zonation_central,
                            name = "Central_hep")
seu <- compute_module_score(seu = seu, features = hep_zonation_markers$hep_zonation_portal,
                            name = "Portal_hep")
seu <- compute_module_score(seu = seu, features = hep_manual_zonation_markers$hep_manual_zon_central,
                            name = "Central_hep_manual")
seu <- compute_module_score(seu = seu, features = hep_manual_zonation_markers$hep_manual_zon_portal,
                            name = "Portal_hep_manual")
seu <- compute_module_score(seu = seu, features = c("NUSAP1", "TOP2A", "SMC4", "TMPO", "NASP", "HELLS"),
                            name = "Cycling")
seu <- cell_cycle_score(seu, s_genes = cc_human_s, g2m_genes = cc_human_g2m)

features <- c("Migratory_hep", "Central_hep", "Portal_hep", "Central_hep_manual", "Portal_hep_manual",
              "Cycling", "S.Score", "G2M.Score", "ANXA2", "ASPH", "S100A10", "MET", "HGF", "SRC")
for (f in features) {
  png(paste0(fig_dir, "SFig3_Hep_condition_split_", f, "_score_label.png"),
      width = 19, height = 6, units = "in", res = 1000)
  print(subset_feature_plot(
    seu = seu, subset.by = "condition", feature = f, min.cutoff = 0,
    max.cutoff = "q98", reduction = "umap", slot = "data",
    ncol = NULL, col_pal = NULL, pt.size = 2, stroke = 0.04,
    legend.position = "top", back.pt.size = 0.7, back.pt.alpha = 0.1,
    back.pt.color = "grey", combine = TRUE, order_points_by_value = FALSE) &
      Seurat::NoAxes())
  dev.off()
}

for (f in features) {
  png(paste0(fig_dir, "SFig3_Hep_condition_split_", f, "_score.png"),
      width = 19, height = 6, units = "in", res = 1000)
  print(subset_feature_plot(
    seu = seu, subset.by = "condition", feature = f, min.cutoff = 0,
    max.cutoff = "q98", reduction = "umap", slot = "data",
    ncol = NULL, col_pal = NULL, pt.size = 2, stroke = 0.04,
    legend.position = "none", back.pt.size = 0.7, back.pt.alpha = 0.1,
    back.pt.color = "grey", combine = TRUE, order_points_by_value = FALSE) &
      ggplot2::ggtitle(NULL) & Seurat::NoAxes())
  dev.off()
}


for (f in features) {
  png(paste0(fig_dir, "SFig3_Hep_expr_", f, "_label.png"), width = 8, height = 6, units = "in", res = 1000)
  print(feature_plot_tailored(
    seu = seu, feature = f, col_pal = "RdYlBu", legend.position = "top", pt.size = 1.5,
    pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04, min.cutoff = 0, order_points_by_value = FALSE) + NoAxes())
  dev.off()

  png(paste0(fig_dir, "SFig3_Hep_expr_", f, "_label_ordered.png"), width = 8, height = 6, units = "in", res = 1000)
  print(feature_plot_tailored(
    seu = seu, feature = f, col_pal = "RdYlBu", legend.position = "top", pt.size = 1.5,
    pt.alpha = 1, pt.shape = 21, pt.stroke = 0.04, min.cutoff = 0, order_points_by_value = TRUE) + NoAxes())
  dev.off()
}
