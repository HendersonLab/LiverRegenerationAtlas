##-----------------------------------------------------------------------
# Perform data integration on all cells after QC
##-----------------------------------------------------------------------

# Load packages
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(dplyr))
set.seed(12345)
here::here()

# Source utils scripts
source("../zzz_utils.R")

##-----------------------------------------------------------------------
# Load QCed data
##-----------------------------------------------------------------------
data_dir <- paste0(results_dir, "/00_qc/")
obj <- readRDS(paste0(data_dir, "/seu_qc.rds"))
seu <- obj$seu
opts <- obj$opts
rm(obj)
gc(verbose = FALSE)


##-----------------------------------------------------------------------
# Additional Settings
##-----------------------------------------------------------------------

##-------------------
# Input/Output
io <- list()
io$out_dir <- paste0(results_dir, "/01_integrate/")


##-------------------
# Opts

# Number of PCs, can be a vector: c(30, 50)
opts$npcs <- c(25)
# Top PCs to perform UMAP and clustering, can be a vector: c(30, 40)
opts$ndims <- c(25)
# Clustering resolutions
opts$res <- seq(from = 0.5, to = 0.7, by = 0.1)
# Batch ID to perform integration
opts$batch_id <- "sample"
# QC information to plot
opts$qc_to_plot <- c("nFeature_Spatial", "nCount_Spatial")
# Metadata columns to plot
opts$metadata_to_plot <- c("sample", "condition")


# Test genes that are detected in a minimum fraction of min.pct cells
opts$min.pct <- 0.05
# Test genes that show, on average, at least X-fold difference
# between the two groups of cells.
opts$logfc.threshold <- 0.2
# Only return positive markers
opts$only.pos <- TRUE
# Maximum number of marker genes to plot for each cluster
opts$topn_genes <- 20
# Retain marker genes per cluster if their
# `pct.1 - pct.2 > diff_cluster_pct`, i.e. they show cluster
# specific expression. Set to -Inf, to ignore this additional filtering.
opts$diff_cluster_pct <- 0.1
# Adjusted p-value threshold to consider marker genes per cluster.
opts$pval_adj <- 0.05
# Should we create feature plots of cluster markers?
opts$plot_cluster_markers = FALSE
# Which PCs to remove prior to running harmony (technical effects)
# Can be a vector e.g. c(1, 2, 4). If NULL, all PCs are used as input.
opts$pcs_to_remove <- NULL
# Maximum cutoff values for plotting continuous features, e.g. gene expression
# Gives better plots where colour scale is not driven by a (few) outlier cells.
# Set to NULL to recoved default Seurat plots
opts$max.cutoff <- "q98"
# Filename of the integrated object to be stored on disk
opts$obj_filename <- "seu_harmony"
# If integrated Harmony file exists, should we force re-analysis
# of Harmony, or read object? For computing time efficiency purposes.
opts$force_reanalysis = FALSE


# Scaling factor for normalisation
opts$scale.factor <- 100


# Number of highly variable genes to compute
opts$n_hvgs <- 100
# Set specific seed for reproducibility
opts$seed <- 1
# Discrete colour palette
opts$discrete_col_pal <- SeuratPipe:::discrete_col_pal
# Whether to label the clusters in 'plot_reduction' space.
opts$label <- TRUE
# Sets size of labels.
opts$label.size <- 6
# Adjust point size for plotting.
opts$pt.size <- 1.7
# Figure resolution in ppi
opts$fig.res = 200
# Crop image?
opts$crop <- TRUE
opts$cont_alpha <- c(0.1, 1)
opts$pt.size.factor <- 1.2
opts$plot_spatial_markers <- FALSE


##-------------------
# Load modules
# Group all modules in names list
opts$modules_group <- list(lineage = lineage_markers,
                           hep_zonation = hep_zonation_markers,
                           hep_migratory = hep_migratory_livreg_markers)


##-----------------------------------------------------------------------
# Analysis pipeline
##-----------------------------------------------------------------------

# Perform data integration using Harmony
seu <- run_harmony_pipeline(
  seu_obj = seu,
  out_dir = io$out_dir,
  batch_id = opts$batch_id,
  npcs = opts$npcs,
  ndims = opts$ndims,
  res = opts$res,
  modules_group = opts$modules_group,
  metadata_to_plot = opts$metadata_to_plot,
  qc_to_plot = opts$qc_to_plot,
  logfc.threshold = opts$logfc.threshold,
  min.pct = opts$min.pct,
  only.pos = opts$only.pos,
  topn_genes = opts$topn_genes,
  diff_cluster_pct = opts$diff_cluster_pct,
  pval_adj = opts$pval_adj,
  pcs_to_remove = opts$pcs_to_remove,
  obj_filename = opts$obj_filename,
  force_reanalysis = opts$force_reanalysis,
  plot_cluster_markers = opts$plot_cluster_markers,
  max.cutoff = opts$max.cutoff,
  n_hvgs = opts$n_hvgs,
  seed = opts$seed,
  discrete_col_pal = opts$discrete_col_pal,
  label = opts$label,
  label.size = opts$label.size,
  pt.size = opts$pt.size,
  fig.res = opts$fig.res,
  scale.factor = opts$scale.factor,
  crop = opts$crop,
  cont_alpha = opts$cont_alpha,
  pt.size.factor = opts$pt.size.factor,
  plot_spatial_markers = opts$plot_spatial_markers,
  ctrl = 6, nbin = 6,
  heatmap_downsample_cols = 600)
