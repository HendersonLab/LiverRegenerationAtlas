##-----------------------------------------------------------------------
# Robustness analysis of integration on downsampled Hepatocytes after
# cleansing. Any sample with > 800 cells will be downsampled.
##-----------------------------------------------------------------------

# Load packages
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
set.seed(12345)
here::here()

# Source utils scripts
source("../../../zzz_utils.R")

##-----------------------------------------------------------------------
# Load Harmony integrated data
##-----------------------------------------------------------------------
lineage_anno <- "Hepatocytes"
iter_npcs <- 100
iter_dims <- 90
iter_res <- 0.6
iter_dir <- paste0(results_dir, lineage_anno, "/01_integrate/00_iter1/")
obj <- readRDS(paste0(iter_dir, "/seu_harmony_npcs", iter_npcs, ".rds"))
seu <- obj$seu
opts <- obj$opts
rm(obj)
gc(verbose = FALSE)


##-----------------------------------------------------------------------
# Subset to retain cells from specified lineage
##-----------------------------------------------------------------------
opts$lineage_anno <- lineage_anno

io <- list() # Input/Output
io$out_dir <- paste0(results_dir, opts$lineage_anno,
                     "_downsample/01_integrate/robustness/")
if (!dir.exists(io$out_dir)) dir.create(io$out_dir, recursive = TRUE)

# Read cluster annotation
seu_anno <- read.csv(
  file = paste0(iter_dir, "seu_harmony_npcs", iter_npcs, "/dim", iter_dims,
                "/clusters/seu_meta_res", iter_res, "_anno.csv")) |>
  tibble::column_to_rownames(var = "X") |>
  select(anno_broad, anno) |>
  mutate_at(c("anno_broad", "anno"), factor)
# Add metadata information to current Seurat object and filter
seu <- AddMetaData(seu, metadata = seu_anno)
Idents(object = seu) <- "anno"
seu <- subset(seu, idents = c("MigratoryHep", "PortalHep", "CentralHep",
                              "MitoHighHep", "Hep", "Cycling", "APOHighHep"), invert = FALSE)


##-----------------------------------------------------------------------
# Additional steps to remove samples with low cell contribution
##-----------------------------------------------------------------------
# Filter samples with low number of cells. Set to NULL to ignore
opts$min_ncells_sample <- NULL
ncells <- seu@meta.data |> dplyr::count(sample)
colnames(ncells) <- c("Sample", opts$lineage_anno)
write.csv(x = ncells, file = paste0(io$out_dir, "ncells_per_sample.csv"))

if (!is.null(opts$min_ncells_sample)) {
  ncells <- ncells |> dplyr::filter(get(opts$lineage_anno) > opts$min_ncells_sample)
  write.csv(x = ncells, file = paste0(io$out_dir, "ncells_per_sample_retained.csv"))
  samples <- ncells$Sample
  Seurat::Idents(object = seu) <- "sample"
  seu <- subset(seu, idents = as.vector(samples), invert = FALSE)
}


##-----------------------------------------------------------------------
# Perform random subsamplings to show robustness of signature
##-----------------------------------------------------------------------
for (s in 50:54) {
  set.seed(s)

  # Samples with large number of cells
  samples <- seu@meta.data |> dplyr::count(sample) |>
    dplyr::filter(n > 800) |> dplyr::pull(sample)
  # Compute mean number of cells across conditions
  mean_ncells_cond <- seu@meta.data |> dplyr::count(sample, condition) |>
    dplyr::group_by(condition) |>
    dplyr::summarise(mean_ncells = round(mean(n)))
  # Number of cells for APAP condition
  ncells_apap <- mean_ncells_cond |> dplyr::filter(condition == "APAP") |>
    dplyr::pull(mean_ncells)
  # Perform downsampling on specified samples
  to_retain <- seu@meta.data |> tibble::rownames_to_column(var = "cellid") |>
    filter(sample %in% samples) |>
    group_by(sample) |> do(sample_n(., ncells_apap))
  # Concatenate downsampled samples and other samples
  to_retain <- rbind(to_retain, seu@meta.data |>
                       tibble::rownames_to_column(var = "cellid") |>
                       filter(!(sample %in% samples)))

  # Subset the Seurat object
  seu_filt <- subset(seu, cells = to_retain$cellid, invert = FALSE)

  # Store sample information post-downsampling
  ncells <- seu_filt@meta.data |> count(sample)
  colnames(ncells) <- c("Sample", opts$lineage_anno)
  write.csv(x = ncells, file = paste0(io$out_dir, "ncells_per_sample_downsampling.csv"))


  ##-----------------------------------------------------------------------
  # Additional Settings
  ##-----------------------------------------------------------------------

  ##-------------------
  # Opts

  # Number of PCs, can be a vector: c(30, 50)
  opts$npcs <- c(40)
  # Top PCs to perform UMAP and clustering, can be a vector: c(30, 40)
  opts$ndims <- c(30)
  # Clustering resolutions
  opts$res <- seq(from = 0.4, to = 0.5, by = 0.1)
  # Batch ID to perform integration
  opts$batch_id <- "sample"
  # QC information to plot
  opts$qc_to_plot <- c("nFeature_RNA", "percent.mito", "doublet_score")
  # Metadata columns to plot
  opts$meta_to_plot <- c("sample", "condition", "doublet_prediction")


  # Test genes that are detected in a minimum fraction of min.pct cells
  opts$min.pct <- 0.25
  # Test genes that show, on average, at least X-fold difference
  # between the two groups of cells.
  opts$logfc.threshold <- 0.6
  # Only return positive markers
  opts$only.pos <- TRUE
  # Maximum number of marker genes to plot for each cluster
  opts$topn_genes <- 50
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
  opts$pcs_to_remove <- c(1, 4)
  # Maximum cutoff values for plotting continuous features, e.g. gene expression
  # Gives better plots where colour scale is not driven by a (few) outlier cells.
  # Set to NULL to recoved default Seurat plots
  opts$max.cutoff <- "q98"
  # Filename of the integrated object to be stored on disk
  opts$obj_filename <- paste0("seu_harmony", "_rep", s)
  # If integrated Harmony file exists, should we force re-analysis
  # of Harmony, or read object? For computing time efficiency purposes.
  opts$force_reanalysis = FALSE


  # Number of highly variable genes to compute
  opts$n_hvgs <- 3000
  # Set specific seed for reproducibility
  opts$seed <- 1
  # Discrete colour palette
  opts$discrete_col_pal <- SeuratPipe:::discrete_col_pal
  # Whether to label the clusters in 'plot_reduction' space.
  opts$label <- TRUE
  # Sets size of labels.
  opts$label.size <- 7
  # Adjust point size for plotting.
  opts$pt.size <- 2.3
  # Figure resolution in ppi
  opts$fig.res = 200


  ##-------------------
  # Load modules
  # Group all modules in names list
  opts$modules_group <- list(hep_zonation = hep_zonation_markers,
                             hep_migratory = hep_migratory_livreg_markers)


  ##-----------------------------------------------------------------------
  # Analysis pipeline
  ##-----------------------------------------------------------------------

  # Perform data integration using Harmony
  seu_filt <- run_harmony_pipeline(
    seu_obj = seu_filt,
    out_dir = io$out_dir,
    batch_id = opts$batch_id,
    npcs = opts$npcs,
    ndims = opts$ndims,
    res = opts$res,
    modules_group = opts$modules_group,
    metadata_to_plot = opts$meta_to_plot,
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
    fig.res = opts$fig.res)
}
