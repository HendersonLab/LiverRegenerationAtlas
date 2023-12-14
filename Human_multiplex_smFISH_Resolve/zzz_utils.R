
###-------------------------------------------------------------
# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(nichenetr))

###-------------------------------------------------------------
# Global settings
if (grepl("machine-name",Sys.info()['nodename'])) {
  base_dir <- "livreg"
  raw_data_dir <- paste0(base_dir, "/raw/")
  results_dir <- paste0(base_dir, "/analysis/hlr_livreg/")
  analysis_dir <- "livreg/analysis/"
} else{
  stop("Computer not recognised")
}
# Source global utils helper functions
source(paste0(analysis_dir, "/zzz_seurat_utils.R"))


###-------------------------------------------------------------
# Loading marker genes
anno_dir <- paste0(base_dir, "/analysis/zzz_features/")


lineage_markers <- list(
  lin_endo = c("PECAM1", "RSPO3", "PLVAP", "STAB2", "LYVE1"),
  lin_chol = c("KRT19", "EPCAM", "GPM6A", "MUC5B"),
  lin_mps = c("PTPRC", "CD68", "CSF1R"),
  lin_kupffer = c("CD5L", "MARCO", "TIMD4"),
  lin_cc = c("MKI67", "PCNA", "CDK1"),
  lin_vsmc = c("RGS5", "MYH11"),
  lin_hsc = c("HGF", "RELN", "TNC"),
  lin_centralhep = c("CYP2E1", "CYP3A4", "ADH1B"),
  lin_portalhep = c("SDS", "HAL", "KRT18"),
  lin_hep = c("CYP2E1", "CYP3A4", "SDS", "HAL", "HNF4A"),
  lin_migratory = c("AKAP12", "FMNL2", "ITGA2", "CREB5", "ANXA2", "NDRG1", "PRAG1", "LUCAT1", "RELN")
)

saveRDS(object = lineage_markers_final,
        file = paste0(anno_dir, "/hep_human_markers/resolve_lineage_signatures.rds"))

##
# Human cell cycling genes
obj <- readRDS(file = paste0(anno_dir, "/cell_cycle_seurat_markers/cc_human_genes.rds"))
cc_human_s <- obj$cc_human_s
cc_human_g2m <- obj$cc_human_g2m


##
# Zonation markers obtained by SPATA analysis on Visium samples
##
hep_zonation <- read.csv(
  file = paste0(anno_dir, "/hep_human_markers/hep_spata_zonation_markers.csv"))
hep_zonation$zonation <- plyr::mapvalues(hep_zonation$zonation,
                                       from = c("Central", "Portal"),
                                       to = c("hep_zonation_central", "hep_zonation_portal"))
hep_zonation_markers <- list()
for (c in unique(hep_zonation$zonation)) {
  hep_zonation_markers[[c]] <- hep_zonation |> filter(zonation == c) |> dplyr::pull(variables)
}
rm(hep_zonation)


# Manually selected genes for Hep zonation
hep_manual_zonation_markers <-  list(hep_manual_zon_portal = c("HAL","SDS"),
                                     hep_manual_zon_central = c("CYP2E1", "CYP3A4"))


##
# Marker genes for the migratory hepatocyte phenotype identified in the LiverReg project
hep_migratory_livreg_markers <- read.csv(
  file = paste0(anno_dir, "/hep_human_markers/hep_livreg_markers.csv")) |>
  dplyr::filter(cluster == 2) |>
  #dplyr::filter(pct.1 - pct.2 > 0.1) |>
  dplyr::arrange(-avg_log2FC) |>
  top_n(25, avg_log2FC) |> dplyr::pull(gene)
hep_migratory_livreg_markers <- list(hep_migratory = hep_migratory_livreg_markers)



load_slide_seq <- function(coord, assay = 'Spatial') {
  slide.seq <- new(
    Class = 'SlideSeq',
    assay = assay,
    coordinates = coord
  )
  return(slide.seq)
}


run_resolve_qc_pipeline <- function(
    data_dir, sample_meta, sample_meta_filename = NULL, min.cells = 1, min.features = 1,
    meta_colnames = c("condition", "pass_qc"),
    out_dir = NULL, qc_to_plot = c("nFeature_Spatial", "nCount_Spatial"),
    alpha = c(0.1, 0.9), pt.size.factor = 1.1, max.cutoff = "q98", min.cutoff = NA,
    spatial_col_pal = "inferno", crop = TRUE, obj_filename = "seu_qc", stroke = 0.05, ...) {

  # Store all parameters for reproducibility
  opts <- c(as.list(environment()), list(...))

  # SO CMD passes without NOTES
  x = y = dens = plot_dir = seu <- NULL

  # Check and return if metadata object or filename is given as input
  opts$sample_meta <- SeuratPipe:::.get_metadata(sample_meta = sample_meta,
                                                 sample_meta_filename = sample_meta_filename)
  # Create plot directory if it doesn't exist
  if (!is.null(out_dir)) {
    plot_dir <- paste0(out_dir, "/plots/")
    if (!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  }

  # Create Seurat object
  seu <- resolve_to_seurat(
    data_dir = data_dir, sample_meta = opts$sample_meta,
    meta_colnames = meta_colnames, min.cells = min.cells, min.features = min.features)

  # In case we have one sample, we make the Seurat object a list
  if (!is.list(seu)) {
    # Extract sample name
    sample <- opts$sample_meta$sample
    seu <- list(seu)
    names(seu) <- sample
  }
  # Instantiate default assay by looking at DefaultAssay of first sample
  assay <- Seurat::DefaultAssay(object = seu[[1]])

  # Plot QC prior to filtering
  if (!is.null(plot_dir) && !is.null(qc_to_plot)) {
    plot_dim <- SeuratPipe:::.plot_dims(feat_len = length(qc_to_plot))
    pdf(paste0(plot_dir, "vln_preqc.pdf"), width = 10, height = 6,
        useDingbats = FALSE)
    for (s in names(seu)) {
      print(Seurat::VlnPlot(seu[[s]], features = qc_to_plot,
                            ncol = plot_dim$ncols, pt.size = 0) &
              ggplot2::geom_jitter(height = 0, size = 0.1, show.legend = FALSE, alpha = 0.1))
    }
    dev.off()

    spat_plot_dim <- SeuratPipe:::.spatial_plot_dims(feat_len = length(qc_to_plot))
    pdf(paste0(plot_dir, "spatial_preqc.pdf"), width = spat_plot_dim$width,
        height = spat_plot_dim$height, useDingbats = FALSE)
    for (s in names(seu)) {
      print(spatial_feature_plot(
        seu[[s]], features = qc_to_plot, alpha = alpha,
        pt.size.factor = pt.size.factor, ncol = spat_plot_dim$ncols,
        max.cutoff = max.cutoff, min.cutoff = min.cutoff, crop = crop,
        col_pal = spatial_col_pal, legend.position = "top", title = s, ...))
    }
    dev.off()



    # Get all combinations of QCs to plot
    combs <- utils::combn(qc_to_plot, 2)
    plot_dim <- SeuratPipe:::.plot_dims(feat_len = NCOL(combs))

    # Create scatter plots for feature-to-feature relationships
    pdf(paste0(plot_dir, "scatter_preqc.pdf"), width = plot_dim$width,
        height = plot_dim$height, useDingbats = FALSE)
    for (s in names(seu)) {
      gg_list <- scatter_meta_plot(seu = seu[[s]], features = qc_to_plot)
      plot(patchwork::wrap_plots(gg_list, ncol = plot_dim$ncols) +
             patchwork::plot_annotation(title = s, theme = ggplot2::theme(
               plot.title = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                  size = 20))))
    }
    dev.off()
  }

  # Pre-QC summary
  if (!is.null(plot_dir)) {
    preqc_summary <- data.frame(sample = character(), cells = numeric(),
                                median_nGenes = numeric(),
                                median_nCounts = numeric(),
                                prop_filt_nGenes = numeric())
    for (s in names(seu)) {
      tmp <- seu[[s]][[]]
      preqc_summary <- rbind(preqc_summary, data.frame(
        sample = s, cells = NROW(tmp),
        median_nGenes = round(median(tmp[[paste0("nFeature_", assay)]])),
        median_nCounts = round(median(tmp[[paste0("nCount_", assay)]])),
        prop_filt_nGenes = round(
          sum(tmp[[paste0("nFeature_", assay)]] <= min.features) / NROW(tmp), 2)
      ))
    }
    write.csv(preqc_summary, file = paste0(plot_dir, "preqc_sample_summary.csv"))
  }
#
#   # Perform actual QC filtering
#   seu <- qc_filter_seurat_object(seu = seu, nfeat_thresh = nfeat_thresh,
#                                  mito_thresh = mito_thresh)

  if (!is.null(plot_dir)) {
    # Plot QC after filtering
    if (!is.null(plot_dir) && !is.null(qc_to_plot)) {
      plot_dim <- SeuratPipe:::.plot_dims(feat_len = length(qc_to_plot))
      pdf(paste0(plot_dir, "vln_qc.pdf"), width = 10, height = 6,
          useDingbats = FALSE)
      for (s in names(seu)) {
        print(Seurat::VlnPlot(seu[[s]], features = qc_to_plot,
                              ncol = spat_plot_dim$ncols, pt.size = 0) &
                ggplot2::geom_jitter(height = 0, size = 0.1, show.legend = FALSE, alpha = 0.1))
      }
      dev.off()

      spat_plot_dim <- SeuratPipe:::.spatial_plot_dims(feat_len = length(qc_to_plot))
      pdf(paste0(plot_dir, "spatial_qc.pdf"), width = spat_plot_dim$width,
          height = spat_plot_dim$height, useDingbats = FALSE)
      for (s in names(seu)) {
        print(spatial_feature_plot(seu[[s]], features = qc_to_plot,
                                   alpha = alpha, pt.size.factor = pt.size.factor,
                                   ncol = spat_plot_dim$ncols, max.cutoff = max.cutoff,
                                   min.cutoff = min.cutoff,
                                   crop = crop, col_pal = spatial_col_pal,
                                   legend.position = "top", title = s, ...))
      }
      dev.off()
    }

    # Post-QC summary
    qc_summary <- data.frame(sample = character(), cells = numeric(),
                             median_nGenes = numeric(),
                             median_nCounts = numeric())
    for (s in names(seu)) {
      tmp <- seu[[s]][[]]
      qc_summary <- rbind(qc_summary, data.frame(
        sample = s, cells = NROW(tmp),
        median_nGenes = round(median(tmp[[paste0("nFeature_", assay)]])),
        median_nCounts = round(median(tmp[[paste0("nCount_", assay)]]))
      ))
    }
    write.csv(qc_summary, file = paste0(plot_dir, "qc_sample_summary.csv"))
  }

  # If we have single sample, remove it from list so we return Seurat object
  if (length(seu) == 1) { seu <- seu[[1]] }

  # Save Seurat object and opts
  if (!is.null(out_dir)) {
    if (obj_filename == "") { obj_filename <- "seu_qc.rds" }
    saveRDS(object = list(seu = seu, opts = opts),
            file = paste0(out_dir,"/",obj_filename,".rds"))
  }
  # Return QCed Seurat object
  return(seu)
}

# Function to convert Resolve output to Seurat object
resolve_to_seurat <- function(
    data_dir, sample_meta, meta_colnames = c("condition", "pass_qc"),
    min.cells = 5, min.features = 5) {
  seu <- list()
  for (i in 1:NROW(sample_meta)) {
    s <- sample_meta$sample[i]
    print(s)

    df <- read.delim(paste0(data_dir, sample_meta$path[i], "/count_mat.csv"), sep = ",") |>
      tibble::column_to_rownames("X") |> t()
    seu[[s]] <- Seurat::CreateSeuratObject(
      counts = df, min.cells = min.cells, min.features = min.features, project = s, assay = "Spatial")
    # Rename cells to distinguish across datasets
    #seu[[s]] <- Seurat::RenameCells(object = seu[[s]], add.cell.id = s)

    # Add required sample metadata information
    seu[[s]]$orig.ident <- as.factor(s)
    seu[[s]]$sample <- as.factor(s)
    seu[[s]]$condition <- as.factor(sample_meta$condition[i])
    Seurat::Idents(seu[[s]]) <- "sample"

    # Add user defined sample metadata information
    for (m in meta_colnames) {
      seu[[s]][[m]] <- as.factor(sample_meta[[m]][i])
    }

    # Add tissue coordinates as if they were from SlideSeq
    coord <- read.delim(paste0(data_dir, sample_meta$path[i], "/cellpos_mat.csv"), sep = ",") |>
      dplyr::mutate(cells = X)
    cell_ids <- data.frame(cells = colnames(seu[[s]]))
    coord <- dplyr::left_join(x = cell_ids, y = coord, by = "cells") |>
      tibble::column_to_rownames("X") |>
      dplyr::select(x, y, cells)
    seu[[s]]@images[[s]] <- load_slide_seq(coord = coord, assay = "Spatial")
    seu[[s]]@images[[s]]@key <- "image_"
  }

  # We have a single sample, so removing from list
  if (length(seu) == 1) { seu <- seu[[1]] }
  # Return generated Seurat object or list of Seurat objects
  return(seu)
}



rescale_data <- function(dim_dt, min.cutoff = NA, max.cutoff = NA,
                         scale_zero_one = TRUE, prob_quantile = 0.1) {
  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(dim_dt[, 3]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff
  )
  max.cutoff <- mapply(
    FUN = function(cutoff) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(dim_dt[, 3]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff
  )

  # Apply cutoffs
  data.feature <- as.vector(x = dim_dt[, 3])
  min.use <- Seurat::SetQuantile(cutoff = min.cutoff, data.feature)
  max.use <- Seurat::SetQuantile(cutoff = max.cutoff, data.feature)
  data.feature[data.feature < min.use] <- min.use
  data.feature[data.feature > max.use] <- max.use

  if (scale_zero_one) {
    data.feature <- scales::rescale(data.feature, to = c(0, 1))
  }
  dim_dt$y <- -dim_dt$y
  dim_dt[, 3] <- data.feature

  idx <- which(dim_dt[, 3] > quantile(dim_dt[,3], probs = prob_quantile))
  return(dim_dt[idx, ])
}



subset_spatial <- function(x, sample_id) {
  samples <- unique(x$sample)
  Idents(x) <- "sample"
  x <- subset(x, idents = sample_id)
  for (s in samples) {
    if (s != sample_id) {
      x@images[[s]] <- NULL
    }
  }
  return(x)
}


roi_to_polygon <- function(seu, roi, sample_id) {
  if (!is.character(sample_id)) {
    stop("Stopping 'sample' should be a character")
  }
  ids <- names(roi)
  foo <- lapply(names(roi), function(x) {
    coord <- as.data.frame(roi[[x]]$coords)
    coord$y <- -coord$y
    coord <- cbind(x, coord)
    colnames(coord) <- c("id", "x", "y")
    return(coord)
  })
  values <- data.frame(
    id = ids
  )
  obj <- merge(values, do.call("rbind", foo), by = c("id"))
  seu@misc$polygons[[sample_id]] <- obj
  return(seu)
}



smfish_spatial_plot <- function(
    seu, s, pops, cols, alpha = c(0.3, 0.9), min.cutoff = 0,
    max.cutoff = "q98", prob_quantile = 0.5, scale_zero_one = TRUE) {
  Idents(seu) <- "sample"
  seu_sub <- subset(seu, idents = s)
  coord <- seu_sub@images[[s]]@coordinates[, c("x", "y")]

  # Plot all cell segmentations
  all_cells <- cbind(seu_sub@misc$polygons[[s]], value = 1)
  gg <- ggplot(mapping = aes(x, y)) +
    geom_polygon(data = all_cells, colour = "grey40", fill = NA,
                 size = 0.05, aes(group = id))

  for (i in 1:length(pops)) {
    df <- data.frame(X = seu_sub[[pops[i]]][[1]])
    colnames(df) <- names(pops[i])
    df <- cbind(coord, df)
    df <- rescale_data(
      dim_dt = df, min.cutoff = min.cutoff, max.cutoff = max.cutoff,
      scale_zero_one = scale_zero_one, prob_quantile = prob_quantile) |>
      tibble::rownames_to_column("id") |> dplyr::select(c("id", names(pops[i])))
    df <- inner_join(x = seu_sub@misc$polygons[[s]], y = df, by = "id")

    min_lim <- min(df[, names(pops[i])], na.rm = TRUE)
    max_lim <- max(df[, names(pops[i])], na.rm = TRUE)

    gg <- gg + ggnewscale::new_scale_fill() +
      geom_polygon(data = df, aes_string(fill = names(pops[i]), alpha = names(pops[i]),
                                         group = "id")) +
      scale_fill_gradient(names(pops[i]), low = cols[i], high = cols[i],
                          limits = c(min_lim, max_lim))
  }
  gg <- gg + scale_alpha(range = alpha) + ggdark::dark_theme_gray() +
    theme(plot.background = element_rect(fill = "grey10"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 5),
          legend.key.size = ggplot2::unit(0.1,"line"),
          legend.key.height = ggplot2::unit(0.15, 'cm'),
          legend.key.width = ggplot2::unit(0.1, "cm")) + Seurat::NoAxes()
  return(gg)
}


# Preprocess function to extract info from Seurat object and
# create required input data for STM
seu2stm <- function(seu) {
  # Extract count matrix and transform to required format for stm
  X <- ceiling(seu@assays[[DefaultAssay(seu)]]@counts)
  if (is(X, "dgCMatrix")) {
    X <- as(X, "dgTMatrix")
  }
  # convert to data frame: convert to 1-based indexing
  docs <- data.frame(Feature = X@i + 1, Cell = X@j + 1, total_reads = X@x) |>
    dplyr::group_by(Cell) |>
    dplyr::group_split(.keep = FALSE)

  # Document structure as required by stm
  docs <- lapply(docs, FUN = function(x) {
    x <- t(x)
    rownames(x) <- NULL
    return(x)
  })
  # Vocabulary
  vocab <- rownames(X)
  # Metadata information
  meta <- seu@meta.data
  return(list(docs = docs, vocab = vocab, meta = meta))
}


run_stm <- function(seu, ntopics, topic_prevalence = 0.05, topic_names_prefix = "STM") {
  stm_input <- seu2stm(seu = seu)
  stm_input <- prepDocuments(stm_input$docs, stm_input$vocab,
                              stm_input$meta, lower.thresh = 5)
  # Fit the structural topic model
  fit <- stm::stm(documents = stm_input$documents, vocab = stm_input$vocab, K = ntopics,
                  max.em.its = 100, init.type = "Spectral", seed = 1)
  fit$meta <- stm_input$meta
  fit$s <- colSums(ceiling(seu@assays[[DefaultAssay(seu)]]@counts))
  fit$L <- fit$theta
  fit$Ln <- fit$L
  fit$F <- t(exp(fit$beta$logbeta[[1]]))
  fit$Fn <- fit$F

  # Filter topics
  idx <- which(colSums(fit$L) / sum(fit$L) > topic_prevalence)
  tmp <- fit$L[, idx]
  fit$L <- tmp / rowSums(tmp)
  colnames(fit$L) <- paste0(topic_names_prefix, 1:NCOL(fit$L))
  rownames(fit$L) <- colnames(seu)
  fit$F <- fit$F[, idx]
  colnames(fit$F) <- paste0(topic_names_prefix, 1:NCOL(fit$F))
  rownames(fit$F) <- stm_input$vocab

  # Add topic signatures as metadata
  seu <- Seurat::AddMetaData(seu, metadata = as.data.frame(fit$L))
  return(list(seu = seu, fit = fit))
}


run_fast_topic <- function(seu, ntopics, topic_prevalence = 0.05, topic_names_prefix = "FT") {
  # Extract count matrix and transpose it
  counts <- t(seu@assays[[DefaultAssay(seu)]]@counts)
  # Fit the topic model
  fit <- fastTopics::fit_topic_model(counts, k = ntopics)

  # Filter topics
  idx <- which(colSums(fit$L) / sum(fit$L) > topic_prevalence)
  tmp <- fit$L[, idx]
  fit$L <- tmp / rowSums(tmp)
  colnames(fit$L) <- paste0(topic_names_prefix, 1:NCOL(fit$L))
  fit$F <- fit$F[, idx]
  colnames(fit$F) <- paste0(topic_names_prefix, 1:NCOL(fit$F))

  # Add topic signatures as metadata
  seu <- Seurat::AddMetaData(seu, metadata = as.data.frame(fit$L))
  return(list(seu = seu, fit = fit))
}





extract_marker_genes <- function(fit, topn_genes = 15) {
  markers <- fit$F |> dplyr::as_tibble(rownames = "gene") |>
    tidyr::pivot_longer(cols = -(gene), names_to = "cluster", values_to = "prob") |>
    mutate(cluster = as.factor(cluster))
  # Filtering and arranging marker genes
  markers <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(n = topn_genes, order_by = prob) |>
    # Remove duplicates
    dplyr::group_by(gene) |>
    dplyr::arrange(prob, .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |> dplyr::ungroup() |>
    # Arrange by cluster and probability
    dplyr::group_by(cluster) |>
    dplyr::arrange(-prob, .by_group = TRUE) |>
    dplyr::pull(gene)
  df <- fit$F[which(rownames(fit$F) %in% markers), ]
  df <- df / rowSums(df)
  return(df)
}
