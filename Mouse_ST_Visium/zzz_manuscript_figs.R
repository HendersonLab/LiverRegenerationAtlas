
###-------------------------------------------------------------
# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(SPATA2))
suppressPackageStartupMessages(library(ggplot2))

# Source utils scripts
source("zzz_utils.R")
source("spata/zzz_spata_utils.R")

###---------------------------------------------------------------------------
# Global settings
###---------------------------------------------------------------------------

# Global colour pallete
cp <- SeuratPipe:::discrete_col_pal
condition_col_pal <- c("lightskyblue3", "indianred3")

# Cycling phase threshold
cycl_thresh <- 0.1



# Output directory for all figures
fig_dir <- paste0(results_dir, "/zzz_manuscript_figs/")
if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive = TRUE)}


## Load QCed Visium samples and normalize
qc_dir <- paste0(results_dir, "/00_qc/")
obj <- readRDS(paste0(qc_dir, "seu_qc.rds"))
seu <- obj$seu
rm(obj, qc_dir)
seu <- lapply(X = seu, FUN = function(x) {
  x <- SeuratPipe::lognormalize_and_pca(x, npcs = NULL)
  return(x)
})

samples <- names(seu)



###---------------------------------------------------------------------------
# Expression of certain marker genes and gene signatures for specific samples
###---------------------------------------------------------------------------
# Modules to perform analysis
modules <- list(Migratory = unname(hep_migratory_livreg_top25_markers$hep_migratory_top25),
                MyoFB = c("Acta2","Col1a1","Col3a1"),
                Hepatocytes = c("Ttr","Tf","Hp","Cyp2e1", "Hal"))
# Compute module scores
for (s in samples) {
  # Compute cycling score
  seu[[s]] <- cell_cycle_score(seu = seu[[s]], s_genes = cc_mouse_s,
                               g2m_genes = cc_mouse_g2m, thresh = cycl_thresh)
  for (m in names(modules)) {
    seu[[s]] <- compute_module_score(seu[[s]], features = modules[[m]], name = m)
    if (m == "Hepatocytes") {
      # Compute normalised values for Heps
      q99 <- quantile(seu[[s]]@meta.data$Hepatocytes, .99)
      seu[[s]]@meta.data$Hepatocytes_norm <- seu[[s]]@meta.data$Hepatocytes / q99
      # Ensure values are between (0, 1)
      seu[[s]]@meta.data$Hepatocytes_norm[which(seu[[s]]@meta.data$Hepatocytes_norm > 1)] <- 1
      seu[[s]]@meta.data$Hepatocytes_norm[which(seu[[s]]@meta.data$Hepatocytes_norm < 0)] <- 0
    }
  }
}

# Merge data to obtain common min and max threshold values across samples
seu_tmp <- seu[samples]
seu_tmp <- merge(x = seu_tmp[[1]], y = seu_tmp[2:length(seu_tmp)])

features <- c("Anxa2", "Migratory", "MyoFB")
probs <- list(Anxa2 = c(0.05, 0.98), Migratory = c(0.05, 0.98),
              MyoFB = c(0.05, 0.98), Hepatocytes = c(0.05, 0.98),
              Hepatocytes_norm = c(0.05, 0.98))

for (f in features) {
  qq <- quantile(Seurat::FetchData(object = seu_tmp, vars = f, slot = "data")[, 1],
                 probs = probs[[f]])
  if (f %in% c("Migratory", "MyoFB", "Hepatocytes", "Hepatocytes_norm")) {
    qq[1] <- 0
  }
  for (s in samples) {
    pdf(paste0(fig_dir, "SFig4_mouse_visium_score_", f, "_", s, ".pdf"), width = 5.5, height = 6)
    print(spatial_feature_plot(seu[[s]], features = f, max.cutoff = qq[2],
                               min.cutoff = qq[1], alpha = c(0, 0.9)) &
            viridis::scale_fill_viridis(option = "inferno", limits = qq))
    dev.off()
  }
}


###---------------------------------------------------------------------------
# SPATA zonation analysis
###---------------------------------------------------------------------------
spata_dir <- paste0(results_dir, "/spata/")
hm_colors <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(100)
sample <- "HEA_S4"
traj <- "zon1"
pattern <- list(Central = c("One peak (reversed)"),
                Portal = c("One peak"))

sp <- SPATA2::loadSpataObject(directory_spata = paste0(spata_dir, "/data/", sample, "/", sample, ".RDS"))
# ... load marker genes
top_traj <- read.csv(file = paste0(spata_dir, "/hea_gene_zonation/spata_zonation_genes.csv"))
top_traj$zonation <- factor(top_traj$zonation, levels = c("Central", "Portal"))

# Plot trajectories on tissue
pdf(paste0(fig_dir, "/SFig4_mouse_visium_spata_Cyp2e1.pdf"), width = 7.5, height = 8)
print(plotTrajectory(object = sp, trajectory_name = traj, color_by = "Cyp2e1", pt_alpha2 = 1,
                     pt_alpha = 0.25, pt_size = 1.75) + legendTop())
dev.off()

pdf(paste0(fig_dir, "/SFig4_mouse_visium_spata_Hal.pdf"), width = 7.5, height = 8)
print(plotTrajectory(object = sp, trajectory_name = traj, color_by = "Hal", pt_alpha2 = 1,
                     pt_alpha = 0.25, pt_size = 1.75) + legendTop())
dev.off()


# SPATA heatmap trajectory plot
pdf(paste0(fig_dir, "SFig4_mouse_visium_spata_zonation_heatmap_genes_label.pdf"), width = 7, height = 6)
print(plot_trajectory_heatmap(
  object = sp, trajectory_name = traj, variables = top_traj$variables,
  arrange_rows = "zonation", colors = rev(hm_colors), show_rownames = TRUE,
  split_columns = FALSE, smooth_span = 0.7, cluster_rows = FALSE,
  treeheight_row = 0, atdf = top_traj))
dev.off()

pdf(paste0(fig_dir, "SFig4_mouse_visium_spata_zon_heatmap_genes_label.pdf"), width = 7, height = 6)
print(plot_trajectory_heatmap(
  object = sp, trajectory_name = traj, variables = top_traj$variables,
  arrange_rows = "zonation", colors = rev(hm_colors), show_rownames = TRUE,
  split_columns = FALSE, smooth_span = 0.7, cluster_rows = FALSE,
  treeheight_row = 0, atdf = top_traj))
dev.off()

pdf(paste0(fig_dir, "SFig4_mouse_visium_spata_zon_heatmap_geneset.pdf"), width = 7, height = 7)
print(plot_trajectory_heatmap(
  object = sp, trajectory_name = traj, variables = top_traj$variables,
  arrange_rows = "zonation", colors = rev(hm_colors), show_rownames = FALSE,
  split_columns = FALSE, smooth_span = 0.7, cluster_rows = FALSE,
  treeheight_row = 0, atdf = top_traj))
dev.off()


##--------------------------------------------------------
# Zonation module scores applied on Visium slides
##--------------------------------------------------------

# Iterate over each sample
feats <- list("hep_zonation_central" = "Central",
              "hep_zonation_portal" = "Portal")
for (s in names(seu)) {
  for (i in names(feats)){
    seu[[s]] <- compute_module_score(seu[[s]], features = hep_zonation_markers[[i]], name = feats[[i]])
  }
  # Plot zonation signature
  pdf(paste0(fig_dir, "SFig4_mouse_visium_zon_central_score_", s, ".pdf"), width = 5.5, height = 6)
  print(spatial_feature_plot(seu[[s]], features = "Central",
                            max.cutoff = "q95", min.cutoff = "q05", alpha = c(0, 0.9)))
  dev.off()
  pdf(paste0(fig_dir, "SFig4_mouse_visium_zon_portal_score_", s, ".pdf"), width = 5.5, height = 6)
  print(spatial_feature_plot(seu[[s]], features = "Portal",
                            max.cutoff = "q95", min.cutoff = "q05", alpha = c(0, 0.9)))
  dev.off()
}
rm(s, i)


##--------------------------------------------------------
# Loss of zonation analysis
##--------------------------------------------------------
df <- list()
for (s in names(seu)) {
  for (i in names(feats)){
    seu[[s]] <- compute_module_score(seu[[s]], features = hep_zonation_markers[[i]], name = feats[[i]])
  }
  # Extract metadata information
  meta <- seu[[s]]@meta.data
  # Perform scaling to (0, 1)
  meta$Central <- scales::rescale(meta$Central, to = c(0,1))
  meta$Portal <- scales::rescale(meta$Portal, to = c(0,1))
  meta$Zonation <- ((meta$Central) / (meta$Central + meta$Portal + 0.01) )
  #apap_seu[[s]]@meta.data$Zonation <- meta$Zonation
  df[[s]] <- data.table::data.table(sample = s, value = meta$Zonation)
  x <- df[[s]]$value
  idx <- x < 1e-2
  if (sum(idx > 0)) {
    x[idx] <- 1e-2 + stats::runif(n = sum(idx), min = 0, max = 1e-4)
  }
  idx <- x > 1 - 1e-2
  if (sum(idx > 0)) {
    x[idx] <- 1 - 1e-2 - stats::runif(n = sum(idx), min = 0, max = 1e-4)
  }
  df[[s]]$value <- x
}
# Join all samples
df_plot <- data.table::rbindlist(df) |> as.data.frame()


for (s in unique(df_plot$sample)) {
  # Specific sample
  foo <- df_plot %>% dplyr::filter(sample == s)

  if (s %in% c("HEA_S1", "HEA_S2", "HEA_S3", "HEA_S4")) {
    col = "#492050"
  } else if (s %in% c("apap24H_S1", "apap24H_S2")) {
    col = "#E5CEEB"
  } else if (s %in% c("apap36H_S1", "apap36H_S2")) {
    col = "#E4EDE4"
  } else{
    col = "#96C597"
  }
  # Histogram plot
  pdf(paste0(fig_dir, "SFig4_mouse_visium_zon_specificity_", s, ".pdf"), width = 6, height = 2)
  print(ggplot(foo, aes(x = value)) +
         geom_histogram(aes(y=..count../sum(..count..)), bins = 50, fill = col,
                        alpha = 0.85, color="grey60") +
         ylim(c(0, 0.13)) + xlim(c(0,1)) + xlab(NULL) +
         ylab(NULL) + theme_classic())
  dev.off()
}

