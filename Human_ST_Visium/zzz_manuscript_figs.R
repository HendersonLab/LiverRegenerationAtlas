
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
# Samples to perform analysis
samples <- c("xx1", "xx2")

# Output directory for all figures
fig_dir <- paste0(results_dir, "/zzz_manuscript_figs/")
if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive = TRUE)}


## Load QCed Visium samples and normalize
qc_dir <- paste0(results_dir, "/00_qc/")
obj <- readRDS(paste0(qc_dir, "seu_qc.rds"))
seu <- obj$seu[samples]
rm(obj, qc_dir)
seu <- lapply(X = seu, FUN = function(x) {
  x <- SeuratPipe::lognormalize_and_pca(x, npcs = NULL)
  return(x)
})



###---------------------------------------------------------------------------
# Expression of certain marker genes and gene signatures for specific samples
###---------------------------------------------------------------------------
# Modules to perform analysis
modules <- list(Migratory = hep_migratory_livreg_top25_markers$hep_migratory_top25,
                MyoFB = c("ACTA2","COL1A1","COL1A2","COL3A1"),
                Hepatocytes = c("TTR","TF","HP","CYP2A6","CYP2E1","CYP3A4", "HAL"))
# Compute module scores
for (s in samples) {
  # Compute cycling score
  seu[[s]] <- cell_cycle_score(seu = seu[[s]], s_genes = cc_human_s,
                               g2m_genes = cc_human_g2m, thresh = cycl_thresh)
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

features <- c("Migratory", "MyoFB", "Hepatocytes", "Hepatocytes_norm")
probs <- list(Migratory = c(0.05, 0.98),
              MyoFB = c(0.05, 0.98), Hepatocytes = c(0.05, 0.98),
              Hepatocytes_norm = c(0.05, 0.98))

for (f in features) {
  qq <- quantile(Seurat::FetchData(object = seu_tmp, vars = f, slot = "data")[, 1],
                 probs = probs[[f]])
  qq[1] <- 0
  for (s in samples) {
    pdf(paste0(fig_dir, "SFig3_visium_score_", f, "_", s, ".pdf"), width = 5.5, height = 6)
    print(spatial_feature_plot(seu[[s]], features = f, max.cutoff = qq[2],
                              min.cutoff = qq[1], alpha = c(0, 0.9)) &
            viridis::scale_fill_viridis(option = "inferno", limits = qq))
    dev.off()
  }
}


## Hepatocyte enriched spots cycling analysis
for (s in samples) {
  idx <- seu[[s]]@meta.data$MyoFB < 0.4 & seu[[s]]@meta.data$Hepatocytes_norm > 0.45
  x <- subset(seu[[s]], cells = colnames(seu[[s]])[idx])
  Idents(x) <- x@meta.data$sample

  phase <- c("SG2Mmax.Score")
  for (p in phase){
    pdf(paste0(fig_dir, "SFig1_visium_cc_score_SG2Mmax_", s, ".pdf"), width = 5.5, height = 6)
    print(spatial_feature_plot(x, features = "SG2Mmax.Score",
                              max.cutoff = 0.1, min.cutoff = 0, alpha = c(0, 0.9)) &
            viridis::scale_fill_viridis(option = "inferno", limits = c(0, 0.1)))
    dev.off()
  }
}


# SPATA apap analysis
###---------------------------------------------------------------------------
spata_dir <- paste0(results_dir, "/spata/")
hm_colors <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(100)
sample <- "apapxx1"
traj <- "vn1"

sp <- SPATA2::loadSpataObject(directory_spata = paste0(spata_dir, "/data/", sample, "/", sample, ".RDS"))
# ... load marker genes
top_traj <- read.csv(file = paste0(spata_dir, "/apap/spata_apap_geneset_", sample, ".csv"))
top_traj$region <- factor(top_traj$region, levels = c("Viable", "Perinecrotic", "Necrotic"))

# Plot trajectories on tissue
pdf(paste0(fig_dir, "/Fig1_visium_spata_COL3A1.pdf"), width = 7.5, height = 8)
print(plotTrajectory(object = sp, trajectory_name = traj, color_by = "COL3A1", pt_alpha2 = 1,
               pt_alpha = 0.25, pt_size = 1.75) + legendTop())
dev.off()

# SPATA heatmap trajectory plot
pdf(paste0(fig_dir, "Fig1_visium_spata_apap_heatmap_geneset_label.pdf"), width = 15.5, height = 20)
print(plot_trajectory_heatmap(
  object = sp, trajectory_name = traj, variables = top_traj$variables,
  arrange_rows = "region", colors = rev(hm_colors), show_rownames = TRUE,
  split_columns = FALSE, smooth_span = 0.7, cluster_rows = FALSE,
  treeheight_row = 0, atdf = top_traj))
dev.off()

pdf(paste0(fig_dir, "Fig1_visium_spata_apap_heatmap_geneset.pdf"), width = 5, height = 8)
print(plot_trajectory_heatmap(
  object = sp, trajectory_name = traj, variables = top_traj$variables,
  arrange_rows = "region", colors = rev(hm_colors), show_rownames = FALSE,
  split_columns = FALSE, smooth_span = 0.7, cluster_rows = FALSE,
  treeheight_row = 0, atdf = top_traj))
dev.off()



###---------------------------------------------------------------------------
# SPATA zonation analysis
###---------------------------------------------------------------------------
spata_dir <- paste0(results_dir, "/spata/")
hm_colors <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(100)
sample <- "heaxx1"
traj <- "zon1"
pattern <- list(Central = c("One peak (reversed)"), Portal = c("One peak"))

sp <- SPATA2::loadSpataObject(directory_spata = paste0(spata_dir, "/data/", sample, "/", sample, ".RDS"))

# ... load spatial score trends
df <- read.csv(file = paste0(spata_dir, "/data/", sample, "/", sample , "_geneset_", traj, ".csv"))
top_traj <- NULL
for (p in names(pattern)) {
  thresh <- ifelse(p == "Portal", 2.5, 2)
  res <- filter_trajectory_trends(
    df = df, auc_thresh = thresh, pvalue_thresh = 1,
    wald_stat_thresh = 0,
    topn = NULL, trends = pattern[[p]], variables_only = FALSE) |>
    dplyr::mutate(zonation = p)
  if (is.null(top_traj)) {
    top_traj <- res
  } else {
    top_traj <- rbind(top_traj, res)
  }
}

# Plot trajectories on tissue
pdf(paste0(fig_dir, "Fig1_visium_spata_CYP3A4.pdf"), width = 5.5, height = 6)
print(plotTrajectory(object = sp, trajectory_name = traj, color_by = "CYP3A4", pt_alpha2 = 1,
                    pt_alpha = 0.25) + legendTop())
dev.off()
pdf(paste0(fig_dir, "Fig1_visium_spata_HAL.pdf"), width = 5.5, height = 6)
print(plotTrajectory(object = sp, trajectory_name = traj, color_by = "HAL", pt_alpha2 = 1,
               pt_alpha = 0.25) + legendTop())
dev.off()


# SPATA heatmap trajectory plot
pdf(paste0(fig_dir, "Fig1_visium_spata_zon_heatmap_geneset_label.pdf"), width = 22.5, height = 30)
print(plot_trajectory_heatmap(
  object = sp, trajectory_name = traj, variables = top_traj$variables,
  arrange_rows = "zonation", colors = rev(hm_colors), show_rownames = TRUE,
  split_columns = TRUE, smooth_span = 0.7, cluster_rows = FALSE,
  treeheight_row = 0, atdf = top_traj))
dev.off()

pdf(paste0(fig_dir, "Fig1_visium_spata_zon_heatmap_geneset.pdf"), width = 8, height = 7)
print(plot_trajectory_heatmap(
  object = sp, trajectory_name = traj, variables = top_traj$variables,
  arrange_rows = "zonation", colors = rev(hm_colors), show_rownames = FALSE,
  split_columns = TRUE, smooth_span = 0.7, cluster_rows = FALSE,
  treeheight_row = 0, atdf = top_traj))
dev.off()


pdf(paste0(fig_dir, "Fig1_visium_spata_hep_central_genes.pdf"), width = 5, height = 5)
print(plotTrajectoryGenes(
  object = sp, trajectory_name = traj, clrp = NULL, clrp_adjust = "black",
  genes = c("CYP2E1", "ADH1B", "ADH1A", "ALAS1", "CXCL2"), binwidth = 5, smooth_method = "loess",
  smooth_span = 0.8, smooth_se = FALSE, display_facets = FALSE, nrow = 5) + legendTop())
dev.off()


pdf(paste0(fig_dir, "Fig1_visium_spata_hep_portal_genes.pdf"), width = 5, height = 5)
print(plotTrajectoryGenes(
  object = sp, trajectory_name = traj, clrp = NULL, clrp_adjust = "black",
  genes = c("SDS", "A2M", "CRP", "CYP2A6", "FGA"), binwidth = 5, smooth_method = "loess",
  smooth_span = 0.8, smooth_se = FALSE, display_facets = FALSE, nrow = 5) + legendTop())
dev.off()


##--------------------------------------------------------
# Zonation module scores applied on Visium slides
##--------------------------------------------------------

# Iterate over each sample
feats <- list("hep_zonation_central" = "Central",
              "hep_zonation_portal" = "Portal")
samples <- c("headxx1", "apapxx2")
for (s in samples) {
  for (i in names(feats)){
    seu[[s]] <- compute_module_score(seu[[s]], features = hep_zonation_markers[[i]], name = feats[[i]])
  }
  # Plot zonation signature
  pdf(paste0(fig_dir, "Fig1_visium_zon_central_score_", s, ".pdf"), width = 5.5, height = 6)
  print(spatial_feature_plot(seu[[s]], features = "Central",
                            max.cutoff = "q95", min.cutoff = "q05", alpha = c(0, 0.9)))
  dev.off()
  pdf(paste0(fig_dir, "Fig1_visium_zon_portal_score_", s, ".pdf"), width = 5.5, height = 6)
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
df_plot <- data.table::rbindlist(df) |> as.data.frame() |>
  dplyr::filter(sample %in% c("heaxx1", "apapxx1"))

pdf(paste0(fig_dir, "Fig1_visium_zon_specificity.pdf"), width = 5, height = 4)
print(ggplot(df_plot, aes(x = value, fill = sample, colour = sample)) +
  geom_histogram(aes(y=..density..), bins = 50, alpha = 0.4, position = "identity") + #..count../sum(..count..))
  scale_color_manual(values = condition_col_pal) +
  scale_fill_manual(values = condition_col_pal) +
  xlab("Zonation specificity score") +
  ylab(NULL) + theme_classic() + theme(legend.position = c(0.85, 0.88)))
dev.off()

