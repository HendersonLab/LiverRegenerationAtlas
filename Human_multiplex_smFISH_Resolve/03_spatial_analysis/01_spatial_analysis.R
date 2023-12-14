
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(RImageJROI))
set.seed(12345)
here::here()

# Source utils scripts
source("../zzz_utils.R")


##-------------------------------------------------------
# Settings
##-------------------------------------------------------
##
# Load Resolve data
npcs <- 25
data_dir <- paste0(results_dir, "/02_integrate/")
obj <- readRDS(paste0(data_dir, "/seu_harmony_npcs", npcs, ".rds"))
seu <- obj$seu
opts <- obj$opts
rm(obj)
gc(verbose = FALSE)

plot_dir <- paste0(results_dir, "/03_spatial_analysis/")
dir.create(plot_dir, recursive = TRUE)

# Read Seurat ROI object
if (!file.exists(paste0(plot_dir, "/seu_harmony_npcs", npcs, ".rds"))) {
  # Add ROIs
  for (s in levels(seu$sample)){
    print(s)
    roi <- read.ijzip(paste0(raw_data_dir, "/resolve/", s, "/roiset.zip"))
    seu <- roi_to_polygon(seu = seu, roi = roi, sample_id = s)
  }
  saveRDS(object = list(seu = seu, opts = opts),
          file = paste0(plot_dir, "/seu_harmony_npcs", npcs, ".rds"))
} else {
  obj <- readRDS(paste0(plot_dir, "/seu_harmony_npcs", npcs, ".rds"))
  seu <- obj$seu
  opts <- obj$opts
  rm(obj)
  gc(verbose = FALSE)
}


##
# Additional settings
opts$npcs <- npcs
opts$dims <- 25 # Top PC dimensions to perform UMAP and clustering
opts$res <- 0.5 # Clustering resolution
samples <- levels(seu$sample)
s_width <- c(15, 13, 13, 13, 11, 13, 11, 7)
s_height <- c(9, 10, 10, 9, 12, 10, 12, 14)
names(s_width) = names(s_height) <- samples
col_pal <- SeuratPipe:::discrete_col_pal
# Output and plot directory
int_dir <- paste0(data_dir, "/seu_harmony_npcs", npcs, "/dim", opts$dims, "/")
plot_dir <- paste0(results_dir, "/03_spatial_analysis/")
dir.create(plot_dir, recursive = TRUE)

# Load modules
opts$modules_group <- list(lineage = lineage_markers,
                           hep_zonation = hep_zonation_markers,
                           hep_migratory = list(hep_migratory = lineage_markers$lin_migratory))

##-------------------------------------------------------
# Add metadata and lineage annotation
##-------------------------------------------------------
# Read UMAP and cluster information
seu <- add_umap_embedding(
  seu = seu, embedding = paste0(int_dir, "umap_embedding.csv"))
seu_cl <- read.csv(
  file = paste0(int_dir, "/clusters/seu_meta_res", opts$res, "_anno.csv")) |>
  tibble::column_to_rownames(var = "X") |>
  dplyr::select(seurat_clusters, anno)
# Add to Seurat metadata
seu <- AddMetaData(seu, metadata = seu_cl, col.name = c("seurat_clusters", "anno"))
seu$anno <- as.factor(seu$anno)
seu$seurat_clusters <- as.factor(seu$seurat_clusters)

##-------------------------------------------------------
# Defining Migratory Hepatocytes
# To show the Migratory signature only from Hepatocyte population,
# we obtain the labels of Hepatocytes and set the signature for the remaining NPCs to 0.
##-------------------------------------------------------
# Compute signatures
seu <- module_score_analysis(
  seu = seu, modules_group = opts$modules_group, ctrl = 10, nbin = 9)
seu$hep_migratory[seu@meta.data$anno != "Hepatocytes"] <- 0
seu$lin_cc[seu@meta.data$anno != "Hepatocytes"] <- 0


modules <- opts$modules_group$hep_zonation
features <- modules[[1]][modules[[1]] %in% rownames(seu)]
features <- modules[[2]][modules[[2]] %in% rownames(seu)]

modules <- opts$modules_group$hep_migratory
features <- modules[[1]][modules[[1]] %in% rownames(seu)]

##-------------------------------------------------------
# Loss of zonation analysis
##-------------------------------------------------------
df <- list()
modules <- list(Central = c("CYP2E1", "CYP3A4"),
                Portal = c("HAL", "SDS", "ASS1"))
for (s in samples) {
  seu_sub <- subset_spatial(x = seu, sample = s)
  for (i in names(modules)){
    seu_sub <- compute_module_score(seu_sub, features = modules[[i]], name = i, ctrl = 6, nbin = 6)
  }
  # Extract metadata information
  meta <- seu_sub@meta.data
  # Perform scaling to (0, 1)
  meta$Central <- scales::rescale(meta$Central, to = c(0,1))
  meta$Portal <- scales::rescale(meta$Portal, to = c(0,1))
  meta$Zonation <- ((meta$Central) / (meta$Central + meta$Portal + 0.01) )
  
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

df_plot$sample <- factor(df_plot$sample)

out_dir <- paste0(plot_dir, "/hep_zonation/")
if (!dir.exists(out_dir)) {dir.create(out_dir)}
pdf(file = paste0(out_dir, "/zonation_specificity.pdf"),
    width = 9, height = 6.5)
ggplot(df_plot, aes(x = value, fill = sample)) +
  geom_histogram(aes(y = ..density..), bins = 80, alpha = 0.6, color = "grey40") +
  scale_fill_manual(values = c("royalblue1", "royalblue4", "red1", "red3")) +
  facet_grid(sample ~ ., scales = "free") +
  xlab("Zonation specificity score") +
  ylab("Zonation specificity score") +
  theme_classic() +
  theme(legend.position = "none")
dev.off()


##-------------------------------------------------------
# Spatial analysis of hepatocyte zonation
##-------------------------------------------------------
pops <- c("hep_zonation_portal", "hep_zonation_central")
cols <- col_pal[c(1, 7)]
cols <- c("red1", "springgreen3")
names(pops) = names(cols) <- c("Portal", "Central")
out_dir <- paste0(plot_dir, "/hep_zonation/")
dir.create(out_dir, recursive = TRUE)
for (s in samples) {
  pdf(file = paste0(out_dir, "/spat_zonation_", s, ".pdf"),
      width = s_width[s], height = s_height[s])
  print(smfish_spatial_plot(
    seu = seu, s = s, pops = pops, cols = cols, alpha = c(0.3, 0.9),
    min.cutoff = 0.15, max.cutoff = "q98", scale_zero_one = TRUE, prob_quantile = 0.6)
  )
  dev.off()
}

###-------------------------------------------------
# Spatial analysis of major cell types
###-------------------------------------------------
pops <- c("hep_zonation_portal", "hep_zonation_central", "lin_endo", "lin_mps",
          "lin_chol", "lin_vsmc", "lin_hsc", "lin_kupffer")
cols <- col_pal[c(24, 24, 13, 31, 20, 33, 5, 34)]
names(pops) = names(cols) <- c("Hepatocytes", "Hepatocytes", "Endothelia", "MPs",
                               "Cholangiocytes", "VSMC", "HSC", "Kupffer")
out_dir <- paste0(plot_dir, "/spatial_cell_types/")
dir.create(out_dir, recursive = TRUE)
for (s in samples) {
  pdf(file = paste0(out_dir, "/spat_celltypes_", s, ".pdf"),
      width = s_width[s], height = s_height[s])
  print(smfish_spatial_plot(
    seu = seu, s = s, pops = pops, cols = cols, alpha = c(0.3, 0.9),
    min.cutoff = 0.15, max.cutoff = "q98", scale_zero_one = TRUE, prob_quantile = 0.6)
  )
  dev.off()
}


###-------------------------------------------------
# Spatial analysis of cell states
###-------------------------------------------------
pops <- c("hep_zonation_portal", "hep_zonation_central", "lin_mps",
          "lin_hsc", "hep_migratory", "lin_vsmc")
# cols <- col_pal[c(14, 14, 32, 13, 19, 33)]
cols <- c(col_pal[c(14, 14, 32)], "peachpuff1", col_pal[c(19, 33)])
names(pops) = names(cols) <- c("Hepatocytes", "Hepatocytes", "MPs",
                               "HSCs", "Migratory", "VSMC")
out_dir <- paste0(plot_dir, "/spatial_cell_states/")
dir.create(out_dir, recursive = TRUE)
for (s in samples) {
  pdf(file = paste0(out_dir, "/spat_migratory_", s, ".pdf"),
      width = s_width[s], height = s_height[s])
  print(smfish_spatial_plot(
    seu = seu, s = s, pops = pops, cols = cols, alpha = c(0.5, 0.9),
    min.cutoff = 0.35, max.cutoff = "q98", scale_zero_one = TRUE, prob_quantile = 0.6)
  )
  dev.off()
}



###-------------------------------------------------
# Spatial analysis of cycling cells
###-------------------------------------------------
pops <- c("hep_zonation_portal", "hep_zonation_central", "lin_mps",
          "lin_hsc", "lin_vsmc", "lin_cc")
cols <- c(col_pal[c(14, 14, 32)], "peachpuff1", col_pal[c(33, 35)])
# pops <- c("hep_zonation_portal", "hep_zonation_central", "lin_mps",
#           "hep_migratory", "lin_cc")
# cols <- col_pal[c(14, 14, 32, 19, 35)]
names(pops) = names(cols) <- c("Hepatocytes", "Hepatocytes", "MPs",
                               "HSCs", "VSMC", "Cycling")
out_dir <- paste0(plot_dir, "/spatial_cell_cycling/")
dir.create(out_dir, recursive = TRUE)
for (s in samples) {
  pdf(file = paste0(out_dir, "/spat_cycling_", s, ".pdf"),
      width = s_width[s], height = s_height[s])
  print(smfish_spatial_plot(
    seu = seu, s = s, pops = pops, cols = cols, alpha = c(0.5, 0.9),
    min.cutoff = 0.35, max.cutoff = "q98", scale_zero_one = TRUE, prob_quantile = 0.6)
  )
  dev.off()
}



