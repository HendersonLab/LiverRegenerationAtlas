library(Seurat)
library(SeuratPipe)
library(dplyr)
library(tibble)

destination <- "/Analysis/Mouse_Nuclei_apap"

# Read pre-cleansed data
iter_npcs <- 80
iter_dims <- 40
iter_res <- 0.4
obj <- readRDS(paste0(destination, "/mesenchyme_raw/ml_apap_mes_harmony_npcs", iter_npcs, ".rds"))
ml_apap_mes <- obj$seu
opts <- obj$opts
rm(obj)

# Read cluster annotation
ml_apap_mes_anno <- read.csv(
  file = paste0(destination, "/mesenchyme_raw/ml_apap_mes_harmony_npcs", iter_npcs, "/dim", iter_dims,
                "/clusters/seu_meta_res", iter_res, ".csv")) %>%
  tibble::column_to_rownames(var = "X") %>%
  select(seurat_clusters) %>%
  mutate_at(c("seurat_clusters"), factor)

# Add metadata information to current Seurat object and filter
ml_apap_mes <- AddMetaData(ml_apap_mes, metadata = ml_apap_mes_anno)
Idents(object = ml_apap_mes) <- "seurat_clusters"
ml_apap_mes <- subset(ml_apap_mes, idents = c("1","3","5","6","8"), invert = TRUE)


# Perform data integration using Harmony
ml_apap_mes <- run_harmony_pipeline(
  seu_obj = ml_apap_mes,
  out_dir = "/Analysis/Mouse_Nuclei_apap/mesenchyme/",
  batch_id = "sample",
  npcs = c(10,20,30,40,50),
  ndims = c(5,10,20,30,40,50),
  res = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0,2.0),
  modules_group = list(),
  metadata_to_plot = c("sample", "condition", "doublet_prediction", "timepoint"),
  qc_to_plot = c("nFeature_RNA", "percent.mito", "doublet_score"),
  logfc.threshold = 0.7,
  min.pct = 0.25,
  only.pos = TRUE,
  topn_genes = 10,
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
