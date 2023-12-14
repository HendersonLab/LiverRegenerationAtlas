library(Seurat)
library(SeuratPipe)
library(dplyr)

destination <- "/Analysis/Mouse_Nuclei_apap"

ml_apap <- readRDS(paste0(destination,"/seu_harmony_npcs80.rds"))
annotation <- read.csv(paste0(destination,"/seu_meta_res0.4_anno.csv"))

ml_apap_mps <- ml_apap$seu
ml_apap_mps$annotation <- annotation$anno
ml_apap_mps <- SetIdent(ml_apap_mps, value="annotation")
ml_apap_mps <- subset(ml_apap_mps, idents=c("MPs"))
rm(ml_apap)
rm(annotation)


# Perform data integration using Harmony
ml_apap_mps <- run_harmony_pipeline(
  seu_obj = ml_apap_mps,
  out_dir = "/Analysis/Mouse_Nuclei_apap/mps_raw/",
  batch_id = "sample",
  npcs = c(10,20,30,40,50,80),
  ndims = c(10,20,30,40,50,80),
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
  obj_filename = "ml_apap_mps_harmony",
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