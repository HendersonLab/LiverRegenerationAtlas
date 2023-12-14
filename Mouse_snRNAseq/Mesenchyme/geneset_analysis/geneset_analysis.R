library(dplyr)

here::here()
source("../../zzz_utils.R")

## Settings
quick_analysis <- TRUE
lineage_anno <- "Mesenchyme"
npcs <- 20 # Number of PCs
dims <- 10 # Top PC dimensions to perform UMAP and clustering
res <- 0.8 # Clustering resolution
pcs_to_remove <- NULL
npcs_dir <- npcs
if(!is.null(pcs_to_remove)) {
  npcs_dir <- paste0(npcs_dir, "_r", paste(pcs_to_remove, collapse = ""))
}

out_dir <- paste0(results_dir, lineage_anno, "/geneset_analysis/", "npcs", npcs_dir,
                  "_dim", dims, "_res", res, "/")

markers_file <- paste0(results_dir, lineage_anno, "/01_integrate/00_iter2/",
                       "/seu_harmony_npcs", npcs_dir, "/dim", dims,
                            "/clusters/seu_markers_res", res, ".csv")


tmp <- geneset_analysis(
  gene_markers = markers_file, out_dir = out_dir, species = "mouse",
  topn_genes = NULL, cluster_pval_adj = 0.05, quick_analysis = quick_analysis,
  ontologies = c("BP"), gs_pval_adj = 0.05, gs_qval = 0.1)

