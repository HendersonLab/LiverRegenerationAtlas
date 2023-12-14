
###-------------------------------------------------------------
# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(nichenetr))

###-------------------------------------------------------------
# Global settings
if (grepl("machine-name",Sys.info()['nodename'])) {
  base_dir <- "livreg/"
  raw_data_dir <- paste0(base_dir, "/raw/")
  results_dir <- paste0(base_dir, "/analysis/mlv_livreg/")
  analysis_dir <- "livreg_analysis/"
} else{
  stop("Computer not recognised")
}
# Source global utils helper functions
source(paste0(analysis_dir, "/zzz_seurat_utils.R"))


###-------------------------------------------------------------
# Loading marker genes
anno_dir <- paste0(base_dir, "/analysis/zzz_features/")


##
# Human cell cycling genes
obj <- readRDS(file = paste0(anno_dir, "/cell_cycle_seurat_markers/cc_mouse_genes.rds"))
cc_mouse_s <- obj$cc_mouse_s
cc_mouse_g2m <- obj$cc_mouse_g2m


##
# Zonation markers obtained by SPATA analysis on Human Visium samples
##
hep_zonation <- read.csv(
  file = paste0(anno_dir, "/hep_human_markers/hep_spata_zonation_markers.csv"))
hep_zonation$zonation <- plyr::mapvalues(hep_zonation$zonation,
                                         from = c("Central", "Portal"),
                                         to = c("hep_zonation_central", "hep_zonation_portal"))
hep_zonation_markers <- list()
for (c in unique(hep_zonation$zonation)) {
  genes <- hep_zonation |> filter(zonation == c) |> dplyr::pull(variables)
  to_mouse_symbol <- c(nichenetr::convert_human_to_mouse_symbols(genes, version = 2))
  to_mouse_symbol <- to_mouse_symbol[!is.na(to_mouse_symbol)]
  hep_zonation_markers[[c]] <- to_mouse_symbol
}
rm(hep_zonation, genes, to_mouse_symbol)

# Manually selected genes for Hep zonation
hep_manual_zonation_markers <-  list(hep_manual_zon_portal = c("Hal", "Sds", "Cyp2f2", "Pck1"),
                                     hep_manual_zon_central = c("Cyp2e1", "Cyp1a2", "Glul", "Gulo"))

##
# Marker genes for the migratory hepatocyte phenotype identified in the LiverReg project
hep_migratory_livreg_top25_markers <- read.csv(
  file = paste0(anno_dir, "/hep_human_markers/hep_livreg_markers.csv")) |>
  dplyr::filter(cluster == 2) |>
  #dplyr::filter(pct.1 - pct.2 > 0.1) |>
  dplyr::arrange(-avg_log2FC) |>
  top_n(25, avg_log2FC) |> dplyr::pull(gene)
to_mouse_symbol <- c(nichenetr::convert_human_to_mouse_symbols(hep_migratory_livreg_top25_markers, version = 2))
to_mouse_symbol <- to_mouse_symbol[!is.na(to_mouse_symbol)]
hep_migratory_livreg_top25_markers <- list(hep_migratory_top25 = to_mouse_symbol)
rm(to_mouse_symbol)
