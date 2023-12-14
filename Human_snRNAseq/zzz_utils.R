
###-------------------------------------------------------------
# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(nichenetr))

###-------------------------------------------------------------
# Global settings
if (grepl("machine-name",Sys.info()['nodename'])) {
  base_dir <- "/liver/"
  raw_data_dir <- paste0(base_dir, "/raw/")
  results_dir <- paste0(base_dir, "/analysis/hln_livreg/")
  analysis_dir <- "/liver/analysis/"
} else{
  stop("Computer not recognised")
}
# Source global utils helper functions
source(paste0(analysis_dir, "/zzz_seurat_utils.R"))


###-------------------------------------------------------------
# Loading marker genes
anno_dir <- paste0(base_dir, "/analysis/zzz_features/")

# Store lineage markers
saveRDS(object = human_liver_lineage_markers,
        file = paste0(anno_dir, "/hep_human_markers/lineage_signatures.rds"))

##
# Human cell cycling genes
obj <- readRDS(file = paste0(anno_dir, "/cell_cycle_seurat_markers/cc_human_genes.rds"))
cc_human_s <- obj$cc_human_s
cc_human_g2m <- obj$cc_human_g2m
rm(obj)


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

##
# Zonation markers obtained by SPATA analysis on Mouse Visium samples
##
hep_mouse_zonation <- read.csv(
  file = paste0(anno_dir, "/hep_mouse_markers/hep_spata_zonation_markers.csv"))
hep_mouse_zonation$zonation <- plyr::mapvalues(hep_mouse_zonation$zonation,
                                               from = c("Central", "Portal"),
                                               to = c("hep_mouse_zonation_central", "hep_mouse_zonation_portal"))
hep_mouse_zonation_markers <- list()
for (c in unique(hep_mouse_zonation$zonation)) {
  genes <- hep_mouse_zonation |> filter(zonation == c) |> dplyr::pull(variables)
  to_human_symbol <- nichenetr::convert_mouse_to_human_symbols(genes)
  to_human_symbol <- to_human_symbol[!is.na(to_human_symbol)]
  hep_mouse_zonation_markers[[c]] <- to_human_symbol
}
rm(hep_mouse_zonation, genes, to_human_symbol)


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

