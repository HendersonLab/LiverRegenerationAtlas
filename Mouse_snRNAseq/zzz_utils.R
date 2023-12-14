
###-------------------------------------------------------------
# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(nichenetr))

###-------------------------------------------------------------
# Global settings
if (grepl("machine-name",Sys.info()['nodename'])) {
  base_dir <- "livreg/"
  raw_data_dir <- paste0(base_dir, "/raw/")
  results_dir <- paste0(base_dir, "/analysis/mln_livreg/")
  analysis_dir <- "livreg_analysis/"
} else{
  stop("Computer not recognised")
}
# Source global utils helper functions
source(paste0(analysis_dir, "/zzz_seurat_utils.R"))


###-------------------------------------------------------------
# Loading marker genes
anno_dir <- paste0(base_dir, "/analysis/zzz_features/")

# Store lineage markers
saveRDS(object = mouse_liver_lineage_markers,
        file = paste0(anno_dir, "/hep_mouse_markers/lineage_signatures.rds"))

##
# Human cell cycling genes
obj <- readRDS(file = paste0(anno_dir, "/cell_cycle_seurat_markers/cc_mouse_genes.rds"))
cc_mouse_s <- obj$cc_mouse_s
cc_mouse_g2m <- obj$cc_mouse_g2m


##
# Zonation markers obtained by SPATA analysis on Human Visium samples
##
hep_human_zonation <- read.csv(
  file = paste0(anno_dir, "/hep_human_markers/hep_spata_zonation_markers.csv"))
hep_human_zonation$zonation <- plyr::mapvalues(hep_human_zonation$zonation,
                                               from = c("Central", "Portal"),
                                               to = c("hep_human_zonation_central", "hep_human_zonation_portal"))
hep_human_zonation_markers <- list()
for (c in unique(hep_human_zonation$zonation)) {
  genes <- hep_human_zonation |> filter(zonation == c) |> dplyr::pull(variables)
  to_mouse_symbol <- c(nichenetr::convert_human_to_mouse_symbols(genes, version = 2))
  to_mouse_symbol <- to_mouse_symbol[!is.na(to_mouse_symbol)]
  hep_human_zonation_markers[[c]] <- to_mouse_symbol
}
rm(hep_human_zonation, genes, to_mouse_symbol)

##
# Zonation markers obtained by SPATA analysis on Human Visium samples
##
hep_mouse_zonation <- read.csv(
  file = paste0(anno_dir, "/hep_mouse_markers/hep_spata_zonation_markers.csv"))
hep_mouse_zonation$zonation <- plyr::mapvalues(hep_mouse_zonation$zonation,
                                               from = c("Central", "Portal"),
                                               to = c("hep_zonation_central", "hep_zonation_portal"))
hep_mouse_zonation_markers <- list()
for (c in unique(hep_mouse_zonation$zonation)) {
  hep_mouse_zonation_markers[[c]] <- hep_mouse_zonation |> filter(zonation == c) |> dplyr::pull(variables)
}
rm(hep_mouse_zonation)


# Manually selected genes for Hep zonation
hep_manual_zonation_markers <-  list(hep_manual_zon_central = c("Cyp2e1", "Cyp1a2", "Glul", "Gulo"),
                                     hep_manual_zon_portal = c("Hal", "Sds", "Cyp2f2", "Pck1"))



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

