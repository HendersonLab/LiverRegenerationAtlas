##-----------------------------------------------------------------------
# QC analysis pipeline
##-----------------------------------------------------------------------

# Load packages
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(dplyr))
set.seed(12345)
here::here()

# Source utils scripts
source("../zzz_utils.R")

##-----------------------------------------------------------------------
# Settings
##-----------------------------------------------------------------------

##-------------------
# Input/Output
io <- list()
io$data_dir <- raw_data_dir
io$meta_file <- "meta_mlv_livreg.csv"
io$out_dir <- paste0(results_dir, "/00_qc/")
io$tenx_dir <- "/outs/"


##-------------------
# Opts
opts <- list()
# Threshold for nFeatures
opts$nfeat_thresh <- 800
# Threshold for mitochondrial content
opts$mito_thresh <- 20
# Columns in metadata file to store in Seurat object
opts$meta_colnames <- c("condition", "pass_qc", "slide_sn",
                        "slide_area", "technology")
# QCs to plot (present in Seurat meta.data)
opts$qc_to_plot <- c("nFeature_Spatial", "nCount_Spatial", "percent.mito")
# Opacity of spots. Vector specifying the min and max range of values (between 0 and 1).
opts$alpha <- c(0.1, 0.9)
# Size factor for spot size to plot on tissue
opts$pt.size.factor <- 1.1
# Crop the imagge?
opts$crop <- FALSE
# Max cutoff for continuous values
opts$max.cutoff <- "q98"
# Filename to store QCed object
opts$obj_filename <- "seu_qc"


##-------------------
# Read and filter sample metadata file
sample_meta <- read.csv(file = paste0(io$data_dir, io$meta_file)) |>
  dplyr::filter(condition %in% c("Healthy", "apap")) |>
  dplyr::filter(pass_qc == "TRUE") |>
  dplyr::arrange(sample)


##-----------------------------------------------------------------------
# Analysis pipeline
##-----------------------------------------------------------------------
seu <- run_spatial_qc_pipeline(
  data_dir = io$data_dir,
  sample_meta = sample_meta,
  sample_meta_filename = NULL,
  nfeat_thresh = opts$nfeat_thresh,
  mito_thresh = opts$mito_thresh,
  meta_colnames = opts$meta_colnames,
  out_dir = io$out_dir,
  qc_to_plot = opts$qc_to_plot,
  alpha = opts$alpha,
  pt.size.factor = opts$pt.size.factor,
  max.cutoff = opts$max.cutoff,
  crop = opts$crop,
  tenx_dir = io$tenx_dir,
  obj_filename = opts$obj_filename)
