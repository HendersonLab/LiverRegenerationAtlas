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
io$meta_file <- "meta_hlr_livreg.csv"
io$out_dir <- paste0(results_dir, "/00_qc/")
io$plot_dir <- paste0(io$out_dir, "/plots/")
if(!dir.exists(io$plot_dir)) {dir.create(io$plot_dir, recursive = TRUE)}


##-------------------
# Opts
opts <- list()
# Columns in metadata file to store in Seurat object
opts$meta_colnames <- c("sample", "condition", "pass_qc", "radius", "section")
# QCs to plot (present in Seurat meta.data)
opts$qc_to_plot <- c("nFeature_Spatial", "nCount_Spatial")
# Filter genes that are expressed in less than `min.cells`.
opts$min.cells <- 3
# Filter cells that express less than `min.features`.
opts$min.features <- 6
# Filename to store QCed object
opts$obj_filename <- "seu_qc"


##-------------------
# Read and filter sample metadata file
opts$sample_meta <- read.csv(file = paste0(io$data_dir, io$meta_file)) |>
  dplyr::filter(pass_qc == "TRUE") |>
  dplyr::arrange(sample)


##-----------------------------------------------------------------------
# Analysis pipeline
##-----------------------------------------------------------------------
seu <- run_resolve_qc_pipeline(
  data_dir = io$data_dir,
  sample_meta = opts$sample_meta,
  sample_meta_filename = NULL,
  min.cells = opts$min.cells,
  min.features = opts$min.features,
  meta_colnames = opts$meta_colnames,
  out_dir = io$out_dir,
  qc_to_plot = opts$qc_to_plot,
  obj_filename = opts$obj_filename,
  stroke = 0.05)
