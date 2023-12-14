##-----------------------------------------------------------------------
# Create spatial trajectories
##-----------------------------------------------------------------------

# Load packages
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SPATA2))
set.seed(12345)
here::here()

# Source utils scripts
source("../zzz_utils.R")

##-----------------------------------------------------------------------
# Load QCed data
##-----------------------------------------------------------------------
data_dir <- paste0(results_dir, "/00_qc/")
obj <- readRDS(paste0(data_dir, "/seu_qc.rds"))
seu_obj <- obj$seu
opts <- obj$opts
rm(obj)
gc(verbose = FALSE)



##-------------------
# Input/Output

# Define the sample we want to work with
sample <- "HEA_S4"
sample_old <- "HEA4"
io <- list()
io$out_dir <- paste0(results_dir, "/spata/data/", sample, "/")
if (!dir.exists(io$out_dir)) {dir.create(io$out_dir, recursive = TRUE)}

##-------------------
# Preprocessing

# Number of PCA dimensions
npcs <- 40
# Number of HVGs
n_hvgs <- 7000

# Keep only specific sample
seu <- seu_obj[[sample]]
seu <- subset(seu, features = rownames(seu)[which(Matrix::rowSums(seu[["Spatial"]]@counts) > 30) ])
seu <- lognormalize_and_pca(seu = seu, npcs = npcs, n_hvgs = n_hvgs)

##-------------------
# Transform to SPATA object

# Transform to SPATA object
sp <- transformSeuratToSpata(seurat_object = seu, sample_name = sample, method = "spatial")
sp <- adjustGeneSetDf(sp, limit = 20)
# Define default directory to store object
sp <- adjustDirectoryInstructions(object = sp, to = "spata_object",
                                  directory_new = paste0(io$out_dir, sample, ".RDS"))
# Adjust default instructions
sp <- adjustDefaultInstructions(
  sp, display_image = TRUE, display_residuals = FALSE, method_gs = "mean",
  pt_clrp = "npg", clrp = "npg", pt_size = 1.4)


##-------------------
# Draw trajectories

# open interactive application
# sp <- createTrajectories(object = sp)

sp_old <- loadSpataObject(directory_spata = paste0(base_dir, "/analysis/mlv_apap/16_spata/data/",
                                                   sample_old, "/", sample_old, ".RDS"))
sp@trajectories[[sample]] <- sp_old@trajectories[[sample_old]]

for (z in names(sp@trajectories[[sample]])) {
  sp@trajectories[[sample]][[z]]@sample <- sample
  sp@trajectories[[sample]][[z]]@compiled_trajectory_df$sample <- sample

  sp@trajectories[[sample]][[z]]@compiled_trajectory_df$barcodes <- gsub(sample_old, sample, sp@trajectories[[sample]][[z]]@compiled_trajectory_df$barcodes)
}


# Store SPATA object
saveSpataObject(object = sp)
