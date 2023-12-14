##-----------------------------------------------------------------------
# Analyse spatial trajectory trends, i.e. identify genes (or genesets)
# that follow specific trends along the spatial trajectories
##-----------------------------------------------------------------------

# Load packages
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SPATA2))
set.seed(12345)
here::here()

# Source utils scripts
source("../zzz_utils.R")
source("zzz_spata_utils.R")

# Samples to perform analysis
samples <- c("xx_sample1", "xx_sample2")

##-----------------------------------------------------------------------
# Load SPATA objects
##-----------------------------------------------------------------------
data_dir <- paste0(results_dir, "/spata/data/")
for (s in samples) {
  sp <- loadSpataObject(directory_spata = paste0(data_dir, s, "/", s, ".RDS"))

  for (tr in getTrajectoryNames(sp)) {
    hvgs <- getGenes(sp)
    tradeseq_test <- run_tradeSeq(spata_obj = sp, traj_name = tr,
                                  sample = s, genes = hvgs, nknots = 4)
    # Assess trajectory trends
    spata_test <- assess_trajectory_trends(
      object = sp, trajectory_name = tr, variables = hvgs, binwidth = 5,
      whole_sample = FALSE, verbose = TRUE, of_sample = NA, method_gs = "mean")

    # Combine tests together
    traj_test <- dplyr::inner_join(spata_test, tradeseq_test,
                              by = c("variables" = "variables"))
      # mutate(joint_score = 1/(auc/10) + waldStat) |> arrange(-joint_score)
    write.csv(x = traj_test, file = paste0(data_dir, s, "/", s, "_gene_", tr, ".csv"))
  }
}
