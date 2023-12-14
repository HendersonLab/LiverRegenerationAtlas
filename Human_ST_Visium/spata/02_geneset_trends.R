##-----------------------------------------------------------------------
# Analyse spatial trajectory trends, i.e. identify genesets
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
samples <- c("xx_sample1", "xx_sample2", "xx_sample3")

##-----------------------------------------------------------------------
# Load SPATA objects abnd perform analysis
##-----------------------------------------------------------------------
data_dir <- paste0(results_dir, "/spata/data/")

for (s in samples) {
  sp <- loadSpataObject(directory_spata = paste0(data_dir, s, "/", s, ".RDS"))
  # Keep genesets that have high proportion of genes in our gex matrix
  sp <- adjustGeneSetDf(sp, limit = 20)

  for (tr in getTrajectoryNames(sp)) {
    # Get highly variable genes
    hvgs <- getGenes(sp)
    # Get genesets for BP.GO
    genesets <- getGeneSets(sp, of_class = "BP.GO")
    tradeseq_test <- run_tradeSeq_geneset(
      spata_obj = sp, traj_name = tr, sample = s,
      genesets = genesets, genes = hvgs, filter_gs = 0.2, nknots = 4)
    write.csv(x = tradeseq_test, file = paste0(data_dir, s, "/", s, "_tradeseq_geneset_", tr, ".csv"))

    # Assess trajectory trends
    spata_test <- assess_trajectory_trends(
      object = sp, trajectory_name = tr, variables = genesets, binwidth = 5,
      whole_sample = FALSE, verbose = TRUE, of_sample = NA, method_gs = "mean")
    write.csv(x = spata_test, file = paste0(data_dir, s, "/", s, "_spata_geneset_", tr, ".csv"))

    # Combine tests together
    traj_test <- dplyr::inner_join(spata_test, tradeseq_test,
                              by = c("variables" = "variables"))
      # mutate(joint_score = 1/(auc/10) + waldStat) |> arrange(-joint_score)
    write.csv(x = traj_test, file = paste0(data_dir, s, "/", s, "_geneset_", tr, ".csv"))
  }
}
