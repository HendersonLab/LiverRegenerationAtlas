# Function for running tradeSeq on a SPATA object to identify genes
# that show correlation with spatial trajectory.
run_tradeSeq <- function(spata_obj, traj_name, sample, genes = NULL, nknots = 4) {
  # Extract trajectory information
  traj_df <- SPATA2::getTrajectoryObject(object = spata_obj, trajectory_name = traj_name,
                                         of_sample = sample)@compiled_trajectory_df |>
    dplyr::mutate(order = 1:n()) |>
    dplyr::group_by(barcodes) |>
    dplyr::slice_head(n = 1) |> dplyr::ungroup() |>
    dplyr::arrange(order) |>
    dplyr::mutate(pseudotime = 1:n()) |>
    dplyr::mutate(spot_weight = 1)

  # Subset to retain the counts from specific spots/cells
  if (is.null(genes)) {
    gex <- spata_obj@data[[sample]]$counts[, traj_df$barcodes]
  } else {
    gex <- spata_obj@data[[sample]]$counts[genes, traj_df$barcodes]
  }

  # Get pseudotime and spot weights
  pseudotime <- traj_df |> select(barcodes, pseudotime) %>% tibble::column_to_rownames(var = "barcodes")
  spot_weight <- traj_df |> select(barcodes, spot_weight) %>% tibble::column_to_rownames(var = "barcodes")

  # Run tradeSeq to identify genes that correlate with spatial trajectory
  sce <- tradeSeq::fitGAM(counts = gex, pseudotime = pseudotime,
                          cellWeights = spot_weight, nknots = nknots, sce = TRUE, verbose = FALSE)
  # Perform the actual test
  tradeseq_test <- tradeSeq::associationTest(sce) |> dplyr::arrange(pvalue) |>
    tibble::rownames_to_column(var = "variables") |> as_tibble()
  return(tradeseq_test)
}


# Function for running tradeSeq on a SPATA object to identify genes
# that show correlation with spatial trajectory.
run_tradeSeq_geneset <- function(spata_obj, traj_name, sample, genesets = NULL, genes = NULL,
                                 nknots = 4, filter_gs = 0.25) {
  # Extract trajectory information
  traj_df <- SPATA2::getTrajectoryObject(object = spata_obj, trajectory_name = traj_name,
                                         of_sample = sample)@compiled_trajectory_df |>
    dplyr::mutate(order = 1:n()) |>
    dplyr::group_by(barcodes) |>
    dplyr::slice_head(n = 1) |> dplyr::ungroup() |>
    dplyr::arrange(order) |>
    dplyr::mutate(pseudotime = 1:n()) |>
    dplyr::mutate(spot_weight = 1)

  # Gene expression matrix
  gex <- spata_obj@data[[sample]]$counts

  # Subset to retain the gex counts from specific spots/cells
  if (is.null(genes)) {
    gex <- gex[, traj_df$barcodes]
  } else {
    gex <- gex[genes, traj_df$barcodes]
  }

  # # Remove lowly expressed genes
  # gex <- gex[rowSums(gex) > 40, ]

  # Get mapping of genes to genesets
  geneset_df <- SPATA2::getGeneSetDf(object = spata_obj)

  # Initiate joined dataframe
  joined_df <- SPATA2::getCoordsDf(spata_obj) |>
    dplyr::filter(barcodes %in% traj_df$barcodes) %>%
    dplyr::select(barcodes)

  for(i in seq_along(genesets)) {
    # get genes of gene set
    genes_df <- dplyr::filter(geneset_df, ont %in% genesets[i])
    n_genes <- base::nrow(genes_df)
    # get genes found in expression matrix
    genes <- dplyr::filter(genes_df, gene %in% base::rownames(gex)) |>
      dplyr::pull(gene)
    n_found_genes <- base::length(genes)

    # calculate percentage of genes found
    p_found_genes <- base::round(n_found_genes/n_genes, digits = 2)
    # make sure that percentage is equal to or higher than the threshold
    if (p_found_genes >= filter_gs) {
      geneset_vls <- ceiling(Matrix::colMeans(gex[genes, ])) |>
        base::as.data.frame() |>
        tibble::rownames_to_column(var = "barcodes")
      colnames(geneset_vls) <- c("barcodes", genesets[i])
      # gradually add gene-set columns to joined_df
      joined_df <- dplyr::left_join(x = joined_df, y = geneset_vls, by = "barcodes")
    }
  }

  # Remove column and transpose matrix
  gs_expr <- joined_df |> tibble::column_to_rownames(var = "barcodes")
  gs_expr <- t(gs_expr)

  # Get pseudotime and spot weights
  pseudotime <- traj_df |> select(barcodes, pseudotime) %>% tibble::column_to_rownames(var = "barcodes")
  spot_weight <- traj_df |> select(barcodes, spot_weight) %>% tibble::column_to_rownames(var = "barcodes")

  # Run tradeSeq to identify genes that correlate with spatial trajectory
  sce <- tradeSeq::fitGAM(counts = gs_expr, pseudotime = pseudotime,
                          cellWeights = spot_weight, nknots = nknots, sce = TRUE, verbose = FALSE)
  # Perform the actual test
  tradeseq_test <- tradeSeq::associationTest(sce) |> dplyr::arrange(pvalue) |>
    tibble::rownames_to_column(var = "variables") |> as_tibble()
  return(tradeseq_test)
}


## Local implementation to manually provide method for assessing trajectory trends
assess_trajectory_trends <- function(
    object, trajectory_name, variables, binwidth = 5, whole_sample = FALSE,
    verbose = TRUE, of_sample = NA, method_gs = "mean") {
  SPATA2::check_object(object)
  of_sample <- SPATA2::check_sample(object, of_sample, desired_length = 1)
  SPATA2:::check_trajectory(object, trajectory_name, of_sample)
  stdf <- SPATA2::getTrajectoryDf(object = object, trajectory_name = trajectory_name,
                          of_sample = of_sample, variables = variables, binwidth = binwidth,
                          verbose = verbose, method_gs = method_gs)
  rtdf <- SPATA2:::hlpr_rank_trajectory_trends(stdf = stdf, verbose = verbose)
  atdf <- SPATA2:::hlpr_assess_trajectory_trends(rtdf = rtdf, verbose = verbose)
  # Store old AUC value
  atdf$default_auc <- atdf$auc
  # Compute relative auc: nomralise by largest AUC value
  atdf$auc <- (atdf$auc / quantile(atdf$auc, .99, na.rm = TRUE)) * 10
  return(atdf)
}

filter_trajectory_trends <- function(
    df, auc_thresh = 4, pvalue_thresh = 0.2, wald_stat_thresh = 0,
    topn = 5, trends = "all", variables_only = TRUE) {
  SPATA2:::check_atdf(df)
  all_patterns <- dplyr::pull(df, var = "pattern") %>% base::unique()
  trajectory_patterns <- c(all_patterns, trajectory_patterns)
  if (base::all(trends == "all")) {
    trends <- trajectory_patterns
  }
  confuns::is_vec(x = trends, mode = "character", "trends")
  trends <- confuns::check_vector(input = trends, against = trajectory_patterns,
                                  verbose = TRUE, ref.input = "argument 'trends'",
                                  ref.against = "known trajectory trends")
  if (base::isTRUE(variables_only)) {
    res <- SPATA2::hlpr_filter_trend(atdf = df, limit = topn, poi = trends)
  }
  else {
    res <- dplyr::filter(.data = df, pattern %in% trends) |>
      dplyr::filter(auc <= auc_thresh) |>
      dplyr::filter(pvalue < pvalue_thresh) |>
      dplyr::filter(waldStat > wald_stat_thresh) |>
      dplyr::group_by(variables) %>%
      dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
      dplyr::group_by(pattern) %>% dplyr::arrange(auc, .by_group = TRUE)
  }
  base::return(res)
}


# Define helper functions
plot_trajectory_heatmap <- function(
    object, trajectory_name, variables, binwidth = 5, arrange_rows = "none",
    colors = NULL, method_gs = NULL, show_rownames = NULL, show_colnames = NULL,
    split_columns = NULL, smooth_span = NULL, verbose = NULL, cluster_rows = FALSE,
    of_sample = NA, atdf = NULL, ...) {
  SPATA2:::hlpr_assign_arguments(object)
  SPATA2:::check_trajectory_binwidth(binwidth)
  confuns::are_values(c("method_gs", "arrange_rows"), mode = "character")
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)
  SPATA2:::check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)
  SPATA2::check_method(method_gs = method_gs)
  variables <- SPATA2::check_variables(variables = variables, all_gene_sets = getGeneSets(object),
                                       all_genes = getGenes(object), max_slots = 1)
  var_type <- "variables"
  smooth <- TRUE
  trajectory_object <- SPATA2::getTrajectoryObject(object = object,
                                                   trajectory_name = trajectory_name, of_sample = of_sample)
  stdf <- SPATA2::hlpr_summarize_trajectory_df(object = object, ctdf = trajectory_object@compiled_trajectory_df,
                                               variables = variables[[1]], binwidth = binwidth,
                                               method_gs = method_gs, verbose = verbose) %>% dplyr::ungroup()
  wide_tdf <- dplyr::group_by(.data = stdf, {{var_type}}) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    dplyr::ungroup() %>% tidyr::pivot_wider(id_cols = dplyr::all_of(var_type),
                                            names_from = c("trajectory_part", "trajectory_order"),
                                            names_sep = "_", values_from = "values")
  n_parts <- base::length(base::unique(trajectory_object@compiled_trajectory_df$trajectory_part))
  if (base::isTRUE(split_columns) && n_parts > 1) {
    gaps_col <- dplyr::select(.data = stdf, trajectory_part,
                              trajectory_part_order) %>% dplyr::distinct() %>%
      dplyr::group_by(trajectory_part) %>% dplyr::summarise(count = dplyr::n()) %>%
      dplyr::mutate(positions = base::cumsum(count) * 10) %>%
      dplyr::pull(positions) %>% base::as.numeric()
  }
  else {
    gaps_col <- NULL
  }
  mtr <- base::as.matrix(dplyr::select(.data = wide_tdf, -{{var_type}}))
  base::rownames(mtr) <- dplyr::pull(.data = wide_tdf, var_type)
  keep <- base::apply(mtr, MARGIN = 1, FUN = function(x) {dplyr::n_distinct(x) != 1})
  n_discarded <- base::sum(!keep)
  if (base::isTRUE(smooth) && n_discarded != 0) {
    discarded <- base::rownames(mtr)[!keep]
    discarded_ref <- stringr::str_c(discarded, collapse = ", ")
    mtr <- mtr[keep, ]
    base::warning(glue::glue("Discarded {n_discarded} variables due to uniform expression.
                             (Can not smooth uniform values.): '{discarded_ref}'"))
  }
  mtr_smoothed <- matrix(0, nrow = nrow(mtr), ncol = ncol(mtr) * 10)
  ## Andreas: Keep only first 50 characters
  base::rownames(mtr_smoothed) <- stringr::str_sub(base::rownames(mtr), 1, 50)
  atdf$variables <- stringr::str_sub(atdf$variables, 1, 50)
  if (base::isTRUE(smooth)) {
    confuns::give_feedback(msg = glue::glue("Smoothing values with smoothing
                                            span: {smooth_span}."), verbose = verbose)
    for (i in 1:base::nrow(mtr)) {
      x <- 1:base::ncol(mtr)
      values <- base::as.numeric(mtr[i, ])
      y <- (values - base::min(values))/(base::max(values) - base::min(values))
      model <- stats::loess(formula = y ~ x, span = smooth_span)
      mtr_smoothed[i, ] <- stats::predict(model, seq(1,base::max(x), length.out = base::ncol(mtr) * 10))
    }
  }
  if (base::all(arrange_rows == "maxima") | base::all(arrange_rows == "minima")) {
    mtr_smoothed <- confuns::arrange_rows(df = base::as.data.frame(mtr_smoothed),
                                          according.to = arrange_rows, verbose = verbose) %>% base::as.matrix()
  } else if (base::all(arrange_rows == "zonation")){
    # Split by zonation
    atdf_list <- atdf %>% group_by(zonation) %>% group_split()
    df <- base::as.data.frame(mtr_smoothed)
    tmp <- list()
    for (z in 1:length(atdf_list)) {
      tmp[[z]] <- confuns::arrange_rows(df = df[atdf_list[[z]]$variables,], according.to = "maxima", verbose = verbose)
    }
    mtr_smoothed <- do.call(rbind, tmp) %>% base::as.matrix()
  } else if (base::all(arrange_rows == "region")){
    # Split by region
    atdf_list <- atdf %>% group_by(region) %>% group_split()
    df <- base::as.data.frame(mtr_smoothed)
    tmp <- list()
    for (z in 1:length(atdf_list)) {
      tmp[[z]] <- confuns::arrange_rows(df = df[atdf_list[[z]]$variables,], according.to = "maxima", verbose = verbose)
    }
    mtr_smoothed <- do.call(rbind, tmp) %>% base::as.matrix()
  }
  pheatmap::pheatmap(mat = mtr_smoothed, cluster_cols = FALSE, cluster_rows = cluster_rows,
                     color = colors, gaps_col = gaps_col[1:(base::length(gaps_col) - 1)],
                     show_colnames = show_colnames, show_rownames = show_rownames, ...)
}
