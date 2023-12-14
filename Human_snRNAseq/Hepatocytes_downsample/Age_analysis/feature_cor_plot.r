feature_cor_plot <- function(seu, features, reduction) {
  PC = QC = cor <- NULL
  # Heatmap showing correlation of metadata with Principal Components
  if (reduction %in% names(seu@reductions)) {
    if (reduction == "pca"){
      var1 <- "PC"
      var2 <- rlang::sym("PC")
    }
    if (reduction == "harmony"){
      var1 <- "harmony"
      var2 <- rlang::sym("harmony")
    }
    df_corr <- abs(stats::cor(seu@reductions[[reduction]]@cell.embeddings,
                              seu[[features]])) |>
      dplyr::as_tibble(rownames = var1) |>
      tidyr::pivot_longer(cols = -c(var1), names_to = "QC", values_to = "cor")
    df_corr[,1] <- factor(
      df_corr |> dplyr::pull(var1), levels = paste0(var1,"_", seq(1, NCOL(seu@reductions[[reduction]]))))
    # Create heatmap
    gg <- ggplot2::ggplot(df_corr, ggplot2::aes(x = !!var2, y = QC, fill = cor)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_distiller(name = "abs(cor)", type="div", palette="RdYlBu", limits = c(0,1)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                         hjust = 1, size = 6))
    return(gg)
  } else {
    return(ggplot2::ggplot())
  }
}
