
###-------------------------------------------------------------
# Load packages
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))


###-------------------------------------------------------------
# Known lineage marker genes
human_liver_lineage_markers <- list(
  lin_endo = c("PECAM1", "CDH5", "ERG", "ENG", "CLEC4G", "ICAM2", "KDR", "STAB1", "CD34"),
  lin_mesenchyme = c("PDGFRB", "PDGFRA", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "DES", "DCN"),
  lin_hepatocyte = c("ALB", "TTR", "TF", "HP", "HNF4A", "CYP2A6"),
  lin_cholangiocyte = c("EPCAM", "KRT19", "CD24"),
  lin_mp = c("CD68", "ITGAM", "ITGAX", "HLA-DRA", "CSF1R", "CD14"),
  lin_tcell = c("CD3D", "CD3E", "CD3G", "CD8A", "CD4"),
  lin_bcell = c("CD79A", "CD79B", "CD19", "MS4A1"),
  lin_ilc = c("KLRF1", "KLRC1", "GZMA", "GZMB", "NKG7"),
  lin_pdc = c("LILRA4", "CLEC4C", "GZMB"),
  lin_mast = c("KIT", "TPSAB1", "TPSB2"),
  lin_neuronal = c("NRXN1", "NRXN3"),
  lin_plasma = c("CD79A", "JCHAIN", "IGHA2"),
  lin_rbc = c("HBA1", "HBA2", "HBB"),
  lin_cycling = c("MKI67", "CCNA2", "CCNB2", "STMN1")
)

mouse_liver_lineage_markers <- list(
  lin_endo = c("Pecam1", "Icam2", "Cdh5", "Kdr", "Eng", "Flt4", "Stab2"),
  lin_mesenchyme = c("Pdgfrb", "Col1a1", "Col1a2", "Col3a1", "Des", "Acta2", "Dcn"),
  lin_hepatocyte = c("Alb","Ttr","Tf","Hp","Hnf4a", "Apoa1", "Cyp2e1"),
  lin_cholangiocyte = c("Epcam","Cdh1","Krt8","Krt19"),
  lin_mesothelia = c("Gpm6a", "Msln", "Pdpn"),
  lin_leucocyte = c("Ptprc", "Lyz2", "Gzma", "Cd3g", "Siglech", "Cd79a", "Csf1r"),
  lin_tcell = c("Cd3g", "Cd3e", "Cd3d"),
  lin_bcell = c("Cd79a", "Cd79b", "Cd19"),
  lin_ilc = c("Gzma", "Gzmb"),
  lin_pdc = c("Siglech"),
  lin_cdc = c("Cd24a", "Clec10a"),
  lin_mp = c("Lyz2", "Csf1r", "Cd68", "Cd5l"),
  lin_neutrophil = c("Cxcr2"),
  lin_cycling = c("Mki67", "Top2a", "Ccnb2")
)


cell_cycle_score <- function(seu, s_genes, g2m_genes, thresh = 0.1) {
  # Compute S and G2M score
  seu <- Seurat::CellCycleScoring(seu, s.features = s_genes,
                                  g2m.features = g2m_genes, set.ident = FALSE)
  # Compute maximum score of S or G2M (idea of a cycling cell irrespective of stage)
  seu$SG2Mmax.Score <- apply(seu@meta.data[, c("S.Score", "G2M.Score")], 1, max)
  # Average cycling score of a cell (not the best approach since the
  # scores are competing with each other)
  seu$SG2M.Score <- (seu$S.Score + seu$G2M.Score)/2

  seu$Phase <- factor(seu$Phase,
                      levels = c("G1","S","G2M"),
                      labels = c("G1","S","G2M"),
                      ordered = TRUE)
  seu$Phase_inferred <- "G1"
  seu$Phase_inferred[seu$S.Score > 0.1] <- "S"
  seu$Phase_inferred[seu$G2M.Score > 0.1] <- "G2M"
  seu$Phase_inferred <- factor(seu$Phase_inferred,
                               levels = c("G1", "S", "G2M"),
                               labels = c("G1", "S", "G2M"),
                               ordered = TRUE)
  return(seu)
}

geneset_analysis <- function(
    gene_markers, out_dir, species, topn_genes = NULL, cluster_pval_adj = 0.05,
    quick_analysis = FALSE, ontologies = c("CC", "MF", "BP"), gs_pval_adj = 0.05, gs_qval = 0.1) {

  # Which species do we perform geneset analysis
  if (species == "human") {
    orgdb <- "org.Hs.eg.db"
    keggdb <- "hsa"
    pathwaydb <- "human"
  } else if (species == "mouse") {
    orgdb <- "org.Mm.eg.db"
    keggdb <- "mmu"
    pathwaydb <- "mouse"
  } else {
    stop("Unrecognised 'species' argument.")
  }

  # Search higher levels of Ontology only in quick_analysis mode
  if (quick_analysis) {
    levs <- seq(4, 6)
  } else {
    levs <- seq(1, 6)
  }

  if (is.null(ontologies)) {
    ontologies <- c("CC", "MF", "BP")
  }

  if (is.character(gene_markers)) {
    # Load marker genes from each cluster, e.g. output of Seurat::FindAllMarkers
    markers <- read.csv(gene_markers, head = TRUE) |>
      dplyr::filter(p_val_adj < cluster_pval_adj)
    # In case we want to retain topn marker genes per cluster
    if (!is.null(topn_genes)) {
      markers <- markers |> dplyr::group_by(cluster) |>
        dplyr::slice_max(n = topn_genes, order_by = avg_log2FC)
    }
  } else {
    markers <- gene_markers
  }

  # Create a list for each geneset
  geneset <- list()
  for (i in unique(markers$cluster)) {
    geneset[[paste0("Cluster", i)]] <- markers[markers$cluster == i, ] |> tibble::remove_rownames()
  }

  # Add entrez ID column
  for (m in names(geneset)) {
    tmp <- clusterProfiler::bitr(geneset[[m]]$gene, fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = orgdb, drop = FALSE) |>
      dplyr::mutate(order = 1:n()) |>
      dplyr::group_by(SYMBOL) |>
      dplyr::slice_head(n = 1) |> dplyr::ungroup() |>
      dplyr::arrange(order)
    geneset[[m]]$entrez <- tmp$ENTREZID
    geneset[[m]] <- na.omit(geneset[[m]]) %>% tibble::remove_rownames()
  }

  ### Gene Ontology Analyses --------------------------------------------------------------------
  # MF = Molecular Function, BP = Biological Process, CC = Cellular Compartment

  ##
  # 1. Check for overlaps between each gene and the gene set associated
  #    with each GO term at a given semantic level.
  ##
  if (!dir.exists(paste0(out_dir, "/go_membership"))) {
    dir.create(paste0(out_dir, "/go_membership"), recursive = TRUE)
  }
  if (!quick_analysis) {
    lapply(X = names(geneset), FUN = function(set) {
      level = 3
      while (level < 7) {
        for (d in ontologies) {
          go_group <- clusterProfiler::groupGO(
            gene = geneset[[set]]$gene, keyType = "SYMBOL", OrgDb = orgdb,
            ont = d, level = level, readable = FALSE)
          go_group <- dplyr::filter(go_group@result, Count >= 1)
          go_group <- go_group[order(-go_group$Count), ]
          write.csv(go_group, paste0(out_dir, "/go_membership/",
                                     set, "_GO_", d, "_group_level", level, ".csv"))
        }
        level = level + 1
      }
    })
  }

  ##
  # 2. GO term enrichment analysis: identifies statistical enrichment
  #    of genes belonging to GO terms.
  ##
  out_dir_enr <- paste0(out_dir, "/go_enrichment/")
  if (!dir.exists(paste0(out_dir_enr, "/outs/"))) {
    dir.create(paste0(out_dir_enr, "/outs/"), recursive = TRUE)
  }
  for (set in names(geneset)) {
    genes_fc <- geneset[[set]]$avg_log2FC
    names(genes_fc) <- geneset[[set]]$gene

    for (d in ontologies) {
      go_enr <- clusterProfiler::enrichGO(
        gene = geneset[[set]]$gene, OrgDb = orgdb, keyType = "SYMBOL", ont = d,
        pAdjustMethod = "BH", pvalueCutoff = gs_pval_adj, qvalueCutoff = gs_qval)

      for (l in levs) {
        go_enr_level <- clusterProfiler::gofilter(go_enr, level = l)
        if (length(go_enr_level@result$geneID) > 0) {
          write.csv(go_enr_level@result,
                    paste0(out_dir_enr, "/outs/", set, "_GO_", d, "_enrichment_level", l, ".csv"))

          # filter for terms which are statistical significant and have more than 2 genes
          tmp <- dplyr::filter(go_enr_level, Count > 2)
          if (length(tmp@result$geneID) > 1 && length(tmp@result$p.adjust) > 1 &&
              tmp@result$p.adjust[[1]] < gs_pval_adj) {
            png(paste0(out_dir_enr, "/", set, "_GO_", d, "_level", l, "_ConceptNetwork.png"),
                width = 10, height = 7, res = 300, units = "in")
            p1 <- enrichplot::cnetplot(tmp, foldChange = genes_fc, colorEdge = TRUE) +
              ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
              ggplot2::ggtitle(paste0("GO Enrichment ", d, ":semantic level ", l))
            print(p1)
            dev.off()

            png(paste0(out_dir_enr, "/", set, "_GO_", d, "_level", l, "_DotPlot.png"),
                width = 13, height = 13, res = 300, units = "in")
            p2 <- enrichplot::dotplot(tmp, showCategory = 20) +
              ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
              ggplot2::theme_classic() +
              ggplot2::ggtitle(paste0("GO Enrichment ", d, ":semantic level ", l))
            print(p2)
            dev.off()

            png(paste0(out_dir_enr, "/", set, "_GO_", d, "_level", l, "_Heatmap.png"),
                width = 15, height = 8, res = 300, units = "in")
            p3 <- enrichplot::heatplot(tmp, foldChange = genes_fc) +
              ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
              ggplot2::ggtitle(paste0("GO Enrichment ", d, ":semantic level ", l))
            print(p3)
            dev.off()
          }
        }
      }
    }
  }


  ### Cluster Comparison ----------------------------------------------------------------------------
  if (!dir.exists(paste0(out_dir, "/cluster_comparison"))) {
    dir.create(paste0(out_dir, "/cluster_comparison"), recursive = TRUE)
  }
  # Keep only a list of genes (entrez IDs) for each cluster
  genes <- list()
  for (set in names(geneset)) {
    genes[[set]] <- geneset[[set]]$entrez
  }

  # Compute enriched genesets
  for (d in ontologies) {
    ego <- clusterProfiler::compareCluster(
      geneCluster = genes, fun = "enrichGO", OrgDb = orgdb, ont = d,
      pvalueCutoff = gs_pval_adj, qvalueCutoff = gs_qval)
    write.csv(ego@compareClusterResult, paste0(out_dir, "/cluster_comparison/enrichGO_", d, ".csv" ))

    png(paste0(out_dir, "/cluster_comparison/enrichGO_", d, ".png"),
        width = 20, height = 25, res = 300, units = "in")
    p <- enrichplot::dotplot(ego) +
      ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
      ggplot2::theme_classic()
    print(p)
    dev.off()
  }

  ekegg <- clusterProfiler::compareCluster(
    geneCluster = genes, fun = "enrichKEGG", organism = keggdb,
    pvalueCutoff = gs_pval_adj, qvalueCutoff = gs_qval)
  write.csv(ekegg@compareClusterResult, paste0(out_dir, "/cluster_comparison/enrichKEGG.csv" ))

  png(paste0(out_dir, "/cluster_comparison/enrichKEGG.png"),
      width = 20, height = 25, res = 300, units = "in")
  p <- enrichplot::dotplot(ekegg) +
    ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
    ggplot2::theme_classic()
  print(p)
  dev.off()
}
