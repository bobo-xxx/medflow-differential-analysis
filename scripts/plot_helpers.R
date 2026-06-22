# plot_helpers.R — Visualization functions for differential-analysis
#
# Provides volcano plot, heatmap, Venn diagram, and chromosome location plot.
# Adapted from original/scripts/diff_volcano.R, diff_heatmap.R, diff_venn.R, diff_locate.R.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

#' Generate volcano plot
#'
#' @param dif Data frame with DE results (columns: gene_id, logFC, p-value column)
#' @param p_name Character, name of p-value column ("Pvalue" or "Padj")
#' @param p_value Numeric, p-value threshold for horizontal line
#' @param logfc_cutoff Numeric, logFC threshold for vertical lines
#' @param top Integer, number of top genes to label (used if gene is NULL/None)
#' @param gene Character, comma-separated gene IDs to label, or NULL/None to auto-label top N
#' @param outfile Character, output PDF file path
#' @return NULL (side effect: writes PDF)
plot_volcano <- function(dif, p_name, p_value, logfc_cutoff, top, gene, outfile) {
  fc_name <- "logFC"
  gene_col <- colnames(dif)[1]

  # Build group classification
  dif_up <- dif %>%
    filter(!!sym(fc_name) > logfc_cutoff & !!sym(p_name) < p_value)
  dif_down <- dif %>%
    filter(!!sym(fc_name) < -logfc_cutoff & !!sym(p_name) < p_value)
  dif_plot <- dif %>% mutate(group = case_when(
    dif[[gene_col]] %in% dif_up[[gene_col]] ~ "Up",
    dif[[gene_col]] %in% dif_down[[gene_col]] ~ "Down",
    TRUE ~ "Not"
  ))

  dif_plot$logP <- -log10(dif_plot[[p_name]])
  dif_plot$change <- factor(dif_plot$group, levels = c("Down", "Not", "Up"))
  dif_plot <- dif_plot[order(abs(dif_plot[[fc_name]]), decreasing = TRUE), ]

  # Determine genes to label
  if (is.null(gene) || gene == "NULL" || gene == "None" || gene == "") {
    if (is.null(top) || top == "NULL" || top == "None") {
      label_genes <- NULL
    } else {
      label_genes <- head(dif_plot[[gene_col]], as.integer(top))
    }
  } else {
    label_genes <- strsplit(as.character(gene), ",")[[1]]
  }

  lab_y <- paste0("-Log10(", p_name, ")")
  x_lim <- max(abs(dif_plot[[fc_name]])) * 1.1

  p <- ggplot(dif_plot, aes(.data[[fc_name]], .data$logP, color = change)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(x = "LogFC", y = lab_y, color = "Significance") +
    geom_hline(yintercept = -log10(p_value), lty = 2) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), lty = 2) +
    scale_x_continuous(limits = c(-x_lim, x_lim)) +
    scale_color_manual(values = c(
      "Down" = "#4DBBD5",
      "Not" = "grey",
      "Up" = "#E64B35"
    )) +
    ggtitle("Volcano") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text = element_text(size = 12)
    ) +
    theme(panel.grid = element_blank())

  if (!is.null(label_genes) && length(label_genes) > 0) {
    dif_plot$label <- mapply(function(x) {
      ifelse(x %in% label_genes, x, "")
    }, dif_plot[[gene_col]])

    p <- p + ggrepel::geom_text_repel(
      data = dif_plot, aes(label = label),
      segment.alpha = 0.4,
      box.padding = 0.5,
      force = 1,
      max.overlaps = Inf,
      min.segment.length = 0.25,
      show.legend = FALSE
    )
  }

  ggsave(file = outfile, p, width = 8, height = 5)
  report_info(sprintf("Volcano plot saved to %s", outfile))
}

#' Generate heatmap of top DEGs
#'
#' @param mat Expression matrix (genes x samples), full dataset
#' @param map Data frame with sample-to-group mapping
#' @param rdegs Data frame of filtered DEGs with gene_id and logFC columns
#' @param top Integer, number of top up/down genes for heatmap display
#' @param color_heat Character, comma-separated heatmap colors (e.g., "blue,white,red")
#' @param outfile Character, output PDF file path
#' @return NULL (side effect: writes PDF)
plot_heatmap <- function(mat, map, rdegs, top, color_heat, outfile) {
  suppressPackageStartupMessages(library(pheatmap))

  map[[2]] <- factor(map[[2]], levels = unique(map[[2]]))

  # Select top N up and top N down genes for heatmap
  if (nrow(rdegs) > top * 2) {
    rdegs_sorted <- rdegs[order(rdegs[["logFC"]]), ]
    top_mat <- rdegs_sorted[c(1:top, (nrow(rdegs_sorted) - top + 1):nrow(rdegs_sorted)), ]
    deg_heatmap <- top_mat[[1]]  # gene_id column
  } else {
    deg_heatmap <- rdegs[[1]]
  }

  dat_heatmap <- mat[deg_heatmap, map[[1]], drop = FALSE]

  color_vec <- strsplit(as.character(color_heat), ",")[[1]]

  pdf(NULL)  # suppress pheatmap implicit pdf call
  p_heatmap <- pheatmap(dat_heatmap,
    scale = "row",
    annotation_col = data.frame(
      Group = map[[2]],
      row.names = map[[1]]
    ),
    annotation_colors = list(Group = setNames(
      c("#E64B35", "#4DBBD5")[seq_len(length(unique(map[[2]])))],
      unique(map[[2]])
    )),
    color = colorRampPalette(color_vec)(50),
    breaks = c(seq(-3, 3, length = 50)),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    labels_col = "",
    border_color = NA,
    main = "Heatmap"
  )
  dev.off()

  ggsave(file = outfile, p_heatmap, width = 8, height = 8)
  report_info(sprintf("Heatmap saved to %s", outfile))
}

#' Generate Venn diagram of DEGs vs related gene set
#'
#' @param dif Data frame of DEGs with group column (Up/Down) and gene_id column
#' @param rgs Character vector of related gene set identifiers
#' @param pheno_abbr Character, phenotype abbreviation for Venn label
#' @param color_panel Character, comma-separated colors for Venn sets
#' @param outfile Character, output PDF file path
#' @return NULL (side effect: writes PDF)
plot_venn <- function(dif, rgs, pheno_abbr, color_panel, outfile) {
  if (!requireNamespace("ggvenn", quietly = TRUE)) {
    report_info("ggvenn not installed, skipping Venn diagram")
    return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(ggvenn))

  gene_col <- colnames(dif)[1]
  deg_genes <- dif[dif[["group"]] == "Up" | dif[["group"]] == "Down", ][[gene_col]]
  colors <- strsplit(as.character(color_panel), ",")[[1]]

  p_venn <- ggvenn(
    setNames(
      list(deg_genes, rgs),
      c("DEGs", pheno_abbr)
    ),
    c("DEGs", pheno_abbr),
    show_percentage = FALSE,
    fill_alpha = 0.5,
    stroke_color = NA,
    fill_color = color_map(colors, c("DEGs", pheno_abbr))
  )

  ggsave(file = outfile, p_venn, width = 8, height = 5)
  report_info(sprintf("Venn diagram saved to %s", outfile))
}

#' Generate chromosome location plot via RCircos
#'
#' @param gene_list Character vector of gene identifiers to plot
#' @param locate Character, path to chromosome annotation CSV (columns: Gene, Chromosome, Start, End)
#' @param tax_id Character, NCBI taxonomy ID ("9606" for human, "10090" for mouse)
#' @param outfile Character, output PDF file path
#' @return NULL (side effect: writes PDF)
plot_locate <- function(gene_list, locate, tax_id, outfile) {
  if (!requireNamespace("RCircos", quietly = TRUE)) {
    report_info("RCircos not installed, skipping chromosome location plot")
    return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(RCircos))

  chr <- read.csv(locate)
  gene <- gene_list

  chr_gene <- chr[which(chr$Gene %in% gene), ]

  if (nrow(chr_gene) == 0) {
    report_info("No genes matched in chromosome annotation, skipping location plot")
    return(invisible(NULL))
  }

  pdf(file = outfile, width = 8, height = 8)

  if (tax_id == "9606") {
    data(UCSC.HG38.Human.CytoBandIdeogram)
    cyto_info <- UCSC.HG38.Human.CytoBandIdeogram
  } else if (tax_id == "10090") {
    data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
    cyto_info <- UCSC.Mouse.GRCm38.CytoBandIdeogram
  } else {
    dev.off()
    report_exception_ndjson(
      "E802_UNSUPPORTED_TAXID", "data_mismatch", "halt",
      sprintf("Unsupported tax_id: %s. Supported: 9606 (human), 10090 (mouse)", tax_id),
      exit_code = 1
    )
    return(invisible(NULL))
  }

  RCircos.Set.Core.Components(cyto_info)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  RCircos.Gene.Connector.Plot(chr_gene, track.num = 1, side = "in")
  RCircos.Gene.Name.Plot(chr_gene, name.col = 4, track.num = 2, side = "in")
  dev.off()

  report_info(sprintf("Chromosome location plot saved to %s", outfile))
}
