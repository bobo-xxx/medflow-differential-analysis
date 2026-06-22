# diff_methods.R — Differential expression analysis methods for deg-analysis
#
# Provides four DE analysis functions adapted from original/scripts/diff.R.
# Each function takes an expression matrix and sample group map, returns
# a data.frame with columns: logFC, Pvalue, Padj, [stat].
# Gene identifiers are preserved as the first column (gene_id).

suppressPackageStartupMessages(library(dplyr))

#' DESeq2 differential expression (for count data)
#'
#' @param df Expression count matrix (genes x samples), raw counts
#' @param group Data frame with sample-to-group mapping (col 1: sample, col 2: group)
#' @return Data frame with columns: gene_id, logFC, Pvalue, Padj
diff_deseq2 <- function(df, group) {
  suppressPackageStartupMessages(library(DESeq2))

  col_data <- data.frame(
    row.names = colnames(df),
    group_list = group[[2]]
  )
  dds <- DESeqDataSetFromMatrix(
    countData = round(df),
    colData = col_data,
    design = ~group_list
  )
  dds2 <- DESeq(dds)
  res <- results(dds2)
  dif <- res %>%
    as.data.frame() %>%
    dplyr::rename(
      logFC = "log2FoldChange",
      Pvalue = "pvalue",
      Padj = "padj"
    ) %>%
    na.omit() %>%
    dplyr::arrange(Pvalue)

  dif$gene_id <- rownames(dif)
  rownames(dif) <- NULL
  dif <- dif[, c("gene_id", "logFC", "Pvalue", "Padj")]
  dif
}

#' limma differential expression (for normalized microarray or RNA-seq data)
#'
#' @param df Expression matrix (genes x samples)
#' @param map Data frame with sample-to-group mapping (col 1: sample, col 2: group)
#' @return Data frame with columns: gene_id, logFC, Pvalue, Padj
diff_limma <- function(df, map) {
  pdf(NULL)  # Suppress limma implicit plotSA call
  suppressPackageStartupMessages(library(limma))

  treat_name <- levels(map[[2]])[2]
  con_name <- levels(map[[2]])[1]

  design <- model.matrix(~ 0 + map[[2]])
  rownames(design) <- map[[1]]
  colnames(design) <- levels(map[[2]])

  fit <- lmFit(df, design)
  cont_matrix <- makeContrasts(
    contrasts = paste0(treat_name, "-", con_name),
    levels = design
  )
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2)
  plotSA(fit2)

  dif <- topTable(fit2, coef = 1, n = Inf)
  dif <- dif %>%
    as.data.frame() %>%
    dplyr::rename(Pvalue = "P.Value", Padj = "adj.P.Val") %>%
    na.omit() %>%
    dplyr::arrange(Pvalue)

  dev.off()

  dif$gene_id <- rownames(dif)
  rownames(dif) <- NULL
  # Keep logFC, Pvalue, Padj, gene_id
  dif <- dif[, c("gene_id", "logFC", "Pvalue", "Padj")]
  dif
}

#' edgeR differential expression (for count data)
#'
#' @param df Expression count matrix (genes x samples), raw counts
#' @param map Data frame with sample-to-group mapping
#' @param norm Normalization method: "TMM", "RLE", "upperquartile", "none"
#' @param model Model fitting method: "glmFit" or "glmQLFit"
#' @return Data frame with columns: gene_id, logFC, Pvalue, Padj
diff_edger <- function(df, map, norm = "TMM", model = "glmFit") {
  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(statmod))

  dgelist <- DGEList(counts = df, group = map[[2]])
  keep <- rowSums(cpm(dgelist) > 1) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist_norm <- calcNormFactors(dgelist, method = norm)
  design <- model.matrix(~ map[[2]])
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE)

  func <- get(model)
  fit <- func(dge, design, robust = TRUE)
  if (model == "glmQLFit") {
    lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
  } else {
    lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
  }

  dif <- lrt %>%
    as.data.frame() %>%
    dplyr::rename(Pvalue = "PValue", Padj = "FDR") %>%
    na.omit() %>%
    dplyr::arrange(Pvalue)

  dif$gene_id <- rownames(dif)
  rownames(dif) <- NULL
  dif <- dif[, c("gene_id", "logFC", "Pvalue", "Padj")]
  dif
}

#' Row-wise fold change calculation
#'
#' @param row Numeric vector of expression values for one gene
#' @param treat Character vector of treatment sample column names
#' @param control Character vector of control sample column names
#' @return log2 fold change (numeric)
logfc_row <- function(row, treat, control) {
  mean_treatment <- mean(row[treat])
  mean_control <- mean(row[control])
  fc <- mean_treatment / mean_control
  if (fc <= 0 || is.na(fc)) return(0)
  log2(fc)
}

#' Row-wise statistical test
#'
#' @param row Numeric vector of expression values
#' @param group Factor vector of group assignments (same length as row)
#' @param func Test function (t.test or wilcox.test)
#' @return Numeric vector c(p.value, statistic)
test_row <- function(row, group, func = t.test) {
  df <- data.frame(x = row, Group = group)
  res <- tryCatch(
    {
      result <- func(x ~ Group, data = df)
      c(result$p.value, result$statistic)
    },
    error = function(e) c(1, 0)
  )
  res
}

#' Statistical test differential expression (t-test or Wilcoxon)
#'
#' @param df Expression matrix (genes x samples)
#' @param map Data frame with sample-to-group mapping
#' @param stat Statistical test: "t" for t-test, "wilcox" for Wilcoxon
#' @return Data frame with columns: gene_id, logFC, Pvalue, Padj, stat
diff_stat <- function(df, map, stat) {
  if (!stat %in% c("t", "wilcox")) {
    stop(paste0("Unsupported statistical test: ", stat,
                ". Use 't' or 'wilcox'."))
  }

  treat_name <- levels(map[[2]])[2]
  con_name <- levels(map[[2]])[1]
  treat <- map[map[[2]] == treat_name, ][[1]]
  control <- map[map[[2]] == con_name, ][[1]]

  res <- data.frame(logFC = apply(df, 1, logfc_row,
    treat = treat, control = control
  ))

  test_func <- if (stat == "t") t.test else wilcox.test
  df_test <- apply(df, 1, test_row, group = map[[2]], func = test_func)
  df_test <- t(df_test)
  colnames(df_test) <- c("Pvalue", "stat")
  res <- cbind(res, df_test)
  res$Padj <- p.adjust(res$Pvalue, method = "BH")

  res <- res %>%
    na.omit() %>%
    dplyr::arrange(Pvalue)

  res$gene_id <- rownames(res)
  rownames(res) <- NULL
  res <- res[, c("gene_id", "logFC", "Pvalue", "Padj", "stat")]
  res
}
