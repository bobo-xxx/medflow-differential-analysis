# filter_helpers.R — Proportion check, cutoff testing, DEG filtering for deg-analysis
#
# All functions report exceptions via report_exception_ndjson() when thresholds
# are violated, using the structured exception codes from the design doc.

suppressPackageStartupMessages(library(dplyr))

#' Check sample proportion ratio
#'
#' Ensures the ratio of larger group to smaller group does not exceed 10:1
#' unless force_imbalanced is TRUE. Emits B1_PROPORTION exception on violation.
#'
#' @param map Data frame with sample-to-group mapping (col 2 = group)
#' @param force_imbalanced Logical, if TRUE bypass the check
#' @return Logical TRUE if proportion is acceptable (or forced)
proportion_check <- function(map, force_imbalanced = FALSE) {
  counts <- table(map[[2]])
  if (length(counts) != 2) {
    report_exception_ndjson(
      "B9_SAMPLE_MISMATCH", "data_corrupt", "halt",
      sprintf("Expected exactly 2 groups, found %d: %s",
              length(counts), paste(names(counts), collapse = ", ")),
      exit_code = 1
    )
    return(FALSE)
  }

  c_num <- as.numeric(counts[1])
  t_num <- as.numeric(counts[2])
  ratio <- max(c_num, t_num) / min(c_num, t_num)

  if (ratio > 10 && !force_imbalanced) {
    g1 <- names(counts)[1]
    g2 <- names(counts)[2]
    report_exception_ndjson(
      "B1_PROPORTION", "data_insufficient", "skip_with_warning",
      sprintf("Sample proportion %s (%d samples) vs %s (%d samples), ratio %.1f:1 exceeds 10:1 limit",
              g1, c_num, g2, t_num, ratio),
      exit_code = 1
    )
    return(FALSE)
  }
  TRUE
}

#' Test sensitivity of DEG count at different logFC cutoffs
#'
#' @param dif Data frame with DE results (must contain fc_name and p_name columns)
#' @param fc_name Character, name of fold-change column (e.g., "logFC")
#' @param p_name Character, name of p-value column (e.g., "Padj" or "Pvalue")
#' @param logfc_test Numeric vector of logFC cutoffs to test
#' @param p_value Numeric, p-value threshold
#' @return Named integer vector with DEG counts at each cutoff
test_cutoff <- function(dif, fc_name, p_name, logfc_test, p_value = 0.05) {
  test <- c()
  for (value in logfc_test) {
    test <- append(
      test,
      sum((abs(dif[[fc_name]]) > value) & (dif[[p_name]] < p_value))
    )
  }
  names(test) <- logfc_test
  return(test)
}

#' Filter DEGs by logFC and p-value cutoffs
#'
#' Classifies genes as Up, Down, or Not. Emits B6_NO_DEGS or B2_FEW_DEGS
#' if results are empty or below cutoff.
#'
#' @param dif Data frame with DE results (must contain fc_name and p_name columns,
#'   and gene identifiers in column 1)
#' @param fc_name Character, name of fold-change column (e.g., "logFC")
#' @param p_name Character, name of p-value column (e.g., "Padj")
#' @param logfc_cutoff Numeric, absolute logFC threshold
#' @param p_value Numeric, p-value threshold
#' @param cutoff Integer, minimum required DEG count
#' @param rgs Optional character vector of related gene set identifiers for intersection
#' @return List with elements: degs (DEG data frame), rdegs (related DEGs),
#'   dif_grouped (full data frame with group column), up (count), down (count), total (count)
filter_degs <- function(dif, fc_name, p_name, logfc_cutoff, p_value, cutoff, rgs = NULL) {
  dif_up <- dif %>%
    filter(!!sym(fc_name) > logfc_cutoff & !!sym(p_name) < p_value)
  dif_down <- dif %>%
    filter(!!sym(fc_name) < -logfc_cutoff & !!sym(p_name) < p_value)

  gene_col <- colnames(dif)[1]
  dif_grouped <- dif %>% mutate(group = case_when(
    dif[[gene_col]] %in% dif_up[[gene_col]] ~ "Up",
    dif[[gene_col]] %in% dif_down[[gene_col]] ~ "Down",
    TRUE ~ "Not"
  ))

  degs <- dif_grouped[dif_grouped[["group"]] == "Up" | dif_grouped[["group"]] == "Down", ]
  n_degs <- nrow(degs)

  if (n_degs == 0) {
    report_exception_ndjson(
      "B6_NO_DEGS", "data_insufficient", "skip_with_warning",
      sprintf("No differentially expressed genes found at |logFC|>%.2f, %s<%.4f",
              logfc_cutoff, p_name, p_value),
      exit_code = 1
    )
    return(list(degs = degs, rdegs = data.frame(), dif_grouped = dif_grouped,
                up = 0, down = 0, total = 0))
  }

  if (n_degs < cutoff) {
    report_exception_ndjson(
      "B2_FEW_DEGS", "data_insufficient", "skip_with_warning",
      sprintf("DEG count %d below cutoff %d at |logFC|>%.2f, %s<%.4f",
              n_degs, cutoff, logfc_cutoff, p_name, p_value),
      exit_code = 1
    )
  }

  n_up <- nrow(dif_up)
  n_down <- nrow(dif_down)

  if (!is.null(rgs) && length(rgs) > 0) {
    genes <- intersect(degs[[gene_col]], rgs)
    rdegs <- degs[degs[[gene_col]] %in% genes, ]
    if (nrow(rdegs) < cutoff) {
      report_exception_ndjson(
        "B7_RGS_INTERSECTION", "data_insufficient", "skip_with_warning",
        sprintf("Intersection of DEGs with related gene set (%d genes) below cutoff %d",
                nrow(rdegs), cutoff),
        exit_code = 1
      )
    }
  } else {
    rdegs <- degs
  }

  list(degs = degs, rdegs = rdegs, dif_grouped = dif_grouped,
       up = n_up, down = n_down, total = n_degs)
}
