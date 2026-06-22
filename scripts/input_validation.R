#!/usr/bin/env Rscript
#
# input_validation.R — Pre-flight input validation for differential-analysis
#
# Standalone executable. Validates inputs before running the DE pipeline.
# Also sourceable by main.R for the validate-input subcommand.
#
# Usage:
#   Rscript input_validation.R --mat expr.csv --map groups.csv
#   Rscript input_validation.R --mat expr.csv --map groups.csv --force-imbalanced
#
# Exit codes:
#   0 — all checks passed
#   1 — validation failed (stderr has reason)
#   2 — usage/argument error

# Determine script directory when running standalone.
# Falls back to "." when sourced (caller is responsible for sourcing deps).
args_with_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(args_with_file_arg) > 0) {
  script_dir <- tryCatch(
    dirname(normalizePath(sub("^--file=", "", args_with_file_arg))),
    warning = function(w) ".",
    error = function(e) "."
  )
} else {
  # When sourced (e.g. from tests or main.R), assume dependencies are already loaded.
  # If they are not, the caller must source report.R, exceptions.R, and
  # filter_helpers.R before sourcing this file.
  script_dir <- "."
}

if (script_dir != "." || !exists("report_error")) {
  source(file.path(script_dir, "report.R"))
  source(file.path(script_dir, "exceptions.R"))
  source(file.path(script_dir, "filter_helpers.R"))
}

# -------------------------------------------------------------------
# Validation function
# -------------------------------------------------------------------

#' Validate inputs for DEG analysis
#'
#' Checks file existence, CSV readability, required columns, sample
#' ID matching, exactly 2 groups, and sample proportion ratio <= 10:1.
#' Emits structured exceptions via report_exception_ndjson for
#' proportion and group-count violations (via filter_helpers::proportion_check).
#'
#' @param opts List with mat (file path), map (file path), force_imbalanced (logical)
#' @return List with valid (logical), reason (character), and optionally n_genes, n_samples
validate_input <- function(opts) {
  # ---- File existence checks ----
  if (is.null(opts$mat) || !file.exists(opts$mat)) {
    return(list(valid = FALSE,
      reason = sprintf("Expression matrix file not found: %s",
                       if (is.null(opts$mat)) "(null)" else opts$mat)))
  }
  if (is.null(opts$map) || !file.exists(opts$map)) {
    return(list(valid = FALSE,
      reason = sprintf("Sample map file not found: %s",
                       if (is.null(opts$map)) "(null)" else opts$map)))
  }

  # ---- Read expression matrix ----
  mat <- tryCatch(
    read.csv(opts$mat, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    return(list(valid = FALSE,
      reason = sprintf("Cannot read expression matrix as CSV: %s", opts$mat)))
  }

  # ---- Check for gene identifier column (first column) ----
  if (ncol(mat) < 2) {
    return(list(valid = FALSE,
      reason = sprintf("Expression matrix has only %d column(s). Need gene ID column + samples.",
                       ncol(mat))))
  }

  # ---- Check for empty matrix ----
  if (nrow(mat) == 0) {
    return(list(valid = FALSE,
      reason = "Expression matrix has 0 rows (empty matrix)"))
  }

  # ---- Read sample group map ----
  map <- tryCatch(
    read.csv(opts$map, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(map)) {
    return(list(valid = FALSE,
      reason = sprintf("Cannot read sample map as CSV: %s", opts$map)))
  }

  if (ncol(map) < 2) {
    return(list(valid = FALSE,
      reason = sprintf("Sample map has only %d column(s). Need sample ID + group columns.",
                       ncol(map))))
  }

  # ---- Check sample columns in mat match sample IDs in map ----
  mat_samples <- colnames(mat)[-1]  # first col is gene ID
  map_samples <- map[[1]]
  missing_in_mat <- setdiff(map_samples, mat_samples)
  missing_in_map <- setdiff(mat_samples, map_samples)

  if (length(missing_in_mat) > 0) {
    return(list(valid = FALSE,
      reason = sprintf("Sample(s) in map but not in matrix: %s",
                       paste(missing_in_mat, collapse = ", "))))
  }
  if (length(missing_in_map) > 0) {
    return(list(valid = FALSE,
      reason = sprintf("Sample(s) in matrix but not in map: %s",
                       paste(missing_in_map, collapse = ", "))))
  }

  # ---- Group structure: exactly 2 groups ----
  groups <- unique(map[[2]])
  if (length(groups) != 2) {
    return(list(valid = FALSE,
      reason = sprintf("Expected exactly 2 groups, found %d: %s",
                       length(groups), paste(groups, collapse = ", "))))
  }

  # ---- Make factor for proportion_check ----
  map[[2]] <- factor(map[[2]], levels = groups)

  # ---- Proportion check ----
  force_imb <- isTRUE(opts$force_imbalanced) ||
               identical(opts$force_imbalanced, "TRUE") ||
               identical(opts$force_imbalanced, "true")
  if (!proportion_check(map, force_imb)) {
    return(list(valid = FALSE,
      reason = "Sample proportion check failed (ratio > 10:1, use --force-imbalanced to override)"))
  }

  list(valid = TRUE, reason = "OK",
       n_genes = nrow(mat), n_samples = length(map_samples),
       groups = paste(groups, collapse = " vs "))
}

# -------------------------------------------------------------------
# CLI argument parsing
# -------------------------------------------------------------------

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) {
    cat("Usage: Rscript input_validation.R --mat <expr.csv> --map <groups.csv> [--force-imbalanced]\n")
    cat("\nOptions:\n")
    cat("  --mat FILE           Expression matrix CSV (genes x samples)\n")
    cat("  --map FILE           Sample-to-group mapping CSV\n")
    cat("  --force-imbalanced   Override 10:1 proportion check\n")
    quit(status = 2)
  }

  opts <- list(mat = NULL, map = NULL, force_imbalanced = FALSE)
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--mat") {
      i <- i + 1; if (i <= length(args)) opts$mat <- args[i]
    } else if (key == "--map") {
      i <- i + 1; if (i <= length(args)) opts$map <- args[i]
    } else if (key == "--force-imbalanced") {
      opts$force_imbalanced <- TRUE
    } else if (startsWith(key, "--mat=")) {
      opts$mat <- sub("^--mat=", "", key)
    } else if (startsWith(key, "--map=")) {
      opts$map <- sub("^--map=", "", key)
    } else if (startsWith(key, "--force-imbalanced=")) {
      opts$force_imbalanced <- TRUE
    } else {
      cat(sprintf("Error: unknown option: %s\n", key), file = stderr())
      quit(status = 2)
    }
    i <- i + 1
  }

  if (is.null(opts$mat)) {
    cat("Error: --mat is required\n", file = stderr())
    quit(status = 2)
  }
  if (is.null(opts$map)) {
    cat("Error: --map is required\n", file = stderr())
    quit(status = 2)
  }
  opts
}

# -------------------------------------------------------------------
# Main dispatch (standalone mode)
# -------------------------------------------------------------------

if (sys.nframe() == 0) {
  opts <- parse_args()
  result <- validate_input(opts)
  if (!result$valid) {
    cat(sprintf("Validation failed: %s\n", result$reason), file = stderr())
    quit(status = 1)
  }
  cat(sprintf("OK: %d genes x %d samples, groups: %s\n",
              result$n_genes, result$n_samples, result$groups))
  quit(status = 0)
}
