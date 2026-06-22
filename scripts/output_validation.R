#!/usr/bin/env Rscript
#
# output_validation.R — Post-hoc output validation for deg-analysis
#
# Standalone executable. Validates outputs after running the DE pipeline.
# Also sourceable by main.R for the validate-output subcommand.
#
# Checks: Diffanalysis.csv existence and required columns (gene_id, logFC,
# Pvalue, Padj), DEGs.csv existence and group column, Volcano.pdf and
# Heatmap.pdf existence and non-zero size.
#
# Usage:
#   Rscript output_validation.R --outdir ./output
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
  # If they are not, the caller must source report.R and exceptions.R before
  # sourcing this file.
  script_dir <- "."
}

if (script_dir != "." || !exists("report_error")) {
  source(file.path(script_dir, "report.R"))
  source(file.path(script_dir, "exceptions.R"))
}

# -------------------------------------------------------------------
# Validation function
# -------------------------------------------------------------------

#' Validate outputs from DEG analysis
#'
#' Checks output CSV column presence (gene_id, logFC, Pvalue, Padj),
#' non-empty results, and plot file existence with non-zero size.
#'
#' @param opts List with outdir (output directory path)
#' @return List with valid (logical), reason (character), file_info (named list)
validate_output <- function(opts) {
  outdir <- opts$outdir
  if (is.null(outdir) || !dir.exists(outdir)) {
    return(list(valid = FALSE,
      reason = sprintf("Output directory not found or not specified: %s",
                       if (is.null(outdir)) "(null)" else outdir)))
  }

  # ---- Check Diffanalysis.csv ----
  diff_csv <- file.path(outdir, "Diffanalysis.csv")
  if (!file.exists(diff_csv)) {
    return(list(valid = FALSE,
      reason = sprintf("Missing Diffanalysis.csv in %s", outdir)))
  }

  dif <- tryCatch(
    read.csv(diff_csv, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(dif)) {
    return(list(valid = FALSE,
      reason = sprintf("Cannot read Diffanalysis.csv: %s", diff_csv)))
  }

  required_cols <- c("gene_id", "logFC", "Pvalue", "Padj")
  missing_cols <- setdiff(required_cols, colnames(dif))
  if (length(missing_cols) > 0) {
    return(list(valid = FALSE,
      reason = sprintf("Diffanalysis.csv missing required columns: %s",
                       paste(missing_cols, collapse = ", "))))
  }

  if (nrow(dif) == 0) {
    return(list(valid = FALSE,
      reason = "Diffanalysis.csv has 0 rows (no genes in result)"))
  }

  diff_rows <- nrow(dif)
  diff_cols <- ncol(dif)

  # ---- Check DEGs.csv ----
  degs_csv <- file.path(outdir, "DEGs.csv")
  if (!file.exists(degs_csv)) {
    return(list(valid = FALSE,
      reason = sprintf("Missing DEGs.csv in %s", outdir)))
  }

  degs <- tryCatch(
    read.csv(degs_csv, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(degs)) {
    return(list(valid = FALSE,
      reason = sprintf("Cannot read DEGs.csv: %s", degs_csv)))
  }

  if (!"group" %in% colnames(degs)) {
    return(list(valid = FALSE,
      reason = "DEGs.csv missing 'group' column"))
  }

  degs_rows <- nrow(degs)

  # ---- Check plot files ----
  required_plots <- c("Volcano.pdf", "Heatmap.pdf")
  for (pf in required_plots) {
    ppath <- file.path(outdir, pf)
    if (!file.exists(ppath)) {
      return(list(valid = FALSE,
        reason = sprintf("Missing plot file: %s", pf)))
    }
    if (file.info(ppath)$size == 0) {
      return(list(valid = FALSE,
        reason = sprintf("Plot file is empty: %s", pf)))
    }
  }

  file_info <- list(
    diff_csv_rows = diff_rows,
    diff_csv_cols = diff_cols,
    degs_csv_rows = degs_rows
  )

  list(valid = TRUE, reason = "OK", file_info = file_info)
}

# -------------------------------------------------------------------
# CLI argument parsing
# -------------------------------------------------------------------

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) {
    cat("Usage: Rscript output_validation.R --outdir <path>\n")
    cat("\nOptions:\n")
    cat("  --outdir DIR   Output directory containing DEG results\n")
    quit(status = 2)
  }

  opts <- list(outdir = NULL)
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--outdir") {
      i <- i + 1; if (i <= length(args)) opts$outdir <- args[i]
    } else if (startsWith(key, "--outdir=")) {
      opts$outdir <- sub("^--outdir=", "", key)
    } else {
      cat(sprintf("Error: unknown option: %s\n", key), file = stderr())
      quit(status = 2)
    }
    i <- i + 1
  }

  if (is.null(opts$outdir)) {
    cat("Error: --outdir is required\n", file = stderr())
    quit(status = 2)
  }
  opts
}

# -------------------------------------------------------------------
# Main dispatch (standalone mode)
# -------------------------------------------------------------------

if (sys.nframe() == 0) {
  opts <- parse_args()
  result <- validate_output(opts)
  if (!result$valid) {
    cat(sprintf("Validation failed: %s\n", result$reason), file = stderr())
    quit(status = 1)
  }
  info <- paste(vapply(names(result$file_info), function(f) {
    sprintf("%s=%s", f, result$file_info[[f]])
  }, ""), collapse = ", ")
  cat(sprintf("OK: %s\n", info))
  quit(status = 0)
}
