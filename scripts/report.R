# report.R — NDJSON reporting helpers for deg-analysis
#
# All output to stdout is valid NDJSON. Each function writes one line.
# Use report_info() for progress, report_result() for final output,
# and report_error() for terminal failures.

# Exception accumulator (filled by report_exception_ndjson during run)
.exceptions <- list()

#' Write an info-level NDJSON message to stdout
#'
#' @param msg Character string with progress message
#' @param ... Additional named fields to include in the JSON object
report_info <- function(msg, ...) {
  extra <- list(...)
  obj <- c(list(level = "info", msg = msg), extra)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}

#' Write a result-level NDJSON message to stdout
#'
#' @param status Character status code ("success", "error", etc.)
#' @param files List of file info (each with path, and optionally rows/cols)
#' @param metadata List of metadata (method, version, n_degs, up, down)
#' @param ... Additional named fields
report_result <- function(status, files = list(), metadata = list(), ...) {
  extra <- list(...)
  obj <- c(list(level = "result", status = status,
                files = files, metadata = metadata), extra)
  if (length(files) == 0) obj$files <- NULL
  if (length(metadata) == 0) obj$metadata <- NULL
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}

#' Write run provenance file for agent inspection
#'
#' Creates .run_result.json in the output directory with full provenance:
#' parameters, output metadata, exceptions, file list, timestamps.
#'
#' @param out_dir Output directory
#' @param result Result list from do_run()
#' @param params Parameter list
#' @param exit_code Integer exit code
#' @param times Character vector of c(started_at, finished_at) ISO8601 timestamps
write_run_result <- function(out_dir, result, params, exit_code, times) {
  output <- list()
  if (!is.null(result$metadata) && length(result$metadata) > 0) {
    for (name in names(result$metadata)) {
      output[[name]] <- result$metadata[[name]]
    }
  }

  files <- list()
  for (fname in c("Diffanalysis.csv", "DEGs.csv", "Volcano.pdf",
                   "Heatmap.pdf", "Venn.pdf", "Chromosome_location.pdf")) {
    fpath <- file.path(out_dir, fname)
    if (file.exists(fpath)) {
      finfo <- list(path = fpath)
      if (endsWith(fname, ".csv")) {
        df <- tryCatch(read.csv(fpath), error = function(e) NULL)
        if (!is.null(df)) {
          finfo$rows <- nrow(df)
          finfo$cols <- ncol(df)
        }
      }
      files <- c(files, list(finfo))
    }
  }

  clean_params <- list()
  for (k in names(params)) {
    if (!is.null(params[[k]]) && k != "outdir") {
      clean_params[[k]] <- params[[k]]
    }
  }

  obj <- list(
    node         = "deg-analysis",
    subcommand   = if (is.null(params$subcommand)) "unknown" else params$subcommand,
    status       = if (is.null(result$status)) "unknown" else result$status,
    exit_code    = exit_code,
    started_at   = times[1],
    finished_at  = times[2],
    parameters   = clean_params,
    output       = output,
    exceptions   = if (length(.exceptions) > 0) .exceptions else list(),
    files        = files
  )

  json_path <- file.path(out_dir, ".run_result.json")
  writeLines(jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE), json_path)
}

#' Write an error-level NDJSON message and exit
#'
#' @param msg Error message
#' @param exit_code Integer exit code (default: 1)
report_error <- function(msg, exit_code = 1) {
  obj <- list(level = "error", msg = msg)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
  quit(status = exit_code)
}
