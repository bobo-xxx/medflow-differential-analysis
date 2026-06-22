# exceptions.R — Structured exception handling for deg-analysis
#
# Categories: B=Data, W=Write, E=Environment
# All exceptions flow through report_exception_ndjson() for machine-readable output.

# Exception accumulator for .run_result.json (also initialized in report.R)
if (!exists(".exceptions")) .exceptions <- list()

#' Check runtime environment
#'
#' Verifies required packages are installed.
#'
#' @return List with status ("ok", "error"), missing packages, and details
check_environment <- function() {
  required <- c("DESeq2", "limma", "edgeR", "statmod",
                "ggplot2", "ggrepel", "pheatmap",
                "dplyr", "data.table", "jsonlite",
                "yaml", "filelock")
  missing  <- character(0)
  warnings <- character(0)

  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }

  if (length(missing) > 0) {
    return(list(
      status = "error",
      missing = missing,
      msg = paste("Missing required packages:", paste(missing, collapse = ", "))
    ))
  }

  # Optional packages
  for (pkg in c("ggvenn", "RCircos")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warnings <- c(warnings, paste("Optional package not available:", pkg))
    }
  }

  list(status = "ok", msg = "Environment OK", warnings = warnings)
}

#' Report an exception as structured NDJSON
#'
#' @param code Exception code (e.g., "B1_PROPORTION")
#' @param nature Exception nature (data_insufficient, data_corrupt, data_mismatch, env_bug)
#' @param action Response action (halt, skip_with_warning, escalate)
#' @param msg Human-readable message
#' @param exit_code Integer exit code (default: 1)
#' @param dry_run Logical, if TRUE do not quit (for testing)
report_exception_ndjson <- function(code, nature, action, msg,
                                     exit_code = 1, dry_run = FALSE) {
  level <- switch(action,
    skip_with_warning = "decision",
    escalate          = "exception",
    "exception"
  )

  obj <- list(
    level  = level,
    code   = code,
    nature = nature,
    action = action,
    msg    = msg
  )
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")

  # Accumulate for .run_result.json
  .exceptions <<- c(.exceptions, list(obj))

  if (action == "halt" && !dry_run) {
    quit(status = exit_code)
  }
}
