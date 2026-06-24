#!/usr/bin/env Rscript
#
# main.R — Single entry point for differential-analysis node
#
# Usage:
#   Rscript scripts/main.R run --mat expr.csv --map groups.csv --outdir ./output
#   Rscript scripts/main.R validate-input --mat expr.csv --map groups.csv
#   Rscript scripts/main.R validate-output --outdir ./output
#
# The first positional argument is the subcommand.
# All parameters declared in SKILL.md frontmatter are accepted.
# Output is NDJSON to stdout.

# Resolve script directory for relative sourcing.
# When run standalone (Rscript main.R), commandArgs includes --file=main.R.
# When sourced (e.g., from tests via -e), --file= is absent and we fall back to "."
# — the caller must have sourced dependencies first.
args_with_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(args_with_file_arg) > 0) {
  script_dir <- tryCatch(
    dirname(normalizePath(sub("^--file=", "", args_with_file_arg))),
    warning = function(w) ".",
    error = function(e) "."
  )
} else {
  script_dir <- "."
}

# Source modules in dependency order (only when not already loaded)
if (script_dir != "." || !exists("report_error")) {
  source(file.path(script_dir, "report.R"))
  source(file.path(script_dir, "exceptions.R"))
  source(file.path(script_dir, "io_helpers.R"))
  source(file.path(script_dir, "diff_methods.R"))
  source(file.path(script_dir, "filter_helpers.R"))
  source(file.path(script_dir, "plot_helpers.R"))
  source(file.path(script_dir, "input_validation.R"))
  source(file.path(script_dir, "output_validation.R"))
}

# -------------------------------------------------------------------
# Argument parsing
# -------------------------------------------------------------------

#' Parse command-line arguments
#'
#' Accepts --key=value and --key value forms.
#'
#' @param args Character vector of CLI args (default: commandArgs)
#' @return Named list of parsed values
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) {
    cat("Usage: Rscript scripts/main.R <subcommand> [options]\n")
    cat("\nSubcommands:\n")
    cat("  run               Run differential expression analysis\n")
    cat("  validate-input    Pre-flight input validation\n")
    cat("  validate-output   Post-hoc output validation\n")
    cat("\nRun options:\n")
    cat("  --mat FILE           Expression matrix CSV (genes x samples)\n")
    cat("  --map FILE           Sample-to-group mapping CSV\n")
    cat("  --method METHOD      DE method: deseq2, limma, edgeR, t, wilcox (default: auto-detect)\n")
    cat("  --p-set P            P-value type: p or padj (default: padj)\n")
    cat("  --pvalue FLOAT       P-value threshold (default: 0.05)\n")
    cat("  --logfc-cutoff FLOAT Log2 fold-change cutoff (default: 1.0)\n")
    cat("  --cutoff INT         Minimum DEG count (default: 10)\n")
    cat("  --norm METHOD        edgeR normalization: TMM, RLE, upperquartile, none (default: TMM)\n")
    cat("  --model METHOD       edgeR model: glmFit, glmQLFit (default: glmFit)\n")
    cat("  --top INT            Top N genes for heatmap/volcano (default: 20)\n")
    cat("  --force-imbalanced   Override 10:1 proportion check\n")
    cat("  --rgs FILE           Related gene set for Venn diagram\n")
    cat("  --locate FILE        Chromosome annotation for location plot\n")
    cat("  --tax-id STRING      NCBI taxonomy ID (default: 9606)\n")
    cat("  --pheno-abbr STRING  Phenotype abbreviation for Venn label\n")
    cat("  --gene STRING        Comma-separated genes to label on volcano\n")
    cat("  --color-heat STRING  Heatmap color palette (default: blue,white,red)\n")
    cat("  --color-panel STRING Comma-separated colors for Venn sets\n")
    cat("  --outdir DIR         Output directory (default: .)\n")
    quit(status = 1)
  }

  subcommand <- args[1]
  valid_subcommands <- c("run", "validate-input", "validate-output")
  if (!subcommand %in% valid_subcommands) {
    report_error(sprintf("Unknown subcommand '%s'. Valid: %s",
      subcommand, paste(valid_subcommands, collapse = ", ")))
  }

  opts <- list(
    subcommand        = subcommand,
    mat               = NULL,
    map               = NULL,
    method            = "deseq2",
    method_explicit   = FALSE,
    p_set             = "padj",
    pvalue            = 0.05,
    logfc_cutoff      = 1.0,
    cutoff            = 10,
    norm              = "TMM",
    model             = "glmFit",
    top               = 20,
    force_imbalanced  = FALSE,
    rgs               = NULL,
    locate            = NULL,
    tax_id            = "9606",
    pheno_abbr        = NULL,
    gene              = NULL,
    color_heat        = "blue,white,red",
    color_panel       = NULL,
    outdir            = "."
  )

  # Map of --key -> opts$key for parsing
  param_map <- list(
    "--mat" = "mat", "--map" = "map",
    "--method" = "method", "--p-set" = "p_set",
    "--pvalue" = "pvalue", "--logfc-cutoff" = "logfc_cutoff",
    "--cutoff" = "cutoff", "--norm" = "norm",
    "--model" = "model", "--top" = "top",
    "--rgs" = "rgs", "--locate" = "locate",
    "--tax-id" = "tax_id", "--pheno-abbr" = "pheno_abbr",
    "--gene" = "gene", "--color-heat" = "color_heat",
    "--color-panel" = "color_panel", "--outdir" = "outdir"
  )

  remaining <- args[-1]
  i <- 1
  while (i <= length(remaining)) {
    key <- remaining[i]

    if (key == "--force-imbalanced") {
      opts$force_imbalanced <- TRUE
      i <- i + 1
      next
    }

    found <- FALSE
    for (flag in names(param_map)) {
      opt_name <- param_map[[flag]]

      # --key=value form
      if (startsWith(key, paste0(flag, "="))) {
        val <- sub(paste0("^", flag, "="), "", key)
        opts[[opt_name]] <- type_convert(val, opt_name)
        if (flag == "--method") opts$method_explicit <- TRUE
        found <- TRUE
        break
      }

      # --key value form
      if (key == flag) {
        i <- i + 1
        if (i <= length(remaining)) {
          opts[[opt_name]] <- type_convert(remaining[i], opt_name)
          if (flag == "--method") opts$method_explicit <- TRUE
        }
        found <- TRUE
        break
      }
    }

    if (!found) {
      report_error(sprintf("Unknown option: %s", key))
    }
    i <- i + 1
  }

  # Type conversions for numeric/integer params
  opts$pvalue <- as.numeric(opts$pvalue)
  opts$logfc_cutoff <- as.numeric(opts$logfc_cutoff)
  opts$cutoff <- as.integer(opts$cutoff)
  opts$top <- as.integer(opts$top)

  # Validate required args per subcommand
  if (opts$subcommand == "run") {
    if (is.null(opts$mat)) report_error("run subcommand requires --mat")
    if (is.null(opts$map)) report_error("run subcommand requires --map")
  }
  if (opts$subcommand == "validate-input") {
    if (is.null(opts$mat)) report_error("validate-input subcommand requires --mat")
    if (is.null(opts$map)) report_error("validate-input subcommand requires --map")
  }
  if (opts$subcommand == "validate-output") {
    # outdir defaults to "." so always present
    NULL
  }

  # Validate choice parameters
  valid_methods <- c("deseq2", "limma", "edgeR", "t", "wilcox")
  if (!opts$method %in% valid_methods) {
    report_error(sprintf("Invalid --method '%s'. Valid: %s",
      opts$method, paste(valid_methods, collapse = ", ")))
  }
  valid_p_sets <- c("p", "padj")
  if (!opts$p_set %in% valid_p_sets) {
    report_error(sprintf("Invalid --p-set '%s'. Valid: %s",
      opts$p_set, paste(valid_p_sets, collapse = ", ")))
  }

  return(opts)
}

#' Convert CLI string to appropriate R type
#'
#' @param val Character value from CLI
#' @param opt_name Parameter name for context
#' @return Converted value
type_convert <- function(val, opt_name) {
  if (opt_name %in% c("pvalue", "logfc_cutoff")) {
    return(as.numeric(val))
  }
  if (opt_name %in% c("cutoff", "top")) {
    return(as.integer(val))
  }
  return(val)
}

# -------------------------------------------------------------------
# Subcommand: run
# -------------------------------------------------------------------

#' Detect write error type and emit appropriate W-code exception
#'
#' Inspects the error message for disk-full or permission-denied
#' patterns and emits W001_DISK_FULL or W002_PERM_DENIED accordingly.
#' Calls quit(1) after emission (action: halt).
#'
#' @param err_msg Character error message from tryCatch
#' @param file_path Character path that was being written to
emit_write_exception <- function(err_msg, file_path) {
  if (grepl("disk full|no space|quota exceeded|cannot allocate",
            err_msg, ignore.case = TRUE)) {
    report_exception_ndjson(
      "W001_DISK_FULL", "resource", "halt",
      sprintf("Disk full or no space left on device while writing %s: %s",
              file_path, err_msg),
      exit_code = 1
    )
  } else if (grepl("permission denied|access denied|read.only",
                   err_msg, ignore.case = TRUE)) {
    report_exception_ndjson(
      "W002_PERM_DENIED", "resource", "halt",
      sprintf("Permission denied while writing %s: %s",
              file_path, err_msg),
      exit_code = 1
    )
  } else {
    # Re-throw unrecognized write errors
    stop(err_msg)
  }
}

#' Run the full DEG analysis pipeline
#'
#' Pipeline: load data -> proportion_check -> diff_analysis ->
#' filter_degs -> plot_volcano -> plot_heatmap ->
#' [plot_venn] -> [plot_locate] -> write outputs -> report_result
#'
#' @param opts Named list of parsed arguments
do_run <- function(opts) {
  report_info("Starting DEG analysis pipeline",
              method = opts$method, p_set = opts$p_set,
              pvalue = opts$pvalue, logfc_cutoff = opts$logfc_cutoff)

  # Load data
  report_info(sprintf("Loading expression matrix: %s", opts$mat))
  mat <- data.table::fread(opts$mat, data.table = FALSE)
  gene_ids <- mat[[1]]
  rownames(mat) <- gene_ids
  mat <- mat[, -1, drop = FALSE]

  map <- data.table::fread(opts$map, data.table = FALSE, header = TRUE)
  rownames(map) <- map[[1]]
  map[[2]] <- factor(map[[2]], levels = unique(map[[2]]))

  treat_name <- levels(map[[2]])[2]
  con_name <- levels(map[[2]])[1]
  report_info(sprintf("Groups: %s (n=%d) vs %s (n=%d)",
    treat_name, sum(map[[2]] == treat_name),
    con_name, sum(map[[2]] == con_name)))

  # Proportion check
  report_info("Checking sample proportion...")
  if (!proportion_check(map, opts$force_imbalanced)) {
    # Determine failure reason and emit the appropriate exception
    counts <- table(map[[2]])
    if (length(counts) != 2) {
      report_exception_ndjson(
        "B9_SAMPLE_MISMATCH", "data_corrupt", "halt",
        sprintf("Expected exactly 2 groups, found %d: %s",
                length(counts), paste(names(counts), collapse = ", ")),
        exit_code = 1
      )
      return(invisible(list(status = "error",
        msg = "Expected exactly 2 groups")))
    }
    c_num <- as.numeric(counts[1])
    t_num <- as.numeric(counts[2])
    ratio <- max(c_num, t_num) / min(c_num, t_num)
    g1 <- names(counts)[1]
    g2 <- names(counts)[2]
    report_exception_ndjson(
      "B1_PROPORTION", "data_insufficient", "skip_with_warning",
      sprintf("Sample proportion %s (%d samples) vs %s (%d samples), ratio %.1f:1 exceeds 10:1 limit",
              g1, c_num, g2, t_num, ratio),
      exit_code = 1
    )
    return(invisible(list(status = "error",
      msg = "Sample proportion check failed")))
  }

  # Subset matrix to map samples
  mat <- mat[, map[[1]], drop = FALSE]

  # Auto-detect data type and adjust default method if needed
  is_count_like <- all(mat >= 0, na.rm = TRUE) &&
    all(abs(mat - round(mat)) < 1e-6, na.rm = TRUE) &&
    max(mat, na.rm = TRUE) > 50

  if (!opts$method_explicit && opts$method %in% c("deseq2", "edgeR") &&
      !is_count_like) {
    report_info(sprintf(
      "Data appears to be normalized/log-transformed (non-integer or low-range). Auto-switching method from '%s' to 'limma'.",
      opts$method))
    opts$method <- "limma"
  }

  # Method consistency check: count-based methods need integer data
  if (opts$method %in% c("deseq2", "edgeR")) {
    if (!all(mat == round(mat), na.rm = TRUE)) {
      report_exception_ndjson(
        "B5_METHOD_MISMATCH", "data_mismatch", "halt",
        sprintf("Method '%s' requires integer count data. Matrix contains non-integer values.", opts$method),
        exit_code = 1
      )
      return(invisible(list(status = "error",
        msg = "Method-data mismatch: count-based method requires integer data")))
    }
  }

  # Run DE analysis
  report_info(sprintf("Running %s differential expression...", opts$method))
  dif <- tryCatch(
    switch(opts$method,
      deseq2 = diff_deseq2(mat, map),
      limma  = diff_limma(mat, map),
      edgeR  = diff_edger(mat, map, norm = opts$norm, model = opts$model),
      t      = diff_stat(mat, map, "t"),
      wilcox = diff_stat(mat, map, "wilcox")
    ),
    error = function(e) {
      msg <- sprintf("DE analysis failed (%s): %s", opts$method, e$message)
      report_exception_ndjson("B5_METHOD_MISMATCH", "data_mismatch", "halt",
        msg, exit_code = 1)
      return(NULL)
    }
  )
  if (is.null(dif)) {
    write_run_result(opts$outdir,
      list(status = "error", msg = "DE analysis failed"),
      opts, 1,
      c(started_at, format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")))
    return(invisible(list(status = "error",
      msg = "DE analysis failed")))
  }

  # Determine p-value column name
  p_name <- ifelse(opts$p_set == "p", "Pvalue", "Padj")
  fc_name <- "logFC"

  # Sensitivity test
  logfc_test <- seq(4, 0, -0.5)
  test <- test_cutoff(dif, fc_name, p_name, logfc_test, opts$pvalue)
  report_info("Cutoff sensitivity test",
              cutoffs = as.list(test))

  # Filter DEGs
  rgs_genes <- NULL
  if (!is.null(opts$rgs) && opts$rgs != "None" && file.exists(opts$rgs)) {
    rgs_df <- data.table::fread(opts$rgs, data.table = FALSE)
    rgs_genes <- rgs_df[[1]]
    report_info(sprintf("Related gene set loaded: %d genes", length(rgs_genes)))
  }

  filt <- filter_degs(dif, fc_name, p_name,
    opts$logfc_cutoff, opts$pvalue, opts$cutoff, rgs_genes)

  report_info(sprintf("%d DEGs: %d up, %d down (|logFC|>%.2f, %s<%.4f)",
    filt$total, filt$up, filt$down,
    opts$logfc_cutoff, p_name, opts$pvalue))

  if (filt$total == 0) {
    # B6_NO_DEGS already emitted by filter_degs
    return(invisible(list(status = "error",
      msg = "No differentially expressed genes found")))
  }

  # Create output directory
  if (!dir.exists(opts$outdir)) {
    dir.create(opts$outdir, recursive = TRUE)
  }

  # Write Diffanalysis.csv (full results)
  diff_path <- file.path(opts$outdir, "Diffanalysis.csv")
  tryCatch(
    write.csv(dif, diff_path, row.names = FALSE),
    error = function(e) {
      emit_write_exception(e$message, diff_path)
    }
  )
  report_info(sprintf("Full DE results written to %s (%d rows)",
    basename(diff_path), nrow(dif)))

  # Write DEGs.csv (filtered with group column)
  degs_path <- file.path(opts$outdir, "DEGs.csv")
  tryCatch(
    write.csv(filt$degs, degs_path, row.names = FALSE),
    error = function(e) {
      emit_write_exception(e$message, degs_path)
    }
  )
  report_info(sprintf("Filtered DEGs written to %s (%d rows)",
    basename(degs_path), nrow(filt$degs)))

  # Volcano plot
  volcano_path <- file.path(opts$outdir, "Volcano.pdf")
  tryCatch(
    plot_volcano(dif, p_name, opts$pvalue, opts$logfc_cutoff,
      opts$top, opts$gene, volcano_path),
    error = function(e) {
      emit_write_exception(e$message, volcano_path)
    }
  )

  # Heatmap
  heatmap_path <- file.path(opts$outdir, "Heatmap.pdf")
  tryCatch(
    plot_heatmap(mat, map, filt$rdegs, opts$top, opts$color_heat, heatmap_path),
    error = function(e) {
      emit_write_exception(e$message, heatmap_path)
    }
  )

  # Optional: Venn diagram
  if (!is.null(opts$rgs) && opts$rgs != "None" && file.exists(opts$rgs)) {
    if (!is.null(opts$pheno_abbr) && !is.null(opts$color_panel)) {
      venn_path <- file.path(opts$outdir, "Venn.pdf")
      tryCatch(
        plot_venn(filt$degs, rgs_genes, opts$pheno_abbr,
          opts$color_panel, venn_path),
        error = function(e) {
          emit_write_exception(e$message, venn_path)
        }
      )
    } else {
      report_info("Skipping Venn diagram: --pheno-abbr and --color-panel required")
    }
  }

  # Optional: Chromosome location plot
  if (!is.null(opts$locate) && opts$locate != "None" && file.exists(opts$locate)) {
    locate_path <- file.path(opts$outdir, "Chromosome_location.pdf")
    tryCatch(
      plot_locate(filt$degs[[1]], opts$locate, opts$tax_id, locate_path),
      error = function(e) {
        emit_write_exception(e$message, locate_path)
      }
    )
  }

  # Build output file list
  outfiles <- list(
    list(path = diff_path, rows = nrow(dif), cols = ncol(dif)),
    list(path = degs_path, rows = nrow(filt$degs), cols = ncol(filt$degs)),
    list(path = volcano_path),
    list(path = heatmap_path)
  )
  venn_pdf <- file.path(opts$outdir, "Venn.pdf")
  if (file.exists(venn_pdf)) outfiles <- c(outfiles, list(list(path = venn_pdf)))
  loc_pdf <- file.path(opts$outdir, "Chromosome_location.pdf")
  if (file.exists(loc_pdf)) outfiles <- c(outfiles, list(list(path = loc_pdf)))

  # Get method version
  pkg_name <- switch(opts$method,
    deseq2 = "DESeq2", limma = "limma", edgeR = "edgeR",
    t = "stats", wilcox = "stats")
  pkg_ver <- as.character(packageVersion(pkg_name))

  # Report result
  report_result(status = "success", files = outfiles,
    metadata = list(
      method = opts$method,
      version = pkg_ver,
      n_degs = filt$total,
      up = filt$up,
      down = filt$down
    ))

  list(status = "success",
    metadata = list(
      method = opts$method,
      version = pkg_ver,
      n_degs = filt$total,
      up = filt$up,
      down = filt$down
    ))
}

# -------------------------------------------------------------------
# Subcommand: validate-input
# -------------------------------------------------------------------

do_validate_input <- function(opts) {
  result <- validate_input(opts)
  if (!result$valid) {
    report_exception_ndjson(
      if (grepl("not found|Cannot read", result$reason)) "B3_MISSING_INPUT"
      else if (grepl("column|empty|Expected exactly 2 groups", result$reason)) "B4_INVALID_COLUMNS"
      else if (grepl("proportion", result$reason)) "B1_PROPORTION"
      else if (grepl("Sample.*in (map|matrix)", result$reason)) "B9_SAMPLE_MISMATCH"
      else "B9_SAMPLE_MISMATCH",
      "data_corrupt", "halt",
      result$reason
    )
    return(invisible(list(status = "error", msg = result$reason)))
  }
  report_info(sprintf("Input validation passed: %d genes, %d samples, groups: %s",
    result$n_genes, result$n_samples, result$groups))
  list(status = "ok")
}

# -------------------------------------------------------------------
# Subcommand: validate-output
# -------------------------------------------------------------------

do_validate_output <- function(opts) {
  result <- validate_output(opts)
  if (!result$valid) {
    cat(sprintf("Output validation failed: %s\n", result$reason), file = stderr())
    quit(status = 1)
  }
  info <- paste(vapply(names(result$file_info), function(f) {
    sprintf("%s=%s", f, result$file_info[[f]])
  }, ""), collapse = ", ")
  cat(sprintf("OK: %s\n", info))
  quit(status = 0)
}

# -------------------------------------------------------------------
# Main dispatch
# -------------------------------------------------------------------

main <- function() {
  started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")

  # Reset exception accumulator to prevent state leakage between consecutive calls
  assign(".exceptions", list(), envir = .GlobalEnv)

  opts <- parse_args()

  # Environment check (after arg parsing — help text shouldn't require packages)
  env <- check_environment()
  if (env$status == "error") {
    report_exception_ndjson("E801_ENV_PKG", "env_bug", "halt", env$msg, exit_code = 3)
    quit(status = 3)
  }
  for (w in env$warnings) {
    report_info(sprintf("Env warning: %s", w))
  }

  result <- switch(opts$subcommand,
    "run"              = do_run(opts),
    "validate-input"   = do_validate_input(opts),
    "validate-output"  = do_validate_output(opts)
  )

  finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")

  if (opts$subcommand == "run") {
    exit_code <- if (is.null(result$status) || result$status == "success") 0 else 1
    write_run_result(opts$outdir, result, opts, exit_code,
      c(started_at, finished_at))
    if (exit_code != 0) {
      quit(status = exit_code)
    }
  }
}

if (sys.nframe() == 0) {
  main()
}
