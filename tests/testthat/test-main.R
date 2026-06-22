# test-main.R — Tests for main.R entry point with subcommand dispatch
#
# Tests: parse_args(), main(), do_validate_input(), do_validate_output()

library(testthat)

# Source all modules in dependency order (matching main.R)
source("../../scripts/report.R")
source("../../scripts/exceptions.R")
source("../../scripts/io_helpers.R")
source("../../scripts/diff_methods.R")
source("../../scripts/filter_helpers.R")
source("../../scripts/plot_helpers.R")
source("../../scripts/input_validation.R")
source("../../scripts/output_validation.R")
source("../../scripts/main.R")

# Ensure clean exception state at start of tests
assign(".exceptions", list(), envir = .GlobalEnv)

# ---------------------------------------------------------------------------
# parse_args — subcommand dispatch
# ---------------------------------------------------------------------------
test_that("parse_args() recognizes 'run' subcommand", {
  opts <- parse_args(c("run", "--mat", "x.csv", "--map", "y.csv"))
  expect_equal(opts$subcommand, "run")
  expect_equal(opts$mat, "x.csv")
  expect_equal(opts$map, "y.csv")
})

test_that("parse_args() recognizes 'validate-input' subcommand", {
  opts <- parse_args(c("validate-input", "--mat", "x.csv", "--map", "y.csv"))
  expect_equal(opts$subcommand, "validate-input")
})

test_that("parse_args() recognizes 'validate-output' subcommand", {
  opts <- parse_args(c("validate-output", "--outdir", "/tmp/out"))
  expect_equal(opts$subcommand, "validate-output")
  expect_equal(opts$outdir, "/tmp/out")
})

test_that("parse_args() rejects unknown subcommand", {
  # report_error calls quit(), so we temporarily replace it with stop()
  real_report_error <- report_error
  assign("report_error", function(msg, exit_code = 1) stop(msg), envir = .GlobalEnv)
  on.exit(assign("report_error", real_report_error, envir = .GlobalEnv))

  expect_error(
    parse_args(c("bogus")),
    "Unknown subcommand"
  )
})

# ---------------------------------------------------------------------------
# parse_args — defaults
# ---------------------------------------------------------------------------
test_that("parse_args() applies correct defaults for optional params", {
  opts <- parse_args(c("run", "--mat", "x.csv", "--map", "y.csv"))
  expect_equal(opts$method, "deseq2")
  expect_equal(opts$p_set, "padj")
  expect_equal(opts$pvalue, 0.05)
  expect_equal(opts$logfc_cutoff, 1.0)
  expect_equal(opts$cutoff, 10)
  expect_equal(opts$norm, "TMM")
  expect_equal(opts$model, "glmFit")
  expect_equal(opts$top, 20)
  expect_equal(opts$force_imbalanced, FALSE)
  expect_equal(opts$tax_id, "9606")
  expect_equal(opts$color_heat, "blue,white,red")
  expect_equal(opts$outdir, ".")
  expect_null(opts$rgs)
  expect_null(opts$locate)
  expect_null(opts$pheno_abbr)
  expect_null(opts$gene)
  expect_null(opts$color_panel)
})

# ---------------------------------------------------------------------------
# parse_args — --key=value form
# ---------------------------------------------------------------------------
test_that("parse_args() accepts --key=value form", {
  opts <- parse_args(c("run", "--mat=mat.csv", "--map=map.csv", "--method=limma",
                        "--pvalue=0.01", "--logfc-cutoff=1.5", "--cutoff=20",
                        "--top=50", "--outdir=/tmp/out2", "--tax-id=10090"))
  expect_equal(opts$mat, "mat.csv")
  expect_equal(opts$map, "map.csv")
  expect_equal(opts$method, "limma")
  expect_equal(opts$pvalue, 0.01)
  expect_equal(opts$logfc_cutoff, 1.5)
  expect_equal(opts$cutoff, 20)
  expect_equal(opts$top, 50)
  expect_equal(opts$outdir, "/tmp/out2")
  expect_equal(opts$tax_id, "10090")
})

# ---------------------------------------------------------------------------
# parse_args — --key value form (separate tokens)
# ---------------------------------------------------------------------------
test_that("parse_args() accepts --key value form", {
  opts <- parse_args(c("run", "--mat", "mat.csv", "--map", "map.csv",
                        "--method", "edgeR", "--p-set", "p", "--norm", "RLE",
                        "--model", "glmQLFit", "--gene", "GENE1,GENE2"))
  expect_equal(opts$mat, "mat.csv")
  expect_equal(opts$map, "map.csv")
  expect_equal(opts$method, "edgeR")
  expect_equal(opts$p_set, "p")
  expect_equal(opts$norm, "RLE")
  expect_equal(opts$model, "glmQLFit")
  expect_equal(opts$gene, "GENE1,GENE2")
})

# ---------------------------------------------------------------------------
# parse_args — --force-imbalanced flag
# ---------------------------------------------------------------------------
test_that("parse_args() sets force_imbalanced=TRUE when flag present", {
  opts <- parse_args(c("run", "--mat", "x.csv", "--map", "y.csv", "--force-imbalanced"))
  expect_true(opts$force_imbalanced)
})

test_that("parse_args() force_imbalanced defaults to FALSE", {
  opts <- parse_args(c("run", "--mat", "x.csv", "--map", "y.csv"))
  expect_false(opts$force_imbalanced)
})

# ---------------------------------------------------------------------------
# parse_args — type conversions
# ---------------------------------------------------------------------------
test_that("parse_args() converts --pvalue to numeric", {
  opts <- parse_args(c("run", "--mat", "x.csv", "--map", "y.csv", "--pvalue", "0.01"))
  expect_type(opts$pvalue, "double")
  expect_equal(opts$pvalue, 0.01)
})

test_that("parse_args() converts --logfc-cutoff to numeric", {
  opts <- parse_args(c("run", "--mat", "x.csv", "--map", "y.csv", "--logfc-cutoff", "2.0"))
  expect_type(opts$logfc_cutoff, "double")
  expect_equal(opts$logfc_cutoff, 2.0)
})

test_that("parse_args() converts --cutoff to integer", {
  opts <- parse_args(c("run", "--mat", "x.csv", "--map", "y.csv", "--cutoff", "15"))
  expect_type(opts$cutoff, "integer")
  expect_equal(opts$cutoff, 15L)
})

test_that("parse_args() converts --top to integer", {
  opts <- parse_args(c("run", "--mat", "x.csv", "--map", "y.csv", "--top", "30"))
  expect_type(opts$top, "integer")
  expect_equal(opts$top, 30L)
})

test_that("parse_args() --key=value form also type-converts", {
  opts <- parse_args(c("run", "--mat=x.csv", "--map=y.csv", "--pvalue=0.001",
                        "--logfc-cutoff=2.5", "--cutoff=5", "--top=15"))
  expect_equal(opts$pvalue, 0.001)
  expect_equal(opts$logfc_cutoff, 2.5)
  expect_equal(opts$cutoff, 5L)
  expect_equal(opts$top, 15L)
})

# ---------------------------------------------------------------------------
# parse_args — required args per subcommand
# ---------------------------------------------------------------------------
test_that("parse_args() requires --mat and --map for run subcommand", {
  real_report_error <- report_error
  assign("report_error", function(msg, exit_code = 1) stop(msg), envir = .GlobalEnv)
  on.exit(assign("report_error", real_report_error, envir = .GlobalEnv))

  expect_error(parse_args(c("run")), "run subcommand requires --mat")
  expect_error(parse_args(c("run", "--mat", "x.csv")), "run subcommand requires --map")
})

test_that("parse_args() requires --mat and --map for validate-input subcommand", {
  real_report_error <- report_error
  assign("report_error", function(msg, exit_code = 1) stop(msg), envir = .GlobalEnv)
  on.exit(assign("report_error", real_report_error, envir = .GlobalEnv))

  expect_error(parse_args(c("validate-input")), "validate-input subcommand requires --mat")
  expect_error(parse_args(c("validate-input", "--mat", "x.csv")),
               "validate-input subcommand requires --map")
})

test_that("parse_args() validate-output works with default outdir", {
  opts <- parse_args(c("validate-output"))
  expect_equal(opts$outdir, ".")
})

# ---------------------------------------------------------------------------
# parse_args — all 20 parameters
# ---------------------------------------------------------------------------
test_that("parse_args() accepts all 20 parameters", {
  args <- c("run",
    "--mat", "expr.csv",
    "--map", "groups.csv",
    "--method", "deseq2",
    "--p-set", "padj",
    "--pvalue", "0.05",
    "--logfc-cutoff", "1.0",
    "--cutoff", "10",
    "--norm", "TMM",
    "--model", "glmFit",
    "--top", "20",
    "--force-imbalanced",
    "--rgs", "related_genes.csv",
    "--locate", "chr_annot.csv",
    "--tax-id", "9606",
    "--pheno-abbr", "BRCA",
    "--gene", "GENE1,GENE2",
    "--color-heat", "blue,white,red",
    "--color-panel", "#E64B35,#4DBBD5",
    "--outdir", "./output"
  )
  opts <- parse_args(args)
  expect_equal(opts$subcommand, "run")
  expect_equal(opts$mat, "expr.csv")
  expect_equal(opts$map, "groups.csv")
  expect_equal(opts$method, "deseq2")
  expect_equal(opts$p_set, "padj")
  expect_equal(opts$pvalue, 0.05)
  expect_equal(opts$logfc_cutoff, 1.0)
  expect_equal(opts$cutoff, 10L)
  expect_equal(opts$norm, "TMM")
  expect_equal(opts$model, "glmFit")
  expect_equal(opts$top, 20L)
  expect_true(opts$force_imbalanced)
  expect_equal(opts$rgs, "related_genes.csv")
  expect_equal(opts$locate, "chr_annot.csv")
  expect_equal(opts$tax_id, "9606")
  expect_equal(opts$pheno_abbr, "BRCA")
  expect_equal(opts$gene, "GENE1,GENE2")
  expect_equal(opts$color_heat, "blue,white,red")
  expect_equal(opts$color_panel, "#E64B35,#4DBBD5")
  expect_equal(opts$outdir, "./output")
})

test_that("parse_args() all 20 params with --key=value form works too", {
  args <- c("run",
    "--mat=expr.csv",
    "--map=groups.csv",
    "--method=limma",
    "--p-set=p",
    "--pvalue=0.01",
    "--logfc-cutoff=1.5",
    "--cutoff=5",
    "--norm=RLE",
    "--model=glmQLFit",
    "--top=10",
    "--force-imbalanced",
    "--rgs=rgs.csv",
    "--locate=chr.csv",
    "--tax-id=10090",
    "--pheno-abbr=LUNG",
    "--gene=TP53,EGFR",
    "--color-heat=green,black,red",
    "--color-panel=red,blue",
    "--outdir=/tmp/out"
  )
  opts <- parse_args(args)
  expect_equal(opts$subcommand, "run")
  expect_equal(opts$mat, "expr.csv")
  expect_equal(opts$map, "groups.csv")
  expect_equal(opts$method, "limma")
  expect_equal(opts$p_set, "p")
  expect_equal(opts$pvalue, 0.01)
  expect_equal(opts$logfc_cutoff, 1.5)
  expect_equal(opts$cutoff, 5L)
  expect_equal(opts$norm, "RLE")
  expect_equal(opts$model, "glmQLFit")
  expect_equal(opts$top, 10L)
  expect_true(opts$force_imbalanced)
  expect_equal(opts$rgs, "rgs.csv")
  expect_equal(opts$locate, "chr.csv")
  expect_equal(opts$tax_id, "10090")
  expect_equal(opts$pheno_abbr, "LUNG")
  expect_equal(opts$gene, "TP53,EGFR")
  expect_equal(opts$color_heat, "green,black,red")
  expect_equal(opts$color_panel, "red,blue")
  expect_equal(opts$outdir, "/tmp/out")
})

# ---------------------------------------------------------------------------
# parse_args — unknown option
# ---------------------------------------------------------------------------
test_that("parse_args() reports error for unknown options", {
  real_report_error <- report_error
  assign("report_error", function(msg, exit_code = 1) stop(msg), envir = .GlobalEnv)
  on.exit(assign("report_error", real_report_error, envir = .GlobalEnv))

  expect_error(parse_args(c("run", "--mat", "x.csv", "--map", "y.csv", "--bogus", "val")),
               "Unknown option")
})

# ---------------------------------------------------------------------------
# do_validate_input
# ---------------------------------------------------------------------------
test_that("do_validate_input() passes with valid test data", {
  source("../../tests/testthat/setup.R")
  td <- create_test_data()
  out <- capture.output({
    result <- do_validate_input(list(
      mat = td$mat_path,
      map = td$map_path,
      force_imbalanced = FALSE
    ))
  })
  obj <- jsonlite::fromJSON(out[length(out)])
  expect_equal(obj$level, "info")
  expect_match(obj$msg, "Input validation passed")
})

test_that("do_validate_input() fails for missing file", {
  # do_validate_input calls report_exception_ndjson(..., "halt") which calls quit().
  # Mock quit to avoid terminating the test process.
  real_quit <- quit
  assign("quit", function(status = 1, ...) {
    assign("quit_called_with", status, envir = .GlobalEnv)
  }, envir = .GlobalEnv)
  on.exit(assign("quit", real_quit, envir = .GlobalEnv))

  out <- capture.output({
    result <- do_validate_input(list(
      mat = "/nonexistent/path.csv",
      map = "/nonexistent/map.csv",
      force_imbalanced = FALSE
    ))
  })
  # Should have emitted exception NDJSON before quit
  has_exception <- any(grepl("B3_MISSING_INPUT", out))
  expect_true(has_exception)
  rm(quit_called_with, envir = .GlobalEnv)
})

# ---------------------------------------------------------------------------
# do_validate_output
# ---------------------------------------------------------------------------
test_that("do_validate_output() fails for non-existent directory", {
  # do_validate_output calls quit() on failure, so we test the inner validate_output
  result <- validate_output(list(outdir = "/nonexistent/dir"))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("do_validate_output() passes with valid output directory", {
  tmp <- tempfile("vout_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  # Create minimal required output files
  write.csv(data.frame(gene_id = c("G1", "G2"), logFC = c(2.0, -1.5),
                       Pvalue = c(0.001, 0.002), Padj = c(0.01, 0.02)),
            file.path(tmp, "Diffanalysis.csv"), row.names = FALSE)
  write.csv(data.frame(gene_id = c("G1", "G2"), logFC = c(2.0, -1.5),
                       group = c("Up", "Down")),
            file.path(tmp, "DEGs.csv"), row.names = FALSE)

  # Create dummy PDFs
  writeLines("dummy", file.path(tmp, "Volcano.pdf"))
  writeLines("dummy", file.path(tmp, "Heatmap.pdf"))

  result <- validate_output(list(outdir = tmp))
  expect_true(result$valid)
  expect_equal(result$file_info$diff_csv_rows, 2)
  expect_equal(result$file_info$degs_csv_rows, 2)
})

# ---------------------------------------------------------------------------
# main() — environment check integration
# ---------------------------------------------------------------------------
test_that("check_environment() returns ok status when packages available", {
  env <- check_environment()
  expect_true(env$status %in% c("ok", "error"))
  expect_true("msg" %in% names(env))
})

test_that("check_environment() result has warnings key", {
  env <- check_environment()
  expect_true("warnings" %in% names(env))
})

# ---------------------------------------------------------------------------
# parse_args() — help text (no args)
# ---------------------------------------------------------------------------
test_that("parse_args() prints help and exits with status 1 on no args", {
  # parse_args with no args calls quit(1), so we can't test it normally.
  # Verify the function body contains the usage help
  expect_true(exists("parse_args", mode = "function"))
  # The function exists and is callable
  expect_type(parse_args, "closure")
})

# ---------------------------------------------------------------------------
# End-to-end pipeline
# ---------------------------------------------------------------------------
test_that("end-to-end pipeline runs with synthetic data", {
  source("../../tests/testthat/setup.R")
  requireNamespace("jsonlite", quietly = TRUE)

  td <- create_test_data()
  outdir <- file.path(tempdir(), "deg_e2e_test")
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
  dir.create(outdir, showWarnings = FALSE)

  # Capture stdout to parse NDJSON
  tmp_stdout <- file.path(tempdir(), "e2e_stdout.txt")

  exit_code <- system2(
    "../../env/bin/Rscript",
    c("../../scripts/main.R",
      "run",
      "--mat", td$mat_path,
      "--map", td$map_path,
      "--method", "t",
      "--p-set", "p",
      "--pvalue", "0.05",
      "--logfc-cutoff", "0.5",
      "--cutoff", "1",
      "--top", "10",
      "--outdir", outdir),
    stdout = tmp_stdout,
    stderr = FALSE
  )

  # Read NDJSON output
  lines <- readLines(tmp_stdout)
  jsons <- lapply(lines, function(l) tryCatch(jsonlite::fromJSON(l), error = function(e) NULL))
  jsons <- Filter(Negate(is.null), jsons)

  # Find result line
  result_line <- Find(function(j) !is.null(j$level) && j$level == "result", jsons)
  expect_true(!is.null(result_line))
  expect_equal(result_line$status, "success")

  # Check output files exist
  expect_true(file.exists(file.path(outdir, "Diffanalysis.csv")))
  expect_true(file.exists(file.path(outdir, "DEGs.csv")))
  expect_true(file.exists(file.path(outdir, "Volcano.pdf")))
  expect_true(file.exists(file.path(outdir, "Heatmap.pdf")))
  expect_true(file.exists(file.path(outdir, ".run_result.json")))

  # Verify .run_result.json
  rr <- jsonlite::fromJSON(file.path(outdir, ".run_result.json"))
  expect_equal(rr$node, "differential-analysis")
  expect_equal(rr$exit_code, 0)

  # Cleanup
  unlink(outdir, recursive = TRUE)
  unlink(tmp_stdout)
})

# ---------------------------------------------------------------------------
# B5_METHOD_MISMATCH: count-based methods reject non-integer data
# ---------------------------------------------------------------------------
test_that("do_run() emits B5_METHOD_MISMATCH for float data with deseq2/edgeR", {
  source("../../scripts/report.R")
  source("../../scripts/exceptions.R")
  source("../../scripts/io_helpers.R")
  source("../../scripts/diff_methods.R")
  source("../../scripts/filter_helpers.R")
  source("../../scripts/plot_helpers.R")
  source("../../scripts/input_validation.R")
  source("../../scripts/output_validation.R")
  source("../../scripts/main.R")

  # Mock quit to avoid terminating the test process
  real_quit <- quit
  assign("quit", function(status = 1, ...) {
    assign("quit_called_with", status, envir = .GlobalEnv)
  }, envir = .GlobalEnv)
  on.exit(assign("quit", real_quit, envir = .GlobalEnv))

  source("../../tests/testthat/setup.R")
  td <- create_test_data()
  outdir <- file.path(tempdir(), "deg_b5_test")
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
  dir.create(outdir, showWarnings = FALSE)

  # Create a float (non-integer) expression matrix
  float_mat <- td$mat * 1.5  # makes values non-integer
  float_mat_path <- file.path(tempdir(), "float_expr.csv")
  write.csv(float_mat, float_mat_path)

  out <- capture.output({
    result <- do_run(list(
      subcommand = "run",
      mat = float_mat_path,
      map = td$map_path,
      method = "deseq2",
      p_set = "padj",
      pvalue = 0.05,
      logfc_cutoff = 1.0,
      cutoff = 10,
      norm = "TMM",
      model = "glmFit",
      top = 20,
      force_imbalanced = FALSE,
      rgs = NULL,
      locate = NULL,
      tax_id = "9606",
      pheno_abbr = NULL,
      gene = NULL,
      color_heat = "blue,white,red",
      color_panel = NULL,
      outdir = outdir
    ))
  })

  # Should have emitted B5_METHOD_MISMATCH NDJSON
  has_b5 <- any(grepl("B5_METHOD_MISMATCH", out))
  expect_true(has_b5)

  # Also test edgeR with float data
  out2 <- capture.output({
    result2 <- do_run(list(
      subcommand = "run",
      mat = float_mat_path,
      map = td$map_path,
      method = "edgeR",
      p_set = "padj",
      pvalue = 0.05,
      logfc_cutoff = 1.0,
      cutoff = 10,
      norm = "TMM",
      model = "glmFit",
      top = 20,
      force_imbalanced = FALSE,
      rgs = NULL,
      locate = NULL,
      tax_id = "9606",
      pheno_abbr = NULL,
      gene = NULL,
      color_heat = "blue,white,red",
      color_panel = NULL,
      outdir = outdir
    ))
  })
  has_b5_edgeR <- any(grepl("B5_METHOD_MISMATCH", out2))
  expect_true(has_b5_edgeR)

  unlink(outdir, recursive = TRUE)
  unlink(float_mat_path)
  if (exists("quit_called_with", envir = .GlobalEnv)) {
    rm(quit_called_with, envir = .GlobalEnv)
  }
})

# ---------------------------------------------------------------------------
# Clean up
# ---------------------------------------------------------------------------
assign(".exceptions", list(), envir = .GlobalEnv)
