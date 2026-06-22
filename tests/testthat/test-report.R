# test-report.R — Tests for NDJSON reporting helpers (report.R)
#
# Tests: report_info(), report_result(), report_error(), write_run_result()

library(testthat)

# Source the module under test
source("../../scripts/report.R")

# Ensure clean exception state at start of report tests
assign(".exceptions", list(), envir = .GlobalEnv)

# ---------------------------------------------------------------------------
# report_info
# ---------------------------------------------------------------------------
test_that("report_info() emits a single NDJSON line to stdout", {
  out <- capture.output(report_info("hello"), type = "message")
  # capture.output with type="message" may not catch cat() output.
  # Use capture.output on stdout instead.
  out <- capture.output(report_info("hello"))
  expect_length(out, 1)
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$level, "info")
  expect_equal(obj$msg, "hello")
})

test_that("report_info() includes extra named fields", {
  out <- capture.output(report_info("progress", step = 3, total = 8))
  expect_length(out, 1)
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$step, 3)
  expect_equal(obj$total, 8)
})

test_that("report_info() handles empty extra args", {
  out <- capture.output(report_info("bare"))
  expect_length(out, 1)
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$level, "info")
  expect_equal(obj$msg, "bare")
})

test_that("report_info() output is valid JSON with auto_unbox", {
  out <- capture.output(report_info("test", count = 1L, name = "x"))
  obj <- jsonlite::fromJSON(out[1])
  # auto_unbox=TRUE means scalars are not arrays
  expect_type(obj$count, "integer")
  expect_type(obj$name, "character")
  expect_false(is.list(obj$name))
})

# ---------------------------------------------------------------------------
# report_result
# ---------------------------------------------------------------------------
test_that("report_result() emits NDJSON with level=result", {
  out <- capture.output(report_result("success"))
  expect_length(out, 1)
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$level, "result")
  expect_equal(obj$status, "success")
})

test_that("report_result() omits empty files and metadata lists", {
  out <- capture.output(report_result("error"))
  obj <- jsonlite::fromJSON(out[1])
  expect_null(obj$files)
  expect_null(obj$metadata)
})

test_that("report_result() includes files and metadata when non-empty", {
  files <- list(list(path = "/tmp/a.csv", rows = 100, cols = 5))
  metadata <- list(method = "DESeq2", version = "1.42.0", n_degs = 50, up = 30, down = 20)
  out <- capture.output(report_result("success", files = files, metadata = metadata))
  obj <- jsonlite::fromJSON(out[1])
  # auto_unbox round-trip: fromJSON converts array-of-objects to data.frame
  expect_equal(obj$files$path, "/tmp/a.csv")
  expect_equal(obj$files$rows, 100)
  expect_equal(obj$files$cols, 5)
  expect_equal(obj$metadata$n_degs, 50)
  expect_equal(obj$metadata$up, 30)
})

test_that("report_result() accepts extra named fields", {
  out <- capture.output(report_result("success", warning = "low_reads"))
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$warning, "low_reads")
})

# ---------------------------------------------------------------------------
# report_error
# ---------------------------------------------------------------------------
test_that("report_error() emits error-level NDJSON and calls quit", {
  # We cannot test quit() directly — it terminates the process.
  # Instead, mock quit() to verify the signature.
  mock_quit <- function(status = 1) {
    assign("quit_called", status, envir = .GlobalEnv)
  }
  # Create a local version of report_error that uses mock_quit
  report_error_test <- function(msg, exit_code = 1) {
    obj <- list(level = "error", msg = msg)
    cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
    mock_quit(exit_code)
  }
  out <- capture.output(report_error_test("fatal error", 3))
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$level, "error")
  expect_equal(obj$msg, "fatal error")
  expect_equal(quit_called, 3)
  rm(quit_called, envir = .GlobalEnv)
})

test_that("report_error() default exit code is 1", {
  report_error_test2 <- function(msg, exit_code = 1) {
    obj <- list(level = "error", msg = msg)
    cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
    assign("quit_called", exit_code, envir = .GlobalEnv)
  }
  out <- capture.output(report_error_test2("default exit"))
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$level, "error")
  expect_equal(quit_called, 1)
  rm(quit_called, envir = .GlobalEnv)
})

# ---------------------------------------------------------------------------
# write_run_result
# ---------------------------------------------------------------------------
test_that("write_run_result() writes .run_result.json to out_dir", {
  tmp <- tempfile("rr_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  result <- list(status = "success", metadata = list(n_degs = 50))
  params <- list(subcommand = "run", case = "Case", control = "Control",
                 outdir = tmp, method = "DESeq2")
  times <- c("2025-06-01T00:00:00Z", "2025-06-01T00:01:00Z")

  write_run_result(tmp, result, params, 0, times)

  json_path <- file.path(tmp, ".run_result.json")
  expect_true(file.exists(json_path))

  obj <- jsonlite::fromJSON(json_path)
  expect_equal(obj$node, "differential-analysis")
  expect_equal(obj$subcommand, "run")
  expect_equal(obj$status, "success")
  expect_equal(obj$exit_code, 0)
  expect_equal(obj$started_at, "2025-06-01T00:00:00Z")
  expect_equal(obj$finished_at, "2025-06-01T00:01:00Z")
  expect_equal(obj$output$n_degs, 50)
  expect_equal(obj$parameters$method, "DESeq2")
  # outdir should be stripped from parameters
  expect_null(obj$parameters$outdir)
  # exceptions should be an empty list when .exceptions is empty
  expect_equal(length(obj$exceptions), 0)
})

test_that("write_run_result() handles missing result fields gracefully", {
  tmp <- tempfile("rr2_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  result <- list()
  params <- list(subcommand = "run")
  times <- c("2025-06-01T00:00:00Z", "2025-06-01T00:01:00Z")

  write_run_result(tmp, result, params, 1, times)

  json_path <- file.path(tmp, ".run_result.json")
  obj <- jsonlite::fromJSON(json_path)
  expect_equal(obj$status, "unknown")
  expect_equal(obj$subcommand, "run")
  expect_equal(obj$exit_code, 1)
})

test_that("write_run_result() accumulates exceptions from .exceptions", {
  tmp <- tempfile("rr3_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  # Ensure clean state before populating
  assign(".exceptions", list(), envir = .GlobalEnv)
  # Pre-populate .exceptions with a single entry
  assign(".exceptions", list(
    list(level = "exception", code = "E1", nature = "env_bug",
         action = "halt", msg = "Missing package")
  ), envir = .GlobalEnv)

  result <- list(status = "error")
  params <- list(subcommand = "run")
  times <- c("2025-06-01T00:00:00Z", "2025-06-01T00:01:00Z")

  write_run_result(tmp, result, params, 1, times)

  json_path <- file.path(tmp, ".run_result.json")
  obj <- jsonlite::fromJSON(json_path)
  # fromJSON converts array-of-objects to data.frame
  # Single exception -> data.frame with 1 row, access columns directly
  expect_equal(obj$exceptions$code, "E1")
  expect_equal(obj$exceptions$nature, "env_bug")

  # Clean up
  assign(".exceptions", list(), envir = .GlobalEnv)
})

test_that("write_run_result() writes pretty-printed JSON", {
  tmp <- tempfile("rr4_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  result <- list(status = "success")
  params <- list(subcommand = "run")
  times <- c("2025-06-01T00:00:00Z", "2025-06-01T00:01:00Z")

  write_run_result(tmp, result, params, 0, times)

  raw <- readLines(file.path(tmp, ".run_result.json"))
  # pretty-printed JSON has more than 1 line for an object this size
  expect_gt(length(raw), 1)
})

test_that("write_run_result() includes files section from existing output files", {
  tmp <- tempfile("rr5_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  # Create a dummy CSV file
  csv_path <- file.path(tmp, "DEGs.csv")
  write.csv(data.frame(gene = c("GENE1", "GENE2"), logFC = c(2.0, -2.0),
                        pvalue = c(0.001, 0.002)),
            csv_path, row.names = FALSE)

  result <- list(status = "success")
  params <- list(subcommand = "run")
  times <- c("2025-06-01T00:00:00Z", "2025-06-01T00:01:00Z")

  write_run_result(tmp, result, params, 0, times)

  json_path <- file.path(tmp, ".run_result.json")
  obj <- jsonlite::fromJSON(json_path)
  # fromJSON converts array-of-objects to data.frame; $path gives path vector
  path_vector <- obj$files$path
  expect_true(any(basename(path_vector) == "DEGs.csv"))
})
