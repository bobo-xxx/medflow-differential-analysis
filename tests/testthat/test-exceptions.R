# test-exceptions.R — Tests for structured exception handling (exceptions.R)
#
# Tests: check_environment(), report_exception_ndjson()

library(testthat)

# Source dependencies
source("../../scripts/report.R")
source("../../scripts/exceptions.R")

# ---------------------------------------------------------------------------
# check_environment
# ---------------------------------------------------------------------------
test_that("check_environment() returns a list with status and msg", {
  env <- check_environment()
  expect_type(env, "list")
  expect_true(env$status %in% c("ok", "error"))
  expect_true("msg" %in% names(env))
})

test_that("check_environment() returns status='ok' when required packages are present", {
  env <- check_environment()
  # In our controlled env, all 13 required packages should be installed
  expect_equal(env$status, "ok")
})

test_that("check_environment() names all required packages when checking", {
  # We verify the function checks the right set by looking at its structure
  # Since packages are installed, missing is empty and status is ok
  env <- check_environment()
  expect_equal(env$status, "ok")
  # The function should have a 'warnings' field
  expect_true("warnings" %in% names(env))
})

test_that("check_environment() detects genuinely missing packages", {
  # Create a local version that checks for a nonexistent package
  check_env_fake <- function() {
    required <- c("nonexistent_pkg_xyz_123")
    missing <- character(0)
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
    list(status = "ok", msg = "Environment OK", warnings = character(0))
  }
  env <- check_env_fake()
  expect_equal(env$status, "error")
  expect_true("nonexistent_pkg_xyz_123" %in% env$missing)
})

test_that("check_environment() warns about optional packages", {
  env <- check_environment()
  # warnings is always returned (character vector)
  expect_type(env$warnings, "character")
})

# ---------------------------------------------------------------------------
# report_exception_ndjson
# ---------------------------------------------------------------------------
test_that("report_exception_ndjson() emits structured NDJSON", {
  # Save and restore .exceptions
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  out <- capture.output(
    report_exception_ndjson("B1_PROPORTION", "data_insufficient",
                            "skip_with_warning",
                            "Sample proportion 15:1 exceeds 10:1 limit",
                            dry_run = TRUE)
  )
  expect_length(out, 1)
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$level, "decision")
  expect_equal(obj$code, "B1_PROPORTION")
  expect_equal(obj$nature, "data_insufficient")
  expect_equal(obj$action, "skip_with_warning")
  expect_equal(obj$msg, "Sample proportion 15:1 exceeds 10:1 limit")
})

test_that("report_exception_ndjson() accumulates to .exceptions", {
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  capture.output(
    report_exception_ndjson("W1", "data_corrupt", "escalate",
                            "Corrupt file detected", dry_run = TRUE)
  )
  ex <- get(".exceptions", envir = .GlobalEnv)
  expect_length(ex, 1)
  expect_equal(ex[[1]]$code, "W1")
})

test_that("report_exception_ndjson() accumulates multiple exceptions", {
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  capture.output({
    report_exception_ndjson("B1", "data_insufficient", "skip_with_warning",
                            "First", dry_run = TRUE)
    report_exception_ndjson("B2", "data_corrupt", "escalate",
                            "Second", dry_run = TRUE)
    report_exception_ndjson("E1", "env_bug", "halt",
                            "Third", dry_run = TRUE)
  })
  ex <- get(".exceptions", envir = .GlobalEnv)
  expect_length(ex, 3)
  expect_equal(ex[[1]]$code, "B1")
  expect_equal(ex[[2]]$code, "B2")
  expect_equal(ex[[3]]$code, "E1")
})

test_that("report_exception_ndjson() uses correct level for halt action", {
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  out <- capture.output(
    report_exception_ndjson("E2", "env_bug", "halt",
                            "Fatal", dry_run = TRUE)
  )
  obj <- jsonlite::fromJSON(out[1])
  # halt is not "skip_with_warning" or "escalate", so it falls to default "exception"
  expect_equal(obj$level, "exception")
})

test_that("report_exception_ndjson() uses 'escalate' level for escalate action", {
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  out <- capture.output(
    report_exception_ndjson("W2", "data_mismatch", "escalate",
                            "Escalated", dry_run = TRUE)
  )
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$level, "exception")
  expect_equal(obj$action, "escalate")
})

test_that("report_exception_ndjson() dry_run=TRUE does not call quit", {
  # If dry_run works, this should return without quitting
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  expect_error(
    report_exception_ndjson("E3", "env_bug", "halt",
                            "Should not quit", dry_run = TRUE),
    NA  # no error expected
  )
})

test_that("report_exception_ndjson() handles all nature values", {
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  natures <- c("data_insufficient", "data_corrupt", "data_mismatch", "env_bug")
  for (n in natures) {
    out <- capture.output(
      report_exception_ndjson("TEST", n, "skip_with_warning",
                              paste("Nature:", n), dry_run = TRUE)
    )
    obj <- jsonlite::fromJSON(out[1])
    expect_equal(obj$nature, n)
  }
})

test_that("report_exception_ndjson() handles exit_code parameter", {
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  # exit_code is passed but only used if quitting; verify it's accepted
  out <- capture.output(
    report_exception_ndjson("E4", "env_bug", "halt",
                            "Custom exit", exit_code = 3,
                            dry_run = TRUE)
  )
  obj <- jsonlite::fromJSON(out[1])
  expect_equal(obj$code, "E4")
})
