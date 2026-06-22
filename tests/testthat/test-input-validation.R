# test-input-validation.R — Tests for input validation (input_validation.R)
#
# Tests: validate_input(), standalone CLI behavior

library(testthat)

# Source dependencies
source("../../scripts/report.R")
source("../../scripts/exceptions.R")
source("../../scripts/filter_helpers.R")
source("../../scripts/input_validation.R")

# ---------------------------------------------------------------------------
# Synthetic test data helpers (each creates isolated subdir under tempdir)
# ---------------------------------------------------------------------------

make_valid_mat <- function(tmpdir = tempfile("vmat_")) {
  dir.create(tmpdir, showWarnings = FALSE)
  genes <- paste0("GENE", 1:10)
  samples <- c("CASE1", "CASE2", "CASE3", "CTRL1", "CTRL2", "CTRL3")
  mat <- matrix(rpois(10 * 6, lambda = 100), nrow = 10, ncol = 6,
                dimnames = list(genes, samples))
  path <- file.path(tmpdir, "test_expr.csv")
  write.csv(mat, path)
  path
}

make_valid_map <- function(tmpdir = tempfile("vmap_")) {
  dir.create(tmpdir, showWarnings = FALSE)
  map <- data.frame(
    sample = c("CASE1", "CASE2", "CASE3", "CTRL1", "CTRL2", "CTRL3"),
    group = c(rep("Case", 3), rep("Control", 3)),
    stringsAsFactors = FALSE
  )
  path <- file.path(tmpdir, "test_map.csv")
  write.csv(map, path, row.names = FALSE)
  path
}

make_imbalanced_map <- function() {
  tmpdir <- tempfile("vimb_")
  dir.create(tmpdir, showWarnings = FALSE)
  map <- data.frame(
    sample = c("CASE1", paste0("CTRL", 1:11)),
    group = c("Case", rep("Control", 11)),
    stringsAsFactors = FALSE
  )
  path <- file.path(tmpdir, "test_map_imbalanced.csv")
  write.csv(map, path, row.names = FALSE)

  # Need a matching expression matrix with 12 samples
  genes <- paste0("GENE", 1:10)
  samples <- map$sample
  mat <- matrix(rpois(10 * 12, lambda = 100), nrow = 10, ncol = 12,
                dimnames = list(genes, samples))
  mat_path <- file.path(tmpdir, "test_expr_imbalanced.csv")
  write.csv(mat, mat_path)

  list(mat_path = mat_path, map_path = path)
}

# ---------------------------------------------------------------------------
# validate_input() function tests
# ---------------------------------------------------------------------------

test_that("validate_input() returns valid=TRUE for correct inputs", {
  old_ex <- get(".exceptions", envir = .GlobalEnv)
  on.exit(assign(".exceptions", old_ex, envir = .GlobalEnv))
  assign(".exceptions", list(), envir = .GlobalEnv)

  mat_path <- make_valid_mat()
  map_path <- make_valid_map()
  result <- validate_input(list(mat = mat_path, map = map_path))
  expect_true(result$valid)
  expect_equal(result$reason, "OK")
  expect_equal(result$n_genes, 10)
  expect_equal(result$n_samples, 6)
  expect_equal(result$groups, "Case vs Control")
})

test_that("validate_input() returns valid=FALSE when mat file does not exist", {
  result <- validate_input(list(mat = "/nonexistent/mat.csv", map = NULL))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("validate_input() returns valid=FALSE when mat is NULL", {
  result <- validate_input(list(mat = NULL, map = NULL))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("validate_input() returns valid=FALSE when map file does not exist", {
  mat_path <- make_valid_mat()
  result <- validate_input(list(mat = mat_path, map = "/nonexistent/map.csv"))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("validate_input() returns valid=FALSE when map is NULL", {
  mat_path <- make_valid_mat()
  result <- validate_input(list(mat = mat_path, map = NULL))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("validate_input() returns valid=FALSE for non-CSV mat file", {
  td <- tempfile("ncs_"); dir.create(td)
  tmp <- file.path(td, "bad.bin")
  writeBin(as.raw(c(0x00, 0xFF, 0x00, 0xFF)), tmp)
  map_path <- make_valid_map()
  result <- validate_input(list(mat = tmp, map = map_path))
  expect_false(result$valid)
  expect_match(result$reason, "Cannot read")
})

test_that("validate_input() returns valid=FALSE for single-column mat", {
  td <- tempfile("onec_"); dir.create(td)
  map_path <- make_valid_map()
  tmp <- file.path(td, "onecol.csv")
  write.csv(data.frame(gene_id = paste0("G", 1:5)), tmp, row.names = FALSE)
  result <- validate_input(list(mat = tmp, map = map_path))
  expect_false(result$valid)
  expect_match(result$reason, "only 1 column")
})

test_that("validate_input() returns valid=FALSE for empty mat (0 rows)", {
  td <- tempfile("empt_"); dir.create(td)
  map_path <- make_valid_map()
  tmp <- file.path(td, "empty.csv")
  df <- data.frame(gene_id = character(0), S1 = character(0))
  write.csv(df, tmp, row.names = FALSE)
  result <- validate_input(list(mat = tmp, map = map_path))
  expect_false(result$valid)
  expect_match(result$reason, "0 rows|empty")
})

test_that("validate_input() returns valid=FALSE for non-CSV map file", {
  td <- tempfile("ncs2_"); dir.create(td)
  mat_path <- make_valid_mat()
  tmp <- file.path(td, "badmap.bin")
  writeBin(as.raw(c(0x00, 0xFF, 0x00, 0xFF)), tmp)
  result <- validate_input(list(mat = mat_path, map = tmp))
  expect_false(result$valid)
  expect_match(result$reason, "Cannot read")
})

test_that("validate_input() returns valid=FALSE for single-column map", {
  td <- tempfile("onem_"); dir.create(td)
  mat_path <- make_valid_mat()
  tmp <- file.path(td, "onecol_map.csv")
  write.csv(data.frame(sample = c("CASE1", "CTRL1")), tmp, row.names = FALSE)
  result <- validate_input(list(mat = mat_path, map = tmp))
  expect_false(result$valid)
  expect_match(result$reason, "only 1 column")
})

test_that("validate_input() returns valid=FALSE when map samples missing from mat", {
  td <- tempfile("msmm_"); dir.create(td)
  mat_path <- make_valid_mat()  # has CASE1,CASE2,CASE3,CTRL1,CTRL2,CTRL3
  tmp <- file.path(td, "missing_map.csv")
  map <- data.frame(
    sample = c("CASE1", "EXTRA_SAMPLE"),
    group = c("Case", "Control"),
    stringsAsFactors = FALSE
  )
  write.csv(map, tmp, row.names = FALSE)
  result <- validate_input(list(mat = mat_path, map = tmp))
  expect_false(result$valid)
  expect_match(result$reason, "Sample.*in map but not in matrix")
})

test_that("validate_input() returns valid=FALSE when mat has extra samples not in map", {
  td <- tempfile("esnm_"); dir.create(td)
  genes <- paste0("GENE", 1:10)
  samples <- c("CASE1", "CASE2", "UNKNOWN_SAMPLE")
  mat <- matrix(rpois(10 * 3, lambda = 100), nrow = 10, ncol = 3,
                dimnames = list(genes, samples))
  mat_path <- file.path(td, "extra_mat.csv")
  write.csv(mat, mat_path)

  map_path <- file.path(td, "extra_map.csv")
  map <- data.frame(
    sample = c("CASE1", "CASE2"),
    group = c("Case", "Control"),
    stringsAsFactors = FALSE
  )
  write.csv(map, map_path, row.names = FALSE)
  result <- validate_input(list(mat = mat_path, map = map_path))
  expect_false(result$valid)
  expect_match(result$reason, "Sample.*in matrix but not in map")
})

test_that("validate_input() returns valid=FALSE when groups != 2", {
  td <- tempfile("grp3_"); dir.create(td)
  genes <- paste0("GENE", 1:10)
  samples <- c("A", "B", "C")
  mat <- matrix(rpois(10 * 3, lambda = 100), nrow = 10, ncol = 3,
                dimnames = list(genes, samples))
  mat_path <- file.path(td, "three_grp_mat.csv")
  write.csv(mat, mat_path)

  map_path <- file.path(td, "three_grp_map.csv")
  map <- data.frame(
    sample = c("A", "B", "C"),
    group = c("Group1", "Group2", "Group3"),
    stringsAsFactors = FALSE
  )
  write.csv(map, map_path, row.names = FALSE)
  result <- validate_input(list(mat = mat_path, map = map_path))
  expect_false(result$valid)
  expect_match(result$reason, "Expected exactly 2 groups")
})

test_that("validate_input() returns valid=FALSE for imbalanced proportion", {
  data <- make_imbalanced_map()
  result <- validate_input(list(mat = data$mat_path, map = data$map_path))
  expect_false(result$valid)
  expect_match(result$reason, "proportion|ratio")
})

test_that("validate_input() accepts imbalanced proportion with force_imbalanced=TRUE", {
  data <- make_imbalanced_map()
  result <- validate_input(list(
    mat = data$mat_path, map = data$map_path, force_imbalanced = TRUE
  ))
  expect_true(result$valid)
})

test_that("validate_input() accepts imbalanced proportion with force_imbalanced=TRUE string", {
  data <- make_imbalanced_map()
  result <- validate_input(list(
    mat = data$mat_path, map = data$map_path, force_imbalanced = "TRUE"
  ))
  expect_true(result$valid)
})

test_that("validate_input() returns valid=FALSE for just 1 group", {
  td <- tempfile("grp1_"); dir.create(td)
  genes <- paste0("GENE", 1:10)
  samples <- c("A", "B", "C")
  mat <- matrix(rpois(10 * 3, lambda = 100), nrow = 10, ncol = 3,
                dimnames = list(genes, samples))
  mat_path <- file.path(td, "one_grp_mat.csv")
  write.csv(mat, mat_path)

  map_path <- file.path(td, "one_grp_map.csv")
  map <- data.frame(
    sample = c("A", "B", "C"),
    group = c("Only", "Only", "Only"),
    stringsAsFactors = FALSE
  )
  write.csv(map, map_path, row.names = FALSE)
  result <- validate_input(list(mat = mat_path, map = map_path))
  expect_false(result$valid)
  expect_match(result$reason, "Expected exactly 2 groups")
})

test_that("validate_input() handles gene column name correctly (non-'gene_id' first col)", {
  td <- tempfile("symb_"); dir.create(td)
  samples <- c("CASE1", "CASE2", "CTRL1", "CTRL2")
  mat <- data.frame(
    Symbol = paste0("GENE", 1:10),
    CASE1 = rpois(10, 100),
    CASE2 = rpois(10, 100),
    CTRL1 = rpois(10, 100),
    CTRL2 = rpois(10, 100)
  )
  mat_path <- file.path(td, "symb_mat.csv")
  write.csv(mat, mat_path, row.names = FALSE)

  map_path <- file.path(td, "symb_map.csv")
  map <- data.frame(
    sample = c("CASE1", "CASE2", "CTRL1", "CTRL2"),
    group = c("Case", "Case", "Control", "Control"),
    stringsAsFactors = FALSE
  )
  write.csv(map, map_path, row.names = FALSE)
  result <- validate_input(list(mat = mat_path, map = map_path))
  expect_true(result$valid)
  expect_equal(result$n_genes, 10)
  expect_equal(result$n_samples, 4)
})

# ---------------------------------------------------------------------------
# validate_input() is sourceable (function exists and is callable)
# ---------------------------------------------------------------------------

test_that("validate_input exists as a function", {
  expect_true(exists("validate_input"))
  expect_type(validate_input, "closure")
})
