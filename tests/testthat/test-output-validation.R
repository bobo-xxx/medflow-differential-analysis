# test-output-validation.R — Tests for output validation (output_validation.R)
#
# Tests: validate_output(), standalone CLI behavior

library(testthat)

# Source dependencies
source("../../scripts/report.R")
source("../../scripts/exceptions.R")
source("../../scripts/output_validation.R")

# ---------------------------------------------------------------------------
# Synthetic output directory helpers
# ---------------------------------------------------------------------------

make_valid_output <- function() {
  outdir <- tempfile("ovalid_")
  dir.create(outdir, showWarnings = FALSE)

  # Diffanalysis.csv with required columns
  dif <- data.frame(
    gene_id = paste0("GENE", 1:50),
    logFC = rnorm(50, 0, 1),
    AveExpr = runif(50, 0, 10),
    t = rnorm(50),
    Pvalue = runif(50, 0, 0.05),
    Padj = p.adjust(runif(50, 0, 0.05), method = "BH"),
    B = rnorm(50),
    stringsAsFactors = FALSE
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)

  # DEGs.csv with group column
  degs <- dif[dif$Padj < 0.05, ]
  degs$group <- ifelse(degs$logFC > 0, "Up", "Down")
  write.csv(degs, file.path(outdir, "DEGs.csv"), row.names = FALSE)

  # Volcano.pdf (non-zero size)
  pdf(file.path(outdir, "Volcano.pdf"), width = 8, height = 6)
  plot(1:10, main = "Volcano Test")
  dev.off()

  # Heatmap.pdf (non-zero size)
  pdf(file.path(outdir, "Heatmap.pdf"), width = 8, height = 6)
  plot(1:10, main = "Heatmap Test")
  dev.off()

  outdir
}

make_empty_output_dir <- function() {
  outdir <- tempfile("oempty_")
  dir.create(outdir, showWarnings = FALSE)
  outdir
}

# ---------------------------------------------------------------------------
# validate_output() function tests
# ---------------------------------------------------------------------------

test_that("validate_output() returns valid=TRUE for complete outputs", {
  outdir <- make_valid_output()
  result <- validate_output(list(outdir = outdir))
  expect_true(result$valid)
  expect_equal(result$reason, "OK")
  expect_true(is.list(result$file_info))
  expect_true(result$file_info$diff_csv_rows > 0)
  expect_true(result$file_info$degs_csv_rows > 0)
})

test_that("validate_output() returns valid=FALSE when outdir is NULL", {
  result <- validate_output(list(outdir = NULL))
  expect_false(result$valid)
  expect_match(result$reason, "not found|not specified")
})

test_that("validate_output() returns valid=FALSE when outdir does not exist", {
  result <- validate_output(list(outdir = "/nonexistent/output/dir"))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("validate_output() returns valid=FALSE when Diffanalysis.csv is missing", {
  outdir <- make_empty_output_dir()
  # Create DEGs.csv and plots but no Diffanalysis.csv
  write.csv(data.frame(gene_id = "G1", logFC = 1.0, group = "Up"),
            file.path(outdir, "DEGs.csv"), row.names = FALSE)
  pdf(file.path(outdir, "Volcano.pdf")); plot(1); dev.off()
  pdf(file.path(outdir, "Heatmap.pdf")); plot(1); dev.off()
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "Missing Diffanalysis.csv")
})

test_that("validate_output() returns valid=FALSE when Diffanalysis.csv is unreadable", {
  outdir <- make_empty_output_dir()
  writeLines(c("garbage", "not,csv,format"), file.path(outdir, "Diffanalysis.csv"))
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "Cannot read Diffanalysis.csv")
})

test_that("validate_output() returns valid=FALSE for Diffanalysis.csv missing gene_id", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(logFC = c(1.0, -0.5), Pvalue = c(0.01, 0.02), Padj = c(0.05, 0.10))
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "missing required columns")
  expect_match(result$reason, "gene_id")
})

test_that("validate_output() returns valid=FALSE for Diffanalysis.csv missing logFC", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(gene_id = c("G1", "G2"), Pvalue = c(0.01, 0.02), Padj = c(0.05, 0.10))
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "missing required columns")
  expect_match(result$reason, "logFC")
})

test_that("validate_output() returns valid=FALSE for Diffanalysis.csv missing Pvalue", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(gene_id = c("G1", "G2"), logFC = c(1.0, -0.5), Padj = c(0.05, 0.10))
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "missing required columns")
  expect_match(result$reason, "Pvalue")
})

test_that("validate_output() returns valid=FALSE for Diffanalysis.csv missing Padj", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(gene_id = c("G1", "G2"), logFC = c(1.0, -0.5), Pvalue = c(0.01, 0.02))
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "missing required columns")
  expect_match(result$reason, "Padj")
})

test_that("validate_output() returns valid=FALSE for empty Diffanalysis.csv (0 rows)", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(
    gene_id = character(0), logFC = numeric(0),
    Pvalue = numeric(0), Padj = numeric(0)
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  # Also create supporting files
  write.csv(data.frame(gene_id = character(0), logFC = numeric(0), group = character(0)),
            file.path(outdir, "DEGs.csv"), row.names = FALSE)
  pdf(file.path(outdir, "Volcano.pdf")); plot(1); dev.off()
  pdf(file.path(outdir, "Heatmap.pdf")); plot(1); dev.off()
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "0 rows")
})

test_that("validate_output() returns valid=FALSE when DEGs.csv is missing", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(
    gene_id = paste0("G", 1:5), logFC = rnorm(5),
    Pvalue = runif(5, 0, 0.05), Padj = runif(5, 0, 0.05)
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  pdf(file.path(outdir, "Volcano.pdf")); plot(1); dev.off()
  pdf(file.path(outdir, "Heatmap.pdf")); plot(1); dev.off()
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "Missing DEGs.csv")
})

test_that("validate_output() returns valid=FALSE when DEGs.csv is unreadable", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(
    gene_id = paste0("G", 1:5), logFC = rnorm(5),
    Pvalue = runif(5, 0, 0.05), Padj = runif(5, 0, 0.05)
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  writeBin(as.raw(c(0x00, 0xFF, 0x00, 0xFF)), file.path(outdir, "DEGs.csv"))
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "Cannot read DEGs.csv")
})

test_that("validate_output() returns valid=FALSE when DEGs.csv missing group column", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(
    gene_id = paste0("G", 1:5), logFC = rnorm(5),
    Pvalue = runif(5, 0, 0.05), Padj = runif(5, 0, 0.05)
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  degs <- data.frame(gene_id = paste0("G", 1:3), logFC = rnorm(3))
  write.csv(degs, file.path(outdir, "DEGs.csv"), row.names = FALSE)
  pdf(file.path(outdir, "Volcano.pdf")); plot(1); dev.off()
  pdf(file.path(outdir, "Heatmap.pdf")); plot(1); dev.off()
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "missing 'group'")
})

test_that("validate_output() returns valid=FALSE when Volcano.pdf is missing", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(
    gene_id = paste0("G", 1:5), logFC = rnorm(5),
    Pvalue = runif(5, 0, 0.05), Padj = runif(5, 0, 0.05)
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  degs <- data.frame(gene_id = "G1", logFC = 1.0, group = "Up")
  write.csv(degs, file.path(outdir, "DEGs.csv"), row.names = FALSE)
  pdf(file.path(outdir, "Heatmap.pdf")); plot(1); dev.off()
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "Missing plot file: Volcano.pdf")
})

test_that("validate_output() returns valid=FALSE when Heatmap.pdf is missing", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(
    gene_id = paste0("G", 1:5), logFC = rnorm(5),
    Pvalue = runif(5, 0, 0.05), Padj = runif(5, 0, 0.05)
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  degs <- data.frame(gene_id = "G1", logFC = 1.0, group = "Up")
  write.csv(degs, file.path(outdir, "DEGs.csv"), row.names = FALSE)
  pdf(file.path(outdir, "Volcano.pdf")); plot(1); dev.off()
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "Missing plot file: Heatmap.pdf")
})

test_that("validate_output() returns valid=FALSE when plot file is empty (0 bytes)", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(
    gene_id = paste0("G", 1:5), logFC = rnorm(5),
    Pvalue = runif(5, 0, 0.05), Padj = runif(5, 0, 0.05)
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  degs <- data.frame(gene_id = "G1", logFC = 1.0, group = "Up")
  write.csv(degs, file.path(outdir, "DEGs.csv"), row.names = FALSE)
  pdf(file.path(outdir, "Volcano.pdf")); plot(1); dev.off()
  # Create an empty Heatmap.pdf
  writeLines(character(0), file.path(outdir, "Heatmap.pdf"))
  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "empty")
})

test_that("validate_output() file_info contains expected keys", {
  outdir <- make_valid_output()
  result <- validate_output(list(outdir = outdir))
  expect_true(result$valid)
  expect_true("diff_csv_rows" %in% names(result$file_info))
  expect_true("diff_csv_cols" %in% names(result$file_info))
  expect_true("degs_csv_rows" %in% names(result$file_info))
})

test_that("validate_output() handles extra optional columns in diff CSV", {
  outdir <- make_empty_output_dir()
  dif <- data.frame(
    gene_id = paste0("G", 1:5),
    logFC = rnorm(5),
    Pvalue = runif(5, 0, 0.05),
    Padj = runif(5, 0, 0.05),
    extra_col1 = 1:5,
    extra_col2 = LETTERS[1:5],
    stringsAsFactors = FALSE
  )
  write.csv(dif, file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  degs <- data.frame(gene_id = "G1", logFC = 1.0, group = "Up")
  write.csv(degs, file.path(outdir, "DEGs.csv"), row.names = FALSE)
  pdf(file.path(outdir, "Volcano.pdf")); plot(1); dev.off()
  pdf(file.path(outdir, "Heatmap.pdf")); plot(1); dev.off()
  result <- validate_output(list(outdir = outdir))
  expect_true(result$valid)
  expect_equal(result$file_info$diff_csv_cols, 6)
})

# ---------------------------------------------------------------------------
# validate_output() is sourceable (function exists and is callable)
# ---------------------------------------------------------------------------

test_that("validate_output exists as a function", {
  expect_true(exists("validate_output"))
  expect_type(validate_output, "closure")
})
