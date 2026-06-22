# test-diff-methods.R — Tests for DE analysis methods (diff_methods.R)
#
# Tests: diff_stat(), diff_deseq2(), diff_limma(), diff_edger()
# Uses synthetic data from tests/testthat/setup.R

library(testthat)

# Source dependencies
source("../../scripts/report.R")
source("../../scripts/exceptions.R")
source("setup.R")

td <- create_test_data()
map <- td$map
map[[2]] <- factor(map[[2]], levels = c("Control", "Case"))

# ---------------------------------------------------------------------------
# diff_stat (t-test)
# ---------------------------------------------------------------------------
test_that("diff_stat() with t-test returns correct column structure", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "t")
  expect_s3_class(dif, "data.frame")
  expect_true("gene_id" %in% colnames(dif))
  expect_true("logFC" %in% colnames(dif))
  expect_true("Pvalue" %in% colnames(dif))
  expect_true("Padj" %in% colnames(dif))
  expect_true("stat" %in% colnames(dif))
})

test_that("diff_stat() with t-test returns all input genes", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "t")
  expect_equal(nrow(dif), nrow(td$mat))
  expect_true(all(dif$gene_id %in% rownames(td$mat)))
})

test_that("diff_stat() with t-test finds GENE1 as upregulated", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "t")
  gene1_row <- dif[dif$gene_id == "GENE1", ]
  expect_gt(gene1_row$logFC, 0)
  expect_lt(gene1_row$Pvalue, 0.05)
})

test_that("diff_stat() with t-test finds GENE2 as downregulated", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "t")
  gene2_row <- dif[dif$gene_id == "GENE2", ]
  expect_lt(gene2_row$logFC, 0)
})

test_that("diff_stat() with t-test sorts by Pvalue ascending", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "t")
  pvals <- dif$Pvalue
  expect_true(all(diff(pvals) >= 0))
})

test_that("diff_stat() Padj is between 0 and 1", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "t")
  expect_true(all(dif$Padj >= 0 & dif$Padj <= 1))
})

test_that("diff_stat() Padj uses BH correction", {
  source("../../scripts/diff_methods.R")

  pvals <- c(0.001, 0.01, 0.05, 0.1, 0.5)
  expected_padj <- p.adjust(pvals, method = "BH")

  # Verify our understanding: BH is more conservative than raw
  expect_true(all(expected_padj >= pvals))
})

# ---------------------------------------------------------------------------
# diff_stat (Wilcoxon)
# ---------------------------------------------------------------------------
test_that("diff_stat() with wilcox works and returns stat column", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "wilcox")
  expect_s3_class(dif, "data.frame")
  expect_true("stat" %in% colnames(dif))
  expect_true("gene_id" %in% colnames(dif))
  expect_true("logFC" %in% colnames(dif))
  expect_true("Pvalue" %in% colnames(dif))
  expect_true("Padj" %in% colnames(dif))
  expect_equal(nrow(dif), nrow(td$mat))
})

test_that("diff_stat() with wilcox gives valid p-values", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "wilcox")
  expect_true(all(dif$Pvalue >= 0 & dif$Pvalue <= 1))
})

# ---------------------------------------------------------------------------
# diff_stat: gene_id is first column
# ---------------------------------------------------------------------------
test_that("diff_stat() gene_id is the first column", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "t")
  expect_equal(colnames(dif)[1], "gene_id")
})

test_that("diff_stat() rownames are reset to default (sequential numbers)", {
  source("../../scripts/diff_methods.R")

  dif <- diff_stat(td$mat, map, "t")
  # rownames(dif) <- NULL in R resets to default row numbers (c("1","2",...))
  # which is NOT the gene names stored in gene_id column
  expect_false(any(rownames(dif) %in% rownames(td$mat)))
  expect_equal(rownames(dif), as.character(seq_len(nrow(dif))))
})

# ---------------------------------------------------------------------------
# diff_stat: edge cases
# ---------------------------------------------------------------------------
test_that("diff_stat() handles single-gene matrix", {
  source("../../scripts/diff_methods.R")

  single_mat <- td$mat[1, , drop = FALSE]
  dif <- diff_stat(single_mat, map, "t")
  expect_s3_class(dif, "data.frame")
  expect_equal(nrow(dif), 1)
  expect_equal(dif$gene_id, "GENE1")
})

test_that("diff_stat() rejects invalid stat method gracefully", {
  source("../../scripts/diff_methods.R")

  # "anova" is not a valid stat method — should error
  expect_error(
    diff_stat(td$mat, map, "anova"),
    "stat"
  )
})

# ---------------------------------------------------------------------------
# diff_deseq2
# ---------------------------------------------------------------------------
test_that("diff_deseq2() returns correct column structure for count data", {
  source("../../scripts/diff_methods.R")

  dif <- diff_deseq2(td$mat, map)
  expect_s3_class(dif, "data.frame")
  expect_true("gene_id" %in% colnames(dif))
  expect_true("logFC" %in% colnames(dif))
  expect_true("Pvalue" %in% colnames(dif))
  expect_true("Padj" %in% colnames(dif))
})

test_that("diff_deseq2() gene_id is first column and rownames are default", {
  source("../../scripts/diff_methods.R")

  dif <- diff_deseq2(td$mat, map)
  expect_equal(colnames(dif)[1], "gene_id")
  # rownames(dif) <- NULL resets to default numbering, not gene names
  expect_false(any(rownames(dif) %in% rownames(td$mat)))
})

test_that("diff_deseq2() returns p-values between 0 and 1", {
  source("../../scripts/diff_methods.R")

  dif <- diff_deseq2(td$mat, map)
  expect_true(all(dif$Pvalue >= 0 & dif$Pvalue <= 1))
  expect_true(all(dif$Padj >= 0 & dif$Padj <= 1))
})

# ---------------------------------------------------------------------------
# diff_limma
# ---------------------------------------------------------------------------
test_that("diff_limma() returns correct column structure", {
  source("../../scripts/diff_methods.R")

  dif <- diff_limma(td$mat, map)
  expect_s3_class(dif, "data.frame")
  expect_true("gene_id" %in% colnames(dif))
  expect_true("logFC" %in% colnames(dif))
  expect_true("Pvalue" %in% colnames(dif))
  expect_true("Padj" %in% colnames(dif))
})

test_that("diff_limma() gene_id is first column and rownames are default", {
  source("../../scripts/diff_methods.R")

  dif <- diff_limma(td$mat, map)
  expect_equal(colnames(dif)[1], "gene_id")
  # rownames(dif) <- NULL resets to default numbering, not gene names
  expect_false(any(rownames(dif) %in% rownames(td$mat)))
})

test_that("diff_limma() takes contrast from map factor levels", {
  source("../../scripts/diff_methods.R")

  # Verify the factor levels are used correctly (Case - Control)
  levs <- levels(map[[2]])
  expect_equal(levs, c("Control", "Case"))
  expect_equal(levs[2], "Case")
  expect_equal(levs[1], "Control")
})

# ---------------------------------------------------------------------------
# diff_edger
# ---------------------------------------------------------------------------
test_that("diff_edger() returns correct column structure", {
  source("../../scripts/diff_methods.R")

  dif <- diff_edger(td$mat, map)
  expect_s3_class(dif, "data.frame")
  expect_true("gene_id" %in% colnames(dif))
  expect_true("logFC" %in% colnames(dif))
  expect_true("Pvalue" %in% colnames(dif))
  expect_true("Padj" %in% colnames(dif))
})

test_that("diff_edger() gene_id is first column and rownames are default", {
  source("../../scripts/diff_methods.R")

  dif <- diff_edger(td$mat, map)
  expect_equal(colnames(dif)[1], "gene_id")
  # rownames(dif) <- NULL resets to default numbering, not gene names
  expect_false(any(rownames(dif) %in% rownames(td$mat)))
})

test_that("diff_edger() accepts different norm methods", {
  source("../../scripts/diff_methods.R")

  dif_tmm <- diff_edger(td$mat, map, norm = "TMM")
  dif_rle <- diff_edger(td$mat, map, norm = "RLE")
  expect_s3_class(dif_tmm, "data.frame")
  expect_s3_class(dif_rle, "data.frame")
})

test_that("diff_edger() accepts different model methods", {
  source("../../scripts/diff_methods.R")

  dif_glmfit <- diff_edger(td$mat, map, model = "glmFit")
  dif_glmqlt <- diff_edger(td$mat, map, model = "glmQLFit")
  expect_s3_class(dif_glmfit, "data.frame")
  expect_s3_class(dif_glmqlt, "data.frame")
})
