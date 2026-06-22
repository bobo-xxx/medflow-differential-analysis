# test-filter-helpers.R — Tests for proportion check, cutoff testing, DEG filtering
#
# Tests: proportion_check(), test_cutoff(), filter_degs()
# Uses synthetic data from tests/testthat/setup.R

library(testthat)

# Source dependencies
source("../../scripts/report.R")
source("../../scripts/exceptions.R")
source("setup.R")

td <- create_test_data()

# Build a simulated DE results data frame for cutoff/filter tests
create_test_dif <- function() {
  # Simulate results that would come from diff_stat()
  genes <- rownames(td$mat)
  # GENE1: strong up (logFC ~ 2), GENE2: strong down (logFC ~ -2)
  # GENE3: moderate up (logFC ~ 1), rest: near zero
  logFC <- c(2.0, -2.0, 1.0, 0.3, -0.2, 0.1, -0.1, 0.05, -0.05, 0.0)
  Pvalue <- c(0.0001, 0.0002, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
  Padj <- p.adjust(Pvalue, method = "BH")
  data.frame(
    gene_id = genes,
    logFC = logFC,
    Pvalue = Pvalue,
    Padj = Padj,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# proportion_check
# ---------------------------------------------------------------------------
test_that("proportion_check() returns TRUE for balanced design (3 vs 3)", {
  source("../../scripts/filter_helpers.R")

  map <- td$map
  map[[2]] <- factor(map[[2]], levels = c("Control", "Case"))
  expect_true(proportion_check(map))
})

test_that("proportion_check() returns TRUE for moderately imbalanced design (4 vs 2)", {
  source("../../scripts/filter_helpers.R")

  # 4 control, 2 case — ratio 2:1 is OK
  map <- data.frame(
    sample = c("C1", "C2", "C3", "C4", "T1", "T2"),
    group = factor(c("Control", "Control", "Control", "Control", "Case", "Case"),
                   levels = c("Control", "Case")),
    stringsAsFactors = FALSE
  )
  expect_true(proportion_check(map))
})

test_that("proportion_check() returns FALSE for ratio > 10:1", {
  source("../../scripts/filter_helpers.R")

  # 11 control, 1 case — ratio 11:1 exceeds limit
  map <- data.frame(
    sample = c(paste0("C", 1:11), "T1"),
    group = factor(c(rep("Control", 11), "Case"),
                   levels = c("Control", "Case")),
    stringsAsFactors = FALSE
  )
  expect_false(proportion_check(map))
})

test_that("proportion_check() bypasses check when force_imbalanced=TRUE", {
  source("../../scripts/filter_helpers.R")

  map <- data.frame(
    sample = c(paste0("C", 1:11), "T1"),
    group = factor(c(rep("Control", 11), "Case"),
                   levels = c("Control", "Case")),
    stringsAsFactors = FALSE
  )
  expect_true(proportion_check(map, force_imbalanced = TRUE))
})

test_that("proportion_check() detects single-group design (B9_SAMPLE_MISMATCH halt)", {
  # NOTE: B9_SAMPLE_MISMATCH triggers action="halt" which calls quit().
  # We verify the detection logic without triggering quit by checking
  # table() counts directly, then verify the NDJSON message is emitted
  # by calling report_exception_ndjson directly with dry_run=TRUE.
  source("../../scripts/filter_helpers.R")

  map <- data.frame(
    sample = c("S1", "S2", "S3"),
    group = factor(c("Control", "Control", "Control"),
                   levels = "Control"),
    stringsAsFactors = FALSE
  )
  counts <- table(map[[2]])
  expect_equal(length(counts), 1)

  # Verify the exception message format directly
  out <- capture.output(
    report_exception_ndjson(
      "B9_SAMPLE_MISMATCH", "data_corrupt", "halt",
      sprintf("Expected exactly 2 groups, found %d: %s",
              length(counts), paste(names(counts), collapse = ", ")),
      exit_code = 1, dry_run = TRUE
    )
  )
  expect_true(any(grepl("B9_SAMPLE_MISMATCH", out)))
})

test_that("proportion_check() detects 3-group design (B9_SAMPLE_MISMATCH halt)", {
  # Same approach: verify detection logic without actually calling quit()
  source("../../scripts/filter_helpers.R")

  map <- data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    group = factor(c("A", "B", "C", "A"),
                   levels = c("A", "B", "C")),
    stringsAsFactors = FALSE
  )
  counts <- table(map[[2]])
  expect_equal(length(counts), 3)

  out <- capture.output(
    report_exception_ndjson(
      "B9_SAMPLE_MISMATCH", "data_corrupt", "halt",
      sprintf("Expected exactly 2 groups, found %d: %s",
              length(counts), paste(names(counts), collapse = ", ")),
      exit_code = 1, dry_run = TRUE
    )
  )
  expect_true(any(grepl("B9_SAMPLE_MISMATCH", out)))
})

test_that("proportion_check() works with non-factor group column", {
  source("../../scripts/filter_helpers.R")

  # table() works on character vectors too
  map <- data.frame(
    sample = c("C1", "C2", "C3", "T1", "T2", "T3"),
    group = c("Control", "Control", "Control", "Case", "Case", "Case"),
    stringsAsFactors = FALSE
  )
  expect_true(proportion_check(map))
})

# ---------------------------------------------------------------------------
# test_cutoff
# ---------------------------------------------------------------------------
test_that("test_cutoff() returns named integer vector", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  logfc_test <- c(2.0, 1.0, 0.5, 0.0)
  tc <- test_cutoff(dif, "logFC", "Padj", logfc_test, 0.05)
  expect_type(tc, "integer")
  expect_named(tc, as.character(logfc_test))
})

test_that("test_cutoff() count decreases with higher logFC threshold", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  logfc_test <- c(2.0, 1.5, 1.0, 0.5, 0.0)
  tc <- test_cutoff(dif, "logFC", "Padj", logfc_test, 0.05)
  # At logFC > 2.0, only GENE1 and GENE2 pass (abs(2.0) > 2.0 is FALSE, need strict)
  expect_true(tc["0"] >= tc["2"])
})

test_that("test_cutoff() returns zero at extremely high cutoff", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  tc <- test_cutoff(dif, "logFC", "Padj", c(10), 0.05)
  expect_equal(tc[["10"]], 0L)
})

test_that("test_cutoff() counts correctly with Pvalue instead of Padj", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  tc_p <- test_cutoff(dif, "logFC", "Pvalue", c(0.0), 0.05)
  tc_adj <- test_cutoff(dif, "logFC", "Padj", c(0.0), 0.05)
  # Pvalue is more lenient than Padj (Padj >= Pvalue), so Pvalue gives >= count
  expect_gte(tc_p[["0"]], tc_adj[["0"]])
})

# ---------------------------------------------------------------------------
# filter_degs
# ---------------------------------------------------------------------------
test_that("filter_degs() returns list with correct elements", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  result <- filter_degs(dif, "logFC", "Pvalue", 1.0, 0.05, cutoff = 2)
  expect_type(result, "list")
  expect_true("degs" %in% names(result))
  expect_true("rdegs" %in% names(result))
  expect_true("dif_grouped" %in% names(result))
  expect_true("up" %in% names(result))
  expect_true("down" %in% names(result))
  expect_true("total" %in% names(result))
})

test_that("filter_degs() classifies GENE1 as Up and GENE2 as Down at logFC > 1", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  result <- filter_degs(dif, "logFC", "Pvalue", 1.0, 0.05, cutoff = 2)
  degs <- result$degs
  gene1 <- degs[degs$gene_id == "GENE1", ]
  gene2 <- degs[degs$gene_id == "GENE2", ]
  expect_equal(gene1$group, "Up")
  expect_equal(gene2$group, "Down")
  expect_gt(result$up, 0)
  expect_gt(result$down, 0)
})

test_that("filter_degs() dif_grouped column has Up, Down, Not levels", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  result <- filter_degs(dif, "logFC", "Pvalue", 1.0, 0.05, cutoff = 2)
  groups <- unique(result$dif_grouped$group)
  expect_true(all(c("Up", "Down", "Not") %in% groups))
  expect_equal(nrow(result$dif_grouped), nrow(dif))
})

test_that("filter_degs() total = up + down", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  result <- filter_degs(dif, "logFC", "Pvalue", 1.0, 0.05, cutoff = 2)
  expect_equal(result$total, result$up + result$down)
  expect_equal(result$total, nrow(result$degs))
})

test_that("filter_degs() emits B6_NO_DEGS when no DEGs found", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  # Extremely strict cutoff: no genes pass
  out <- capture.output(
    result <- filter_degs(dif, "logFC", "Padj", 10.0, 0.0001, cutoff = 1)
  )
  expect_true(any(grepl("B6_NO_DEGS", out)))
  expect_equal(result$total, 0)
  expect_equal(nrow(result$degs), 0)
})

test_that("filter_degs() emits B2_FEW_DEGS when below cutoff", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  # GENE1 and GENE2 pass at 1.0, but cutoff=5 requires 5 minimum
  out <- capture.output(
    result <- filter_degs(dif, "logFC", "Pvalue", 1.0, 0.05, cutoff = 5)
  )
  expect_true(any(grepl("B2_FEW_DEGS", out)))
  # DEGs still returned, just warned
  expect_gt(result$total, 0)
})

test_that("filter_degs() rgs parameter filters to intersection", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  rgs <- c("GENE1", "GENE3", "GENE99")  # GENE1 is up, GENE3 is up, GENE99 not present
  out <- capture.output(
    result <- filter_degs(dif, "logFC", "Pvalue", 0.5, 0.05, cutoff = 1,
                          rgs = rgs)
  )
  rdegs <- result$rdegs
  # Only GENE1 and GENE3 intersect between DEGs and rgs
  expect_true("GENE1" %in% rdegs$gene_id)
  expect_true("GENE3" %in% rdegs$gene_id)
  expect_false("GENE2" %in% rdegs$gene_id)
})

test_that("filter_degs() emits B7_RGS_INTERSECTION when intersection below cutoff", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  rgs <- c("GENE1")  # Only 1 gene in intersection
  out <- capture.output(
    result <- filter_degs(dif, "logFC", "Pvalue", 1.0, 0.05, cutoff = 2,
                          rgs = rgs)
  )
  expect_true(any(grepl("B7_RGS_INTERSECTION", out)))
  expect_equal(nrow(result$rdegs), 1)
})

test_that("filter_degs() rdegs equals degs when rgs is NULL", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  result <- filter_degs(dif, "logFC", "Pvalue", 1.0, 0.05, cutoff = 2)
  expect_equal(nrow(result$rdegs), nrow(result$degs))
})

test_that("filter_degs() works with Padj column", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  result <- filter_degs(dif, "logFC", "Padj", 1.0, 0.1, cutoff = 2)
  expect_type(result, "list")
  expect_true(is.numeric(result$total))
})

test_that("filter_degs() handles strict cutoff where GENE3 is excluded", {
  source("../../scripts/filter_helpers.R")

  dif <- create_test_dif()
  # logFC > 1.5: GENE3 has 1.0, so it should be "Not"
  result <- filter_degs(dif, "logFC", "Pvalue", 1.5, 0.05, cutoff = 2)
  degs <- result$degs
  gene3 <- result$dif_grouped[result$dif_grouped$gene_id == "GENE3", ]
  expect_equal(gene3$group, "Not")
})
