# setup.R — Synthetic test fixtures for deg-analysis tests
#
# Provides create_test_data() which returns a list with:
#   $mat  — 10-gene x 6-sample count matrix (3 case, 3 control)
#   $map  — sample-to-group mapping CSV path (on disk)
#   $mat_path — expression matrix CSV path (on disk)
#   $expected_up   — gene names expected to be upregulated
#   $expected_down — gene names expected to be downregulated

create_test_data <- function(tmpdir = tempdir()) {
  # 10 genes, 6 samples
  genes <- paste0("GENE", 1:10)
  samples <- c("CASE1", "CASE2", "CASE3", "CTRL1", "CTRL2", "CTRL3")

  set.seed(42)

  # Base counts: Poisson with lambda=100 for all genes
  mat <- matrix(rpois(10 * 6, lambda = 100), nrow = 10, ncol = 6,
                dimnames = list(genes, samples))

  # GENE1: 4x up in case (logFC ~ 2)
  mat["GENE1", 1:3] <- mat["GENE1", 1:3] * 4

  # GENE2: 4x down in case (logFC ~ -2)
  mat["GENE2", 4:6] <- mat["GENE2", 4:6] * 3

  # GENE3: subtle 2x up (logFC ~ 1) — near cutoff
  mat["GENE3", 1:3] <- mat["GENE3", 1:3] * 2

  mat_path <- file.path(tmpdir, "test_expr.csv")
  write.csv(mat, mat_path)

  map <- data.frame(
    sample = samples,
    group = c(rep("Case", 3), rep("Control", 3)),
    stringsAsFactors = FALSE
  )
  map_path <- file.path(tmpdir, "test_map.csv")
  write.csv(map, map_path, row.names = FALSE)

  list(
    mat = mat,
    map = map,
    mat_path = mat_path,
    map_path = map_path,
    expected_up = "GENE1",
    expected_down = "GENE2"
  )
}
