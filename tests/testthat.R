# testthat.R — Test runner for deg-analysis node
# This is not an R package, so we source test files directly.

library(testthat)

# Set working directory to project root for consistent source paths
# Tests reference ../../scripts/ from tests/testthat/

test_dir("testthat", reporter = "summary")
