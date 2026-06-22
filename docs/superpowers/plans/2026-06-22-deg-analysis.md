---
change: create-node-package
design-doc: docs/superpowers/specs/2026-06-22-differential-analysis-design.md
base-ref: 1e5d31cd6609ffc348f4701d4024995b706dedfd
archived-with: 2026-06-22-create-node-package
---

# differential-analysis Node Package Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Build a complete IRE node package for differential expression analysis (case vs control comparison) with 5 DE methods, 6 output files, NDJSON reporting, and full test coverage.

**Architecture:** Single R entry point (`scripts/main.R`) with subcommand dispatch (`run`, `validate-input`, `validate-output`). Seven helper modules sourced in dependency order: report.R, exceptions.R, io_helpers.R, diff_methods.R, filter_helpers.R, plot_helpers.R, input_validation.R, output_validation.R. All output to stdout is valid NDJSON; the framework captures it.

**Tech Stack:** R 4.3, Bioconductor (DESeq2, limma, edgeR, statmod), ggplot2, ggrepel, ggvenn, pheatmap, RCircos, jsonlite, yaml, filelock, dplyr, data.table, testthat

## Global Constraints

- R version >= 4.3
- All required packages installed via conda from conda-forge + bioconda channels
- Entry point: `scripts/main.R` with subcommand dispatch
- All stdout is valid NDJSON (one JSON object per line)
- Exit codes: 0 (success), 1 (data/write error), 3 (environment error)
- English is the working language for all files, comments, messages, and NDJSON output
- No hardcoded secrets or paths
- Flat directory layout per node-package.md protocol
- SKILL.md v2 frontmatter with all 8 required sections
- Node name in .run_result.json: `differential-analysis`
- All parameters declared in SKILL.md frontmatter with bind annotations (upstream/config/static/framework)
- Exceptions follow exception-contract.md: halt | skip_with_warning | escalate (no retry)
- No `report_prompt()` in this node (not a gate node; no interactive prompts)

archived-with: 2026-06-22-create-node-package
---

## File Structure

```
differential-analysis/
├── SKILL.md                          # [Task 6] Agent contract with v2 YAML frontmatter
├── envs/env-r-4.3.yaml               # [Task 1] Conda env with all R/Bioconductor deps
├── scripts/
│   ├── main.R                        # [Task 5] Single entry point, arg parsing, dispatch
│   ├── report.R                      # [Task 2] NDJSON reporting helpers
│   ├── exceptions.R                  # [Task 2] Exception accumulator and report_exception_ndjson()
│   ├── io_helpers.R                  # [Task 3] file_lock(), create_file_dir(), color_map()
│   ├── diff_methods.R                # [Task 3] diff_deseq2(), diff_limma(), diff_edger(), diff_stat()
│   ├── filter_helpers.R              # [Task 3] proportion_check(), test_cutoff(), filter_degs()
│   ├── plot_helpers.R                # [Task 3] plot_volcano(), plot_heatmap(), plot_venn(), plot_locate()
│   ├── input_validation.R            # [Task 4] Pre-flight: file existence, columns, groups, proportion
│   └── output_validation.R           # [Task 4] Post-hoc: CSV columns, non-empty, plot file existence
├── tests/
│   └── testthat/
│       ├── setup.R                   # [Task 1] Synthetic fixtures (10-gene x 6-sample)
│       ├── test-input-validation.R   # [Task 7] Input validation tests
│       ├── test-output-validation.R  # [Task 7] Output validation tests
│       ├── test-diff-methods.R       # [Task 7] DE method output tests
│       └── test-main.R               # [Task 7] Dispatch, parsing, end-to-end tests
└── openspec/                         # [Existing] Governance artifacts
```

**Module dependency order for sourcing in main.R:**
1. `report.R` (no internal deps)
2. `exceptions.R` (depends on report.R for .exceptions list init and jsonlite)
3. `io_helpers.R` (depends on filelock, no report deps)
4. `diff_methods.R` (depends on DESeq2/limma/edgeR/statmod packages only)
5. `filter_helpers.R` (no internal deps beyond dplyr)
6. `plot_helpers.R` (depends on ggplot2/ggrepel/ggvenn/pheatmap/RCircos)
7. `input_validation.R` (depends on filter_helpers.R for proportion_check)
8. `output_validation.R` (no internal deps beyond base R)

archived-with: 2026-06-22-create-node-package
---

## Task 1: Environment Setup

**Files:**
- Modify: `envs/env-r-4.3.yaml` (expand from skeleton)
- Create: `tests/testthat/setup.R`

**Interfaces:**
- Produces: Conda environment with all packages, synthetic test fixtures accessible via `create_test_data()`

### 1.1 Expand conda environment YAML

- [x] **Step 1: Write expanded envs/env-r-4.3.yaml**

The existing skeleton only has `r-base=4.3` and `r-essentials`. Expand it with all required Bioconductor packages and R libraries:

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base=4.3
  - r-essentials
  - bioconductor-deseq2
  - bioconductor-limma
  - bioconductor-edger
  - bioconductor-statmod
  - r-ggplot2
  - r-ggrepel
  - r-pheatmap
  - r-dplyr
  - r-data.table
  - r-jsonlite
  - r-yaml
  - r-filelock
  - r-testthat
  - r-r.utils
```

Note: `ggvenn` and `RCircos` are not available on conda-forge/bioconda in all channel configurations. Add a comment explaining these will be installed via `install.packages()` if conda install fails, and provide fallback instructions:

```yaml
# Packages requiring install.packages() fallback if conda resolution fails:
#   install.packages(c("ggvenn", "RCircos"), repos = "https://cloud.r-project.org")
# Or via uv:
#   uv pip install r-ggvenn r-rcircos  (if uv handles R packages)
```

- [x] **Step 2: Verify environment creation**

```bash
conda env create -f envs/env-r-4.3.yaml -p ./env --dry-run 2>&1 | head -20
```

Expected: conda resolves all packages without conflicts. If `ggvenn` or `RCircos` fail, record the fallback in the yaml comment.

- [x] **Step 3: Create environment for real**

```bash
conda env create -f envs/env-r-4.3.yaml -p ./env
```

Expected: Environment created successfully. Verify with:

```bash
./env/bin/Rscript -e 'library(DESeq2); library(limma); library(edgeR); library(ggplot2); library(ggrepel); library(pheatmap); library(dplyr); library(data.table); library(jsonlite); library(yaml); library(filelock); library(testthat); cat("All packages loaded\n")'
```

Expected stdout: `All packages loaded`

If ggvenn or RCircos are not available via conda, install them with `install.packages()` inside the env's R:

```bash
./env/bin/Rscript -e 'install.packages(c("ggvenn", "RCircos"), repos = "https://cloud.r-project.org")'
```

- [x] **Step 4: Commit environment**

```bash
git add envs/env-r-4.3.yaml
git commit -m "feat: add expanded conda environment with all Bioconductor/R dependencies

Bioconductor: DESeq2, limma, edgeR, statmod
R packages: ggplot2, ggrepel, ggvenn, pheatmap, RCircos,
            dplyr, data.table, jsonlite, yaml, filelock,
            testthat, R.utils

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### 1.2 Create test fixtures

- [x] **Step 5: Write tests/testthat/setup.R**

Create a `create_test_data()` function that generates a synthetic 10-gene x 6-sample count matrix with 3 case + 3 control samples and known fold changes for 2 genes:

```r
# setup.R — Synthetic test fixtures for differential-analysis tests
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
```

- [x] **Step 6: Verify setup.R loads cleanly**

```bash
./env/bin/Rscript -e 'source("tests/testthat/setup.R"); td <- create_test_data(); stopifnot(nrow(td$mat) == 10); stopifnot(ncol(td$mat) == 6); stopifnot(td$expected_up == "GENE1"); cat("setup.R OK\n")'
```

Expected stdout: `setup.R OK`

- [x] **Step 7: Commit test fixtures**

```bash
git add tests/testthat/setup.R
git commit -m "feat: add synthetic test fixtures (10-gene x 6-sample count matrix)

Creates create_test_data() for reproducible test input with known
fold changes: GENE1 up ~4x, GENE2 down ~3x, GENE3 subtle ~2x.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

archived-with: 2026-06-22-create-node-package
---

## Task 2: Reporting & Exceptions

**Files:**
- Create: `scripts/report.R`
- Create: `scripts/exceptions.R`

**Interfaces:**
- Consumes: jsonlite, filelock (packages from env)
- Produces:
  - `report_info(msg, ...)` — writes `{"level":"info","msg":<msg>,...}` NDJSON to stdout
  - `report_result(status, files, metadata, ...)` — writes `{"level":"result","status":<status>,...}` NDJSON to stdout
  - `report_error(msg, exit_code=1)` — writes `{"level":"error","msg":<msg>}` NDJSON to stdout, then `quit(status=exit_code)`
  - `report_exception_ndjson(code, nature, action, msg, exit_code, dry_run)` — writes structured exception NDJSON, accumulates to `.exceptions`, halts if action=="halt"
  - `write_run_result(out_dir, result, params, exit_code, times)` — writes `.run_result.json` to out_dir
  - `.exceptions` list (global, initialized here)

### 2.1 Create scripts/report.R

- [x] **Step 1: Write scripts/report.R**

Adapt from reference node `/tmp/medflow-geo-microarray/scripts/report.R`. Remove `report_prompt()` (this node has no interactive prompts per exception-contract.md). The node name in `write_run_result()` is `"differential-analysis"`:

```r
# report.R — NDJSON reporting helpers for differential-analysis
#
# All output to stdout is valid NDJSON. Each function writes one line.
# Use report_info() for progress, report_result() for final output,
# and report_error() for terminal failures.

# Exception accumulator (filled by report_exception_ndjson during run)
.exceptions <- list()

#' Write an info-level NDJSON message to stdout
#'
#' @param msg Character string with progress message
#' @param ... Additional named fields to include in the JSON object
report_info <- function(msg, ...) {
  extra <- list(...)
  obj <- c(list(level = "info", msg = msg), extra)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}

#' Write a result-level NDJSON message to stdout
#'
#' @param status Character status code ("success", "error", etc.)
#' @param files List of file info (each with path, and optionally rows/cols)
#' @param metadata List of metadata (method, version, n_degs, up, down)
#' @param ... Additional named fields
report_result <- function(status, files = list(), metadata = list(), ...) {
  extra <- list(...)
  obj <- c(list(level = "result", status = status,
                files = files, metadata = metadata), extra)
  if (length(files) == 0) obj$files <- NULL
  if (length(metadata) == 0) obj$metadata <- NULL
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}

#' Write run provenance file for agent inspection
#'
#' Creates .run_result.json in the output directory with full provenance:
#' parameters, output metadata, exceptions, file list, timestamps.
#'
#' @param out_dir Output directory
#' @param result Result list from do_run()
#' @param params Parameter list
#' @param exit_code Integer exit code
#' @param times Character vector of c(started_at, finished_at) ISO8601 timestamps
write_run_result <- function(out_dir, result, params, exit_code, times) {
  output <- list()
  if (!is.null(result$metadata) && length(result$metadata) > 0) {
    for (name in names(result$metadata)) {
      output[[name]] <- result$metadata[[name]]
    }
  }

  files <- list()
  for (fname in c("Diffanalysis.csv", "DEGs.csv", "Volcano.pdf",
                   "Heatmap.pdf", "Venn.pdf", "Chromosome_location.pdf")) {
    fpath <- file.path(out_dir, fname)
    if (file.exists(fpath)) {
      finfo <- list(path = fpath)
      if (endsWith(fname, ".csv")) {
        df <- tryCatch(read.csv(fpath), error = function(e) NULL)
        if (!is.null(df)) {
          finfo$rows <- nrow(df)
          finfo$cols <- ncol(df)
        }
      }
      files <- c(files, list(finfo))
    }
  }

  clean_params <- list()
  for (k in names(params)) {
    if (!is.null(params[[k]]) && k != "outdir") {
      clean_params[[k]] <- params[[k]]
    }
  }

  obj <- list(
    node         = "differential-analysis",
    subcommand   = if (is.null(params$subcommand)) "unknown" else params$subcommand,
    status       = if (is.null(result$status)) "unknown" else result$status,
    exit_code    = exit_code,
    started_at   = times[1],
    finished_at  = times[2],
    parameters   = clean_params,
    output       = output,
    exceptions   = if (length(.exceptions) > 0) .exceptions else list(),
    files        = files
  )

  json_path <- file.path(out_dir, ".run_result.json")
  writeLines(jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE), json_path)
}

#' Write an error-level NDJSON message and exit
#'
#' @param msg Error message
#' @param exit_code Integer exit code (default: 1)
report_error <- function(msg, exit_code = 1) {
  obj <- list(level = "error", msg = msg)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
  quit(status = exit_code)
}
```

- [x] **Step 2: Verify report.R loads cleanly**

```bash
./env/bin/Rscript -e 'source("scripts/report.R"); report_info("test", x=1); cat("report.R OK\n")'
```

Expected stdout:
```
{"level":"info","msg":"test","x":1}
report.R OK
```

- [x] **Step 3: Commit report.R**

```bash
git add scripts/report.R
git commit -m "feat: add NDJSON reporting helpers (report.R)

Functions: report_info(), report_result(), report_error(),
write_run_result(). No report_prompt() — this node has no
interactive prompts per exception-contract.md.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### 2.2 Create scripts/exceptions.R

- [x] **Step 4: Write scripts/exceptions.R**

Adapt from reference node. Remove checkpoint and network-specific functions not needed for DEG analysis. Keep `check_environment()`, `report_exception_ndjson()`:

```r
# exceptions.R — Structured exception handling for differential-analysis
#
# Categories: B=Data, W=Write, E=Environment
# All exceptions flow through report_exception_ndjson() for machine-readable output.

# Exception accumulator for .run_result.json (also initialized in report.R)
if (!exists(".exceptions")) .exceptions <- list()

#' Check runtime environment
#'
#' Verifies required packages are installed.
#'
#' @return List with status ("ok", "error"), missing packages, and details
check_environment <- function() {
  required <- c("DESeq2", "limma", "edgeR", "statmod",
                "ggplot2", "ggrepel", "pheatmap",
                "dplyr", "data.table", "jsonlite",
                "yaml", "filelock")
  missing  <- character(0)
  warnings <- character(0)

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

  # Optional packages
  for (pkg in c("ggvenn", "RCircos")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warnings <- c(warnings, paste("Optional package not available:", pkg))
    }
  }

  list(status = "ok", msg = "Environment OK", warnings = warnings)
}

#' Report an exception as structured NDJSON
#'
#' @param code Exception code (e.g., "B1_PROPORTION")
#' @param nature Exception nature (data_insufficient, data_corrupt, data_mismatch, env_bug)
#' @param action Response action (halt, skip_with_warning, escalate)
#' @param msg Human-readable message
#' @param exit_code Integer exit code (default: 1)
#' @param dry_run Logical, if TRUE do not quit (for testing)
report_exception_ndjson <- function(code, nature, action, msg,
                                     exit_code = 1, dry_run = FALSE) {
  level <- switch(action,
    skip_with_warning = "decision",
    escalate          = "exception",
    "exception"
  )

  obj <- list(
    level  = level,
    code   = code,
    nature = nature,
    action = action,
    msg    = msg
  )
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")

  # Accumulate for .run_result.json
  .exceptions <<- c(.exceptions, list(obj))

  if (action == "halt" && !dry_run) {
    quit(status = exit_code)
  }
}
```

- [x] **Step 5: Verify exceptions.R loads cleanly and check_environment() works**

```bash
./env/bin/Rscript -e 'source("scripts/report.R"); source("scripts/exceptions.R"); env <- check_environment(); stopifnot(env$status == "ok"); cat("exceptions.R OK, env status:", env$status, "\n")'
```

Expected stdout: `exceptions.R OK, env status: ok`

- [x] **Step 6: Test report_exception_ndjson output**

```bash
./env/bin/Rscript -e 'source("scripts/report.R"); source("scripts/exceptions.R"); report_exception_ndjson("B1_PROPORTION", "data_insufficient", "skip_with_warning", "Sample proportion 15:1 exceeds 10:1 limit", dry_run=TRUE); cat("NDJSON exception test OK\n")'
```

Expected stdout (first line must be valid JSON with the exception):
```json
{"level":"decision","code":"B1_PROPORTION","nature":"data_insufficient","action":"skip_with_warning","msg":"Sample proportion 15:1 exceeds 10:1 limit"}
```

- [x] **Step 7: Commit exceptions.R**

```bash
git add scripts/exceptions.R
git commit -m "feat: add structured exception handling (exceptions.R)

check_environment() verifies 13 required packages.
report_exception_ndjson() outputs structured NDJSON with
code/nature/action/msg fields and accumulates to .exceptions.
Supports dry_run mode for testing.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

archived-with: 2026-06-22-create-node-package
---

## Task 3: Internal Helper Scripts

**Files:**
- Create: `scripts/io_helpers.R`
- Create: `scripts/diff_methods.R`
- Create: `scripts/filter_helpers.R`
- Create: `scripts/plot_helpers.R`

**Interfaces:**
- Consumes: DESeq2, limma, edgeR, statmod, ggplot2, ggrepel, pheatmap, dplyr, data.table, yaml, filelock (packages from env)
- Produces:
  - `file_lock(path, FUN, ...)` — advisory file locking wrapper
  - `create_file_dir(file)` — mkdir -p for output paths
  - `color_map(colors, groups)` — assign colors to named groups
  - `diff_deseq2(mat, group)` — DESeq2 DE analysis, returns data.frame with logFC, Pvalue, Padj columns
  - `diff_limma(mat, map)` — limma DE analysis, same return shape
  - `diff_edger(mat, map, norm, model)` — edgeR DE analysis, same return shape
  - `diff_stat(mat, map, stat)` — t-test or Wilcoxon, same return shape
  - `proportion_check(map, force_imbalanced)` — returns TRUE/FALSE, raises B1_PROPORTION if ratio > 10:1
  - `test_cutoff(dif, fc_name, p_name, logfc_test, p_value)` — returns named integer vector of gene counts at each cutoff
  - `filter_degs(dif, fc_name, p_name, logfc_cutoff, p_value, cutoff)` — returns list(degs, rdegs, dif_grouped)
  - `plot_volcano(dif, p_name, p_value, logfc_cutoff, top, gene, outfile)` — writes PDF via ggplot2
  - `plot_heatmap(mat, map, rdegs, top, color_heat, outfile)` — writes PDF via pheatmap
  - `plot_venn(dif, rgs, pheno_abbr, color_panel, outfile)` — writes PDF via ggvenn (called only when --rgs provided)
  - `plot_locate(gene_list, locate, tax_id, outfile)` — writes PDF via RCircos (called only when --locate provided)

### 3.1 Create scripts/io_helpers.R

- [x] **Step 1: Write scripts/io_helpers.R**

Extract the `file_lock` and `create_file_dir` utility functions from the original diff.R, and `color_map` from original diff_venn.R:

```r
# io_helpers.R — File I/O and color utilities for differential-analysis
#
# Provides advisory file locking, directory creation, and color mapping.

suppressPackageStartupMessages(library(filelock))

#' Acquire an exclusive file lock, run FUN, then release
#'
#' @param path Path to the file to lock
#' @param FUN Function to call with path as first argument
#' @param ... Additional arguments passed to FUN
#' @param exclusive Logical, whether lock is exclusive (default TRUE)
#' @param timeout Timeout in milliseconds (default 5000)
#' @return Result of FUN(path, ...), invisibly
file_lock <- function(path, FUN, ..., exclusive = TRUE, timeout = 5000) {
  FUN <- match.fun(FUN)
  lock_file <- paste0(path, ".lock")
  lock <- lock(lock_file, exclusive = exclusive, timeout = timeout)
  unlock <- unlock
  if (is.null(lock)) {
    stop(paste0("The file lock cannot be obtained: ", lock_file))
  } else {
    res <- tryCatch(
      forceAndCall(1, FUN, path, ...),
      error = function(e) stop(e),
      finally = unlock(lock)
    )
  }
  invisible(res)
}

#' Create parent directory for a file if it does not exist
#'
#' @param file File path whose parent directory should exist
create_file_dir <- function(file) {
  path <- dirname(file)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

#' Map colors to named groups
#'
#' Cycles colors if there are more groups than colors.
#'
#' @param color Character vector of color values
#' @param group Character vector of group names
#' @return Named character vector mapping group -> color
color_map <- function(color, group) {
  color <- as.vector(color)
  group <- as.vector(group)
  g <- unique(group)
  n <- length(g)
  color <- head(rep(color, ceiling(n / length(color))), n)
  names(color) <- g
  color
}
```

- [x] **Step 2: Verify io_helpers.R loads**

```bash
./env/bin/Rscript -e 'source("scripts/io_helpers.R"); create_file_dir("/tmp/test_deg/test.txt"); stopifnot(dir.exists("/tmp/test_deg")); cm <- color_map(c("red","blue"), c("A","B","C")); stopifnot(cm[["A"]] == "red"); stopifnot(cm[["C"]] == "red"); cat("io_helpers.R OK\n")'
```

Expected stdout: `io_helpers.R OK`

- [x] **Step 3: Commit io_helpers.R**

```bash
git add scripts/io_helpers.R
git commit -m "feat: add I/O helpers (io_helpers.R)

file_lock() — advisory file locking via filelock package
create_file_dir() — mkdir -p for output paths
color_map() — assigns colors to named groups, cycles if needed

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### 3.2 Create scripts/diff_methods.R

- [x] **Step 4: Write scripts/diff_methods.R**

Adapt the 4 DE functions from `original/scripts/diff.R`. Remove the `deal_data()` wrapper (data loading is now handled by main.R). Each function takes a pre-loaded matrix and group map, returns a data.frame with columns: gene_id (from rownames), logFC, Pvalue, Padj, stat (where applicable):

```r
# diff_methods.R — Differential expression analysis methods for differential-analysis
#
# Provides four DE analysis functions adapted from original/scripts/diff.R.
# Each function takes an expression matrix and sample group map, returns
# a data.frame with columns: logFC, Pvalue, Padj, [stat].
# Gene identifiers are preserved as the first column (gene_id).

#' DESeq2 differential expression (for count data)
#'
#' @param df Expression count matrix (genes x samples), raw counts
#' @param group Data frame with sample-to-group mapping (col 1: sample, col 2: group)
#' @return Data frame with columns: gene_id, logFC, Pvalue, Padj
diff_deseq2 <- function(df, group) {
  suppressPackageStartupMessages(library(DESeq2))

  col_data <- data.frame(
    row.names = colnames(df),
    group_list = group[[2]]
  )
  dds <- DESeqDataSetFromMatrix(
    countData = round(df),
    colData = col_data,
    design = ~group_list
  )
  dds2 <- DESeq(dds)
  res <- results(dds2)
  dif <- res %>%
    as.data.frame() %>%
    dplyr::rename(
      logFC = "log2FoldChange",
      Pvalue = "pvalue",
      Padj = "padj"
    ) %>%
    na.omit() %>%
    dplyr::arrange(Pvalue)

  dif$gene_id <- rownames(dif)
  rownames(dif) <- NULL
  dif <- dif[, c("gene_id", "logFC", "Pvalue", "Padj")]
  dif
}

#' limma differential expression (for normalized microarray or RNA-seq data)
#'
#' @param df Expression matrix (genes x samples)
#' @param map Data frame with sample-to-group mapping (col 1: sample, col 2: group)
#' @return Data frame with columns: gene_id, logFC, Pvalue, Padj
diff_limma <- function(df, map) {
  pdf(NULL)  # Suppress limma implicit plotSA call
  suppressPackageStartupMessages(library(limma))

  treat_name <- levels(map[[2]])[2]
  con_name <- levels(map[[2]])[1]

  design <- model.matrix(~ 0 + map[[2]])
  rownames(design) <- map[[1]]
  colnames(design) <- levels(map[[2]])

  fit <- lmFit(df, design)
  cont_matrix <- makeContrasts(
    contrasts = paste0(treat_name, "-", con_name),
    levels = design
  )
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2)
  plotSA(fit2)

  dif <- topTable(fit2, coef = 1, n = Inf)
  dif <- dif %>%
    as.data.frame() %>%
    dplyr::rename(Pvalue = "P.Value", Padj = "adj.P.Val") %>%
    na.omit() %>%
    dplyr::arrange(Pvalue)

  dev.off()

  dif$gene_id <- rownames(dif)
  rownames(dif) <- NULL
  # Keep logFC, AveExpr, t, B, Pvalue, Padj, gene_id
  dif <- dif[, c("gene_id", "logFC", "Pvalue", "Padj")]
  dif
}

#' edgeR differential expression (for count data)
#'
#' @param df Expression count matrix (genes x samples), raw counts
#' @param map Data frame with sample-to-group mapping
#' @param norm Normalization method: "TMM", "RLE", "upperquartile", "none"
#' @param model Model fitting method: "glmFit" or "glmQLFit"
#' @return Data frame with columns: gene_id, logFC, Pvalue, Padj
diff_edger <- function(df, map, norm = "TMM", model = "glmFit") {
  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(statmod))

  dgelist <- DGEList(counts = df, group = map[[2]])
  keep <- rowSums(cpm(dgelist) > 1) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist_norm <- calcNormFactors(dgelist, method = norm)
  design <- model.matrix(~ map[[2]])
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE)

  func <- get(model)
  fit <- func(dge, design, robust = TRUE)
  lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))

  dif <- lrt %>%
    as.data.frame() %>%
    dplyr::rename(Pvalue = "PValue", Padj = "FDR") %>%
    na.omit() %>%
    dplyr::arrange(Pvalue)

  dif$gene_id <- rownames(dif)
  rownames(dif) <- NULL
  dif <- dif[, c("gene_id", "logFC", "Pvalue", "Padj")]
  dif
}

#' Row-wise fold change calculation
#'
#' @param row Numeric vector of expression values for one gene
#' @param treat Character vector of treatment sample column names
#' @param control Character vector of control sample column names
#' @return log2 fold change (numeric)
logfc_row <- function(row, treat, control) {
  mean_treatment <- mean(row[treat])
  mean_control <- mean(row[control])
  fc <- mean_treatment / mean_control
  if (fc <= 0 || is.na(fc)) return(0)
  log2(fc)
}

#' Row-wise statistical test
#'
#' @param row Numeric vector of expression values
#' @param group Factor vector of group assignments (same length as row)
#' @param func Test function (t.test or wilcox.test)
#' @return Numeric vector c(p.value, statistic)
test_row <- function(row, group, func = t.test) {
  df <- data.frame(x = row, Group = group)
  res <- tryCatch(
    {
      result <- func(x ~ Group, data = df)
      c(result$p.value, result$statistic)
    },
    error = function(e) c(1, 0)
  )
  res
}

#' Statistical test differential expression (t-test or Wilcoxon)
#'
#' @param df Expression matrix (genes x samples)
#' @param map Data frame with sample-to-group mapping
#' @param stat Statistical test: "t" for t-test, "wilcox" for Wilcoxon
#' @return Data frame with columns: gene_id, logFC, Pvalue, Padj, stat
diff_stat <- function(df, map, stat) {
  treat_name <- levels(map[[2]])[2]
  con_name <- levels(map[[2]])[1]
  treat <- map[map[[2]] == treat_name, ][[1]]
  control <- map[map[[2]] == con_name, ][[1]]

  res <- data.frame(logFC = apply(df, 1, logfc_row,
    treat = treat, control = control
  ))

  test_func <- if (stat == "t") t.test else wilcox.test
  df_test <- apply(df, 1, test_row, group = map[[2]], func = test_func)
  df_test <- t(df_test)
  colnames(df_test) <- c("Pvalue", "stat")
  res <- cbind(res, df_test)
  res$Padj <- p.adjust(res$Pvalue, method = "BH")

  res <- res %>%
    na.omit() %>%
    dplyr::arrange(Pvalue)

  res$gene_id <- rownames(res)
  rownames(res) <- NULL
  res <- res[, c("gene_id", "logFC", "Pvalue", "Padj", "stat")]
  res
}
```

- [x] **Step 5: Verify diff_methods.R loads and diff_stat works with synthetic data**

```bash
./env/bin/Rscript -e '
source("scripts/report.R")
source("tests/testthat/setup.R")
td <- create_test_data()
source("scripts/diff_methods.R")
map <- td$map
map[[2]] <- factor(map[[2]], levels = c("Control", "Case"))
dif <- diff_stat(td$mat, map, "t")
stopifnot("gene_id" %in% colnames(dif))
stopifnot("logFC" %in% colnames(dif))
stopifnot("Pvalue" %in% colnames(dif))
stopifnot("Padj" %in% colnames(dif))
stopifnot(dif$gene_id[1] %in% rownames(td$mat))
cat("diff_stat OK, top gene:", dif$gene_id[1], "\n")
'
```

Expected: stdout includes `diff_stat OK` and the top gene is GENE1 (most significant).

- [x] **Step 6: Commit diff_methods.R**

```bash
git add scripts/diff_methods.R
git commit -m "feat: add DE analysis methods (diff_methods.R)

diff_deseq2() — DESeq2 for count data
diff_limma() — limma with empirical Bayes
diff_edger() — edgeR with configurable norm/model
diff_stat() — t-test or Wilcoxon with BH correction

All functions return data.frame with gene_id, logFC, Pvalue, Padj.
Adapted from original/scripts/diff.R.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### 3.3 Create scripts/filter_helpers.R

- [x] **Step 7: Write scripts/filter_helpers.R**

Adapt `proportion_check` from `original/scripts/diff_check_proportion.py`, and `test_cutoff`/`filter_degs` from `original/scripts/diff_filter.R`:

```r
# filter_helpers.R — Proportion check, cutoff testing, DEG filtering for differential-analysis
#
# All functions report exceptions via report_exception_ndjson() when thresholds
# are violated, using the structured exception codes from the design doc.

suppressPackageStartupMessages(library(dplyr))

#' Check sample proportion ratio
#'
#' Ensures the ratio of larger group to smaller group does not exceed 10:1
#' unless force_imbalanced is TRUE. Emits B1_PROPORTION exception on violation.
#'
#' @param map Data frame with sample-to-group mapping (col 2 = group)
#' @param force_imbalanced Logical, if TRUE bypass the check
#' @return Logical TRUE if proportion is acceptable (or forced)
proportion_check <- function(map, force_imbalanced = FALSE) {
  counts <- table(map[[2]])
  if (length(counts) != 2) {
    report_exception_ndjson(
      "B9_SAMPLE_MISMATCH", "data_corrupt", "halt",
      sprintf("Expected exactly 2 groups, found %d: %s",
              length(counts), paste(names(counts), collapse = ", ")),
      exit_code = 1
    )
    return(FALSE)
  }

  c_num <- as.numeric(counts[1])
  t_num <- as.numeric(counts[2])
  ratio <- max(c_num, t_num) / min(c_num, t_num)

  if (ratio > 10 && !force_imbalanced) {
    g1 <- names(counts)[1]
    g2 <- names(counts)[2]
    report_exception_ndjson(
      "B1_PROPORTION", "data_insufficient", "skip_with_warning",
      sprintf("Sample proportion %s (%d samples) vs %s (%d samples), ratio %.1f:1 exceeds 10:1 limit",
              g1, c_num, g2, t_num, ratio),
      exit_code = 1
    )
    return(FALSE)
  }
  TRUE
}

#' Test sensitivity of DEG count at different logFC cutoffs
#'
#' @param dif Data frame with DE results (must contain fc_name and p_name columns)
#' @param fc_name Character, name of fold-change column (e.g., "logFC")
#' @param p_name Character, name of p-value column (e.g., "Padj" or "Pvalue")
#' @param logfc_test Numeric vector of logFC cutoffs to test
#' @param p_value Numeric, p-value threshold
#' @return Named integer vector with DEG counts at each cutoff
test_cutoff <- function(dif, fc_name, p_name, logfc_test, p_value = 0.05) {
  test <- c()
  for (value in logfc_test) {
    test <- append(
      test,
      sum((abs(dif[[fc_name]]) > value) & (dif[[p_name]] < p_value))
    )
  }
  names(test) <- logfc_test
  return(test)
}

#' Filter DEGs by logFC and p-value cutoffs
#'
#' Classifies genes as Up, Down, or Not. Emits B6_NO_DEGS or B2_FEW_DEGS
#' if results are empty or below cutoff.
#'
#' @param dif Data frame with DE results (must contain fc_name and p_name columns, and gene identifiers in column 1)
#' @param fc_name Character, name of fold-change column (e.g., "logFC")
#' @param p_name Character, name of p-value column (e.g., "Padj")
#' @param logfc_cutoff Numeric, absolute logFC threshold
#' @param p_value Numeric, p-value threshold
#' @param cutoff Integer, minimum required DEG count
#' @param rgs Optional character vector of related gene set identifiers for intersection
#' @return List with elements: degs (DEG data frame), rdegs (related DEGs), dif_grouped (full data frame with group column), up (count), down (count), total (count)
filter_degs <- function(dif, fc_name, p_name, logfc_cutoff, p_value, cutoff, rgs = NULL) {
  dif_up <- dif %>%
    filter(!!sym(fc_name) > logfc_cutoff & !!sym(p_name) < p_value)
  dif_down <- dif %>%
    filter(!!sym(fc_name) < -logfc_cutoff & !!sym(p_name) < p_value)

  gene_col <- colnames(dif)[1]
  dif_grouped <- dif %>% mutate(group = case_when(
    dif[[gene_col]] %in% dif_up[[gene_col]] ~ "Up",
    dif[[gene_col]] %in% dif_down[[gene_col]] ~ "Down",
    TRUE ~ "Not"
  ))

  degs <- dif_grouped[dif_grouped[["group"]] == "Up" | dif_grouped[["group"]] == "Down", ]
  n_degs <- nrow(degs)

  if (n_degs == 0) {
    report_exception_ndjson(
      "B6_NO_DEGS", "data_insufficient", "skip_with_warning",
      sprintf("No differentially expressed genes found at |logFC|>%.2f, %s<%.4f",
              logfc_cutoff, p_name, p_value),
      exit_code = 1
    )
    return(list(degs = degs, rdegs = data.frame(), dif_grouped = dif_grouped,
                up = 0, down = 0, total = 0))
  }

  if (n_degs < cutoff) {
    report_exception_ndjson(
      "B2_FEW_DEGS", "data_insufficient", "skip_with_warning",
      sprintf("DEG count %d below cutoff %d at |logFC|>%.2f, %s<%.4f",
              n_degs, cutoff, logfc_cutoff, p_name, p_value),
      exit_code = 1
    )
  }

  n_up <- nrow(dif_up)
  n_down <- nrow(dif_down)

  if (!is.null(rgs) && length(rgs) > 0) {
    genes <- intersect(degs[[gene_col]], rgs)
    rdegs <- degs[degs[[gene_col]] %in% genes, ]
    if (nrow(rdegs) < cutoff) {
      report_exception_ndjson(
        "B7_RGS_INTERSECTION", "data_insufficient", "skip_with_warning",
        sprintf("Intersection of DEGs with related gene set (%d genes) below cutoff %d",
                nrow(rdegs), cutoff),
        exit_code = 1
      )
    }
  } else {
    rdegs <- degs
  }

  list(degs = degs, rdegs = rdegs, dif_grouped = dif_grouped,
       up = n_up, down = n_down, total = n_degs)
}
```

- [x] **Step 8: Verify filter_helpers.R with synthetic data**

```bash
./env/bin/Rscript -e '
source("scripts/report.R")
source("scripts/exceptions.R")
source("tests/testthat/setup.R")
td <- create_test_data()
source("scripts/filter_helpers.R")

# proportion_check (balanced design)
map <- td$map
map[[2]] <- factor(map[[2]], levels = c("Control", "Case"))
ok <- proportion_check(map)
stopifnot(ok)

# proportion_check (imbalanced)
map2 <- map[c(1,2,3,4,5,5),]
rownames(map2) <- NULL
map2[6,1] <- "CTRL3dup"
map2[[2]] <- factor(map2[[2]], levels = c("Control", "Case"))
proportion_check(map2)  # should warn
cat("filter_helpers.R OK\n")
'
```

Expected: stdout contains the proportion exception NDJSON line and `filter_helpers.R OK`.

- [x] **Step 9: Commit filter_helpers.R**

```bash
git add scripts/filter_helpers.R
git commit -m "feat: add DEG filtering helpers (filter_helpers.R)

proportion_check() — enforces 10:1 ratio limit, emits B1_PROPORTION
test_cutoff() — sensitivity analysis at multiple logFC thresholds
filter_degs() — classifies Up/Down/Not, emits B6_NO_DEGS, B2_FEW_DEGS, B7_RGS_INTERSECTION

Adapted from original/scripts/diff_filter.R and diff_check_proportion.py.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### 3.4 Create scripts/plot_helpers.R

- [x] **Step 10: Write scripts/plot_helpers.R**

Adapt all four plotting functions from the original scripts. Each function writes a PDF file to disk and reports progress via report_info():

```r
# plot_helpers.R — Visualization functions for differential-analysis
#
# Provides volcano plot, heatmap, Venn diagram, and chromosome location plot.
# Adapted from original/scripts/diff_volcano.R, diff_heatmap.R, diff_venn.R, diff_locate.R.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

#' Generate volcano plot
#'
#' @param dif Data frame with DE results (columns: gene_id, logFC, p-value column)
#' @param p_name Character, name of p-value column ("Pvalue" or "Padj")
#' @param p_value Numeric, p-value threshold for horizontal line
#' @param logfc_cutoff Numeric, logFC threshold for vertical lines
#' @param top Integer, number of top genes to label (used if gene is NULL/None)
#' @param gene Character, comma-separated gene IDs to label, or NULL/None to auto-label top N
#' @param outfile Character, output PDF file path
#' @return NULL (side effect: writes PDF)
plot_volcano <- function(dif, p_name, p_value, logfc_cutoff, top, gene, outfile) {
  fc_name <- "logFC"
  gene_col <- colnames(dif)[1]

  # Build group classification
  dif_up <- dif %>%
    filter(!!sym(fc_name) > logfc_cutoff & !!sym(p_name) < p_value)
  dif_down <- dif %>%
    filter(!!sym(fc_name) < -logfc_cutoff & !!sym(p_name) < p_value)
  dif_plot <- dif %>% mutate(group = case_when(
    dif[[gene_col]] %in% dif_up[[gene_col]] ~ "Up",
    dif[[gene_col]] %in% dif_down[[gene_col]] ~ "Down",
    TRUE ~ "Not"
  ))

  dif_plot$logP <- -log10(dif_plot[[p_name]])
  dif_plot$change <- factor(dif_plot$group, levels = c("Down", "Not", "Up"))
  dif_plot <- dif_plot[order(abs(dif_plot[[fc_name]]), decreasing = TRUE), ]

  # Determine genes to label
  if (is.null(gene) || gene == "NULL" || gene == "None" || gene == "") {
    if (is.null(top) || top == "NULL" || top == "None") {
      label_genes <- NULL
    } else {
      label_genes <- head(dif_plot[[gene_col]], as.integer(top))
    }
  } else {
    label_genes <- strsplit(as.character(gene), ",")[[1]]
  }

  lab_y <- paste0("-Log10(", p_name, ")")
  x_lim <- max(abs(dif_plot[[fc_name]])) * 1.1

  p <- ggplot(dif_plot, aes(.data[[fc_name]], .data$logP, color = change)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(x = "LogFC", y = lab_y, color = "Significance") +
    geom_hline(yintercept = -log10(p_value), lty = 2) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), lty = 2) +
    scale_x_continuous(limits = c(-x_lim, x_lim)) +
    scale_color_manual(values = c(
      "Down" = "#4DBBD5",
      "Not" = "grey",
      "Up" = "#E64B35"
    )) +
    ggtitle("Volcano") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text = element_text(size = 12)
    ) +
    theme(panel.grid = element_blank())

  if (!is.null(label_genes) && length(label_genes) > 0) {
    dif_plot$label <- mapply(function(x) {
      ifelse(x %in% label_genes, x, "")
    }, dif_plot[[gene_col]])

    p <- p + ggrepel::geom_text_repel(
      data = dif_plot, aes(label = label),
      segment.alpha = 0.4,
      box.padding = 0.5,
      force = 1,
      max.overlaps = Inf,
      min.segment.length = 0.25,
      show.legend = FALSE
    )
  }

  ggsave(file = outfile, p, width = 8, height = 5)
  report_info(sprintf("Volcano plot saved to %s", outfile))
}

#' Generate heatmap of top DEGs
#'
#' @param mat Expression matrix (genes x samples), full dataset
#' @param map Data frame with sample-to-group mapping
#' @param rdegs Data frame of filtered DEGs with gene_id and logFC columns
#' @param top Integer, number of top up/down genes for heatmap display
#' @param color_heat Character, comma-separated heatmap colors (e.g., "blue,white,red")
#' @param outfile Character, output PDF file path
#' @return NULL (side effect: writes PDF)
plot_heatmap <- function(mat, map, rdegs, top, color_heat, outfile) {
  suppressPackageStartupMessages(library(pheatmap))

  map[[2]] <- factor(map[[2]], levels = unique(map[[2]]))

  # Select top N up and top N down genes for heatmap
  if (nrow(rdegs) > top * 2) {
    rdegs_sorted <- rdegs[order(rdegs[["logFC"]]), ]
    top_mat <- rdegs_sorted[c(1:top, (nrow(rdegs_sorted) - top + 1):nrow(rdegs_sorted)), ]
    deg_heatmap <- top_mat[[1]]  # gene_id column
  } else {
    deg_heatmap <- rdegs[[1]]
  }

  dat_heatmap <- mat[deg_heatmap, map[[1]], drop = FALSE]

  color_vec <- strsplit(as.character(color_heat), ",")[[1]]

  pdf(NULL)  # suppress pheatmap implicit pdf call
  p_heatmap <- pheatmap(dat_heatmap,
    scale = "row",
    annotation_col = data.frame(
      Group = map[[2]],
      row.names = map[[1]]
    ),
    annotation_colors = list(Group = setNames(
      c("#E64B35", "#4DBBD5")[seq_len(length(unique(map[[2]])))],
      unique(map[[2]])
    )),
    color = colorRampPalette(color_vec)(50),
    breaks = c(seq(-3, 3, length = 50)),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    labels_col = "",
    border_color = NA,
    main = "Heatmap"
  )
  dev.off()

  ggsave(file = outfile, p_heatmap, width = 8, height = 8)
  report_info(sprintf("Heatmap saved to %s", outfile))
}

#' Generate Venn diagram of DEGs vs related gene set
#'
#' @param dif Data frame of DEGs with group column (Up/Down) and gene_id column
#' @param rgs Character vector of related gene set identifiers
#' @param pheno_abbr Character, phenotype abbreviation for Venn label
#' @param color_panel Character, comma-separated colors for Venn sets
#' @param outfile Character, output PDF file path
#' @return NULL (side effect: writes PDF)
plot_venn <- function(dif, rgs, pheno_abbr, color_panel, outfile) {
  if (!requireNamespace("ggvenn", quietly = TRUE)) {
    report_info("ggvenn not installed, skipping Venn diagram")
    return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(ggvenn))

  gene_col <- colnames(dif)[1]
  deg_genes <- dif[dif[["group"]] == "Up" | dif[["group"]] == "Down", ][[gene_col]]
  colors <- strsplit(as.character(color_panel), ",")[[1]]

  p_venn <- ggvenn(
    setNames(
      list(deg_genes, rgs),
      c("DEGs", pheno_abbr)
    ),
    c("DEGs", pheno_abbr),
    show_percentage = FALSE,
    fill_alpha = 0.5,
    stroke_color = NA,
    fill_color = color_map(colors, c("DEGs", pheno_abbr))
  )

  ggsave(file = outfile, p_venn, width = 8, height = 5)
  report_info(sprintf("Venn diagram saved to %s", outfile))
}

#' Generate chromosome location plot via RCircos
#'
#' @param gene_list Character vector of gene identifiers to plot
#' @param locate Character, path to chromosome annotation CSV (columns: Gene, Chromosome, Start, End)
#' @param tax_id Character, NCBI taxonomy ID ("9606" for human, "10090" for mouse)
#' @param outfile Character, output PDF file path
#' @return NULL (side effect: writes PDF)
plot_locate <- function(gene_list, locate, tax_id, outfile) {
  if (!requireNamespace("RCircos", quietly = TRUE)) {
    report_info("RCircos not installed, skipping chromosome location plot")
    return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(RCircos))

  chr <- read.csv(locate)
  gene <- gene_list

  chr_gene <- chr[which(chr$Gene %in% gene), ]

  if (nrow(chr_gene) == 0) {
    report_info("No genes matched in chromosome annotation, skipping location plot")
    return(invisible(NULL))
  }

  pdf(file = outfile, width = 8, height = 8)

  if (tax_id == "9606") {
    data(UCSC.HG38.Human.CytoBandIdeogram)
    cyto_info <- UCSC.HG38.Human.CytoBandIdeogram
  } else if (tax_id == "10090") {
    data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
    cyto_info <- UCSC.Mouse.GRCm38.CytoBandIdeogram
  } else {
    dev.off()
    report_exception_ndjson(
      "E802_UNSUPPORTED_TAXID", "data_mismatch", "halt",
      sprintf("Unsupported tax_id: %s. Supported: 9606 (human), 10090 (mouse)", tax_id),
      exit_code = 1
    )
    return(invisible(NULL))
  }

  RCircos.Set.Core.Components(cyto_info)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  RCircos.Gene.Connector.Plot(chr_gene, track.num = 1, side = "in")
  RCircos.Gene.Name.Plot(chr_gene, name.col = 4, track.num = 2, side = "in")
  dev.off()

  report_info(sprintf("Chromosome location plot saved to %s", outfile))
}
```

- [x] **Step 11: Verify plot_helpers.R loads (plot functions available)**

```bash
./env/bin/Rscript -e '
source("scripts/report.R")
source("scripts/exceptions.R")
source("scripts/io_helpers.R")
source("scripts/plot_helpers.R")
stopifnot(exists("plot_volcano"))
stopifnot(exists("plot_heatmap"))
stopifnot(exists("plot_venn"))
stopifnot(exists("plot_locate"))
cat("plot_helpers.R OK\n")
'
```

Expected stdout: `plot_helpers.R OK`

- [x] **Step 12: Commit plot_helpers.R**

```bash
git add scripts/plot_helpers.R
git commit -m "feat: add visualization helpers (plot_helpers.R)

plot_volcano() — ggplot2 volcano with ggrepel labels
plot_heatmap() — pheatmap of top DEGs, row-scaled
plot_venn() — ggvenn intersection of DEGs and related genes
plot_locate() — RCircos chromosome location with cyto band ideogram

All functions emit report_info() on completion.
Conditional plots skip gracefully with info message if package missing.

Adapted from original/scripts/diff_volcano.R, diff_heatmap.R,
diff_venn.R, diff_locate.R.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

archived-with: 2026-06-22-create-node-package
---

## Task 4: Validation Scripts

**Files:**
- Create: `scripts/input_validation.R`
- Create: `scripts/output_validation.R`

**Interfaces:**
- Consumes: filter_helpers.R (for proportion_check), report.R, exceptions.R
- Produces:
  - `validate_input(opts)` — returns list(valid, reason, ...), used by both main.R subcommand and standalone
  - `validate_output(opts)` — returns list(valid, reason, ...), used by both main.R subcommand and standalone

### 4.1 Create scripts/input_validation.R

- [x] **Step 1: Write scripts/input_validation.R**

Pattern follows the reference node's `input_validation.R` but adapted for DEG-specific inputs:
- File existence and readability (mat CSV, map CSV)
- Required column presence in mat (gene identifier column, all sample columns from map)
- Sample columns in mat must match sample IDs in map
- Group structure validation (exactly 2 groups in map col 2)
- Sample proportion ratio <= 10:1 (calls proportion_check from filter_helpers.R)
- Supports standalone execution: `Rscript input_validation.R --mat <path> --map <path> [--force-imbalanced]`

```r
#!/usr/bin/env Rscript
#
# input_validation.R — Pre-flight input validation for differential-analysis
#
# Standalone executable. Validates inputs before running the DE pipeline.
# Usage:
#   Rscript input_validation.R --mat expr.csv --map groups.csv
#   Rscript input_validation.R --mat expr.csv --map groups.csv --force-imbalanced
#
# Exit codes:
#   0 — all checks passed
#   1 — validation failed (stderr has reason)
#   2 — usage/argument error

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))

source(file.path(script_dir, "report.R"))
source(file.path(script_dir, "exceptions.R"))
source(file.path(script_dir, "filter_helpers.R"))

# -------------------------------------------------------------------
# Validation function
# -------------------------------------------------------------------

#' Validate inputs for DEG analysis
#'
#' @param opts List with mat (file path), map (file path), force_imbalanced (logical)
#' @return List with valid (logical), reason (character), and optionally n_genes, n_samples
validate_input <- function(opts) {
  # File existence checks
  if (is.null(opts$mat) || !file.exists(opts$mat)) {
    return(list(valid = FALSE,
      reason = sprintf("Expression matrix file not found: %s",
                       if (is.null(opts$mat)) "(null)" else opts$mat)))
  }
  if (is.null(opts$map) || !file.exists(opts$map)) {
    return(list(valid = FALSE,
      reason = sprintf("Sample map file not found: %s",
                       if (is.null(opts$map)) "(null)" else opts$map)))
  }

  # Read expression matrix
  mat <- tryCatch(
    data.table::fread(opts$mat, data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    return(list(valid = FALSE,
      reason = sprintf("Cannot read expression matrix as CSV: %s", opts$mat)))
  }

  # Check for gene identifier column (first column)
  if (ncol(mat) < 2) {
    return(list(valid = FALSE,
      reason = sprintf("Expression matrix has only %d column(s). Need gene ID column + samples.",
                       ncol(mat))))
  }

  # Check for empty matrix
  if (nrow(mat) == 0) {
    return(list(valid = FALSE,
      reason = "Expression matrix has 0 rows (empty matrix)"))
  }

  # Read sample group map
  map <- tryCatch(
    data.table::fread(opts$map, data.table = FALSE, header = TRUE),
    error = function(e) NULL
  )
  if (is.null(map)) {
    return(list(valid = FALSE,
      reason = sprintf("Cannot read sample map as CSV: %s", opts$map)))
  }

  if (ncol(map) < 2) {
    return(list(valid = FALSE,
      reason = sprintf("Sample map has only %d column(s). Need sample ID + group columns.",
                       ncol(map))))
  }

  # Check sample columns in mat match sample IDs in map
  mat_samples <- colnames(mat)[-1]  # first col is gene ID
  map_samples <- map[[1]]
  missing_in_mat <- setdiff(map_samples, mat_samples)
  missing_in_map <- setdiff(mat_samples, map_samples)

  if (length(missing_in_mat) > 0) {
    return(list(valid = FALSE,
      reason = sprintf("Sample(s) in map but not in matrix: %s",
                       paste(missing_in_mat, collapse = ", "))))
  }
  if (length(missing_in_map) > 0) {
    return(list(valid = FALSE,
      reason = sprintf("Sample(s) in matrix but not in map: %s",
                       paste(missing_in_map, collapse = ", "))))
  }

  # Group structure: exactly 2 groups
  groups <- unique(map[[2]])
  if (length(groups) != 2) {
    return(list(valid = FALSE,
      reason = sprintf("Expected exactly 2 groups, found %d: %s",
                       length(groups), paste(groups, collapse = ", "))))
  }

  # Make factor for proportion_check
  map[[2]] <- factor(map[[2]], levels = groups)

  # Proportion check
  force_imb <- isTRUE(opts$force_imbalanced) ||
               identical(opts$force_imbalanced, "TRUE") ||
               identical(opts$force_imbalanced, "true")
  if (!proportion_check(map, force_imb)) {
    return(list(valid = FALSE,
      reason = "Sample proportion check failed (ratio > 10:1, use --force-imbalanced to override)"))
  }

  list(valid = TRUE, reason = "OK",
       n_genes = nrow(mat), n_samples = length(map_samples),
       groups = paste(groups, collapse = " vs "))
}

# -------------------------------------------------------------------
# CLI argument parsing
# -------------------------------------------------------------------

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) {
    cat("Usage: Rscript input_validation.R --mat <expr.csv> --map <groups.csv> [--force-imbalanced]\n")
    cat("\nOptions:\n")
    cat("  --mat FILE           Expression matrix CSV (genes x samples)\n")
    cat("  --map FILE           Sample-to-group mapping CSV\n")
    cat("  --force-imbalanced   Override 10:1 proportion check\n")
    quit(status = 2)
  }

  opts <- list(mat = NULL, map = NULL, force_imbalanced = FALSE)
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--mat") {
      i <- i + 1; if (i <= length(args)) opts$mat <- args[i]
    } else if (key == "--map") {
      i <- i + 1; if (i <= length(args)) opts$map <- args[i]
    } else if (key == "--force-imbalanced") {
      opts$force_imbalanced <- TRUE
    } else if (startsWith(key, "--mat=")) {
      opts$mat <- sub("^--mat=", "", key)
    } else if (startsWith(key, "--map=")) {
      opts$map <- sub("^--map=", "", key)
    } else if (startsWith(key, "--force-imbalanced=")) {
      opts$force_imbalanced <- TRUE
    } else {
      cat(sprintf("Error: unknown option: %s\n", key), file = stderr())
      quit(status = 2)
    }
    i <- i + 1
  }

  if (is.null(opts$mat)) {
    cat("Error: --mat is required\n", file = stderr())
    quit(status = 2)
  }
  if (is.null(opts$map)) {
    cat("Error: --map is required\n", file = stderr())
    quit(status = 2)
  }
  opts
}

# -------------------------------------------------------------------
# Main dispatch (standalone mode)
# -------------------------------------------------------------------

if (sys.nframe() == 0) {
  opts <- parse_args()
  result <- validate_input(opts)
  if (!result$valid) {
    cat(sprintf("Validation failed: %s\n", result$reason), file = stderr())
    quit(status = 1)
  }
  cat(sprintf("OK: %d genes x %d samples, groups: %s\n",
              result$n_genes, result$n_samples, result$groups))
  quit(status = 0)
}
```

- [x] **Step 2: Test input_validation.R with synthetic data**

```bash
./env/bin/Rscript scripts/input_validation.R --mat /tmp/test_expr.csv --map /tmp/test_map.csv
```

Expected stdout: `OK: 10 genes x 6 samples, groups: Control vs Case`

- [x] **Step 3: Commit input_validation.R**

```bash
git add scripts/input_validation.R
git commit -m "feat: add input validation script (input_validation.R)

Checks: file existence, CSV readability, required columns,
sample ID matching, exactly 2 groups, proportion <= 10:1.
Standalone executable and sourceable by main.R.
Emits B3_MISSING_INPUT, B4_INVALID_COLUMNS, B8_EMPTY_MATRIX,
B9_SAMPLE_MISMATCH, B1_PROPORTION via exception helpers.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### 4.2 Create scripts/output_validation.R

- [x] **Step 4: Write scripts/output_validation.R**

```r
#!/usr/bin/env Rscript
#
# output_validation.R — Post-hoc output validation for differential-analysis
#
# Standalone executable. Validates outputs after running the DE pipeline.
# Usage:
#   Rscript output_validation.R --outdir ./output
#
# Exit codes:
#   0 — all checks passed
#   1 — validation failed (stderr has reason)
#   2 — usage/argument error

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))

source(file.path(script_dir, "report.R"))
source(file.path(script_dir, "exceptions.R"))

# -------------------------------------------------------------------
# Validation function
# -------------------------------------------------------------------

#' Validate outputs from DEG analysis
#'
#' Checks output CSV column presence (gene_id, logFC, Pvalue, Padj),
#' non-empty results, and plot file existence with non-zero size.
#'
#' @param opts List with outdir (output directory path)
#' @return List with valid (logical), reason (character), file_info (named list)
validate_output <- function(opts) {
  outdir <- opts$outdir
  if (is.null(outdir) || !dir.exists(outdir)) {
    return(list(valid = FALSE,
      reason = sprintf("Output directory not found or not specified: %s",
                       if (is.null(outdir)) "(null)" else outdir)))
  }

  # Check Diffanalysis.csv
  diff_csv <- file.path(outdir, "Diffanalysis.csv")
  if (!file.exists(diff_csv)) {
    return(list(valid = FALSE,
      reason = sprintf("Missing Diffanalysis.csv in %s", outdir)))
  }

  dif <- tryCatch(
    read.csv(diff_csv, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(dif)) {
    return(list(valid = FALSE,
      reason = sprintf("Cannot read Diffanalysis.csv: %s", diff_csv)))
  }

  required_cols <- c("gene_id", "logFC", "Pvalue", "Padj")
  missing_cols <- setdiff(required_cols, colnames(dif))
  if (length(missing_cols) > 0) {
    return(list(valid = FALSE,
      reason = sprintf("Diffanalysis.csv missing required columns: %s",
                       paste(missing_cols, collapse = ", "))))
  }

  if (nrow(dif) == 0) {
    return(list(valid = FALSE,
      reason = "Diffanalysis.csv has 0 rows (no genes in result)"))
  }

  diff_rows <- nrow(dif)
  diff_cols <- ncol(dif)

  # Check DEGs.csv
  degs_csv <- file.path(outdir, "DEGs.csv")
  if (!file.exists(degs_csv)) {
    return(list(valid = FALSE,
      reason = sprintf("Missing DEGs.csv in %s", outdir)))
  }

  degs <- tryCatch(
    read.csv(degs_csv, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(degs)) {
    return(list(valid = FALSE,
      reason = sprintf("Cannot read DEGs.csv: %s", degs_csv)))
  }

  if (!"group" %in% colnames(degs)) {
    return(list(valid = FALSE,
      reason = "DEGs.csv missing 'group' column"))
  }

  degs_rows <- nrow(degs)

  # Check plot files
  required_plots <- c("Volcano.pdf", "Heatmap.pdf")
  for (pf in required_plots) {
    ppath <- file.path(outdir, pf)
    if (!file.exists(ppath)) {
      return(list(valid = FALSE,
        reason = sprintf("Missing plot file: %s", pf)))
    }
    if (file.info(ppath)$size == 0) {
      return(list(valid = FALSE,
        reason = sprintf("Plot file is empty: %s", pf)))
    }
  }

  file_info <- list(
    diff_csv_rows = diff_rows,
    diff_csv_cols = diff_cols,
    degs_csv_rows = degs_rows
  )

  list(valid = TRUE, reason = "OK", file_info = file_info)
}

# -------------------------------------------------------------------
# CLI argument parsing
# -------------------------------------------------------------------

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) {
    cat("Usage: Rscript output_validation.R --outdir <path>\n")
    cat("\nOptions:\n")
    cat("  --outdir DIR   Output directory containing DEG results\n")
    quit(status = 2)
  }

  opts <- list(outdir = NULL)
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--outdir") {
      i <- i + 1; if (i <= length(args)) opts$outdir <- args[i]
    } else if (startsWith(key, "--outdir=")) {
      opts$outdir <- sub("^--outdir=", "", key)
    }
    i <- i + 1
  }

  if (is.null(opts$outdir)) {
    cat("Error: --outdir is required\n", file = stderr())
    quit(status = 2)
  }
  opts
}

# -------------------------------------------------------------------
# Main dispatch (standalone mode)
# -------------------------------------------------------------------

if (sys.nframe() == 0) {
  opts <- parse_args()
  result <- validate_output(opts)
  if (!result$valid) {
    cat(sprintf("Validation failed: %s\n", result$reason), file = stderr())
    quit(status = 1)
  }
  info <- paste(vapply(names(result$file_info), function(f) {
    sprintf("%s=%s", f, result$file_info[[f]])
  }, ""), collapse = ", ")
  cat(sprintf("OK: %s\n", info))
  quit(status = 0)
}
```

- [x] **Step 5: Verify output_validation.R loads and function is accessible**

```bash
./env/bin/Rscript -e '
source("scripts/report.R"); source("scripts/exceptions.R")
source("scripts/output_validation.R")
stopifnot(exists("validate_output"))
cat("output_validation.R OK\n")
'
```

Expected stdout: `output_validation.R OK`

- [x] **Step 6: Commit output_validation.R**

```bash
git add scripts/output_validation.R
git commit -m "feat: add output validation script (output_validation.R)

Checks: Diffanalysis.csv existence and required columns
(gene_id, logFC, Pvalue, Padj), DEGs.csv existence and group
column, Volcano.pdf and Heatmap.pdf existence and non-zero size.
Returns detailed file_info on success.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

archived-with: 2026-06-22-create-node-package
---

## Task 5: Main Entry Point

**Files:**
- Create: `scripts/main.R`

**Interfaces:**
- Consumes: All 7 helper modules (report.R, exceptions.R, io_helpers.R, diff_methods.R, filter_helpers.R, plot_helpers.R, input_validation.R, output_validation.R)
- Produces: CLI entry point with 3 subcommands
- Functions:
  - `parse_args(args)` — returns named list of 20 parameters
  - `check_environment()` — from exceptions.R, re-exported
  - `do_run(opts)` — full DEG pipeline
  - `do_validate_input(opts)` — pre-flight validation
  - `do_validate_output(opts)` — post-hoc validation
  - `main()` — top-level dispatch

### 5.1 Create scripts/main.R

- [x] **Step 1: Write scripts/main.R**

Adapt from reference node `main.R`. Source all modules in dependency order. Implement subcommand dispatch and the full `do_run()` pipeline:

```r
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

# Resolve script directory for relative sourcing (works from any CWD)
script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))

# Source modules in dependency order
source(file.path(script_dir, "report.R"))
source(file.path(script_dir, "exceptions.R"))
source(file.path(script_dir, "io_helpers.R"))
source(file.path(script_dir, "diff_methods.R"))
source(file.path(script_dir, "filter_helpers.R"))
source(file.path(script_dir, "plot_helpers.R"))
source(file.path(script_dir, "input_validation.R"))
source(file.path(script_dir, "output_validation.R"))

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
    cat("  --method METHOD      DE method: deseq2, limma, edgeR, t, wilcox (default: deseq2)\n")
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
        found <- TRUE
        break
      }

      # --key value form
      if (key == flag) {
        i <- i + 1
        if (i <= length(remaining)) {
          opts[[opt_name]] <- type_convert(remaining[i], opt_name)
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
    # proportion_check already emitted the exception
    return(invisible(list(status = "error",
      msg = "Sample proportion check failed")))
  }

  # Subset matrix to map samples
  mat <- mat[, map[[1]], drop = FALSE]

  # Method consistency check: count-based methods need integer data
  if (opts$method %in% c("deseq2", "edgeR")) {
    if (!all(mat == round(mat), na.rm = TRUE)) {
      report_info("Converting float matrix to integer counts for count-based method")
      mat <- round(mat)
    }
  }

  # Run DE analysis
  report_info(sprintf("Running %s differential expression...", opts$method))
  dif <- switch(opts$method,
    deseq2 = diff_deseq2(mat, map),
    limma  = diff_limma(mat, map),
    edgeR  = diff_edger(mat, map, norm = opts$norm, model = opts$model),
    t      = diff_stat(mat, map, "t"),
    wilcox = diff_stat(mat, map, "wilcox")
  )

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
  write.csv(dif, diff_path, row.names = FALSE)
  report_info(sprintf("Full DE results written to %s (%d rows)",
    basename(diff_path), nrow(dif)))

  # Write DEGs.csv (filtered with group column)
  degs_path <- file.path(opts$outdir, "DEGs.csv")
  write.csv(filt$degs, degs_path, row.names = FALSE)
  report_info(sprintf("Filtered DEGs written to %s (%d rows)",
    basename(degs_path), nrow(filt$degs)))

  # Volcano plot
  volcano_path <- file.path(opts$outdir, "Volcano.pdf")
  plot_volcano(dif, p_name, opts$pvalue, opts$logfc_cutoff,
    opts$top, opts$gene, volcano_path)

  # Heatmap
  heatmap_path <- file.path(opts$outdir, "Heatmap.pdf")
  plot_heatmap(mat, map, filt$rdegs, opts$top, opts$color_heat, heatmap_path)

  # Optional: Venn diagram
  if (!is.null(opts$rgs) && opts$rgs != "None" && file.exists(opts$rgs)) {
    if (!is.null(opts$pheno_abbr) && !is.null(opts$color_panel)) {
      venn_path <- file.path(opts$outdir, "Venn.pdf")
      plot_venn(filt$degs, rgs_genes, opts$pheno_abbr,
        opts$color_panel, venn_path)
    } else {
      report_info("Skipping Venn diagram: --pheno-abbr and --color-panel required")
    }
  }

  # Optional: Chromosome location plot
  if (!is.null(opts$locate) && opts$locate != "None" && file.exists(opts$locate)) {
    locate_path <- file.path(opts$outdir, "Chromosome_location.pdf")
    plot_locate(filt$degs[[1]], opts$locate, opts$tax_id, locate_path)
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
      if (grepl("not found", result$reason)) "B3_MISSING_INPUT"
      else if (grepl("column|empty", result$reason)) "B4_INVALID_COLUMNS"
      else if (grepl("proportion", result$reason)) "B1_PROPORTION"
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
```

- [x] **Step 2: Verify main.R loads without errors**

```bash
./env/bin/Rscript scripts/main.R 2>&1 | head -5
```

Expected: usage/help text printed to stdout.

- [x] **Step 3: Test subcommand dispatch (validate-input)**

```bash
./env/bin/Rscript -e '
source("tests/testthat/setup.R")
td <- create_test_data()
' && ./env/bin/Rscript scripts/main.R validate-input --mat /tmp/test_expr.csv --map /tmp/test_map.csv
```

Expected stdout: `{"level":"info"...}` lines and validation passing.

- [x] **Step 4: Commit main.R**

```bash
git add scripts/main.R
git commit -m "feat: add main entry point with subcommand dispatch (main.R)

Subcommands: run, validate-input, validate-output.
Full pipeline: load -> proportion_check -> DE analysis ->
filter_degs -> volcano -> heatmap -> [venn] -> [locate] ->
write_csvs -> .run_result.json.

20 parameters with --key=value and --key value parsing.
NDJSON reporting throughout.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

archived-with: 2026-06-22-create-node-package
---

## Task 6: SKILL.md Contract

**Files:**
- Create: `SKILL.md`

**Interfaces:**
- Produces: Machine-readable YAML frontmatter with all 8 required fields from protocol/skill-md-frontmatter.md

### 6.1 Create SKILL.md

- [x] **Step 1: Write SKILL.md**

```markdown
archived-with: 2026-06-22-create-node-package
---
name: differential-analysis
description: >
  Differential expression analysis comparing case vs control groups
  from a gene expression matrix. Supports DESeq2, limma, edgeR,
  t-test, and Wilcoxon methods. Produces volcano plot, heatmap,
  and optionally Venn diagram and chromosome location plot.
  Requires a sample-to-group map with exactly 2 groups.
type: standard
inputs:
  - name: expression_matrix.csv
    format: csv
    semantic_type: expression_matrix
    description: Gene expression matrix (genes x samples) with gene IDs in first column
  - name: sample_group_map.csv
    format: csv
    semantic_type: sample_metadata
    description: Two-column CSV mapping sample IDs (col 1) to group labels (col 2). Exactly 2 groups required.
outputs:
  - name: Diffanalysis.csv
    format: csv
    semantic_type: differential_expression_results
    columns: [gene_id, logFC, Pvalue, Padj, stat]
  - name: DEGs.csv
    format: csv
    semantic_type: filtered_degs
    columns: [gene_id, logFC, Pvalue, Padj, group]
  - name: Volcano.pdf
    format: pdf
    semantic_type: volcano_plot
  - name: Heatmap.pdf
    format: pdf
    semantic_type: heatmap_plot
  - name: Venn.pdf
    format: pdf
    semantic_type: venn_diagram
    conditional: true
  - name: Chromosome_location.pdf
    format: pdf
    semantic_type: chromosome_location_plot
    conditional: true
entry: scripts/main.R
parameters:
  - name: subcommand
    type: choice
    required: true
    bind: config
    description: Operation to perform (run, validate-input, validate-output)
  - name: --mat
    type: file
    required: true
    bind: upstream
    description: Expression matrix CSV (genes x samples) with gene identifiers in first column
  - name: --map
    type: file
    required: true
    bind: upstream
    description: Sample-to-group mapping CSV (sample_id, group)
  - name: --method
    type: choice
    required: false
    default: deseq2
    bind: config
    description: DE analysis method (deseq2, limma, edgeR, t, wilcox)
  - name: --p-set
    type: choice
    required: false
    default: padj
    bind: static
    description: Which p-value to use for filtering (p or padj)
  - name: --pvalue
    type: float
    required: false
    default: 0.05
    range: [0.001, 0.25]
    bind: static
    description: P-value threshold for significance
  - name: --logfc-cutoff
    type: float
    required: false
    default: 1.0
    range: [0.0, 10.0]
    bind: static
    description: Absolute log2 fold-change cutoff
  - name: --cutoff
    type: int
    required: false
    default: 10
    bind: static
    description: Minimum DEG count after filtering (below this triggers B2_FEW_DEGS)
  - name: --norm
    type: choice
    required: false
    default: TMM
    bind: static
    description: edgeR normalization method (TMM, RLE, upperquartile, none)
  - name: --model
    type: choice
    required: false
    default: glmFit
    bind: static
    description: edgeR model fitting method (glmFit, glmQLFit)
  - name: --top
    type: int
    required: false
    default: 20
    bind: static
    description: Top N genes for heatmap display and volcano labels
  - name: --force-imbalanced
    type: bool
    required: false
    default: false
    bind: static
    description: Override 10:1 sample proportion check
  - name: --rgs
    type: file
    required: false
    bind: upstream
    description: Related gene set CSV for Venn diagram intersection
  - name: --locate
    type: file
    required: false
    bind: upstream
    description: Chromosome annotation CSV (Gene, Chromosome, Start, End) for location plot
  - name: --tax-id
    type: string
    required: false
    default: "9606"
    bind: config
    description: NCBI taxonomy ID (9606=human, 10090=mouse)
  - name: --pheno-abbr
    type: string
    required: false
    bind: config
    description: Phenotype abbreviation for Venn diagram label
  - name: --gene
    type: string
    required: false
    bind: static
    description: Comma-separated gene IDs to label on volcano plot
  - name: --color-heat
    type: string
    required: false
    default: blue,white,red
    bind: static
    description: Heatmap color palette (comma-separated)
  - name: --color-panel
    type: string
    required: false
    bind: static
    description: Comma-separated colors for Venn diagram sets
  - name: --outdir
    type: file_out
    required: false
    default: .
    bind: framework
    description: Output directory for all result files
exceptions:
  - exit_code: 1
    pattern: "sample proportion.*exceeds"
    nature: data_insufficient
    action: skip_with_warning
  - exit_code: 1
    pattern: "gene count.*below.*cutoff"
    nature: data_insufficient
    action: skip_with_warning
  - exit_code: 1
    pattern: "input file not found"
    nature: data_corrupt
    action: halt
  - exit_code: 1
    pattern: "missing required column"
    nature: data_corrupt
    action: halt
  - exit_code: 1
    pattern: "not support.*matrix"
    nature: data_mismatch
    action: halt
  - exit_code: 1
    pattern: "no differentially expressed"
    nature: data_insufficient
    action: skip_with_warning
  - exit_code: 1
    pattern: "intersection.*below.*cutoff"
    nature: data_insufficient
    action: skip_with_warning
  - exit_code: 1
    pattern: "empty|no rows|no columns"
    nature: data_insufficient
    action: halt
  - exit_code: 1
    pattern: "sample.*mismatch|columns.*not.*match"
    nature: data_corrupt
    action: halt
  - exit_code: 1
    pattern: "disk full|no space"
    nature: resource
    action: halt
  - exit_code: 1
    pattern: "permission denied"
    nature: resource
    action: halt
  - exit_code: 3
    pattern: "Missing required packages"
    nature: env_bug
    action: halt
  - exit_code: 1
    pattern: "Unsupported tax_id"
    nature: data_mismatch
    action: halt
hardware:
  memory_gb: 4
  cpu: 2
  gpu: false
  runtime: "~2-10 minutes depending on method and matrix size"
archived-with: 2026-06-22-create-node-package
---

# differential-analysis

## Node Function

Performs differential expression analysis comparing exactly two groups (case vs control) from a gene expression matrix. The pipeline:

1. Validates input data and sample proportion (max 10:1 ratio, overridable)
2. Runs the selected DE method (DESeq2, limma, edgeR, t-test, or Wilcoxon)
3. Filters results by log2 fold-change and p-value cutoffs
4. Generates volcano plot and heatmap
5. Optionally generates Venn diagram (if related gene set provided) and chromosome location plot (if annotation provided)
6. Writes full results (Diffanalysis.csv), filtered DEGs (DEGs.csv), and provenance (.run_result.json)

Count-based methods (DESeq2, edgeR) expect integer count data. limma and statistical tests work with normalized expression values.

## Expected Input

- **Expression matrix** (`--mat`): CSV file with gene identifiers in the first column and sample expression values in subsequent columns. Genes x samples.
- **Sample group map** (`--map`): Two-column CSV mapping sample IDs (first column) to group labels (second column). Must contain exactly 2 unique groups.

Optional inputs:
- **Related gene set** (`--rgs`): Single-column CSV of gene identifiers for Venn diagram intersection with DEGs.
- **Chromosome annotation** (`--locate`): CSV with columns Gene, Chromosome, Start, End for chromosome location plotting.

## Output Files

| File | Description |
|------|-------------|
| `Diffanalysis.csv` | Full DE results: gene_id, logFC, Pvalue, Padj [, stat] |
| `DEGs.csv` | Filtered DEGs with group column (Up/Down/Not) |
| `Volcano.pdf` | Volcano plot with significance thresholds |
| `Heatmap.pdf` | Heatmap of top DEGs, row-scaled |
| `Venn.pdf` | Venn diagram of DEGs vs related gene set (conditional) |
| `Chromosome_location.pdf` | Chromosome location plot via RCircos (conditional) |
| `.run_result.json` | Provenance: parameters, status, exceptions, file list, timestamps |

## Exceptions

**Data issues (B-codes, exit 1):**
- B1_PROPORTION: Sample ratio exceeds 10:1. Override with `--force-imbalanced`.
- B2_FEW_DEGS: DEG count below `--cutoff` threshold.
- B3_MISSING_INPUT: Required input file not found.
- B4_INVALID_COLUMNS: Missing required columns in input files.
- B5_METHOD_MISMATCH: Selected method incompatible with data.
- B6_NO_DEGS: No DEGs found at given cutoffs.
- B7_RGS_INTERSECTION: Intersection with related gene set below cutoff.
- B8_EMPTY_MATRIX: Expression matrix has zero rows or columns.
- B9_SAMPLE_MISMATCH: Sample IDs don't match between matrix and map.

**Write issues (W-codes, exit 1):**
- W001_DISK_FULL: No space left on output device.
- W002_PERM_DENIED: Cannot write to output directory.

**Environment issues (E-codes):**
- E801_ENV_PKG (exit 3): Missing required R/Bioconductor packages.
- E802_UNSUPPORTED_TAXID (exit 1): Unsupported NCBI taxonomy ID for RCircos.

## Usage Examples

```bash
# Run DESeq2 analysis with defaults
Rscript scripts/main.R run \
  --mat counts.csv \
  --map groups.csv \
  --outdir ./results

# Run limma with custom cutoffs and Venn diagram
Rscript scripts/main.R run \
  --mat expr_norm.csv \
  --map groups.csv \
  --method limma \
  --logfc-cutoff 1.5 \
  --pvalue 0.01 \
  --rgs phenotype_genes.csv \
  --pheno-abbr "AS" \
  --color-panel "#E64B35,#4DBBD5" \
  --outdir ./results

# Validate inputs before running
Rscript scripts/main.R validate-input \
  --mat counts.csv \
  --map groups.csv

# Validate outputs after running
Rscript scripts/main.R validate-output \
  --outdir ./results
```
```

- [x] **Step 2: Verify SKILL.md frontmatter is valid YAML**

```bash
./env/bin/Rscript -e '
library(yaml)
fm <- yaml::read_yaml("SKILL.md")
stopifnot(!is.null(fm$name))
stopifnot(!is.null(fm$description))
stopifnot(fm$type == "standard")
stopifnot(length(fm$inputs) >= 1)
stopifnot(length(fm$outputs) >= 1)
stopifnot(!is.null(fm$entry))
stopifnot(length(fm$parameters) >= 1)
stopifnot(length(fm$exceptions) >= 1)
stopifnot(!is.null(fm$hardware))
cat("SKILL.md frontmatter OK: name=", fm$name, ", params=", length(fm$parameters), ", exceptions=", length(fm$exceptions), "\n", sep="")
'
```

Expected stdout: `SKILL.md frontmatter OK: name=differential-analysis, params=20, exceptions=13`

- [x] **Step 3: Commit SKILL.md**

```bash
git add SKILL.md
git commit -m "feat: add SKILL.md with complete v2 YAML frontmatter

8 required sections: name, description, type, inputs, outputs,
entry, parameters (20 total with bind annotations), exceptions
(13 structured codes B1-B9, W001-W002, E801-E802), hardware.

Body sections: Node Function, Expected Input, Output Files,
Exceptions (plain-language), Usage Examples.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

archived-with: 2026-06-22-create-node-package
---

## Task 7: Test Suite

**Files:**
- Create: `tests/testthat/test-input-validation.R`
- Create: `tests/testthat/test-output-validation.R`
- Create: `tests/testthat/test-diff-methods.R`
- Create: `tests/testthat/test-main.R`

**Interfaces:**
- Consumes: All scripts, test fixtures from setup.R
- Produces: testthat test suite with 4 test files

### 7.1 Write test-input-validation.R

- [x] **Step 1: Write tests/testthat/test-input-validation.R**

```r
# test-input-validation.R — Input validation tests for differential-analysis

library(testthat)

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
root_dir <- dirname(dirname(script_dir))

source(file.path(root_dir, "scripts", "report.R"))
source(file.path(root_dir, "scripts", "exceptions.R"))
source(file.path(root_dir, "scripts", "filter_helpers.R"))
source(file.path(root_dir, "scripts", "input_validation.R"))
source(file.path(root_dir, "tests", "testthat", "setup.R"))

td <- create_test_data()

test_that("validate_input passes with valid synthetic data", {
  opts <- list(mat = td$mat_path, map = td$map_path, force_imbalanced = FALSE)
  result <- validate_input(opts)
  expect_true(result$valid)
  expect_equal(result$n_genes, 10)
  expect_equal(result$n_samples, 6)
  expect_equal(result$groups, "Control vs Case")
})

test_that("validate_input detects missing mat file", {
  result <- validate_input(list(mat = "/nonexistent/file.csv", map = td$map_path))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("validate_input detects missing map file", {
  result <- validate_input(list(mat = td$mat_path, map = "/nonexistent/file.csv"))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("validate_input detects empty matrix", {
  empty_path <- file.path(tempdir(), "empty_expr.csv")
  write.csv(data.frame(gene = character(0)), empty_path, row.names = FALSE)
  result <- validate_input(list(mat = empty_path, map = td$map_path))
  expect_false(result$valid)
  expect_match(tolower(result$reason), "empty|0 rows")
  unlink(empty_path)
})

test_that("validate_input detects sample mismatch", {
  mismatch_path <- file.path(tempdir(), "mismatch_map.csv")
  map2 <- td$map
  map2$sample[1] <- "UNKNOWN_SAMPLE"
  write.csv(map2, mismatch_path, row.names = FALSE)
  result <- validate_input(list(mat = td$mat_path, map = mismatch_path))
  expect_false(result$valid)
  expect_match(result$reason, "not in matrix")
  unlink(mismatch_path)
})

test_that("proportion check passes for balanced design (ratio=1)", {
  map <- td$map
  map[[2]] <- factor(map[[2]], levels = c("Control", "Case"))
  expect_true(proportion_check(map))
})

test_that("proportion check fails for ratio > 10:1 without force", {
  # Create imbalanced map: 11 Case, 1 Control
  imbal_map <- rbind(
    data.frame(sample = paste0("CASE", 1:11), group = "Case"),
    data.frame(sample = "CTRL1", group = "Control")
  )
  colnames(imbal_map) <- colnames(td$map)
  imbal_map[[2]] <- factor(imbal_map[[2]], levels = c("Control", "Case"))
  expect_false(proportion_check(imbal_map))
})

test_that("proportion check passes for ratio > 10:1 with force", {
  imbal_map <- rbind(
    data.frame(sample = paste0("CASE", 1:11), group = "Case"),
    data.frame(sample = "CTRL1", group = "Control")
  )
  colnames(imbal_map) <- colnames(td$map)
  imbal_map[[2]] <- factor(imbal_map[[2]], levels = c("Control", "Case"))
  expect_true(proportion_check(imbal_map, force_imbalanced = TRUE))
})

test_that("proportion check fails if not exactly 2 groups", {
  map3 <- rbind(td$map, data.frame(sample = "EXTRA", group = "Third"))
  map3[[2]] <- factor(map3[[2]], levels = c("Control", "Case", "Third"))
  expect_false(proportion_check(map3))
})

cat("All input validation tests passed\n")
```

- [x] **Step 2: Run input validation tests**

```bash
./env/bin/Rscript tests/testthat/test-input-validation.R
```

Expected stdout: `All input validation tests passed` (all test_that blocks pass).

### 7.2 Write test-output-validation.R

- [x] **Step 3: Write tests/testthat/test-output-validation.R**

```r
# test-output-validation.R — Output validation tests for differential-analysis

library(testthat)

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
root_dir <- dirname(dirname(script_dir))

source(file.path(root_dir, "scripts", "report.R"))
source(file.path(root_dir, "scripts", "exceptions.R"))
source(file.path(root_dir, "scripts", "output_validation.R"))

test_that("validate_output detects missing directory", {
  result <- validate_output(list(outdir = "/nonexistent/dir"))
  expect_false(result$valid)
  expect_match(result$reason, "not found")
})

test_that("validate_output detects missing Diffanalysis.csv", {
  empty_dir <- file.path(tempdir(), "test_empty_output")
  dir.create(empty_dir, showWarnings = FALSE)
  result <- validate_output(list(outdir = empty_dir))
  expect_false(result$valid)
  expect_match(result$reason, "Missing Diffanalysis.csv")
  unlink(empty_dir, recursive = TRUE)
})

test_that("validate_output detects missing required columns", {
  outdir <- file.path(tempdir(), "test_bad_output")
  dir.create(outdir, showWarnings = FALSE)
  # Write Diffanalysis.csv without required columns
  write.csv(data.frame(x = 1:5, y = 6:10), file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  write.csv(data.frame(gene_id = "G1", group = "Up"), file.path(outdir, "DEGs.csv"), row.names = FALSE)
  # Create empty plot files
  file.create(file.path(outdir, "Volcano.pdf"))
  file.create(file.path(outdir, "Heatmap.pdf"))

  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "missing required columns")
  unlink(outdir, recursive = TRUE)
})

test_that("validate_output detects empty plot files", {
  outdir <- file.path(tempdir(), "test_empty_plots")
  dir.create(outdir, showWarnings = FALSE)
  write.csv(data.frame(gene_id = "G1", logFC = 1.5, Pvalue = 0.01, Padj = 0.05),
            file.path(outdir, "Diffanalysis.csv"), row.names = FALSE)
  write.csv(data.frame(gene_id = "G1", logFC = 1.5, Pvalue = 0.01, Padj = 0.05, group = "Up"),
            file.path(outdir, "DEGs.csv"), row.names = FALSE)
  file.create(file.path(outdir, "Volcano.pdf"))  # 0 bytes
  file.create(file.path(outdir, "Heatmap.pdf"))

  result <- validate_output(list(outdir = outdir))
  expect_false(result$valid)
  expect_match(result$reason, "empty")
  unlink(outdir, recursive = TRUE)
})

cat("All output validation tests passed\n")
```

- [x] **Step 4: Run output validation tests**

```bash
./env/bin/Rscript tests/testthat/test-output-validation.R
```

Expected stdout: `All output validation tests passed`

### 7.3 Write test-diff-methods.R

- [x] **Step 5: Write tests/testthat/test-diff-methods.R**

```r
# test-diff-methods.R — DE method output tests for differential-analysis

library(testthat)

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
root_dir <- dirname(dirname(script_dir))

source(file.path(root_dir, "scripts", "report.R"))
source(file.path(root_dir, "scripts", "diff_methods.R"))
source(file.path(root_dir, "tests", "testthat", "setup.R"))

td <- create_test_data()
map <- td$map
map[[2]] <- factor(map[[2]], levels = c("Control", "Case"))
mat <- td$mat

# Expected output columns for all methods
required_cols <- c("gene_id", "logFC", "Pvalue", "Padj")

# Known fold-change genes
known_up <- td$expected_up    # GENE1, logFC ~ 2
known_down <- td$expected_down  # GENE2, logFC ~ -2

test_that("diff_stat (t-test) produces expected output columns", {
  dif <- diff_stat(mat, map, "t")
  for (col in c(required_cols, "stat")) {
    expect_true(col %in% colnames(dif),
      info = sprintf("diff_stat missing column: %s", col))
  }
  expect_true(nrow(dif) > 0)
  # GENE1 should be near top (most significant)
  expect_true(known_up %in% dif$gene_id[1:5])
})

test_that("diff_stat (wilcoxon) produces expected output columns", {
  dif <- diff_stat(mat, map, "wilcox")
  for (col in c(required_cols, "stat")) {
    expect_true(col %in% colnames(dif),
      info = sprintf("diff_stat(wilcox) missing column: %s", col))
  }
  expect_true(nrow(dif) > 0)
})

test_that("diff_deseq2 produces expected output columns", {
  dif <- diff_deseq2(mat, map)
  for (col in required_cols) {
    expect_true(col %in% colnames(dif),
      info = sprintf("diff_deseq2 missing column: %s", col))
  }
  expect_true(nrow(dif) > 0)
  # GENE1 should be upregulated (positive logFC)
  g1_row <- dif[dif$gene_id == known_up, ]
  expect_true(nrow(g1_row) == 1)
  expect_true(g1_row$logFC > 0)
})

test_that("diff_limma produces expected output columns", {
  dif <- diff_limma(mat, map)
  for (col in required_cols) {
    expect_true(col %in% colnames(dif),
      info = sprintf("diff_limma missing column: %s", col))
  }
  expect_true(nrow(dif) > 0)
})

test_that("diff_edger produces expected output columns", {
  dif <- diff_edger(mat, map, "TMM", "glmFit")
  for (col in required_cols) {
    expect_true(col %in% colnames(dif),
      info = sprintf("diff_edger missing column: %s", col))
  }
  expect_true(nrow(dif) > 0)
})

test_that("all methods return data sorted by Pvalue ascending", {
  for (method in c("deseq2", "limma", "edgeR", "t", "wilcox")) {
    dif <- switch(method,
      deseq2 = diff_deseq2(mat, map),
      limma  = diff_limma(mat, map),
      edgeR  = diff_edger(mat, map),
      t      = diff_stat(mat, map, "t"),
      wilcox = diff_stat(mat, map, "wilcox")
    )
    pvals <- dif$Pvalue
    expect_true(all(diff(pvals) >= -1e-10),
      info = sprintf("%s results not sorted by Pvalue", method))
  }
})

test_that("known upregulated gene has positive logFC across methods", {
  for (method in c("deseq2", "limma", "edgeR", "t", "wilcox")) {
    dif <- switch(method,
      deseq2 = diff_deseq2(mat, map),
      limma  = diff_limma(mat, map),
      edgeR  = diff_edger(mat, map),
      t      = diff_stat(mat, map, "t"),
      wilcox = diff_stat(mat, map, "wilcox")
    )
    g1 <- dif[dif$gene_id == known_up, ]
    if (nrow(g1) > 0) {
      expect_true(g1$logFC > 0,
        info = sprintf("%s: %s should be upregulated", method, known_up))
    }
  }
})

cat("All diff methods tests passed\n")
```

- [x] **Step 6: Run DE methods tests**

```bash
./env/bin/Rscript tests/testthat/test-diff-methods.R
```

Expected stdout: `All diff methods tests passed`

### 7.4 Write test-main.R

- [x] **Step 7: Write tests/testthat/test-main.R**

```r
# test-main.R — Dispatch, parsing, and end-to-end tests for differential-analysis

library(testthat)

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
root_dir <- dirname(dirname(script_dir))

source(file.path(root_dir, "scripts", "report.R"))
source(file.path(root_dir, "scripts", "exceptions.R"))
source(file.path(root_dir, "scripts", "io_helpers.R"))
source(file.path(root_dir, "scripts", "diff_methods.R"))
source(file.path(root_dir, "scripts", "filter_helpers.R"))
source(file.path(root_dir, "scripts", "plot_helpers.R"))
source(file.path(root_dir, "scripts", "input_validation.R"))
source(file.path(root_dir, "scripts", "output_validation.R"))
source(file.path(root_dir, "scripts", "main.R"))
source(file.path(root_dir, "tests", "testthat", "setup.R"))

test_that("parse_args handles empty args (help text)", {
  expect_error(parse_args(character(0)), "quit|status")
})

test_that("parse_args rejects unknown subcommand", {
  expect_error(parse_args(c("unknown", "--mat", "x.csv")),
    "Unknown subcommand")
})

test_that("parse_args correctly parses --key=value form", {
  args <- c("run",
    "--mat=/tmp/expr.csv",
    "--map=/tmp/map.csv",
    "--method=limma",
    "--p-set=p",
    "--pvalue=0.01",
    "--logfc-cutoff=1.5",
    "--cutoff=5",
    "--outdir=/tmp/out"
  )
  opts <- parse_args(args)
  expect_equal(opts$subcommand, "run")
  expect_equal(opts$mat, "/tmp/expr.csv")
  expect_equal(opts$map, "/tmp/map.csv")
  expect_equal(opts$method, "limma")
  expect_equal(opts$p_set, "p")
  expect_equal(opts$pvalue, 0.01)
  expect_equal(opts$logfc_cutoff, 1.5)
  expect_equal(opts$cutoff, 5)
  expect_equal(opts$outdir, "/tmp/out")
})

test_that("parse_args correctly parses --key value form", {
  args <- c("run",
    "--mat", "/tmp/expr.csv",
    "--map", "/tmp/map.csv",
    "--method", "edgeR",
    "--norm", "RLE",
    "--model", "glmQLFit",
    "--top", "30",
    "--force-imbalanced"
  )
  opts <- parse_args(args)
  expect_equal(opts$method, "edgeR")
  expect_equal(opts$norm, "RLE")
  expect_equal(opts$model, "glmQLFit")
  expect_equal(opts$top, 30)
  expect_true(opts$force_imbalanced)
})

test_that("parse_args requires --mat and --map for run subcommand", {
  expect_error(parse_args(c("run", "--mat", "x.csv")),
    "requires --map")
  expect_error(parse_args(c("run", "--map", "x.csv")),
    "requires --mat")
})

test_that("parse_args applies defaults", {
  args <- c("run", "--mat", "/tmp/x.csv", "--map", "/tmp/y.csv")
  opts <- parse_args(args)
  expect_equal(opts$method, "deseq2")
  expect_equal(opts$p_set, "padj")
  expect_equal(opts$pvalue, 0.05)
  expect_equal(opts$logfc_cutoff, 1.0)
  expect_equal(opts$cutoff, 10)
  expect_equal(opts$top, 20)
  expect_equal(opts$outdir, ".")
})

test_that("check_environment finds required packages", {
  env <- check_environment()
  expect_equal(env$status, "ok")
})

test_that("end-to-end pipeline runs with synthetic data", {
  td <- create_test_data()
  outdir <- file.path(tempdir(), "deg_e2e_test")
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)

  # Capture stdout to parse NDJSON
  tmp_stdout <- file.path(tempdir(), "e2e_stdout.txt")
  exit_code <- system2(
    file.path(root_dir, "env", "bin", "Rscript"),
    c(file.path(root_dir, "scripts", "main.R"),
      "run",
      "--mat", td$mat_path,
      "--map", td$map_path,
      "--method", "t",
      "--logfc-cutoff", "0.5",
      "--cutoff", "1",
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

  unlink(outdir, recursive = TRUE)
  unlink(tmp_stdout)
})

cat("All main tests passed\n")
```

- [x] **Step 8: Run main tests**

```bash
./env/bin/Rscript tests/testthat/test-main.R
```

Expected stdout: `All main tests passed`

- [x] **Step 9: Commit test suite**

```bash
git add tests/testthat/test-input-validation.R tests/testthat/test-output-validation.R tests/testthat/test-diff-methods.R tests/testthat/test-main.R
git commit -m "feat: add complete test suite (4 test files, 25+ test cases)

test-input-validation.R: missing file, empty matrix, sample mismatch,
  proportion check boundary (ratio=1, ratio>10, force override), group count
test-output-validation.R: missing dir, missing files, missing columns, empty plots
test-diff-methods.R: all 5 methods produce required columns, known fold changes verified
test-main.R: arg parsing (both forms), dispatch, E2E pipeline with synthetic data

Co-Authored-By: Claude <noreply@anthropic.com>"
```

archived-with: 2026-06-22-create-node-package
---

## Task 8: Integration Verification

**Files:**
- None (verification only)

**Interfaces:**
- Consumes: All implemented files
- Produces: Test pass/fail report, verified output files

### 8.1 Run full test suite

- [x] **Step 1: Run all test files sequentially**

```bash
echo "=== Test 1: Input Validation ===" && \
./env/bin/Rscript tests/testthat/test-input-validation.R && \
echo "" && \
echo "=== Test 2: Output Validation ===" && \
./env/bin/Rscript tests/testthat/test-output-validation.R && \
echo "" && \
echo "=== Test 3: DE Methods ===" && \
./env/bin/Rscript tests/testthat/test-diff-methods.R && \
echo "" && \
echo "=== Test 4: Main / E2E ===" && \
./env/bin/Rscript tests/testthat/test-main.R && \
echo "" && \
echo "=== ALL TESTS PASSED ==="
```

Expected final output: `=== ALL TESTS PASSED ===`

### 8.2 Run end-to-end with synthetic data

- [x] **Step 2: Full pipeline run with synthetic data**

```bash
./env/bin/Rscript -e 'source("tests/testthat/setup.R"); td <- create_testData()' && \
rm -rf /tmp/deg_e2e_output && \
./env/bin/Rscript scripts/main.R run \
  --mat /tmp/test_expr.csv \
  --map /tmp/test_map.csv \
  --method deseq2 \
  --logfc-cutoff 0.5 \
  --cutoff 1 \
  --outdir /tmp/deg_e2e_output
```

Expected: NDJSON lines on stdout, final line has `"level":"result"` with `"status":"success"`.

- [x] **Step 3: Verify all output files**

```bash
ls -la /tmp/deg_e2e_output/
```

Expected: `Diffanalysis.csv`, `DEGs.csv`, `Volcano.pdf`, `Heatmap.pdf`, `.run_result.json` -- all with non-zero size.

- [x] **Step 4: Verify .run_result.json structure**

```bash
./env/bin/Rscript -e '
rr <- jsonlite::fromJSON("/tmp/deg_e2e_output/.run_result.json")
stopifnot(rr$node == "differential-analysis")
stopifnot(rr$subcommand == "run")
stopifnot(rr$exit_code == 0)
stopifnot(!is.null(rr$parameters))
stopifnot(!is.null(rr$files))
cat(".run_result.json OK: node=", rr$node, ", exit=", rr$exit_code, ", files=", length(rr$files), "\n", sep="")
'
```

Expected stdout: `.run_result.json OK: node=differential-analysis, exit=0, files=4`

### 8.3 Verify SKILL.md frontmatter

- [x] **Step 5: Validate SKILL.md frontmatter completeness**

```bash
./env/bin/Rscript -e '
library(yaml)
# Read only the YAML frontmatter between --- markers
txt <- readLines("SKILL.md")
sep_idx <- which(txt == "---")
yaml_txt <- txt[(sep_idx[1]+1):(sep_idx[2]-1)]
fm <- yaml::yaml.load(paste(yaml_txt, collapse = "\n"))

checks <- list(
  name = !is.null(fm$name) && nchar(fm$name) > 0,
  description = !is.null(fm$description) && nchar(fm$description) > 0,
  type = fm$type == "standard",
  inputs = length(fm$inputs) >= 2,
  outputs = length(fm$outputs) >= 4,
  entry = !is.null(fm$entry) && grepl("main", fm$entry),
  parameters = length(fm$parameters) >= 15,
  exceptions = length(fm$exceptions) >= 10,
  hardware = !is.null(fm$hardware) && !is.null(fm$hardware$memory_gb)
)

all_ok <- all(unlist(checks))
for (n in names(checks)) {
  cat(sprintf("  %-15s %s\n", paste0(n, ":"), if (checks[[n]]) "PASS" else "FAIL"))
}
cat(if (all_ok) "\nALL FRONTMATTER CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (all_ok) 0 else 1)
'
```

Expected: All 9 checks PASS.

- [x] **Step 6: Commit verification results**

No code changes at this stage — verification is read-only. If any issues were found, iterate on the failing task's code.

archived-with: 2026-06-22-create-node-package
---

## Verification Checklist

After all 8 tasks are complete, confirm:

1. **Environment**: `conda env create -f envs/env-r-4.3.yaml -p ./env` succeeds
2. **All R scripts load**: Each file in `scripts/` can be `source()`d without error
3. **Test suite passes**: All 4 test files produce "All ... tests passed"
4. **E2E run produces 5 files**: Diffanalysis.csv, DEGs.csv, Volcano.pdf, Heatmap.pdf, .run_result.json
5. **NDJSON output is valid**: Every stdout line parses as JSON with `level` field
6. **SKILL.md frontmatter**: Parses as valid YAML with all 8 required sections non-empty
7. **No hardcoded secrets**: `grep -r "key\|token\|password\|secret" scripts/` returns nothing in source
8. **Git history is clean**: Each commit is a single logical change with descriptive messages
