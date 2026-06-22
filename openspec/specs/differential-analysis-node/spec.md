# differential-analysis-node Specification

## Purpose
TBD - created by archiving change create-node-package. Update Purpose after archive.
## Requirements
### Requirement: Node package provides SKILL.md contract

The node package SHALL include a `SKILL.md` file at the repository root with valid YAML frontmatter containing all 8 required fields: name, description, type, inputs, outputs, entry, parameters, exceptions, and hardware.

#### Scenario: SKILL.md frontmatter is valid YAML

- **WHEN** the SKILL.md file is parsed by a YAML parser
- **THEN** the frontmatter SHALL contain `name: differential-analysis`, `type: standard`, non-empty `description`, at least two `inputs` entries (`--mat`, `--map`), at least four `outputs` entries (Diffanalysis.csv, DEGs.csv, Volcano.pdf, Heatmap.pdf), `entry: scripts/main.R`, at least 15 `parameters` entries including `subcommand`, `--mat`, `--map`, `--method`, `--outdir`, at least 10 `exceptions` entries with structured codes, and a `hardware` section with `memory_gb: 4`, `cpu: 2`, `gpu: false`, and runtime estimate

#### Scenario: All exception patterns are in English

- **WHEN** the SKILL.md exceptions section is inspected
- **THEN** all `pattern` values SHALL be English substrings

### Requirement: Node accepts subcommand dispatch via main.R

The node SHALL provide `scripts/main.R` as its single entry point. The first positional argument SHALL be the subcommand name. Valid subcommands are `run`, `validate-input`, and `validate-output`.

#### Scenario: Run subcommand executes full pipeline

- **WHEN** `Rscript scripts/main.R run --mat <expr.csv> --map <group.csv> --method deseq2 --outdir <path>` is invoked with valid inputs
- **THEN** the script SHALL execute proportion check, DE analysis, filtering, volcano plot, and heatmap generation
- **AND** exit with code 0

#### Scenario: Validate-input subcommand checks inputs

- **WHEN** `Rscript scripts/main.R validate-input --mat <expr.csv> --map <group.csv>` is invoked
- **THEN** the script SHALL verify file existence, column presence, and group structure
- **AND** exit 0 if valid, non-zero with stderr reason if invalid

#### Scenario: Validate-output subcommand checks outputs

- **WHEN** `Rscript scripts/main.R validate-output --outdir <path>` is invoked
- **THEN** the script SHALL verify output CSV columns, non-empty rows, and plot file existence
- **AND** exit 0 if valid, non-zero with stderr reason if invalid

#### Scenario: Unknown subcommand produces error

- **WHEN** `Rscript scripts/main.R unknown --outdir <path>` is invoked
- **THEN** the script SHALL write an NDJSON error line to stdout and exit with code 1

### Requirement: Node reports progress via NDJSON stdout

The node SHALL write valid JSON to stdout, one object per line. Each line SHALL have a `level` field. Supported levels: `info`, `result`, `exception`.

#### Scenario: Progress lines during execution

- **WHEN** the `run` subcommand is executing
- **THEN** at least one `{"level":"info","msg":"..."}` line SHALL be written to stdout before the result line
- **AND** all info messages SHALL be in English

#### Scenario: Result line on success

- **WHEN** the pipeline completes successfully
- **THEN** the final stdout line SHALL be `{"level":"result","status":"success","files":[...],"metadata":{...}}` listing all output files with paths and key metadata (method, version, DEG counts)

#### Scenario: Exception line on error

- **WHEN** a declared exception condition is encountered
- **THEN** the script SHALL write `{"level":"exception","code":"...","nature":"...","action":"...","msg":"..."}` to stdout

### Requirement: Node persists provenance record

The node SHALL write a `.run_result.json` file to `--outdir` containing parameters, status, exit code, timestamps, output metadata, exceptions, and file list.

#### Scenario: Run result written on completion

- **WHEN** the `run` subcommand completes (success or handled error)
- **THEN** a `.run_result.json` file SHALL exist in `--outdir`
- **AND** the file SHALL contain `node: "differential-analysis"`, `subcommand`, `status`, `exit_code`, `started_at`, `finished_at`, `parameters`, `files`

### Requirement: Node runs environment check before dispatch

The `main()` function SHALL call `check_environment()` before any subcommand dispatch to verify all required R packages are installed.

#### Scenario: Missing package detected

- **WHEN** a required R package is not installed
- **THEN** the script SHALL emit `{"level":"exception","code":"E801_ENV_PKG",...}` and exit with code 3

### Requirement: Node supports five DE methods

The node SHALL support five differential expression methods selectable via the `--method` parameter: `deseq2`, `limma`, `edgeR`, `t`, and `wilcox`.

#### Scenario: DESeq2 method for count data

- **WHEN** `--method deseq2` is specified with an integer count matrix
- **THEN** DESeq2 SHALL be used for DE analysis
- **AND** results SHALL include log2FoldChange, pvalue, and padj columns

#### Scenario: Method/data type mismatch

- **WHEN** DESeq2 or edgeR is specified with a non-count matrix
- **THEN** the script SHALL emit exception `B5_METHOD_MISMATCH` and exit with code 1

#### Scenario: Invalid method produces error

- **WHEN** `--method invalid` is specified
- **THEN** the script SHALL write an error to stdout and exit with code 1

### Requirement: Node provides conda environment definition

The node SHALL include `envs/env-r-4.3.yaml` declaring all R and Bioconductor package dependencies needed to run the pipeline.

#### Scenario: Environment file declares R 4.3 and Bioconductor packages

- **WHEN** the env.yaml is inspected
- **THEN** it SHALL declare `r-base>=4.3` in its dependencies
- **AND** it SHALL include channels: conda-forge and bioconda
- **AND** it SHALL include: bioconductor-deseq2, bioconductor-limma, bioconductor-edger, bioconductor-rcircos, r-yaml, r-filelock, r-dplyr, r-data.table, r-ggplot2, r-ggrepel, r-ggvenn, r-pheatmap, r-jsonlite, r-testthat

### Requirement: Node validates inputs before analysis

The `scripts/input_validation.R` script SHALL verify that input files exist, have required columns, and contain valid data before the analysis proceeds. The proportion check SHALL be a reusable function also callable from `main.R`'s pipeline.

#### Scenario: Missing input file detected

- **WHEN** `--mat` points to a non-existent file
- **THEN** the script SHALL emit exception `B3_MISSING_INPUT` and exit with code 1

#### Scenario: Sample proportion exceeds threshold

- **WHEN** the case/control sample ratio exceeds 10:1 and `--force-imbalanced` is not set
- **THEN** the script SHALL emit exception `B1_PROPORTION` and exit with code 1

#### Scenario: Force-imbalanced override

- **WHEN** the case/control sample ratio exceeds 10:1 and `--force-imbalanced` is true
- **THEN** the pipeline SHALL proceed past the proportion check with an info-level warning

### Requirement: Node validates outputs after analysis

The `scripts/output_validation.R` script SHALL verify that output files have expected columns, non-zero rows, and valid content.

#### Scenario: DEG results have required columns

- **WHEN** output validation runs on the DEG results CSV
- **THEN** the CSV SHALL contain columns: gene identifier, logFC, Pvalue, Padj
- **AND** the file SHALL have at least one data row

### Requirement: Node handles errors with structured exceptions

The node SHALL wrap DE method dispatch in error handling that emits structured NDJSON exceptions on unexpected Bioconductor failures, and SHALL write `.run_result.json` even on DE analysis failure paths.

#### Scenario: DE method throws unexpected error

- **WHEN** a DE method (DESeq2, limma, edgeR, t-test, Wilcoxon) throws an unexpected error during execution
- **THEN** the script SHALL emit a structured NDJSON exception with level `"exception"` and write `.run_result.json` with `status: "error"`
- **AND** exit with code 1

#### Scenario: Consecutive main() calls produce independent results

- **WHEN** `main()` is called twice in the same R session
- **THEN** the second call's `.exceptions` accumulator SHALL contain only exceptions from the second run
- **AND** the second call's `.run_result.json` SHALL not contain exceptions from the first run

### Requirement: Node supports edgeR configurable parameters

The node SHALL expose edgeR normalization method via `--norm` (choices: TMM, RLE, upperquartile, none, default TMM) and model fitting via `--model` (choices: glmFit, glmQLFit, default glmFit).

#### Scenario: edgeR with custom normalization

- **WHEN** `--method edgeR --norm RLE` is specified
- **THEN** edgeR SHALL use `calcNormFactors(method = "RLE")` for normalization

