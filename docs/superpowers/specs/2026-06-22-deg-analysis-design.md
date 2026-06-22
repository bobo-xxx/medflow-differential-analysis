---
comet_change: create-node-package
role: technical-design
canonical_spec: openspec
archived-with: 2026-06-22-create-node-package
status: final
---

# differential-analysis Node Package — Technical Design

## Architecture Overview

```
                    ┌──────────────────────────┐
                    │       scripts/main.R      │
                    │  parse_args() → dispatch  │
                    │  main() {                 │
                    │    check_environment()    │
                    │    switch(subcommand)     │
                    │  }                        │
                    └─────┬────────────────────┘
                          │
          ┌───────────────┼───────────────┐
          ▼               ▼               ▼
    ┌──────────┐   ┌──────────┐   ┌──────────────┐
    │   run    │   │validate- │   │validate-     │
    │          │   │  input   │   │  output      │
    └────┬─────┘   └──────────┘   └──────────────┘
         │
         ▼
┌─────────────────────────────────────────────────┐
│              Pipeline (run subcommand)           │
│                                                 │
│  proportion_check() → diff_analysis()           │
│  → filter_degs() → plot_volcano()               │
│  → plot_heatmap() → [plot_venn()]               │
│  → [plot_locate()] → write_run_result()         │
│                                                 │
│  All steps emit NDJSON via report_info()        │
│  Final output: report_result() + .run_result.json│
└─────────────────────────────────────────────────┘
```

## Module Responsibilities

### `scripts/main.R`
Single entry point. Parses CLI args, dispatches to subcommand handlers.
Sources all helper modules. Runs `check_environment()` before dispatch.

**Subcommands:**
| Subcommand | Handler | Purpose |
|------------|---------|---------|
| `run` | `do_run(opts)` | Full DEG analysis pipeline |
| `validate-input` | `do_validate_input(opts)` | Pre-flight input validation |
| `validate-output` | `do_validate_output(opts)` | Post-hoc output validation |

### `scripts/report.R`
NDJSON reporting to stdout (ephemeral, framework-captured):
- `report_info(msg, ...)` — progress messages
- `report_result(status, files, metadata, ...)` — final output
- `report_error(msg, exit_code)` — terminal failure
- `report_exception_ndjson(code, nature, action, msg)` — structured exception
- `write_run_result(out_dir, result, params, exit_code, times)` — persisted `.run_result.json`

### `scripts/exceptions.R`
Exception accumulator (`.exceptions` list) and `report_exception_ndjson()`.
Structured codes: B=data, W=write, E=environment.

### `scripts/diff_methods.R`
Five DE analysis functions adapted from `original/scripts/diff.R`:
- `diff_deseq2(mat, group)` — DESeq2 for count data
- `diff_limma(mat, group)` — limma with empirical Bayes
- `diff_edger(mat, group, norm, model)` — edgeR with configurable normalization
- `diff_stat(mat, group, stat)` — t-test or Wilcoxon

### `scripts/filter_helpers.R`
- `proportion_check(map, force_imbalanced)` — ratio must be ≤ 10:1 unless overridden
- `test_cutoff(dif, fc_name, p_name, logfc_test, p_value)` — sensitivity analysis
- `filter_degs(dif, fc_name, p_name, logfc_cutoff, p_value, cutoff)` — DEG classification

### `scripts/plot_helpers.R`
- `plot_volcano(dif, p_name, p_value, logfc_cutoff, top, gene, outfile)` — ggplot2
- `plot_heatmap(mat, map, rdegs, top, color_heat, outfile)` — pheatmap
- `plot_venn(dif, rgs, pheno_abbr, color_panel, outfile)` — ggvenn (conditional)
- `plot_locate(gene_list, locate, tax_id, outfile)` — RCircos (conditional)

### `scripts/io_helpers.R`
- `file_lock(path, FUN, ...)` — advisory file locking
- `create_file_dir(file)` — mkdir -p for output paths
- `color_map(colors, groups)` — assign colors to named groups

### `scripts/input_validation.R`
Standalone pre-flight checks, also sourceable by `main.R`:
- File existence and readability
- Required column presence (gene identifier, sample columns match group map)
- Group structure validation (exactly 2 groups)
- Sample proportion ratio ≤ 10:1 (or `--force-imbalanced` set)

### `scripts/output_validation.R`
Post-run checks:
- Output CSV column presence (gene_id, logFC, Pvalue, Padj)
- Non-empty results
- Plot file existence and non-zero size

## Parameter Schema

| Parameter | Type | Default | Range | Bind | Description |
|-----------|------|---------|-------|------|-------------|
| `subcommand` | choice: run, validate-input, validate-output | — | — | config | Operation to perform |
| `--mat` | file_in | — | — | upstream | Expression matrix CSV (genes × samples) |
| `--map` | file_in | — | — | upstream | Sample-to-group mapping CSV |
| `--method` | choice: deseq2, limma, edgeR, t, wilcox | deseq2 | — | config | DE analysis method |
| `--p-set` | choice: p, padj | padj | — | static | Which p-value to use for filtering |
| `--pvalue` | float | 0.05 | [0.001, 0.25] | static | P-value threshold |
| `--logfc-cutoff` | float | 1.0 | [0.0, 10.0] | static | Absolute log2 fold-change cutoff |
| `--cutoff` | int | 10 | — | static | Minimum DEG count after filtering |
| `--norm` | choice: TMM, RLE, upperquartile, none | TMM | — | static | edgeR normalization method |
| `--model` | choice: glmFit, glmQLFit | glmFit | — | static | edgeR model fitting method |
| `--top` | int | 20 | — | static | Top N genes for heatmap/volcano labels |
| `--force-imbalanced` | bool | false | — | static | Override 10:1 proportion check |
| `--rgs` | file_in | — | — | upstream | Related gene set for Venn diagram |
| `--locate` | file_in | — | — | upstream | Chromosome annotation for location plot |
| `--tax-id` | string | 9606 | — | config | NCBI taxonomy ID |
| `--pheno-abbr` | string | — | — | config | Phenotype abbreviation for Venn label |
| `--gene` | string | — | — | static | Comma-separated genes to label on volcano |
| `--color-heat` | string | blue,white,red | — | static | Heatmap color palette |
| `--color-panel` | string | — | — | static | Comma-separated colors for Venn sets |
| `--outdir` | file_out | . | — | framework | Output directory |

## Exception Codes

| Code | Exit | Pattern | Nature | Action |
|------|------|---------|--------|--------|
| B1_PROPORTION | 1 | `sample proportion.*exceeds` | data_insufficient | skip_with_warning |
| B2_FEW_DEGS | 1 | `gene count.*below.*cutoff` | data_insufficient | skip_with_warning |
| B3_MISSING_INPUT | 1 | `input file not found` | data_corrupt | halt |
| B4_INVALID_COLUMNS | 1 | `missing required column` | data_corrupt | halt |
| B5_METHOD_MISMATCH | 1 | `not support.*matrix` | data_mismatch | halt |
| B6_NO_DEGS | 1 | `no differentially expressed` | data_insufficient | skip_with_warning |
| B7_RGS_INTERSECTION | 1 | `intersection.*below.*cutoff` | data_insufficient | skip_with_warning |
| B8_EMPTY_MATRIX | 1 | `empty\|no rows\|no columns` | data_insufficient | halt |
| B9_SAMPLE_MISMATCH | 1 | `sample.*mismatch\|columns.*not.*match` | data_corrupt | halt |
| W001_DISK_FULL | 1 | `disk full\|no space` | resource | halt |
| W002_PERM_DENIED | 1 | `permission denied` | resource | halt |
| E801_ENV_PKG | 3 | `Missing required packages` | env_bug | halt |
| E802_UNSUPPORTED_TAXID | 1 | `Unsupported tax_id` | data_mismatch | halt |

## Output Files

All written to `--outdir`:

| File | Format | Conditional | Description |
|------|--------|-------------|-------------|
| `.run_result.json` | JSON | always | Provenance: params, status, exceptions, file list, timestamps |
| `Diffanalysis.csv` | CSV | always | Full DE results: gene_id, logFC, Pvalue, Padj, stat |
| `DEGs.csv` | CSV | always | Filtered DEGs with group column (Up/Down/Not) |
| `Volcano.pdf` | PDF | always | Volcano plot (ggplot2 + ggrepel) |
| `Heatmap.pdf` | PDF | always | Heatmap of top DEGs (pheatmap, row-scaled) |
| `Venn.pdf` | PDF | `--rgs` provided | Venn diagram (DEGs ∩ phenotype genes) |
| `Chromosome_location.pdf` | PDF | `--locate` provided | Chromosome location plot (RCircos) |

## NDJSON Contract

Progress (stdout, ephemeral):
```json
{"level":"info","msg":"Checking sample proportion (8 vs 8)..."}
{"level":"info","msg":"Running DESeq2 differential expression..."}
{"level":"info","msg":"342 DEGs: 187 up, 155 down (|logFC|>1.0, Padj<0.05)"}
```

Result (stdout, ephemeral):
```json
{"level":"result","status":"success","files":[
  {"path":"Diffanalysis.csv","rows":18420,"cols":5},
  {"path":"Volcano.pdf"},{"path":"Heatmap.pdf"}
],"metadata":{"method":"deseq2","version":"1.42.0","n_degs":342,"up":187,"down":155}}
```

Exception (stdout, ephemeral):
```json
{"level":"exception","code":"B1_PROPORTION","nature":"data_insufficient","action":"skip_with_warning","msg":"Sample proportion 15:1 exceeds 10:1 limit"}
```

## Environment Check

`check_environment()` runs before any subcommand dispatch:
1. Verify all required R packages are installed
2. If missing, emit `report_exception_ndjson("E801_ENV_PKG", ...)` and exit 3
3. Report package versions in info NDJSON lines

## Hardware

```yaml
memory_gb: 4
cpu: 2
gpu: false
runtime: "~2-10 minutes depending on method and matrix size"
```

## Testing Strategy

```
tests/testthat/
├── test-input-validation.R   # File existence, columns, proportion check
├── test-output-validation.R  # Column presence, non-empty results
├── test-diff-methods.R       # All 5 methods produce expected output columns
└── test-main.R               # Arg parsing, dispatch, end-to-end pipeline
```

Synthetic fixtures: 10-gene × 6-sample count matrix (3 case + 3 control) with
known fold changes for 2 genes.

## Open Questions

- Exact conda channel priority for bioconda vs conda-forge — resolve during env testing
- RCircos data availability for non-human/non-mouse species — blocked by Bioconductor data packages
