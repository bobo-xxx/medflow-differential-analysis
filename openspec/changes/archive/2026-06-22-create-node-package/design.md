## Context

The `deg-analysis` node performs differential expression analysis comparing case vs control groups. Seven reference scripts in `original/` implement the pipeline: proportion check → DE analysis (5 methods) → filter → volcano → heatmap → (optional venn + chromosome location). These scripts use CLI args directly and write CSVs/PDFs without structured reporting. They need to be adapted into the IRE node package format with standardized contracts, NDJSON reporting, and a single entry point.

**Constraints:**
- R 4.3 environment via conda (Bioconductor packages from bioconda channel)
- Single entry point at `scripts/main.R`
- English-only for all artifacts, messages, and error patterns
- No hardcoded secrets or paths
- Follow the `medflow-geo-microarray` reference node patterns (report.R, exceptions.R, .run_result.json)

## Goals / Non-Goals

**Goals:**
- Wrap all 7 reference scripts into a single coherent node with `run`, `validate-input`, `validate-output` subcommands
- Define complete SKILL.md v2 frontmatter contract (inputs, outputs, 20 parameters, 13 structured exceptions, hardware)
- Provide conda environment with all R/Bioconductor dependencies
- Implement NDJSON stdout reporting with `report_info()`, `report_result()`, `report_exception_ndjson()`
- Persist `.run_result.json` provenance file per `medflow-geo-microarray` convention
- Provide input and output validation scripts (standalone + as subcommands)
- Write testthat test suite with synthetic fixtures

**Non-Goals:**
- New DE methods or algorithm improvements
- Performance optimization beyond the reference
- New visualization types
- Production test data generation
- Core framework integration (registry.yaml entry)

## Decisions

### 1. Module Architecture

**Decision:** Nine internal modules under `scripts/`, all sourced by `main.R`.

```
scripts/
├── main.R                    # Entry point: arg parsing, env check, subcommand dispatch
├── report.R                  # report_info(), report_result(), report_error(),
│                             #   report_exception_ndjson(), write_run_result()
├── exceptions.R              # Exception accumulator, exception NDJSON emitter
├── diff_methods.R            # diff_deseq2(), diff_limma(), diff_edger(), diff_stat()
├── plot_helpers.R            # plot_volcano(), plot_heatmap(), plot_venn(), plot_locate()
├── filter_helpers.R          # test_cutoff(), filter_degs(), proportion_check()
├── io_helpers.R              # file_lock(), create_file_dir(), color_map()
├── input_validation.R        # Pre-flight checks (file, columns, proportion)
└── output_validation.R       # Post-run checks (columns, non-empty, file existence)
```

### 2. Subcommand Design

**Decision:** Three subcommands on `main.R`: `run` (full pipeline), `validate-input`, `validate-output`.

`validate-input` and `validate-output` are both standalone scripts AND `main.R` subcommands, following the `medflow-geo-microarray` convention.

### 3. Reporting Pattern

**Decision:** NDJSON to stdout (ephemeral, framework-captured) + `.run_result.json` to `--outdir` (persisted).

- `{"level":"info","msg":"..."}` — progress
- `{"level":"result","status":"success","files":[...],"metadata":{...}}` — final output
- `{"level":"exception","code":"B1_PROPORTION","nature":"data_insufficient","action":"skip_with_warning","msg":"..."}` — structured error
- Exceptions reported as NDJSON to stdout (NOT bare stderr), matching the reference node pattern
- `report_error()` for terminal failures — writes error NDJSON then quits

### 4. Parameter Binding Strategy

20 parameters covering DE method selection, filtering thresholds, edgeR-specific options, visualization tuning, and escape hatches. Full schema in the Design Doc at `docs/superpowers/specs/2026-06-22-deg-analysis-design.md`.

Key bindings:
- `--mat`, `--map`, `--rgs`, `--locate` → `upstream` (wired from prior nodes)
- `--method`, `--tax-id`, `--pheno-abbr` → `config` (protocol-level choices)
- `--pvalue`, `--logfc-cutoff`, `--cutoff`, `--norm`, `--model`, `--top`, `--force-imbalanced`, etc. → `static` (tuning knobs with defaults)
- `--outdir` → `framework`

### 5. Exception Model

**Decision:** 13 structured exception codes following the B=data / W=write / E=environment prefix scheme from the reference node.

B1_PROPORTION through B9_SAMPLE_MISMATCH cover data issues. W001_DISK_FULL and W002_PERM_DENIED cover write failures. E801_ENV_PKG (exit 3) and E802_UNSUPPORTED_TAXID cover environment/config issues.

All reported via `report_exception_ndjson(code, nature, action, msg)` to stdout.

### 6. Single-language (R)

**Decision:** Rewrite `diff_check_proportion.py` as an R function in `filter_helpers.R`, sourced by both `main.R` (runtime) and `input_validation.R` (standalone pre-flight).

### 7. Environment Check

**Decision:** `check_environment()` runs before any subcommand dispatch in `main()`. Verifies all required R packages are installed. On missing packages, emits `E801_ENV_PKG` and exits 3.

### 8. env.yaml Location

**Decision:** Place at `envs/env-r-4.3.yaml` per the existing `envs/` directory convention.

## Risks / Trade-offs

- [Bioconductor package availability] Some packages (RCircos) may have version constraints with R 4.3 → Mitigation: declare versions loosely, test resolution during build phase
- [RCircos species support] Only human (9606) and mouse (10090) supported → Mitigation: `E802_UNSUPPORTED_TAXID` catches other species
- [Method/data mismatch] DESeq2/edgeR require counts, limma expects normalized → Mitigation: `B5_METHOD_MISMATCH` guards against wrong combinations
- [Test data] Synthetic fixtures needed → Mitigation: generate minimal 10-gene × 6-sample matrices inline in test setup

## Open Questions

- Exact conda channel priority for bioconda vs conda-forge — resolve during env testing
