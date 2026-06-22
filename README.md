# medflow-deg-analysis

Differential expression analysis (transcriptomics) package in medflow.

## Quick Start

```bash
# Create conda environment
conda env create -f envs/env-r-4.3.yaml -p ./env

# Activate environment
conda activate ./env

# Run DEG analysis
./env/bin/Rscript scripts/main.R run \
  --mat expression_matrix.csv \
  --map sample_group_map.csv \
  --method deseq2 \
  --outdir ./output
```

## Node Package Format

| File | Description |
|------|-------------|
| `SKILL.md` | Agent contract with v2 YAML frontmatter (inputs, outputs, parameters, exceptions, hardware) |
| `envs/env-r-4.3.yaml` | Conda environment with R 4.3 and Bioconductor packages |
| `scripts/main.R` | Single entry point with subcommand dispatch (run, validate-input, validate-output) |
| `scripts/report.R` | NDJSON reporting helpers and .run_result.json provenance |
| `scripts/exceptions.R` | Structured exception handling with B/W/E code prefixes |
| `scripts/diff_methods.R` | DE analysis: DESeq2, limma, edgeR, t-test, Wilcoxon |
| `scripts/filter_helpers.R` | Proportion check, logFC/p-value filtering |
| `scripts/plot_helpers.R` | Volcano plot, heatmap, Venn diagram, chromosome location |
| `scripts/io_helpers.R` | File locking, directory creation, color utilities |
| `scripts/input_validation.R` | Pre-flight input checks |
| `scripts/output_validation.R` | Post-hoc output checks |
| `tests/testthat/` | Test suite (390 tests) |

## Supported Methods

| Method | Data Type | Package |
|--------|-----------|---------|
| DESeq2 | Integer counts | DESeq2 |
| edgeR | Integer counts | edgeR |
| limma | Normalized expression | limma |
| t-test | Any numeric | stats |
| Wilcoxon | Any numeric | stats |

## Reference

Implementation logic adapted from `original/scripts/`:
- `diff.R` — Core DE analysis functions
- `diff_filter.R` — logFC/p-value filtering and DEG classification
- `diff_volcano.R` — Volcano plot via ggplot2
- `diff_heatmap.R` — Heatmap via pheatmap
- `diff_venn.R` — Venn diagram via ggvenn
- `diff_locate.R` — Chromosome location via RCircos
- `diff_check_proportion.py` — Sample proportion validation

## License

MIT
