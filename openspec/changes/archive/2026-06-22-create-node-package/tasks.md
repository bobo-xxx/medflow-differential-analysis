## 1. Environment Setup

- [x] 1.1 Create `envs/env-r-4.3.yaml` with conda-forge + bioconda channels, r-base>=4.3, and all Bioconductor/R package dependencies (DESeq2, limma, edgeR, statmod, ggplot2, ggrepel, ggvenn, pheatmap, RCircos, yaml, filelock, dplyr, data.table, jsonlite, testthat)
- [x] 1.2 Create `tests/testthat/` directory and test infrastructure (setup.R with synthetic 10-gene × 6-sample count/normalized fixtures)

## 2. Reporting & Exceptions

- [x] 2.1 Create `scripts/report.R` with report_info(), report_result(), report_error(), report_exception_ndjson(), write_run_result() following medflow-geo-microarray conventions
- [x] 2.2 Create `scripts/exceptions.R` with exception accumulator (.exceptions list) and structured code definitions (B=data, W=write, E=environment)

## 3. Internal Helper Scripts

- [x] 3.1 Create `scripts/diff_methods.R` adapting DE functions (deseq2, limma, edgeR, t-test, wilcoxon) from `original/scripts/diff.R`
- [x] 3.2 Create `scripts/filter_helpers.R` with proportion_check(), test_cutoff(), filter_degs() adapted from `original/scripts/diff_filter.R` and `original/scripts/diff_check_proportion.py`
- [x] 3.3 Create `scripts/plot_helpers.R` adapting volcano, heatmap, venn, and chromosome location plotting from `original/scripts/diff_volcano.R`, `diff_heatmap.R`, `diff_venn.R`, `diff_locate.R`
- [x] 3.4 Create `scripts/io_helpers.R` with file_lock(), create_file_dir(), color_map() utilities

## 4. Validation Scripts

- [x] 4.1 Create `scripts/input_validation.R` with checks for file existence, required columns (gene identifier, sample columns), group structure (exactly 2 groups), and sample proportion ratio ≤ 10:1 (source proportion_check from filter_helpers.R; accepts --mode flag for standalone vs main.R usage)
- [x] 4.2 Create `scripts/output_validation.R` with checks for output CSV column presence (gene_id, logFC, Pvalue, Padj), non-empty results, and plot file existence

## 5. Main Entry Point

- [x] 5.1 Create `scripts/main.R` with check_environment() (verifies R packages, exits 3 on missing), parse_args() (20 parameters, --key=value and --key value forms), subcommand dispatch: do_run (full pipeline: proportion_check → diff_analysis → filter_degs → plot_volcano → plot_heatmap → optional plot_venn → optional plot_locate → write_run_result), do_validate_input, do_validate_output

## 6. SKILL.md Contract

- [x] 6.1 Create `SKILL.md` with complete v2 YAML frontmatter: name, description, type, inputs (--mat expression_matrix CSV, --map sample_group_map CSV), outputs (Diffanalysis.csv, DEGs.csv, Volcano.pdf, Heatmap.pdf, conditional Venn.pdf and Chromosome_location.pdf), entry, 20 parameters with bind annotations (upstream/config/static/framework), 13 structured exceptions (B1-B9, W001-W002, E801-E802) with code/pattern/nature/action, hardware (4 GB, 2 CPU, no GPU, ~2-10 min)
- [x] 6.2 Write SKILL.md body sections: Node Function, Expected Input, Output Files, Exceptions (plain-language), Usage Examples

## 7. Test Suite

- [x] 7.1 Write `tests/testthat/test-input-validation.R` testing: missing file detection, invalid column handling, proportion check boundary (ratio=1, ratio=12, force-imbalanced override)
- [x] 7.2 Write `tests/testthat/test-output-validation.R` testing: column presence, empty result detection, missing plot file
- [x] 7.3 Write `tests/testthat/test-diff-methods.R` testing: all 5 DE methods produce expected column output with synthetic count/matrix data, method mismatch detection
- [x] 7.4 Write `tests/testthat/test-main.R` testing: subcommand dispatch (run, validate-input, validate-output, unknown), parameter parsing (--key=value, --key value), end-to-end pipeline with synthetic input, environment check

## 8. Integration Verification

- [x] 8.1 Run full test suite and confirm all tests pass
- [x] 8.2 Run `scripts/main.R run` with synthetic data end-to-end and verify all output files + .run_result.json are produced
- [x] 8.3 Verify SKILL.md frontmatter parses as valid YAML with all 8 required fields present and non-empty
