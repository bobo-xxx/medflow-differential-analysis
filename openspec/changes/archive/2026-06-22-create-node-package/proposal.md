## Why

The `deg-analysis` node needs to be packaged as an IRE-compliant node module. Reference DEG analysis scripts exist in `original/` but lack the standardized contract (SKILL.md frontmatter, env.yaml, NDJSON reporting, CLI dispatch, validation scripts, tests) required by the IRE agentic bioinformatics framework. This change wraps the proven reference implementation into the node package format so it can be discovered, invoked, and orchestrated by the framework.

## What Changes

- Create `SKILL.md` with v2 YAML frontmatter declaring the node's input/output contract, parameters, exceptions, and hardware requirements
- Create `envs/env-r-4.3.yaml` declaring conda dependencies for R 4.3 with Bioconductor packages (DESeq2, limma, edgeR, RCircos)
- Create `scripts/main.R` as single entry point with `run` subcommand dispatching the full DEG analysis pipeline
- Create `scripts/input_validation.R` for executable input checks (file existence, column presence, group structure)
- Create `scripts/output_validation.R` for executable output checks (expected columns, non-empty results)
- Create `tests/testthat/` test suite covering the main pipeline and validation scripts
- Wrap the 7 reference scripts (diff.R, diff_filter.R, diff_volcano.R, diff_heatmap.R, diff_venn.R, diff_locate.R, diff_check_proportion.py) as internal modules callable from `main.R`

## Capabilities

### New Capabilities

- `deg-analysis-node`: Complete IRE node package for differential expression analysis supporting DESeq2, limma, edgeR, t-test, and Wilcoxon methods. Includes sample proportion validation, result filtering, volcano plot, heatmap, Venn diagram, and chromosome location visualization.

### Modified Capabilities

None — this is a new node, not modifying existing capabilities.

## Impact

- New files at repo root: `SKILL.md`, `envs/env-r-4.3.yaml`
- New files in `scripts/`: `main.R`, `input_validation.R`, `output_validation.R`
- New files in `tests/testthat/`: test suite for main pipeline and validation
- Internal modules adapted from `original/scripts/` into `scripts/` as helpers
- No changes to existing `protocols/`, `openspec/`, or `CLAUDE.md`
- Depends on conda-forge and bioconda channels for R package availability
