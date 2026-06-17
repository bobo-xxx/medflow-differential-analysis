# Standard Node Developer Guide

Build a new analysis node for the IRE agentic bioinformatics workflow framework.

## What You're Building

A node is an independently versioned, self-contained computational unit. It does one thing (DEG analysis, QC, normalization, ...) and communicates with the framework through a standardized contract. The framework — not your node — handles orchestration, data routing, and provenance.

## Quick Start

```bash
# 1. Create the repo
mkdir my-node && cd my-node && git init

# 2. Create the required files (see templates below)
touch SKILL.md env.yaml scripts/main.R

# 3. Test locally
Rscript scripts/main.R <subcommand> --<param> <value> --outdir ./output

# 4. Register in core
# Add entry to the core's registry.yaml with your repo URL, commit hash, and sha256
```

## Required Files

```
my-node/
├── SKILL.md                    # Agent contract (YAML frontmatter + body)
├── env.yaml                    # Conda/mamba environment
├── scripts/
│   ├── main.R                  # Single entry point
│   ├── input_validation.R      # Recommended: executable input checks
│   └── output_validation.R     # Recommended: executable output checks
├── tests/                      # Test suite
└── references/                 # Optional: reference data, docs
```

## Step 1: Write SKILL.md

The frontmatter is the machine contract. The body is human narrative. Start with this template:

```yaml
---
name: my-analysis-node
description: >
  <One sentence: what it does + when to use it.
  Primary basis for agent orchestration selection.>
type: standard

inputs:
  - name: input_data.csv
    format: csv
    semantic_type: expression_matrix
    description: "Normalized expression matrix (genes × samples)"

outputs:
  - name: results.csv
    format: csv
    semantic_type: deg_results
    columns: [gene, log2FC, pvalue, fdr]

entry: scripts/main.R

parameters:
  - name: subcommand
    type: choice
    choices: [analyze, qc]
    required: true
    bind: config

  - name: --outdir
    type: file_out
    required: true
    bind: framework

  - name: --pvalue-threshold
    type: float
    default: 0.05
    range: [0.001, 0.25]
    bind: static

exceptions:
  - exit_code: 1
    pattern: "insufficient samples"
    nature: data_insufficient
    action: skip_with_warning

  - exit_code: 1
    pattern: "input file corrupt"
    nature: data_corrupt
    action: escalate

hardware:
  memory_gb: 4
  cpu: 2
  gpu: false
  runtime: "~5 min"
---

# Node Function

<Describe what your node does, what algorithm it uses, and its intent.>

# Expected Input

<Plain-language description of the input data. What it means, what format,
what conditions make it valid.>

# Exceptions

<Plain-language: what can go wrong, why, and what the user should check.>
```

### Key Frontmatter Rules

- `semantic_type` is the type system. Pick a descriptive token and use it consistently.
- `bind: upstream` for parameters wired from prior nodes. `bind: config` for values from the protocol config. `bind: static` for tuning knobs with defaults. `bind: framework` for `--outdir`.
- `exceptions` declare only what you CANNOT handle internally. No `retry` — handle transient errors yourself.
- All exception patterns, messages, and documentation must be in English.

## Step 2: Write env.yaml

```yaml
name: my-node
channels:
  - conda-forge
  - bioconda          # only if using R/Bioconductor
dependencies:
  - python=3.13       # or r-base=4.3
  - pandas
  - pip:
    - <pip-only-package>
```

## Step 3: Write scripts/main.R (or main.py)

```r
#!/usr/bin/env Rscript

# Single entry point. First positional arg is the subcommand.
args <- commandArgs(trailingOnly = TRUE)
subcommand <- args[1]
param_args <- args[-1]

# Parse named parameters
parse_arg <- function(args, name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

outdir <- parse_arg(param_args, "--outdir", ".")
pval_threshold <- as.numeric(parse_arg(param_args, "--pvalue-threshold", "0.05"))

# Dispatch
if (subcommand == "analyze") {
  do_analysis(outdir, pval_threshold)
} else if (subcommand == "qc") {
  do_qc(outdir)
} else {
  cat(jsonlite::toJSON(list(level="error",
    msg=paste0("Unknown subcommand: ", subcommand, ". Valid: analyze, qc"))), "\n")
  quit(status = 1)
}
```

### Script Rules

1. **Single entry**: one `main.R` (or `main.py`). Multi-action via subcommand.
2. **Write to `--outdir`**: all outputs go there. Never use relative paths.
3. **NDJSON to stdout**: `{"level":"info","msg":"..."}` for progress. `{"level":"result","status":"...","files":[...]}` as final output.
4. **Errors to stderr**: use `message()` in R or `print(..., file=sys.stderr)` in Python.
5. **Exit 0 on success, non-zero on failure**: match your SKILL.md exception exit codes.
6. **Handle network/transient errors internally**: retry with backoff. Don't delegate to the framework.

## Step 4: Write validation scripts (recommended)

### input_validation.R

```r
# Validate that input files exist and match expectations.
# Exit 0 if valid, non-zero + stderr message if invalid.
```

### output_validation.R

```r
# Validate that output files have expected columns, non-zero rows, etc.
# Exit 0 if valid, non-zero + stderr message if invalid.
```

## Step 5: Write tests

```bash
tests/
└── testthat/
    ├── test-main.R
    ├── test-analyze.R
    └── test-qc.R
```

Run: `Rscript -e 'testthat::test_dir("tests/testthat/")'`

## Step 6: Test Locally with the Core

```bash
# In the core repo:
# 1. Place your node package at nodes/my-node@1.0.0/
cp -r ~/my-node nodes/my-node@1.0.0/

# 2. Verify the core discovers it
npx tsx src/cli.ts list

# 3. Create a test protocol .md that references your node
# 4. Run through the full pipeline
npx tsx test-run.ts <gse-id> runs/test-output
```

## Step 7: Publish

1. Push your node repo to GitHub
2. Tag the commit: `git tag v1.0.0 && git push --tags`
3. Compute the package sha256: `find . -type f ! -path './.git/*' -exec cat {} + | sha256sum`
4. Add an entry to the core's `registry.yaml`:

```yaml
  my-node:
    source: git
    url: https://github.com/<user>/my-node.git
    versions:
      - version: "1.0.0"
        commit: "<git-sha>"
        sha256: "<sha256>"
```

5. Submit a PR to the core repo with the registry update.

## Protocol Reference

Your node must satisfy these protocols:

| Protocol | File | What it covers |
|----------|------|---------------|
| Package Format | `protocols/node-package.md` | Directory structure, flat layout |
| SKILL.md Contract | `protocols/skill-md-frontmatter.md` | v2 frontmatter fields |
| CLI Contract | `protocols/cli-contract.md` | Invocation, NDJSON, exit codes |
| Exception Contract | `protocols/exception-contract.md` | Error patterns, halt/skip/escalate |
| Core Connection | `protocols/core-connection.md` | registry.yaml, ire sync |

## Quick Checklist

- [ ] Flat layout at repo root
- [ ] `SKILL.md` with all 8 frontmatter fields (name, description, type, inputs, outputs, parameters, exceptions, hardware)
- [ ] `env.yaml` with `name`, `channels`, `dependencies`
- [ ] `scripts/main.<ext>` as single entry point
- [ ] NDJSON stdout: `info` + `result` lines
- [ ] `--outdir` parameter with `bind: framework`
- [ ] `exceptions` with `pattern`, `nature`, `action` (no `retry`)
- [ ] No hardcoded secrets
- [ ] All content in English
- [ ] Test suite passing
- [ ] `registry.yaml` entry ready for the core
