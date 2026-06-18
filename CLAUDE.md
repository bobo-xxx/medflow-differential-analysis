# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with this node package.

## Project Identity

This is a node package in the IRE agentic bioinformatics workflow framework. It is a standalone git repository — one repo per node.

**Node:** `deg-analysis`  
**Purpose:** Differential expression analysis — case vs control comparison.  
**Type:** standard node

## Language

**English is the working language.** SKILL.md, env.yaml, code comments, commit messages, error messages, NDJSON reporting, and OpenSpec artifacts are all in English.

## Development Environment

This node has its own conda environment, declared in `envs/env-r-4.3.yaml`.

```bash
conda env create -f envs/env-r-4.3.yaml
conda activate <env-name>

# Python: use uv inside conda
conda install -c conda-forge uv && uv pip install <package>

# R: use conda-forge R
conda install -c conda-forge r-base r-essentials
# R packages: prefer conda-forge, fall back to install.packages()
```

### Rules

1. Use the node's own env from `envs/env-r-4.3.yaml`
2. Python packages via `uv pip install`
3. R packages prefer `conda install -c conda-forge r-<package>`
4. No hardcoded secrets
5. Preferred language: Python > R > shell

## Node Package Reference

The `protocols/` directory contains the authoritative protocol documents:

| Protocol | File | Covers |
|----------|------|--------|
| Package Format | `protocols/node-package.md` | Flat layout, directory structure |
| SKILL.md Contract | `protocols/skill-md-frontmatter.md` | v2 frontmatter fields |
| CLI Contract | `protocols/cli-contract.md` | Invocation, NDJSON, exit codes |
| Exception Contract | `protocols/exception-contract.md` | Error patterns, actions |
| Core Connection | `protocols/core-connection.md` | registry.yaml linkage |
| Developer Guide | `protocols/node-developer-guide.md` | 7-step walkthrough |

## TODO

- [ ] Create `SKILL.md` following the v2 frontmatter protocol
- [ ] Create `env.yaml` with required dependencies
- [ ] Create `scripts/main.<ext>` as single entry point
- [ ] Create `scripts/input_validation.<ext>` and `scripts/output_validation.<ext>`
- [ ] Write tests
- [ ] Register in core's `registry.yaml`
