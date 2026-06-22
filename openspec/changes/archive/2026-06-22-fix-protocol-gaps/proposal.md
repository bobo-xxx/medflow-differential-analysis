## Why

Protocol compliance audit found 2 gaps: `env.yaml` missing at repo root, and exception patterns not written to stderr for framework matching.

## What Changes

1. Add `env.yaml` at repo root (symlink → `envs/env-r-4.3.yaml`)
2. Echo exception pattern to stderr in `report_exception_ndjson()` (`scripts/report.R`)

## Impact

- `env.yaml`: new symlink at root
- `scripts/report.R`: ~1 line added per exception emission
- No API changes, no parameter changes
