## Context

Two protocol gaps identified in compliance audit against `protocols/node-package.md` and `protocols/exception-contract.md`.

## Fixes

### 1. env.yaml at root

Protocol requires `env.yaml` at repo root. Current file at `envs/env-r-4.3.yaml`.
Solution: create symlink `env.yaml → envs/env-r-4.3.yaml`.

### 2. Exception patterns to stderr

Exception contract says framework matches `exceptions[].pattern` against stderr.
Current implementation emits NDJSON to stdout only.
Solution: in `report_exception_ndjson()`, also `cat(paste0(code, ": ", msg), file = stderr())`.
