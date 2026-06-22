# Verification Report: fix-exception-leaks

- **Date:** 2026-06-22
- **Change:** fix-exception-leaks
- **Verify Mode:** light (4 tasks, 1 file)

## Test Results
[ FAIL 0 | WARN 8 | SKIP 0 | PASS 390 ]

## Fix 1: .exceptions reset at main() entry
scripts/main.R line 521: `assign(".exceptions", list(), envir = .GlobalEnv)` added.

## Fix 2: tryCatch around DE dispatch
scripts/main.R lines 326-345: DE method switch() wrapped in tryCatch.
On error: structured NDJSON exception + .run_result.json written + graceful return.

## Final Assessment
Both fixes verified. All tests pass. No regressions.
