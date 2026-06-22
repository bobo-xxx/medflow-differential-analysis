## Why

The deep review identified two critical exception-handling gaps in `scripts/main.R`:
1. Unhandled Bioconductor errors during DE method dispatch crash with raw R output instead of structured NDJSON exceptions
2. The global `.exceptions` accumulator leaks state between `main()` calls in the same R session

## What Changes

- Wrap `switch()` DE method dispatch in `tryCatch` to emit structured NDJSON exception and write `.run_result.json` on failure
- Add `.exceptions` reset at `main()` entry point

## Capabilities

### New Capabilities

None — targeted bug fixes only.

### Modified Capabilities

- `differential-analysis-node`: Two critical exception-handling fixes to `scripts/main.R`

## Impact

- `scripts/main.R`: ~10 lines added (tryCatch block + accumulator reset)
- No API changes, no parameter changes, no new dependencies
