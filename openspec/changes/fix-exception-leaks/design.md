## Context

Two critical findings from comprehensive code review of `scripts/main.R`:
1. **DE dispatch unprotected** (line 326): `switch()` calls to `diff_deseq2()`, `diff_limma()`, `diff_edger()`, `diff_stat()` are not wrapped in error handling. Bioconductor errors (singular design matrix, convergence failure) produce raw R stack traces instead of structured NDJSON exceptions.
2. **Exception state leakage** (line 520): `main()` reuses the global `.exceptions` list without clearing it. Consecutive `main()` calls in the same R session accumulate exceptions across runs.

## Goals / Non-Goals

**Goals:**
- Wrap DE dispatch in `tryCatch` → emit structured NDJSON exception + write `.run_result.json`
- Reset `.exceptions` at `main()` entry

**Non-Goals:**
- No new exception codes (reuse existing patterns)
- No parameter changes
- No new tests beyond verifying the fix

## Decisions

### Fix 1: tryCatch around DE dispatch

Wrap lines 326-332 in `tryCatch`:
```r
dif <- tryCatch(
  switch(opts$method, ...),
  error = function(e) {
    report_exception_ndjson("B5_METHOD_MISMATCH", "data_mismatch", "halt",
      sprintf("DE analysis failed: %s", e$message), exit_code = 1)
    return(NULL)
  }
)
if (is.null(dif)) {
  write_run_result(opts$outdir, list(status = "error", msg = "DE analysis failed"), opts, 1, c(started_at, format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")))
  return(invisible(list(status = "error", msg = "DE analysis failed")))
}
```

### Fix 2: Clear .exceptions at main() entry

Add at line 521 (first line of `main()`):
```r
assign(".exceptions", list(), envir = .GlobalEnv)
```

## Risks / Trade-offs

- None. Both fixes are defensive and backward-compatible.
