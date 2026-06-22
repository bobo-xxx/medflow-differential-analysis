## 1. Fix Exception Leaks

- [ ] 1.1 Add `assign(".exceptions", list(), envir = .GlobalEnv)` at `main()` entry (scripts/main.R line 521)
- [ ] 1.2 Wrap DE method `switch()` dispatch in `tryCatch` with structured NDJSON exception and `.run_result.json` write on error (scripts/main.R lines 326-332)
- [ ] 1.3 Run full test suite and confirm 390+ tests pass
- [ ] 1.4 Run end-to-end pipeline with synthetic data and verify `.run_result.json` isolation between consecutive `main()` calls
