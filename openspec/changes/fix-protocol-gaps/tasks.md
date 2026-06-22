## 1. Fix Protocol Gaps

- [ ] 1.1 Create `env.yaml` symlink at repo root pointing to `envs/env-r-4.3.yaml`
- [ ] 1.2 Add `cat(paste0(code, ": ", msg), file = stderr())` to `report_exception_ndjson()` in `scripts/report.R`
- [ ] 1.3 Run full test suite and confirm 390 tests pass
