## MODIFIED Requirements

### Requirement: Node handles errors with structured exceptions

The node SHALL wrap DE method dispatch in error handling that emits structured NDJSON exceptions on unexpected Bioconductor failures, and SHALL write `.run_result.json` even on DE analysis failure paths.

#### Scenario: DE method throws unexpected error

- **WHEN** a DE method (DESeq2, limma, edgeR, t-test, Wilcoxon) throws an unexpected error during execution
- **THEN** the script SHALL emit a structured NDJSON exception with level `"exception"` and write `.run_result.json` with `status: "error"`
- **AND** exit with code 1

#### Scenario: Consecutive main() calls produce independent results

- **WHEN** `main()` is called twice in the same R session
- **THEN** the second call's `.exceptions` accumulator SHALL contain only exceptions from the second run
- **AND** the second call's `.run_result.json` SHALL not contain exceptions from the first run
