# Protocol: CLI Invocation Contract

Every node must accept command-line arguments matching its SKILL.md frontmatter parameters. The framework constructs the CLI and captures results.

## Entry Point

```
scripts/main.<ext>          # Single entry
scripts/main.R fetch ...    # Multi-action: first positional is the subcommand
scripts/main.R qc ...       # Agent never calls any other file directly
```

## Argument Convention

- Flag parameters: `--name value`
- Positional parameters (no `--` prefix): `value` only (typically `subcommand`)
- Output directory: `--outdir <path>` (set by framework, `bind: framework`)

## Invocation

```bash
# Framework constructs from SKILL.md parameters + bound values:
Rscript scripts/main.R fetch --gse-id GSE100155 --outdir /runs/001/fetch/

# Agent never hardcodes paths. All paths come from bindings.
```

## Stdout: NDJSON

Every line written to stdout must be valid JSON:

```json
{"level": "info", "msg": "Fetching GEO data for GSE100155..."}
{"level": "result", "status": "success_matrix", "files": [...], "metadata": {...}}
```

- `level: "info"` — progress messages. May include arbitrary fields.
- `level: "result"` — final output. Must include `status`. The framework captures this as the execution result.

Gate nodes additionally include `"decision"` or `"metrics"` in the result line.

## Stderr

Used for error diagnostics and exception pattern matching. Must contain the substring declared in SKILL.md `exceptions[].pattern` for the framework to route recovery actions.

## Exit Codes

- `0` — success
- `1` — error (matched against SKILL.md exceptions for recovery)
- Other exit codes as declared in SKILL.md exceptions
