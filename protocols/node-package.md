# Protocol: Node Package Format

Every node module connected to the IRE framework must follow this package structure.

## Directory Layout

```
node-repo/
├── SKILL.md          # Agent contract (YAML frontmatter + human narrative)
├── env.yaml          # Declarative conda/mamba environment
├── env-<variant>.yaml # Optional: variant environment for CI (e.g. env-4.5.yaml)
├── scripts/
│   ├── main.<ext>             # Single entry point (.py preferred, .R secondary)
│   ├── input_validation.<ext>  # Optional: executable input checks
│   ├── output_validation.<ext> # Optional: executable output checks
│   └── ...                    # Internal modules called by main
├── references/       # Optional: static assets, reference data
├── tests/            # Test suite
├── openspec/         # Project governance
├── CLAUDE.md         # Developer guidance
└── .gitignore
```

## Requirements

1. `SKILL.md` must be present at the repo root
2. `env.yaml` must declare `name`, `channels`, `dependencies`
3. `scripts/main.<ext>` must be the single entry point
4. Multi-action nodes use subcommand dispatch: `Rscript main.R fetch --gse-id ...`
5. No hardcoded secrets in any file
6. English is the working language for all files

## Legacy Support

Repos with a `node/` subdirectory (containing `node/SKILL.md`) are supported for backward compatibility. New nodes must use flat layout.

## Connection

Each node lives in its own git repository. The core discovers it via `registry.yaml`.
