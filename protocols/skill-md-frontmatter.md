# Protocol: SKILL.md Frontmatter Contract

The YAML frontmatter of SKILL.md is the machine-readable contract between the node and AI agents. Agents structurally parse the frontmatter — they do not regex body prose.

## Required Fields

```yaml
---
name: <short-kebab-case-name>
description: >
  <function + applicable data/preconditions/typical scenarios.
  Primary basis for Agent orchestration selection.>
type: standard          # standard | bridge
---
```

## Structured Fields

### Input/Output Contract

```yaml
inputs:
  - name: <file-name>
    format: <pickle|csv|tsv|h5ad|...>
    semantic_type: <type-token>
    description: <plain-language>

outputs:
  - name: <file-name>
    format: <csv|tsv|...>
    semantic_type: <type-token>
    columns: [<expected-columns>]    # optional
```

### Parameters

```yaml
entry: scripts/main.<ext>

parameters:
  - name: --<flag>           # or: subcommand (positional, no --)
    type: file|file_out|int|float|bool|string|choice
    required: true|false
    default: <value>         # static params only
    range: [<min>, <max>]    # float/int params, optional
    bind: upstream|config|static|framework
    description: <plain-language>
```

**Bind annotations:**

| `bind` | Meaning |
|--------|---------|
| `upstream` | Wired from a prior node's output |
| `config` | From the protocol's `config:` block |
| `static` | Has a default; no wiring required |
| `framework` | Set by the orchestrator (e.g., `--outdir`) |

### Exceptions

```yaml
exceptions:
  - exit_code: <N>
    pattern: "<stderr-substring>"
    nature: data_insufficient|data_corrupt|data_mismatch|env_bug
    action: halt|skip_with_warning|escalate
```

No `retry` action. Nodes must handle transient errors internally.

### Hardware

```yaml
hardware:
  memory_gb: <N>
  cpu: <N>
  gpu: true|false
  runtime: "<estimate>"
```

## Body Sections

Body sections are human/LLM narrative. Agents read them semantically for context. They do not parse them structurally.

| Section | Purpose |
|---------|---------|
| `# Node Function` | What it does, algorithm choice, intent |
| `# Expected Input` | Plain-language data description |
| `# Exceptions` | Plain-language: what can go wrong, why |
