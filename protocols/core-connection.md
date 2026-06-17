# Protocol: Core Connection via registry.yaml

Every node module lives in its own git repository. The IRE core framework discovers and links to nodes through a declarative `registry.yaml`.

## Registry Entry

The core's `registry.yaml` maps node names to their git repositories:

```yaml
nodes:
  geo-microarray-processing:
    source: git
    url: https://github.com/bobo-xxx/medflow-geo-microarray.git
    versions:
      - version: "1.0.0"
        commit: "<full-git-sha>"
        sha256: "<sha256-of-package-files>"
```

## How ire sync Works

```
ire sync
  → reads registry.yaml
  → for each node, checks if requested version exists in nodes/
  → if missing: git clone <url> into nodes/<name>@<version>/
  → checkout pinned commit
  → verify sha256 integrity
  → verify package root (flat layout; legacy node/ subdirectory supported)
  → done — node is available for execution
```

## Package Root Resolution

1. Check repo root for `SKILL.md` (flat layout — default)
2. If absent, check `node/SKILL.md` (legacy `node/` subdirectory)
3. Extract the package root to `nodes/<name>@<version>/`

## What the Node Repo Must Contain

```
SKILL.md          # Agent contract
env.yaml          # Conda environment
scripts/          # Entry point + internal modules
references/       # Optional static assets
```

## What the Node Repo Must NOT Contain

- `registry.yaml` — belongs to the core only
- Framework engine code
- Other node packages

## Publishing a New Version

1. Commit and push changes to the node repo
2. Tag the commit
3. Compute sha256 of package files
4. Add a new entry under `versions:` in the core's `registry.yaml`
5. The core maintainer verifies and merges

## Local Testing

Place the package directory at `nodes/<name>@<version>/` in a local core checkout. The `NodeRegistry.discover()` method indexes it without `ire sync`.
