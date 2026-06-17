# Protocol: Exception Handling Contract

Nodes declare known failure modes in SKILL.md frontmatter. The framework matches stderr against declared patterns and routes recovery actions.

## Exception Declaration

```yaml
exceptions:
  - exit_code: 1
    pattern: "样本数不足"           # substring matched against stderr
    nature: data_insufficient      # data_insufficient | data_corrupt | data_mismatch | env_bug
    action: skip_with_warning      # halt | skip_with_warning | escalate
```

## Actions

| Action | Framework Response |
|--------|-------------------|
| `halt` | Stop the pipeline. Node cannot produce output. |
| `skip_with_warning` | Mark node as skipped. Record reason. Continue if downstream accepts partial results. |
| `escalate` | Pause the run. Present stderr to human. Wait for decision. |

## No Retry

Nodes must handle transient errors internally (network retries, algorithm fallbacks). If a transient error reaches the framework, the node's script is incomplete. There is no `retry` action.

## How the Framework Matches

```
1. Node exits non-zero
2. Framework captures stderr
3. For each exception in SKILL.md with matching exit_code:
   - If 'pattern' appears as substring in stderr → follow declared action
4. If no exception matches → escalate to human with full stderr
```

## Best Practices

1. Include the error context in stderr so pattern matching works
2. Keep patterns specific enough to not false-match unrelated errors
3. Use distinct exit codes for different error classes
4. Handle everything you can internally — exceptions are for delegating to the agent
