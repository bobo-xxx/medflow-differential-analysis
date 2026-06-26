# Protocol: License Requirement

All node packages distributed under MedFlow MUST use the Apache License 2.0.

## Requirements

### Apache 2.0 License File

Every node repository MUST contain a `LICENSE` file at the root with the full Apache License 2.0 text.

- Source: https://www.apache.org/licenses/LICENSE-2.0.txt
- Copyright line: `2026 MedFlow Contributors`

### Verification

After creating or pushing the repository, verify:

```bash
head -3 LICENSE | grep -q "Apache License" || ( \
  curl -sL https://www.apache.org/licenses/LICENSE-2.0.txt -o LICENSE && \
  sed -i 's/\[yyyy\] \[name of copyright owner\]/2026 MedFlow Contributors/' LICENSE && \
  git add LICENSE && git commit --amend -m "chore: add Apache 2.0 LICENSE" )
```

### README Consistency

The README.md `## License` section MUST state "Apache 2.0" (NOT MIT).

### Comet Verify Checklist

Before archive, these checks MUST pass:

```
□ LICENSE: head -3 LICENSE | grep "Apache License" returns match
□ LICENSE: Copyright is "2026 MedFlow Contributors" (NOT xxx, NOT placeholder)
□ README: License section says "Apache 2.0" (NOT MIT)
```
