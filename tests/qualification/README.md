# Sealed Scientific Qualification Acceptance Tests

This directory contains sealed acceptance test suites that verify mathematical and physical invariants.
The system protects these tests with cryptographic hashes under [EP-5](../../docs/eps/ep-5-specification-driven-agentic-workflow.md).

## Integrity Protection (`.manifest`)

The `.manifest` file records the SHA-256 hash of each qualification test file.
The test harness normalises all line endings to LF before hashing.
This normalisation ensures deterministic verification on Linux, macOS, and Windows.

```bash
# Verify integrity of sealed qualification tests
python3 scripts/verify_traceability.py --verify-manifest tests/qualification/.manifest
```

## Governance & Modification Procedure

Autonomous AI agents and contributors must not loosen numerical tolerances or alter invariants without human approval.

1. **Review Required:** Every change in this directory requires review from maintainers listed in `.github/CODEOWNERS`.
2. **Pull Request Label:** The pull request must include the `qualification-override` label.
3. **Regenerate Manifest:** After intentional updates, regenerate `.manifest`:
   ```bash
   python3 scripts/verify_traceability.py --generate-manifest tests/qualification
   ```
