# Agent Constraints

- **CRITICAL**: You must NEVER modify, create, edit, or delete any files inside the `research/` directory. This directory is strictly read-only for agents.

## Mise Task & Automation Standards

- **Strict Colon-Only Naming**: All `mise` tasks must strictly use colon-separated notation (`<namespace>:<action>`, e.g. `pre:commit`, `pre:push`, `ci:plan`, `ci:verify:pr`, `ci:test:cpp`, `test:affected`). Never use hyphens (`-`) or mixed notation in task names.
- **Contract Encapsulation**: All CI workflow steps (`.github/workflows/*.yml`) and local Git hooks (`.githooks/*`) must be encapsulated behind `mise` targets (`mise run <target>`). Never call raw Python scripts (`scripts/*.py`), CMake, or CTest directly in workflows or git hooks.

## Agent Workflow & Concurrency Rules

- **Non-Blocking Main Thread**: Never block or poll the main agent thread waiting for long-running processes (CTest suites, full builds, remote GitHub Actions CI runs). Always delegate monitoring of long-running CI or background tasks to a lightweight subagent (model tier `flash_lite` or `flash`) or use reactive wakeup.

## Code Commenting & Documentation Rules

- **Unit Tests (`*/tests/*`)**: Do not use Doxygen markup tags (`\f$`, `\f[`, `@brief`, `@param`, `@return`). Use standard line (`//`) or block (`/* ... */`) comments with ASCII/Markdown math.
- **Public Headers (`*.hpp`)**: Use the Dual-Tier Doxygen pattern: keep `@brief` in clean ASCII/code spans for IDE hover tooltips; put complex LaTeX formulas (`\f$ ... \f$`) in the body.
- **Explain "WHY", Not "WHAT"**: Do not paraphrase code syntax in English. Explain physical invariants, symmetries, magic factors, and tolerance derivations.
