# Agent Constraints

- **CRITICAL**: You must NEVER modify, create, edit, or delete any files inside the `research/` directory. This directory is strictly read-only for agents.

## Code Commenting & Documentation Rules

- **Unit Tests (`*/tests/*`)**: Do not use Doxygen markup tags (`\f$`, `\f[`, `@brief`, `@param`, `@return`). Use standard line (`//`) or block (`/* ... */`) comments with ASCII/Markdown math.
- **Public Headers (`*.hpp`)**: Use the Dual-Tier Doxygen pattern: keep `@brief` in clean ASCII/code spans for IDE hover tooltips; put complex LaTeX formulas (`\f$ ... \f$`) in the body.
- **Explain "WHY", Not "WHAT"**: Do not paraphrase code syntax in English. Explain physical invariants, symmetries, magic factors, and tolerance derivations.

