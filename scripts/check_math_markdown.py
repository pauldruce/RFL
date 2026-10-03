#!/usr/bin/env python3
r"""Validates mathematical notation in Markdown files for GitHub and VS Code rendering compatibility.

Enforces rules defined in docs/theory/README.md:
- R1: Standalone display math blocks ($$...$$) on own lines; no ```math code fences.
- R2 / C1: Inline math spans ($...$) must not contain punctuation before subscripts ([})\]|]_)
           when 2 or more underscores exist across the paragraph (preventing CommonMark <em> collisions).
- R3 / C2: No inline math ending in ')' directly followed by ')' (prevents ')$)' delimiter bug).
- R4: No inline math spans inside bold delimiters (**...$...**).
- R5: Headings must not contain multi-subscript math violating R2.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Regex to match Markdown bold spans: **content**
BOLD_SPAN_PATTERN = re.compile(r"(?<!\*)\*\*([^*\n]+?)\*\*(?!\*)")

# Regex to detect inline math ending in ')' followed by ')'
PAREN_DELIM_PATTERN = re.compile(r"\)\$\)")

# Regex to detect ```math fence
MATH_FENCE_PATTERN = re.compile(r"^```math\b", re.MULTILINE)

# Regex to match inline math spans $...$ (not $$...$$ and not escaped \$)
INLINE_MATH_PATTERN = re.compile(r"(?<!\\)(?<!\$)\$(?!\$)(.+?)(?<!\\)(?<!\$)\$")

# Regex to detect subscripts after closing delimiters or punctuation: }, ), ], |
DELIMITER_UNDERSCORE_PATTERN = re.compile(r"[})\]|]_")


def check_markdown_math(file_path: Path) -> list[str]:
    """Validates math notation in a single Markdown file."""
    errors: list[str] = []
    try:
        content = file_path.read_text(encoding="utf-8")
    except Exception as exc:
        return [f"{file_path}: Failed to read file: {exc}"]

    lines = content.splitlines()

    # Rule R1: check for ```math code fence
    for line_no, line in enumerate(lines, start=1):
        if line.strip().startswith("```math"):
            errors.append(
                f"{file_path}:{line_no}: [R1] Use standalone '$$' display blocks instead of '```math' code fences."
            )

    # State tracking: ignore code blocks and display math blocks
    in_code_block = False
    in_display_math = False

    # Group lines into paragraphs with starting line numbers
    paragraphs: list[tuple[int, list[str]]] = []
    current_para: list[str] = []
    para_start_line = 1

    for line_no, line in enumerate(lines, start=1):
        stripped = line.strip()

        # Handle fenced code blocks (``` or ~~~)
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_code_block = not in_code_block
            continue

        if in_code_block:
            continue

        # Handle display math blocks ($$)
        if stripped.startswith("$$"):
            if not in_display_math:
                if stripped.endswith("$$") and len(stripped) > 2:
                    # Single-line display math
                    pass
                else:
                    in_display_math = True
            else:
                in_display_math = False
            continue

        if in_display_math:
            if r"\{" in line or r"\}" in line:
                errors.append(
                    f"{file_path}:{line_no}: [R6] Escaped brace in display math: "
                    "Markdown consumes the backslash before math renders. Use '\\lbrace' and '\\rbrace' instead."
                )
            continue

        # Headings are treated as independent single-line paragraphs
        if stripped.startswith("#"):
            if current_para:
                paragraphs.append((para_start_line, current_para))
                current_para = []
            paragraphs.append((line_no, [line]))
            para_start_line = line_no + 1
            continue

        # Blank line separates paragraphs
        if not stripped:
            if current_para:
                paragraphs.append((para_start_line, current_para))
                current_para = []
            para_start_line = line_no + 1
        else:
            if not current_para:
                para_start_line = line_no
            current_para.append(line)

    if current_para:
        paragraphs.append((para_start_line, current_para))

    # Validate each paragraph
    for start_line, para_lines in paragraphs:
        para_text = " ".join(para_lines)

        # Rule R4: Check for math inside bold asterisks
        for match in BOLD_SPAN_PATTERN.finditer(para_text):
            if "$" in match.group(1):
                errors.append(
                    f"{file_path}:{start_line}: [R4] Inline math inside bold asterisks '**...$math$...**': "
                    f"'{match.group(0)}'. Move math outside asterisks (e.g. '**Label** ($O(z)$):')."
                )

        # Rule R3 / C2: Check for ')$)' collision
        for match in PAREN_DELIM_PATTERN.finditer(para_text):
            # Find context around match
            snippet_start = max(0, match.start() - 15)
            snippet_end = min(len(para_text), match.end() + 15)
            snippet = para_text[snippet_start:snippet_end]
            errors.append(
                f"{file_path}:{start_line}: [R3/C2] Inline math ending in ')' directly followed by ')' ')$)': "
                f"'{snippet}'. Reword to separate delimiters (e.g. 'order $O(z)$:')."
            )

        # Rule R2 / C1: Check for punctuation-underscores when paragraph has >= 2 underscores
        # Remove inline code spans (`...`) first to avoid false positives
        clean_para = re.sub(r"`[^`]+`", "", para_text)
        total_underscores = clean_para.count("_")

        if total_underscores >= 2:
            math_spans = INLINE_MATH_PATTERN.findall(clean_para)
            for span in math_spans:
                delim_match = DELIMITER_UNDERSCORE_PATTERN.search(span)
                if delim_match:
                    errors.append(
                        f"{file_path}:{start_line}: [R2/C1] Delimiter followed by underscore in inline math "
                        f"with multiple underscores in paragraph: '${span}$'. "
                        "Move equation to a '$$' display block or drop subscript in prose."
                    )
                    break

        # Rule R6: Check for escaped braces (\{ or \}) in inline math
        for span in INLINE_MATH_PATTERN.findall(clean_para):
            if r"\{" in span or r"\}" in span:
                errors.append(
                    f"{file_path}:{start_line}: [R6] Escaped brace in inline math '${span}$': "
                    "Markdown consumes the backslash before math renders. Use '\\lbrace' and '\\rbrace' instead."
                )

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify math notation compatibility in Markdown files."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        default=[REPO_ROOT / "docs"],
        help="Files or directories to scan (defaults to docs/)",
    )
    args = parser.parse_args()

    files_to_check: list[Path] = []
    for path in args.paths:
        if path.is_file() and path.suffix == ".md":
            files_to_check.append(path)
        elif path.is_dir():
            files_to_check.extend(path.rglob("*.md"))

    files_to_check.sort()
    all_errors: list[str] = []

    for file_path in files_to_check:
        errs = check_markdown_math(file_path)
        all_errors.extend(errs)

    if all_errors:
        print(f"❌ Found {len(all_errors)} math formatting violation(s):", file=sys.stderr)
        for err in all_errors:
            print(f"  • {err}", file=sys.stderr)
        return 1

    print(f"✅ Verified math formatting across {len(files_to_check)} Markdown file(s). Zero violations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
