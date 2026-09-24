#!/usr/bin/env python3
"""Check the public formula, current defaults and local documentation links."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from urllib.parse import unquote

import yaml

from pace_livestock.config import DEFAULTS

RETIRED_TERMS = (
    "U(G)",
    "U_G",
    "unassigned_mass",
    "missing_geometric",
    "competition_power",
    "PACE.Score",
    "首次实现限定0或1",
    "0 or 1 only",
)
RETIRED_TREES = ("legacy", "workflow", "example_quantified", "config", "example", "reference")


def heading_anchors(text: str) -> set[str]:
    """Return GitHub-style anchors for the ATX headings used in our docs."""
    anchors = set()
    fenced = False
    for line in text.splitlines():
        if line.startswith("```"):
            fenced = not fenced
        if fenced or not re.match(r"^#{1,6} ", line):
            continue
        heading = re.sub(r"^#+ | +#+$", "", line).strip().lower()
        base = re.sub(r"[^\w -]", "", heading).replace(" ", "-")
        anchor, suffix = base, 0
        while anchor in anchors:
            suffix += 1
            anchor = f"{base}-{suffix}"
        anchors.add(anchor)
    return anchors


def check_markdown(path: Path, root: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    relative = path.relative_to(root)
    issues = []
    for term in RETIRED_TERMS:
        if term in text:
            issues.append(f"{relative}: retired model term {term!r}")
    if re.search(r"(?i)(?:default|默认)[^\n]{0,80}(?:\\eta|eta)\s*[=$ ]+1(?!\d)", text):
        issues.append(f"{relative}: incompatible default allocation exponent")
    if "\\[" in text or "\\]" in text:
        issues.append(f"{relative}: use GitHub-supported fenced math blocks")
    if "$$" in text:
        issues.append(f"{relative}: use fenced math blocks to preserve LaTeX escapes on GitHub")
    if re.search(r"\\operatorname\b", text):
        issues.append(f"{relative}: unsupported GitHub math macro operatorname; use mathrm")
    for match in re.finditer(r"!?\[[^\]\n]*\]\(([^)\n]+)\)", text):
        target = match.group(1).strip().strip("<>")
        if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*:", target):
            continue
        filename, _, fragment = unquote(target).partition("#")
        destination = path.parent / filename if filename else path
        if not destination.exists():
            issues.append(f"{relative}: broken local link {target}")
        elif fragment and destination.suffix == ".md":
            if fragment not in heading_anchors(destination.read_text(encoding="utf-8")):
                issues.append(f"{relative}: missing heading in local link {target}")
    return issues


def audit(root: Path) -> list[str]:
    docs = sorted([*root.glob("*.md"), *(root / "docs").rglob("*.md")])
    issues = [problem for path in docs for problem in check_markdown(path, root)]
    formula = (root / "README.md").read_text(encoding="utf-8")
    total = re.search(r"```math\n\\boxed\{.*?\n```", formula, flags=re.S)
    if total is None:
        issues.append("README.md: main equation missing")
    else:
        expected = total.group(0)
        if "\\sum_{e\\in\\mathcal E(G)}" not in expected:
            issues.append("README.md: full-candidate denominator missing")
        if "\\eta" in expected or "B(E,G)" in expected:
            issues.append("README.md: allocation belongs in the experimental extension")
    reference = (root / "docs/ADVANCED.md").read_text(encoding="utf-8")
    block = re.search(
        r"<!-- configuration-defaults:start -->\s*```yaml\n(.*?)\n```"
        r"\s*<!-- configuration-defaults:end -->",
        reference,
        flags=re.S,
    )
    if block is None or yaml.safe_load(block.group(1)) != DEFAULTS:
        issues.append("ADVANCED.md: documented defaults differ from installed configuration")
    for name in RETIRED_TREES:
        tree = root / name
        # Ignore local untracked caches/results left by earlier checkouts.
        if tree.exists() and any(
            p.suffix in (".py", ".md", ".smk", ".yaml", ".yml")
            for p in tree.rglob("*")
            if "results" not in p.parts
        ):
            issues.append(f"{name}: retired code or documentation is present in current tree")
    return issues


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    issues = audit(root)
    if issues:
        print("\n".join(issues), file=sys.stderr)
        return 1
    print("Public formula, math macros, defaults, local links and retired-tree checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
