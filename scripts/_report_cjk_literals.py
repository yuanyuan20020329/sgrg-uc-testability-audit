# -*- coding: utf-8 -*-
"""Review aid: report non-docstring CJK string literals in the generators.

Shipped with the code archive (a reader checking the v1.0.1 language claim may want to
re-run it), but it is a review aid, not a gate: it inspects *source* text and cannot tell a
shipped payload value from a console message. The real gate is
``scripts/audit_data_descriptor.py``, which inspects the payload itself, and
``scripts/test_language_normalisation.py``, which diffs it against the frozen v1.0.0 bytes.
"""
from __future__ import annotations

import ast
import pathlib
import sys


def docstrings(tree: ast.AST) -> set[str]:
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            d = ast.get_docstring(node, clean=False)
            if d:
                out.add(d)
    return out


def main(paths: list[str]) -> int:
    total = 0
    for f in paths:
        src = pathlib.Path(f).read_text(encoding="utf-8")
        tree = ast.parse(src)
        docs = docstrings(tree)
        bad = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                if node.value in docs:
                    continue
                if any("\u4e00" <= c <= "\u9fff" for c in node.value):
                    bad.append((node.lineno, node.value[:70]))
        print(f"{f}: non-docstring CJK literals = {len(bad)}")
        for ln, v in bad:
            print(f"    line {ln}: {v!r}")
        total += len(bad)
    print(f"\nTOTAL non-docstring CJK literals: {total}")
    return 0 if total else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
