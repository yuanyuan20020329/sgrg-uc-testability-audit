# -*- coding: utf-8 -*-
"""Apply the v1.0.1 English language normalisation to the generated result files.

Why not just re-run the generators?
-----------------------------------
``scripts/partD_resfeas_ledger.py`` and ``scripts/partD_resfeas_libwide.py`` each emit a
figure next to their tables, and a Matplotlib PDF embeds a creation timestamp — so
re-running them would move figure bytes for reasons unrelated to this change. v1.0.1 is a
language normalisation, so the figures must not move.

(The generator *sources* are updated to English too, so a future re-run reproduces the
v1.0.1 files byte-for-byte in the fields this script touches. This script exists to bring
the already-generated files in line without re-executing the plotting code.)

How the substitution stays faithful
-----------------------------------
Rather than re-deriving each reason from the row's numbers, this script parses the
*existing* Chinese value with a regex built from the Chinese template and re-emits the
**captured substrings** into the English template. The digits in the English string are
therefore literally the digits from the Chinese string: numeric identity is structural
rather than merely asserted. If a value fails to parse against exactly one known template,
the script aborts instead of guessing.

Line endings
------------
Every file in ``results/`` was written in Windows text mode and is CRLF-terminated (the
JSON files have no trailing newline). Bytes outside the normalised values must not move,
so this script detects the convention per file, rewrites in it, and writes with
``write_bytes`` so no implicit translation can occur.

Usage
-----
    python scripts/i18n_normalise_reason.py            # apply
    python scripts/i18n_normalise_reason.py --check    # report only; exit 1 if pending
    python scripts/i18n_normalise_reason.py --from-frozen   # re-derive from the v1.0.0 payload
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FROZEN = ROOT / "resource_paper" / "frozen_deposit_v1.0.0"
sys.path.insert(0, str(Path(__file__).resolve().parent))

import resfeas_reason_i18n as I18N  # noqa: E402

LEDGER = RESULTS / "partD_resfeas_ledger.tsv"
SUMMARY = RESULTS / "partD_resfeas_libwide_summary.json"
VERDICT = RESULTS / "partD_resfeas_verdict.json"
MAP_OUT = RESULTS / "partD_resfeas_ledger_reason_i18n.tsv"
FROZEN_LEDGER = FROZEN / "partD_resfeas_ledger.tsv"

CJK = re.compile(r"[\u4e00-\u9fff]")

# --------------------------------------------------------------- %-format toolkit
_DIRECTIVE = re.compile(r"%(?:\.\d+)?[dsf]|%%")
_PATTERNS = {"d": r"(-?\d+)", "f": r"(-?\d+\.\d+|nan)", "s": r"(.+)"}


def _directives(template: str) -> list[tuple[str, str]]:
    """Decompose ``template`` into ``('lit', text)`` / ``('d'|'f'|'s'|'%', token)``."""
    out, pos = [], 0
    for m in _DIRECTIVE.finditer(template):
        if m.start() > pos:
            out.append(("lit", template[pos:m.start()]))
        tok = m.group(0)
        out.append(("%" if tok == "%%" else tok[-1], tok))
        pos = m.end()
    if pos < len(template):
        out.append(("lit", template[pos:]))
    return out


def template_to_regex(template: str) -> re.Pattern:
    """Compile any ``%``-template (Chinese or English) into an anchored matcher."""
    parts = []
    for kind, literal in _directives(template):
        if kind == "lit":
            parts.append(re.escape(literal))
        elif kind == "%":
            parts.append("%")
        else:
            parts.append(_PATTERNS[kind])
    return re.compile("^" + "".join(parts) + "$")


def split_zh(template: str, text: str) -> list[str] | None:
    """Return the captured groups, or ``None`` when the text does not match.

    Note the deliberate ``is not None`` distinction: a template with no format
    directive matches with *zero* captures, and an empty list is falsy.
    """
    m = template_to_regex(template).match(text)
    return None if m is None else list(m.groups())


def fill_en(template: str, values: list[str]) -> str:
    """Emit ``template`` with each directive replaced by the matching captured value."""
    out, it = [], iter(values)
    try:
        for kind, literal in _directives(template):
            if kind == "lit":
                out.append(literal)
            elif kind == "%":
                out.append("%")
            else:
                out.append(next(it))
    except StopIteration:
        raise SystemExit(f"[FAIL] template {template!r} expects more values than given")
    return "".join(out)


def translate_reason(text: str) -> str:
    """Translate one ``reason`` cell from the v1.0.0 Chinese to the v1.0.1 English."""
    if not CJK.search(text):
        return text
    components = text.split(I18N.REASON_SEPARATOR)
    out = []
    for comp in components:
        hits = []
        for key, pair in I18N.REASON.items():
            values = split_zh(pair["zh"], comp)
            if values is not None:
                hits.append((key, values))
        if len(hits) != 1:
            raise SystemExit(f"[FAIL] {comp!r} matched {len(hits)} templates "
                             f"({[k for k, _ in hits]}); refusing to guess")
        key, values = hits[0]
        out.append(fill_en(I18N.REASON[key]["en"], values))
    return I18N.REASON_SEPARATOR.join(out)


# ------------------------------------------------------------------- file plumbing
def read_lines(path: Path) -> tuple[list[str], str, bool]:
    """Return (logical lines, newline token, whether the file ends with one)."""
    raw = path.read_text(encoding="utf-8", newline="")
    nl = "\r\n" if "\r\n" in raw else "\n"
    trailing = raw.endswith(nl)
    body = raw[:-len(nl)] if trailing else raw
    return body.split(nl), nl, trailing


def write_lines(path: Path, lines: list[str], nl: str, trailing: bool) -> None:
    text = nl.join(lines) + (nl if trailing else "")
    path.write_bytes(text.encode("utf-8"))


# ------------------------------------------------------------------------ the work
def normalise_ledger(dry: bool, from_frozen: bool = False) -> dict:
    lines, nl, trailing = read_lines(LEDGER)
    header = lines[0].split("\t")
    i_reason = header.index("reason")
    assert header[-1] == "reason", "reason is expected to be the final column"
    assert len(lines) - 1 == 29, f"expected 29 ledger rows, found {len(lines) - 1}"

    frozen, _, _ = read_lines(FROZEN_LEDGER)
    assert len(frozen) == len(lines), "frozen and current ledger differ in length"

    # `--from-frozen` always re-derives the English from the v1.0.0 Chinese instead of from
    # whatever is on disk. Without it, a value that is already English is left alone (the
    # idempotent path); with it, an English value can be regenerated after a template is
    # reworded, which is the only way to change the v1.0.1 wording without hand-editing a
    # generated file. Either way the source of the digits stays the frozen ledger.
    src_lines = frozen if from_frozen else lines

    changed, nonreason_diffs = 0, []
    fresh = [lines[0]]
    for idx, (cur, ref, src) in enumerate(
            zip(lines[1:], frozen[1:], src_lines[1:]), start=1):
        cells, ref_cells = cur.split("\t"), ref.split("\t")
        src_cells = src.split("\t")
        if [c for i, c in enumerate(cells) if i != i_reason] != \
           [c for i, c in enumerate(ref_cells) if i != i_reason]:
            nonreason_diffs.append(idx)
        if [c for i, c in enumerate(cells) if i != i_reason] != \
           [c for i, c in enumerate(src_cells) if i != i_reason]:
            raise SystemExit(f"[FAIL] row {idx}: the file to translate disagrees with the "
                             f"frozen v1.0.0 outside `reason`; refusing to translate")
        before = cells[i_reason]
        after = translate_reason(src_cells[i_reason])
        changed += after != before
        cells[i_reason] = after
        fresh.append("\t".join(cells))

    if nonreason_diffs:
        raise SystemExit(f"[FAIL] non-reason columns differ from v1.0.0 on rows "
                         f"{nonreason_diffs}; this is not a language-only change")
    if not dry:
        write_lines(LEDGER, fresh, nl, trailing)
    return {"file": LEDGER.name, "rows": len(lines) - 1, "changed": changed}


def build_map(dry: bool) -> int:
    """Emit the v1.0.0 -> v1.0.1 row-by-row mapping, always sourced from the frozen copy."""
    frozen, _, _ = read_lines(FROZEN_LEDGER)
    current, _, _ = read_lines(LEDGER)
    head = frozen[0].split("\t")
    i_reason = head.index("reason")
    assert len(frozen) == len(current), "row count changed"

    cols = ["row_index", "source", "resource", "dataset",
            "reason_v1_0_0_zh", "reason_v1_0_1_en"]
    out = ["\t".join(cols)]
    for idx, (ref, cur) in enumerate(zip(frozen[1:], current[1:]), start=1):
        ref_c, cur_c = ref.split("\t"), cur.split("\t")
        for i, name in enumerate(head):
            if i == i_reason:
                continue
            if ref_c[i] != cur_c[i]:
                raise SystemExit(f"[FAIL] row {idx} column {name!r} differs from v1.0.0")
        rec = {"row_index": idx, "source": ref_c[0], "resource": ref_c[1],
               "dataset": ref_c[2], "reason_v1_0_0_zh": ref_c[i_reason],
               "reason_v1_0_1_en": cur_c[i_reason]}
        out.append("\t".join(str(rec[c]) for c in cols))
    if not dry:
        # Same CRLF convention as its siblings in results/
        write_lines(MAP_OUT, out, "\r\n", True)
    return len(out) - 1


def _walk(doc, keypath: tuple, filename: str):
    """Follow ``keypath`` into ``doc``, failing loudly when a key is absent."""
    node = doc
    for k in keypath:
        node = node.get(k) if isinstance(node, dict) else None
        if node is None:
            raise SystemExit(f"[FAIL] {filename}: missing {'.'.join(map(str, keypath))}")
    return node


def _normalise_json(path: Path, specs: list[tuple[tuple[str, ...], dict]], dry: bool,
                    from_frozen: bool = False) -> dict:
    raw = path.read_text(encoding="utf-8", newline="")
    nl = "\r\n" if "\r\n" in raw else "\n"
    trailing = raw.endswith(nl)
    doc = json.loads(raw)
    frozen_doc = None
    if from_frozen:
        fp = FROZEN / path.name
        if not fp.exists():
            raise SystemExit(f"[FAIL] no frozen reference to re-derive from: {fp}")
        frozen_doc = json.loads(fp.read_text(encoding="utf-8"))

    changed = 0
    for keypath, pair in specs:
        node = _walk(doc, keypath[:-1], path.name)
        leaf = keypath[-1]
        if not isinstance(node, dict) or leaf not in node:
            raise SystemExit(f"[FAIL] {path.name}: missing {'.'.join(keypath)}")
        current = node[leaf]
        if current == pair["en"]:
            continue                                     # already normalised
        if from_frozen:
            # Re-derive from the v1.0.0 payload, so a reworded template can be re-applied
            # without hand-editing a generated file. If the frozen copy does not hold the
            # Chinese the module claims, the module and the reference have drifted.
            reference = _walk(frozen_doc, keypath, f"{path.name} (frozen v1.0.0)")
            if reference != pair["zh"]:
                raise SystemExit(
                    f"[FAIL] {path.name}: {'.'.join(keypath)} in the frozen v1.0.0 payload "
                    f"is {str(reference)[:60]!r}, not the Chinese recorded in "
                    f"resfeas_reason_i18n.py; refusing to overwrite")
        elif current != pair["zh"]:
            raise SystemExit(
                f"[FAIL] {path.name}: {'.'.join(keypath)} holds an unexpected value "
                f"({str(current)[:60]!r}); refusing to overwrite")
        node[leaf] = pair["en"]
        changed += 1

    out = json.dumps(doc, indent=2, ensure_ascii=False, default=str)
    out = out.replace("\n", nl) + (nl if trailing else "")
    if not dry:
        path.write_bytes(out.encode("utf-8"))
    return {"file": path.name, "changed": changed}


def normalise_summary(dry: bool, from_frozen: bool = False) -> dict:
    specs = [(("metric_definitions", k), v) for k, v in I18N.METRIC_DEFINITIONS.items()]
    specs.append((("caveat",), I18N.SUMMARY_CAVEAT))
    return _normalise_json(SUMMARY, specs, dry, from_frozen)


def normalise_verdict(dry: bool, from_frozen: bool = False) -> dict:
    return _normalise_json(
        VERDICT, [(("geo_max_pairs_mirna_example",), I18N.VERDICT_MIRNA_EXAMPLE)], dry,
        from_frozen)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="report pending replacements and exit non-zero if any remain")
    ap.add_argument("--from-frozen", action="store_true",
                    help="re-derive every English value from the frozen v1.0.0 payload "
                         "instead of from the current (already-English) file; use after "
                         "rewording a template")
    args = ap.parse_args()

    res = [normalise_ledger(args.check, args.from_frozen),
           normalise_summary(args.check, args.from_frozen),
           normalise_verdict(args.check, args.from_frozen)]
    n_map = build_map(args.check)

    print("language normalisation (Chinese -> English, v1.0.0 -> v1.0.1)")
    print("-" * 72)
    for r in res:
        extra = f"  ({r['rows']} rows)" if "rows" in r else ""
        print(f"  {r['file']:42s} {r['changed']:4d} value(s) replaced{extra}")
    print(f"  {MAP_OUT.name:42s} {n_map:4d} row(s) mapped")

    pending = sum(r["changed"] for r in res)
    if args.check:
        print(f"\nCHECK: {pending} value(s) still in Chinese"
              + (" - nothing to do." if pending == 0 else " - run without --check."))
        return 1 if pending else 0
    print("\nDONE." if pending else "\nDONE (already normalised; no change).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
