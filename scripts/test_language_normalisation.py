# -*- coding: utf-8 -*-
"""Regression test: v1.0.1 must be a *language-only* change against v1.0.0.

v1.0.1 exists for one reason — the free-text values shipped in v1.0.0 were written in
Chinese and the deposit is cited by an English-language data descriptor. That is the whole
scope. This test is the machine-checkable statement of it.

What it proves, file by file:

* the four tabular records keep their column names, column order, row order, row count,
  every numeric value, every boolean and the four gates plus ``USABLE``;
* the only cells that differ are the ones a language normalisation may touch;
* the frozen source snapshot and the three provenance-free tables are byte-identical;
* ``datapackage.json`` keeps its resource list, field order, types and constraints —
  only the version, the two prose descriptions and the example are allowed to move;
* no Han character is left anywhere in the payload.

Run it with no arguments to compare the current ``resource_paper/deposit/`` against the
frozen v1.0.0 reference; pass ``--old``/``--new`` to compare other pairs.

    python scripts/test_language_normalisation.py
    python scripts/test_language_normalisation.py --json-out results/language_normalisation.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

import resfeas_reason_i18n as I18N  # noqa: E402
import i18n_normalise_reason as NORM  # noqa: E402

OLD_DEFAULT = ROOT / "resource_paper" / "frozen_deposit_v1.0.0"
NEW_DEFAULT = ROOT / "resource_paper" / "deposit"

CJK = re.compile(r"[\u4e00-\u9fff]")

# Files that a language normalisation must not touch at all.
FROZEN_FILES = {
    "cxg_coll_detail_snapshot_v1.0.0.zip",
    "partD_resfeas_geo_series.tsv",
    "partD_resfeas_libwide_collections.tsv",
    "partD_resfeas_libwide_datasets.tsv",
}

# JSON files whose changed key paths must stay inside an allow-list.
JSON_ALLOWED = {
    "partD_resfeas_libwide_summary.json": (
        {("metric_definitions", k) for k in I18N.METRIC_DEFINITIONS} | {("caveat",)}
    ),
    "partD_resfeas_verdict.json": {("geo_max_pairs_mirna_example",)},
    "datapackage.json": {
        ("version",), ("created",),
        ("resources", 3, "description"),
        ("resources", 2, "description"),
        ("resources", 3, "schema", "fields", 15, "description"),
    },
}

GATE_COLS = ["G1_mRNA", "G2_UC", "G3_powered_control", "G4_in_vivo", "USABLE"]


# --------------------------------------------------------------------------- helpers
def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def lines_of(p: Path) -> tuple[list[str], str, bool]:
    raw = p.read_text(encoding="utf-8", newline="")
    nl = "\r\n" if "\r\n" in raw else "\n"
    trailing = raw.endswith(nl)
    body = raw[:-len(nl)] if trailing else raw
    return body.split(nl), nl, trailing


def leaves(node, path=()):
    """Flatten a JSON document to ``{path_tuple: value}``."""
    if isinstance(node, dict):
        for k, v in node.items():
            yield from leaves(v, path + (k,))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from leaves(v, path + (i,))
    else:
        yield path, node


def allowed(name: str, path: tuple) -> bool:
    """A path is allowed to change if an allow-listed prefix for *that file* matches."""
    for pre in JSON_ALLOWED.get(name, set()):
        if path[:len(pre)] == pre:
            return True
    return False


class Report:
    def __init__(self) -> None:
        self.rows: list[tuple[bool, str, str]] = []

    def check(self, ok: bool, label: str, detail: str = "") -> bool:
        self.rows.append((bool(ok), label, detail))
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {label}" + (f"  -- {detail}" if detail and not ok else ""))
        return bool(ok)

    @property
    def failed(self) -> int:
        return sum(1 for ok, _, _ in self.rows if not ok)


# ------------------------------------------------------------------------------ test
def run(old: Path, new: Path, rep: Report) -> None:
    # ---- A0. module contract: the templates must be safe to join and to translate
    sep = I18N.REASON_SEPARATOR
    dirty = [k for k, v in I18N.REASON.items() if sep in v["en"]]
    rep.check(not dirty, f"no English reason template contains the join token {sep!r}",
              f"templates that would be ambiguous to split: {dirty}")
    mismatch = []
    for key, pair in I18N.REASON.items():
        n_zh = len(NORM._directives(pair["zh"]))
        n_en = len(NORM._directives(pair["en"]))
        if n_zh != n_en:
            mismatch.append(f"{key}: {n_zh} vs {n_en} directives")
    rep.check(not mismatch, "every reason template keeps its zh/en directive counts",
              f"{mismatch}")

    # ---- A. inventory
    rep.check(old.is_dir(), f"reference payload exists at {old}")
    rep.check(new.is_dir(), f"current payload exists at {new}")
    old_files = sorted(p.name for p in old.iterdir() if p.name != "_FROZEN_REFERENCE.json")
    new_files = sorted(p.name for p in new.iterdir())
    rep.check(old_files == new_files, "file inventory unchanged",
              f"only in old: {set(old_files) - set(new_files)}; "
              f"only in new: {set(new_files) - set(old_files)}")
    rep.check(len(new_files) == 14, f"payload has 14 files (found {len(new_files)})")

    # ---- B. tables that must not move
    for name in sorted(FROZEN_FILES):
        if (old / name).exists() and (new / name).exists():
            same = sha256(old / name) == sha256(new / name)
            rep.check(same, f"{name} is byte-identical")

    # ---- C. the ledger
    if (old / "partD_resfeas_ledger.tsv").exists():
        ledger_checks(old, new, rep)

    # ---- D. JSON documents
    for name in sorted(JSON_ALLOWED):
        o, n = old / name, new / name
        if not (o.exists() and n.exists()):
            continue
        json_checks(name, o, n, rep)

    # ---- E. nothing in Han remains
    for p in sorted(new.iterdir()):
        if p.suffix.lower() in {".tsv", ".json", ".md", ".txt", ".cff"}:
            bad = [i + 1 for i, l in enumerate(
                p.read_text(encoding="utf-8", errors="replace").splitlines()) if CJK.search(l)]
            rep.check(not bad, f"no Han characters left in {p.name}",
                      f"{len(bad)} line(s), e.g. {bad[:5]}")


def ledger_checks(old: Path, new: Path, rep: Report) -> None:
    o_lines, o_nl, o_tr = lines_of(old / "partD_resfeas_ledger.tsv")
    n_lines, n_nl, n_tr = lines_of(new / "partD_resfeas_ledger.tsv")
    rep.check(o_nl == n_nl and o_tr == n_tr, "ledger keeps its line-ending convention",
              f"old {o_nl!r}/{o_tr} vs new {n_nl!r}/{n_tr}")
    rep.check(len(o_lines) == len(n_lines) == 30,
              f"ledger row count unchanged (old {len(o_lines) - 1}, new {len(n_lines) - 1})")

    o_hdr, n_hdr = o_lines[0].split("\t"), n_lines[0].split("\t")
    rep.check(o_hdr == n_hdr, "ledger header (field names and order) unchanged",
              f"{o_hdr} vs {n_hdr}")
    i_reason = n_hdr.index("reason")

    gate_ok = True
    reason_changed = reason_ascii = 0
    for idx, (o, n) in enumerate(zip(o_lines[1:], n_lines[1:]), start=1):
        oc, nc = o.split("\t"), n.split("\t")
        for i, col in enumerate(n_hdr):
            if col == "reason":
                continue
            if oc[i] != nc[i]:
                gate_ok = False
                rep.check(False, f"ledger row {idx} column {col} differs", f"{oc[i]!r} -> {nc[i]!r}")
        if oc[i_reason] != nc[i_reason]:
            reason_changed += 1
        if not CJK.search(nc[i_reason]):
            reason_ascii += 1
        else:
            rep.check(False, f"ledger row {idx} reason is still Han", nc[i_reason][:50])

    rep.check(gate_ok, "every ledger column except `reason` is identical, row by row "
                       "(includes G1_mRNA, G2_UC, G3_powered_control, G4_in_vivo, USABLE)")
    rep.check(reason_changed == 29, f"all 29 `reason` cells changed (found {reason_changed})")
    rep.check(reason_ascii == 29, f"all 29 `reason` cells are ASCII (found {reason_ascii})")

    # every English component must be one a declared template produces
    rxs = [NORM.template_to_regex(v["en"]) for v in I18N.REASON.values()]
    unknown = []
    for line in n_lines[1:]:
        for comp in line.split("\t")[i_reason].split(I18N.REASON_SEPARATOR):
            if not any(r.match(comp) for r in rxs):
                unknown.append(comp)
    rep.check(not unknown, "every `reason` component matches a declared English template",
              f"unmatched: {unknown[:3]}")

    # gate columns must be identical as a set, explicitly
    o_gate = [[l.split("\t")[o_hdr.index(c)] for c in GATE_COLS] for l in o_lines[1:]]
    n_gate = [[l.split("\t")[n_hdr.index(c)] for c in GATE_COLS] for l in n_lines[1:]]
    rep.check(o_gate == n_gate, "gate matrix and USABLE verdict identical (row order included)")


def json_checks(name: str, o: Path, n: Path, rep: Report) -> None:
    old = json.loads(o.read_text(encoding="utf-8"))
    new = json.loads(n.read_text(encoding="utf-8"))
    ol, nl = dict(leaves(old)), dict(leaves(new))

    rep.check(set(ol) == set(nl), f"{name}: key set unchanged",
              f"diff {sorted(set(ol) ^ set(nl))[:4]}")

    moved = [k for k in ol if k in nl and ol[k] != nl[k]]
    bad = [k for k in moved if not allowed(name, k)]
    rep.check(not bad, f"{name}: every changed key path is allow-listed",
              f"unexpected changes at {['.'.join(map(str, k)) for k in bad][:4]}")

    # numeric fidelity: same value and same repr everywhere
    drift = [k for k in ol if isinstance(ol[k], (int, float)) and not isinstance(ol[k], bool)
             and (k not in nl or ol[k] != nl[k] or repr(ol[k]) != repr(nl[k]))]
    rep.check(not drift, f"{name}: every numeric value byte-identical",
              f"{['.'.join(map(str, k)) for k in drift][:4]}")

    if name == "partD_resfeas_libwide_summary.json":
        same = [k for k in ol if isinstance(ol[k], str) and ol[k] == nl[k]]
        rep.check(len(same) >= 1, f"{name}: prose keys replaced, not all text touched")
    if name == "datapackage.json":
        o_sch = [r.get("schema", {}).get("fields") for r in old.get("resources", [])]
        n_sch = [r.get("schema", {}).get("fields") for r in new.get("resources", [])]
        rep.check(len(o_sch) == len(n_sch) == 4, "datapackage.json: four tabular resources")
        for a, b in zip(o_sch, n_sch):
            if a and b:
                rep.check([f["name"] for f in a] == [f["name"] for f in b],
                          "datapackage.json: field names and order unchanged in a resource")
        # constraints and types must be untouched
        def sig(fields):
            return [(f["name"], f.get("type"), json.dumps(f.get("constraints"), sort_keys=True))
                    for f in (fields or [])]
        rep.check([sig(s) for s in o_sch] == [sig(s) for s in n_sch],
                  "datapackage.json: field types and constraints unchanged")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", type=Path, default=OLD_DEFAULT)
    ap.add_argument("--new", type=Path, default=NEW_DEFAULT)
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args()

    print("=" * 74)
    print(f"language-normalisation regression: {args.old.name}  ->  {args.new.name}")
    print("=" * 74)
    rep = Report()
    run(args.old, args.new, rep)

    total, failed = len(rep.rows), rep.failed
    print("-" * 74)
    print(f"checks: {total - failed}/{total} passed" + (f", {failed} FAILED" if failed else ""))
    print("RESULT: " + ("language-only change confirmed by every check"
                        if not failed else "NOT a language-only change"))

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps({
            "old": str(args.old), "new": str(args.new),
            "total": total, "failed": failed,
            "checks": [{"ok": ok, "label": lbl, "detail": d} for ok, lbl, d in rep.rows],
        }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        print(f"report -> {args.json_out}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
