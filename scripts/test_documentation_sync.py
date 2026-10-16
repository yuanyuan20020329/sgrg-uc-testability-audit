# -*- coding: utf-8 -*-
"""Regression test: v1.0.2 must be a *documentation-only* change against v1.0.1.

v1.0.2 exists for one reason — the v1.0.1 ``datapackage.json`` described
``G3_powered_control`` as a *"power gate"*, while the Data Descriptor that cites the
deposit explains the same field as a **preset design-eligibility rule, not a power
calculation**. A reader who downloads the package must not be told a different semantics
than the manuscript. The same release writes the candidate-pairing and control-vocabulary
warnings, and the three per-table UC conventions, into the dictionary and the README.

That is the whole scope. This test is the machine-checkable statement of it.

What it proves, file by file:

* the four tabular records, the two summary JSONs and the frozen snapshot are
  **byte-identical** to v1.0.1 — no value, gate, row or accession moved;
* ``datapackage.json`` keeps its resource list, field order, field names, types and
  constraints — only the version, the creation date and the prose descriptions of the
  fields listed in ``DESCRIPTION_FIELDS`` (plus the top-level blurb) may differ;
* ``README.md`` / ``CHANGELOG.md`` are reissued and must state the new version, and
  ``LICENSE-DATA.md`` keeps its CC BY 4.0 statement;
* ``manifest.json`` is regenerated but keeps the same file inventory;
* no Han character is left anywhere in the payload.

Run it with no arguments to compare the current ``resource_paper/deposit/`` against the
frozen v1.0.1 reference; pass ``--old``/``--new`` to compare other pairs.

    python scripts/test_documentation_sync.py
    python scripts/test_documentation_sync.py --json-out results/documentation_sync.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

OLD_DEFAULT = ROOT / "resource_paper" / "frozen_deposit_v1.0.1"
NEW_DEFAULT = ROOT / "resource_paper" / "deposit"

NEW_VERSION = "1.0.2"

CJK = re.compile(r"[\u4e00-\u9fff]")
DOI = re.compile(r"10\.5281/zenodo\.\d{6,}")

# Files that a documentation-only change must not touch at all.
FROZEN_FILES = {
    "cxg_coll_detail_snapshot_v1.0.0.zip",
    "partD_resfeas_geo_series.tsv",
    "partD_resfeas_libwide_collections.tsv",
    "partD_resfeas_libwide_datasets.tsv",
    "partD_resfeas_ledger.tsv",
    "partD_resfeas_libwide_summary.json",
    "partD_resfeas_verdict.json",
}

# Files that are reissued by the release and may therefore differ freely.
REISSUED = {"README.md", "CHANGELOG.md", "manifest.json", "MANIFEST.sha256"}

# datapackage.json leaves that are allowed to change, addressed by resource + field name
# rather than by position, so a future re-ordering cannot silently widen the allowance.
DESCRIPTION_FIELDS: dict[str, set[str]] = {
    "libwide_datasets": {"diseases"},
    "libwide_collections": {"has_UC", "has_CD"},
    "geo_series": {"diseases", "n_paired_donors", "n_paired_same_site"},
    "gate_ledger": {"G2_UC", "G3_powered_control", "USABLE"},
}
TOP_LEVEL_ALLOWED = {"version", "created", "description"}


# --------------------------------------------------------------------------- helpers
def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


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
    rep.check(old.is_dir(), f"reference payload exists at {old}")
    rep.check(new.is_dir(), f"current payload exists at {new}")

    # ---- A. inventory
    old_files = sorted(p.name for p in old.iterdir() if p.name != "_FROZEN_REFERENCE.json")
    new_files = sorted(p.name for p in new.iterdir())
    rep.check(old_files == new_files, "file inventory unchanged",
              f"only in old: {set(old_files) - set(new_files)}; "
              f"only in new: {set(new_files) - set(old_files)}")
    rep.check(len(new_files) == 14, f"payload has 14 files (found {len(new_files)})")

    # ---- B. files that must not move
    for name in sorted(FROZEN_FILES):
        if (old / name).exists() and (new / name).exists():
            rep.check(sha256(old / name) == sha256(new / name), f"{name} is byte-identical")

    # ---- C. datapackage.json — structure frozen, descriptions allow-listed
    dp_checks(old / "datapackage.json", new / "datapackage.json", rep)

    # ---- D. reissued prose files
    readme = (new / "README.md").read_text(encoding="utf-8")
    changelog = (new / "CHANGELOG.md").read_text(encoding="utf-8")
    lic = (new / "LICENSE-DATA.md").read_text(encoding="utf-8")
    rep.check(f"Version {NEW_VERSION}" in readme,
              f"deposit README states Version {NEW_VERSION}")
    rep.check(f"## [{NEW_VERSION}]" in changelog,
              f"CHANGELOG.md carries a [{NEW_VERSION}] section")
    rep.check("## [1.0.1]" in changelog and "## [1.0.0]" in changelog,
              "CHANGELOG.md keeps the [1.0.1] and [1.0.0] entries (history, not overwritten)")
    rep.check("Known limitations" in readme,
              "README.md carries the new 'Known limitations' warning section")
    for needle in ("design-eligibility gate", "candidate", "normal pressure hydrocephalus"):
        rep.check(needle in readme, f"README.md states the warning keyword {needle!r}")
    rep.check("CC BY 4.0" in lic, "LICENSE-DATA.md keeps its CC BY 4.0 statement")
    rep.check(bool(DOI.search(lic)), "LICENSE-DATA.md carries a resolved dataset DOI",
              "no 10.5281/zenodo.* token found")

    # ---- E. manifest.json — regenerated, inventory preserved
    mo = json.loads((old / "manifest.json").read_text(encoding="utf-8"))
    mn = json.loads((new / "manifest.json").read_text(encoding="utf-8"))
    rep.check(mn.get("version") == NEW_VERSION, "manifest.json version is " + NEW_VERSION,
              str(mn.get("version")))
    names_o = sorted(f["name"] for f in mo.get("files", []))
    names_n = sorted(f["name"] for f in mn.get("files", []))
    rep.check(names_o == names_n, "manifest.json keeps the same file inventory",
              f"only in old: {set(names_o) - set(names_n)}; "
              f"only in new: {set(names_n) - set(names_o)}")

    # ---- F. nothing in Han remains
    for p in sorted(new.iterdir()):
        if p.suffix.lower() in {".tsv", ".json", ".md", ".txt", ".cff"}:
            bad = [i + 1 for i, l in enumerate(
                p.read_text(encoding="utf-8", errors="replace").splitlines()) if CJK.search(l)]
            rep.check(not bad, f"no Han characters left in {p.name}",
                      f"{len(bad)} line(s), e.g. {bad[:5]}")


def dp_checks(o_path: Path, n_path: Path, rep: Report) -> None:
    old = json.loads(o_path.read_text(encoding="utf-8"))
    new = json.loads(n_path.read_text(encoding="utf-8"))

    rep.check(new.get("version") == NEW_VERSION,
              f"datapackage.json version is {NEW_VERSION}", str(new.get("version")))

    # top level: only the allow-listed scalar/blurb keys may differ
    moved_top = sorted(k for k in set(old) | set(new)
                       if old.get(k) != new.get(k) and not isinstance(old.get(k), (list, dict)))
    bad_top = [k for k in moved_top if k not in TOP_LEVEL_ALLOWED]
    rep.check(not bad_top, "datapackage.json: only allow-listed top-level keys changed",
              f"unexpected: {bad_top}")

    o_res, n_res = old.get("resources", []), new.get("resources", [])
    rep.check(len(o_res) == len(n_res) == 4, "datapackage.json: four tabular resources")
    rep.check([r["name"] for r in o_res] == [r["name"] for r in n_res],
              "datapackage.json: resource names and order unchanged")

    for ro, rn in zip(o_res, n_res):
        rname = rn["name"]
        fo = ro.get("schema", {}).get("fields") or []
        fn = rn.get("schema", {}).get("fields") or []
        rep.check([f["name"] for f in fo] == [f["name"] for f in fn],
                  f"{rname}: field names and order unchanged")

        def sig(fields):
            return [(f["name"], f.get("type"), json.dumps(f.get("constraints"), sort_keys=True))
                    for f in fields]
        rep.check(sig(fo) == sig(fn), f"{rname}: field types and constraints unchanged")

        allowed = DESCRIPTION_FIELDS.get(rname, set())
        changed, bad = [], []
        for a, b in zip(fo, fn):
            if a.get("description") != b.get("description"):
                changed.append(a["name"])
                if a["name"] not in allowed:
                    bad.append(a["name"])
        rep.check(not bad,
                  f"{rname}: only allow-listed field descriptions changed",
                  f"unexpected description change(s): {bad}")
        rep.check(set(changed) == allowed,
                  f"{rname}: exactly the intended field descriptions were updated",
                  f"changed {sorted(changed)} vs intended {sorted(allowed)}")

    # numeric fidelity at the schema level
    def all_nums(node, path=()):
        if isinstance(node, dict):
            for k, v in node.items():
                yield from all_nums(v, path + (str(k),))
        elif isinstance(node, list):
            for i, v in enumerate(node):
                yield from all_nums(v, path + (str(i),))
        elif isinstance(node, (int, float)) and not isinstance(node, bool):
            yield path, node
    on, nn = dict(all_nums(old)), dict(all_nums(new))
    drift = [k for k in on if k in nn and repr(on[k]) != repr(nn[k])]
    rep.check(not drift, "datapackage.json: every numeric constraint byte-identical",
              f"{['.'.join(k) for k in drift][:4]}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", type=Path, default=OLD_DEFAULT)
    ap.add_argument("--new", type=Path, default=NEW_DEFAULT)
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args()

    print("=" * 74)
    print(f"documentation-sync regression: {args.old.name}  ->  {args.new.name}")
    print("=" * 74)
    rep = Report()
    run(args.old, args.new, rep)

    total, failed = len(rep.rows), rep.failed
    print("-" * 74)
    print(f"checks: {total - failed}/{total} passed" + (f", {failed} FAILED" if failed else ""))
    print("RESULT: " + ("documentation-only change confirmed by every check"
                        if not failed else "NOT a documentation-only change"))

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
