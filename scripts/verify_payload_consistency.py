# -*- coding: utf-8 -*-
"""Produce the *actual-file* evidence that the v1.0.2 payload is intact and coherent.

The fourth and fifth review rounds both asked for the same thing: statements such as
"14 files, 3,108,031 bytes, 12 SHA-256 entries, the four TSVs are byte-identical to
v1.0.1" must be backed by a check run against the real files, not by prose.  This script
runs that check and writes a durable report, so the claim can be re-verified by anyone
who has the working tree.

What it does, all against ``resource_paper/deposit/`` as it exists on disk:

  A  inventory          file count, per-file bytes and sha256, grand total
  B  MANIFEST.sha256    every listed hash recomputed; the two self-referential exclusions
                        are exactly right
  C  manifest.json      per-file bytes/sha256 agree with the actual files; version,
                        licence and release date
  D  frozen files       the seven data files (four TSVs + two JSONs + the snapshot zip)
                        are byte-identical to the frozen v1.0.1 reference
  E  datapackage.json   version/date/licence, four tabular resources, field names, types,
                        constraints and numeric bounds unchanged, and *exactly* the
                        allow-listed descriptions updated
  F  prose files        README / CHANGELOG / LICENSE content requirements
  G  language           no Han character anywhere in the payload
  H  PENDING            things that are legitimately not final yet, listed rather than
                        failed, because the version DOIs do not exist until the record is
                        published
  I  manuscript         no DOI placeholder, interim sentence present the right number of
                        times, and no DOI cited that disagrees with the DOI registry

Usage:

    python scripts/verify_payload_consistency.py
    python scripts/verify_payload_consistency.py --json-out results/v1p0p2_payload_verification.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

PAYLOAD_DEFAULT = ROOT / "resource_paper" / "deposit"
REFERENCE_DEFAULT = ROOT / "resource_paper" / "frozen_deposit_v1.0.1"
JSON_DEFAULT = ROOT / "results" / "v1p0p2_payload_verification.json"
MD_DEFAULT = ROOT / "results" / "v1p0p2_payload_verification.md"

VERSION = "1.0.2"
PREVIOUS_VERSION = "1.0.1"
RELEASE_DATE = "2026-10-16"
EXPECTED_FILES = 14
EXPECTED_BYTES = 3_108_031
EXPECTED_MANIFEST_ENTRIES = 12

CJK = re.compile(r"[\u4e00-\u9fff]")
DOI = re.compile(r"10\.5281/zenodo\.(\d{6,})")
PLACEHOLDER = re.compile(r"zenodo\.(XXXXXXX|YYYYYYY)")

# Data files a documentation-only release must not touch.
FROZEN_FILES = [
    "partD_resfeas_libwide_datasets.tsv",
    "partD_resfeas_libwide_collections.tsv",
    "partD_resfeas_geo_series.tsv",
    "partD_resfeas_ledger.tsv",
    "partD_resfeas_libwide_summary.json",
    "partD_resfeas_verdict.json",
    "cxg_coll_detail_snapshot_v1.0.0.zip",
]
FOUR_TSV = [f for f in FROZEN_FILES if f.endswith(".tsv")]

# datapackage.json leaves allowed to change, addressed by resource + field name.
DESCRIPTION_FIELDS: dict[str, set[str]] = {
    "libwide_datasets": {"diseases"},
    "libwide_collections": {"has_UC", "has_CD"},
    "geo_series": {"diseases", "n_paired_donors", "n_paired_same_site"},
    "gate_ledger": {"G2_UC", "G3_powered_control", "USABLE"},
}
TOP_LEVEL_ALLOWED = {"version", "created", "description"}

CN_MANUSCRIPT = ROOT / "resource_paper" / "Data_Descriptor_v1_CN.md"
INTERIM_PHRASE = "1.0.2 版已完成修订，待正式归档后以版本 DOI 引用"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_doi_registry(path: Path) -> dict[str, str]:
    """First line = version DOI; '# concept doi:' / '# version:' lines parsed too."""
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines()):
        line = line.strip()
        if not line:
            continue
        if i == 0 and line.startswith("10.5281/"):
            out["version_doi"] = line
        elif line.lower().startswith("# concept doi:"):
            out["concept_doi"] = line.split(":", 1)[1].strip()
        elif line.lower().startswith("# version:"):
            out["version"] = line.split(":", 1)[1].strip()
    return out


class Report:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def check(self, ok: bool, cid: str, label: str, detail: str = "") -> bool:
        ok = bool(ok)
        self.rows.append({"ok": ok, "id": cid, "label": label, "detail": detail})
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {cid}  {label}" + (f"  -- {detail}" if detail and not ok else ""))
        return ok

    def note(self, cid: str, label: str, detail: str = "") -> None:
        self.rows.append({"ok": True, "id": cid, "label": label, "detail": detail,
                          "pending": True})
        print(f"  [ .. ] {cid}  {label}" + (f"  -- {detail}" if detail else ""))

    @property
    def failed(self) -> int:
        return sum(1 for r in self.rows if not r["ok"])

    @property
    def pending(self) -> list[dict]:
        return [r for r in self.rows if r.get("pending")]


def section(title: str) -> None:
    print(f"\n{title}\n" + "-" * 74)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--payload", type=Path, default=PAYLOAD_DEFAULT)
    ap.add_argument("--reference", type=Path, default=REFERENCE_DEFAULT)
    ap.add_argument("--json-out", type=Path, default=JSON_DEFAULT)
    ap.add_argument("--md-out", type=Path, default=MD_DEFAULT)
    args = ap.parse_args()

    payload, ref = args.payload, args.reference
    rep = Report()
    inventory: list[dict] = []

    print("=" * 74)
    print(f"payload consistency verification: {payload}")
    print(f"reference (v{PREVIOUS_VERSION} frozen):    {ref}")
    print("=" * 74)

    rep.check(payload.is_dir(), "A0", "payload directory exists", str(payload))
    rep.check(ref.is_dir(), "A0b", "frozen reference directory exists", str(ref))
    if not payload.is_dir() or not ref.is_dir():
        return 1

    files = sorted(p for p in payload.iterdir() if p.is_file())
    names = [p.name for p in files]

    # ---------------------------------------------------------------- A. inventory
    section("A. inventory — file count, per-file bytes and sha256")
    total = 0
    for p in files:
        b, h = p.stat().st_size, sha256(p)
        total += b
        inventory.append({"name": p.name, "bytes": b, "sha256": h})
    for row in inventory:
        print(f"    {row['bytes']:>9}  {row['sha256'][:16]}  {row['name']}")
    rep.check(len(files) == EXPECTED_FILES, "A1",
              f"payload has exactly {EXPECTED_FILES} files", f"found {len(files)}")
    rep.check(total == EXPECTED_BYTES, "A2",
              f"payload totals {EXPECTED_BYTES:,} bytes", f"found {total:,}")
    print(f"    ---- {len(files)} files, {total:,} bytes")

    # -------------------------------------------------------------- B. MANIFEST.sha256
    section(f"B. MANIFEST.sha256 — every listed hash recomputed")
    mf_path = payload / "MANIFEST.sha256"
    rep.check(mf_path.exists(), "B0", "MANIFEST.sha256 present")
    entries: list[tuple[str, str]] = []
    if mf_path.exists():
        for line in mf_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                entries.append((parts[-1].lstrip("*"), parts[0].lower()))
    rep.check(len(entries) == EXPECTED_MANIFEST_ENTRIES, "B1",
              f"MANIFEST.sha256 covers {EXPECTED_MANIFEST_ENTRIES} payload files",
              f"found {len(entries)}")
    listed = [n for n, _ in entries]
    rep.check(set(listed) <= set(names), "B2",
              "every MANIFEST entry exists in the payload",
              f"missing: {sorted(set(listed) - set(names))}")
    bad: list[str] = []
    for name, declared in entries:
        actual = sha256(payload / name).lower() if (payload / name).exists() else ""
        if actual != declared:
            bad.append(f"{name}: declared {declared[:12]} actual {actual[:12] or '-'}")
    rep.check(not bad, "B3", "every declared sha256 matches the file on disk",
              "; ".join(bad))
    excluded = sorted(set(names) - set(listed))
    rep.check(excluded == ["MANIFEST.sha256", "manifest.json"], "B4",
              "the only unlisted files are the two self-referential manifests",
              f"unlisted: {excluded}")

    # ---------------------------------------------------------------- C. manifest.json
    section("C. manifest.json — version, licence, per-file digests")
    mp = payload / "manifest.json"
    man = json.loads(mp.read_text(encoding="utf-8")) if mp.exists() else {}
    rep.check(man.get("version") == VERSION, "C1", f"manifest.json version is {VERSION}",
              str(man.get("version")))
    rep.check(str(man.get("released", "")) == RELEASE_DATE, "C2",
              f"manifest.json released is {RELEASE_DATE}", str(man.get("released")))
    rep.check(str(man.get("licence", "")).upper() == "CC-BY-4.0", "C3",
              "manifest.json licence is CC-BY-4.0", str(man.get("licence")))
    man_files = {f["name"]: f for f in man.get("files", [])}
    rep.check(sorted(man_files) == sorted(listed), "C4",
              "manifest.json file list equals the MANIFEST.sha256 file list",
              f"only in manifest: {sorted(set(man_files) - set(listed))}; "
              f"only in MANIFEST: {sorted(set(listed) - set(man_files))}")
    inventory_map = {r["name"]: r for r in inventory}
    drift = []
    for name in sorted(set(man_files) & set(inventory_map)):
        m = man_files[name]
        a = inventory_map[name]
        if m.get("sha256") != a["sha256"] or m.get("bytes") != a["bytes"]:
            drift.append(name)
    rep.check(not drift, "C5",
              "manifest.json bytes+sha256 agree with the actual files", f"{drift}")
    snap = man.get("source_snapshot", {})
    rep.check(snap.get("n_files") == 397 and snap.get("frozen") == "2026-10-09", "C6",
              "manifest.json source snapshot is 397 files frozen 2026-10-09",
              f"{snap}")

    # ---------------------------------------------------------------- D. frozen files
    section(f"D. frozen data files — byte-identical to v{PREVIOUS_VERSION}")
    for name in FROZEN_FILES:
        o, n = ref / name, payload / name
        if not (o.exists() and n.exists()):
            rep.check(False, f"D·{name}", f"{name} exists on both sides",
                      f"old={o.exists()} new={n.exists()}")
            continue
        rep.check(sha256(o) == sha256(n), f"D·{name}", f"{name} is byte-identical")
    for name in FOUR_TSV:
        rep.check(sha256(ref / name) == sha256(payload / name), f"D·4TSV·{name}",
                  f"TSV {name} unchanged byte for byte")

    # ------------------------------------------------------------ E. datapackage.json
    section("E. datapackage.json — structure frozen, descriptions allow-listed")
    dp = json.loads((payload / "datapackage.json").read_text(encoding="utf-8"))
    dpo = json.loads((ref / "datapackage.json").read_text(encoding="utf-8"))
    rep.check(dp.get("version") == VERSION, "E1", f"datapackage.json version is {VERSION}",
              str(dp.get("version")))
    rep.check(str(dp.get("created")) == RELEASE_DATE, "E2",
              f"datapackage.json created is {RELEASE_DATE}", str(dp.get("created")))
    rep.check([l.get("name") for l in dp.get("licenses", [])] == ["CC-BY-4.0"], "E3",
              "datapackage.json declares the CC-BY-4.0 licence",
              str(dp.get("licenses")))
    res = dp.get("resources", [])
    rep.check(len(res) == 4, "E4", "four tabular resources", f"found {len(res)}")
    rep.check(all(r.get("profile") == "tabular-data-resource" for r in res), "E5",
              "every resource profile is tabular-data-resource",
              str([(r['name'], r.get('profile')) for r in res]))
    rep.check([r["name"] for r in res] == [r["name"] for r in dpo.get("resources", [])],
              "E6", "resource names and order unchanged")
    num_drift: list[str] = []
    desc_changed: dict[str, list[str]] = {}
    desc_unexpected: list[str] = []
    for ro, rn in zip(dpo.get("resources", []), res):
        rname = rn["name"]
        fo = ro.get("schema", {}).get("fields") or []
        fn = rn.get("schema", {}).get("fields") or []
        rep.check([f["name"] for f in fo] == [f["name"] for f in fn], f"E·{rname}",
                  f"{rname}: field names and order unchanged")
        sig = lambda fs: [(f["name"], f.get("type"),
                           json.dumps(f.get("constraints"), sort_keys=True)) for f in fs]
        rep.check(sig(fo) == sig(fn), f"E·{rname}·sig",
                  f"{rname}: field types and constraints unchanged")
        changed = [a["name"] for a, b in zip(fo, fn)
                   if a.get("description") != b.get("description")]
        allowed = DESCRIPTION_FIELDS.get(rname, set())
        desc_changed[rname] = changed
        extra = [n for n in changed if n not in allowed]
        missed = sorted(allowed - set(changed))
        desc_unexpected += [f"{rname}.{n}" for n in extra]
        rep.check(not extra, f"E·{rname}·allow",
                  f"{rname}: no description changed outside the allow-list", f"{extra}")
        rep.check(not missed, f"E·{rname}·intended",
                  f"{rname}: every intended field description was updated", f"missing {missed}")

        def nums(node, path=()):
            if isinstance(node, dict):
                for k, v in node.items():
                    yield from nums(v, path + (str(k),))
            elif isinstance(node, list):
                for i, v in enumerate(node):
                    yield from nums(v, path + (str(i),))
            elif isinstance(node, (int, float)) and not isinstance(node, bool):
                yield path, node
        on, nn = dict(nums(ro)), dict(nums(rn))
        for k in on:
            if k in nn and repr(on[k]) != repr(nn[k]):
                num_drift.append(f"{rname}." + ".".join(k))
    rep.check(not num_drift, "E7", "every numeric constraint is unchanged", f"{num_drift[:4]}")
    top_moved = sorted(k for k in set(dpo) | set(dp)
                       if dpo.get(k) != dp.get(k) and not isinstance(dpo.get(k), (list, dict)))
    rep.check([k for k in top_moved if k not in TOP_LEVEL_ALLOWED] == [], "E8",
              "only allow-listed top-level keys changed", f"{top_moved}")

    # ------------------------------------------------------------------ F. prose files
    section("F. prose files — README / CHANGELOG / LICENSE")
    readme = (payload / "README.md").read_text(encoding="utf-8")
    changelog = (payload / "CHANGELOG.md").read_text(encoding="utf-8")
    lic = (payload / "LICENSE-DATA.md").read_text(encoding="utf-8")
    rep.check(f"Version {VERSION}" in readme, "F1", f"README states Version {VERSION}")
    rep.check("Known limitations" in readme, "F2",
              "README carries the 'Known limitations' section")
    for needle in ("design-eligibility gate", "candidate", "normal pressure hydrocephalus"):
        rep.check(needle in readme, f"F·{needle}",
                  f"README states the warning keyword {needle!r}")
    for sec in (VERSION, PREVIOUS_VERSION, "1.0.0"):
        rep.check(f"## [{sec}]" in changelog, f"F·{sec}",
                  f"CHANGELOG carries a [{sec}] section")
    rep.check("CC BY 4.0" in lic, "F3", "LICENSE-DATA.md keeps its CC BY 4.0 statement")

    # ----------------------------------------------------------------- G. no Han
    section("G. language — no Han character in the payload")
    han: list[str] = []
    for p in files:
        if p.suffix.lower() in {".tsv", ".json", ".md", ".txt", ".cff"}:
            txt = p.read_text(encoding="utf-8", errors="replace")
            if CJK.search(txt):
                han.append(p.name)
    rep.check(not han, "G1", "no Han character anywhere in the payload", f"{han}")

    # ------------------------------------------------------- H. legitimately pending
    section("H. pending — must be resolved by the publication step, not a defect")
    prev_data = read_doi_registry(ROOT / "resource_paper" / "dataset_doi.txt")
    prev_code = read_doi_registry(ROOT / "resource_paper" / "code_doi.txt")
    embedded = sorted({m.group(0) for p in files
                       if p.suffix.lower() in {".tsv", ".json", ".md", ".txt"}
                       for m in DOI.finditer(p.read_text(encoding="utf-8", errors="replace"))})
    stale = [d for d in embedded
             if d.endswith(prev_data.get("version_doi", "").rsplit(".", 1)[-1])
             or d.endswith(prev_code.get("version_doi", "").rsplit(".", 1)[-1])]
    rep.note("H1", "payload still embeds the previous version DOIs "
                   f"({len(stale)} distinct: {stale})",
             "expected before publication: the v1.0.2 DOI does not exist yet. "
             "publish_dataset.py stage 2 re-runs build_data_deposit.py --data-doi "
             "<reserved DOI>, which rewrites exactly these lines; re-run this script "
             "afterwards and A1/A2/B*/C*/D*/E* must stay green")
    cl_line = next((l.strip() for l in changelog.splitlines()
                    if f"re-released as v{VERSION}" in l), "")
    rep.note("H2", "CHANGELOG v1.0.2 line cites the companion code record",
             cl_line or "(line not found)")
    rep.note("H3", "the code DOI inside the payload is read from code_doi.txt at build "
                   "time, so it becomes correct only after backfill_v1p0p2_dois.py "
                   "--code-doi has run and the payload has been rebuilt")
    rep.note("H4", "reserved dataset DOI in the registry is "
                   f"{prev_data.get('version_doi')} ({prev_data.get('version')}); code "
                   f"{prev_code.get('version_doi')} ({prev_code.get('version')})")

    # -------------------------------------------------------------- I. manuscript
    section("I. manuscript — no placeholder, interim sentence counted, DOIs agree")
    if CN_MANUSCRIPT.exists():
        md = CN_MANUSCRIPT.read_text(encoding="utf-8")
        ph = PLACEHOLDER.findall(md)
        rep.check(not ph, "I1", "no DOI placeholder token anywhere in the CN manuscript",
                  f"{ph}")
        n_interim = md.count(INTERIM_PHRASE)
        published = prev_data.get("version") == VERSION
        if published:
            rep.check(n_interim == 0, "I2",
                      "interim 'pending archival' sentence removed after publication",
                      f"{n_interim} left")
            rep.check(md.count(prev_data.get("version_doi", "\x00")) >= 3, "I3",
                      "the published dataset version DOI is cited",
                      f"{md.count(prev_data.get('version_doi',''))} occurrence(s)")
        else:
            rep.check(n_interim == 5, "I2",
                      "interim 'pending archival' sentence appears 5x "
                      "(abstract / software / data records / data availability / code "
                      "availability)", f"found {n_interim}")
        rep.check(md.count("10.5281/zenodo.23297702") >= 1, "I4",
                  "data concept DOI is cited in the manuscript")
        rep.check(md.count("10.5281/zenodo.23291758") >= 1, "I5",
                  "code concept DOI is cited in the manuscript")
        rep.check(not re.search(r"10\.5281/zenodo\.23299341|10\.5281/zenodo\.23298914", md),
                  "I6", "the v1.0.1 version DOIs are no longer cited in the manuscript")
    else:
        rep.check(False, "I0", "CN manuscript exists", str(CN_MANUSCRIPT))

    # ------------------------------------------------------------------- summary
    total_checks = len([r for r in rep.rows if not r.get("pending")])
    failed = rep.failed
    print("\n" + "=" * 74)
    print(f"checks: {total_checks - failed}/{total_checks} passed"
          + (f", {failed} FAILED" if failed else ""))
    print(f"pending (resolved by publication): {len(rep.pending)}")
    print("RESULT: " + ("payload is internally consistent and ready to publish"
                        if not failed else "NOT ready — see FAILED checks"))

    payload_json = {
        "payload": str(payload), "reference": str(ref),
        "version": VERSION, "expected_files": EXPECTED_FILES,
        "expected_bytes": EXPECTED_BYTES, "actual_files": len(files),
        "actual_bytes": total,
        "inventory": inventory,
        "manifest_entries": len(entries),
        "frozen_files": FROZEN_FILES,
        "description_changes": desc_changed,
        "checks_total": total_checks, "checks_failed": failed,
        "pending": rep.pending,
        "checks": rep.rows,
    }
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload_json, indent=2, ensure_ascii=False)
                                 + "\n", encoding="utf-8", newline="\n")
        print(f"report -> {args.json_out}")

    if args.md_out:
        lines = ["# v1.0.2 数据包载荷一致性核验（实测）", "",
                 f"- 核验对象：`{payload.relative_to(ROOT).as_posix()}/`",
                 f"- 参照（v{PREVIOUS_VERSION} 冻结副本）：`{ref.relative_to(ROOT).as_posix()}/`",
                 f"- 核验时间：{RELEASE_DATE}",
                 f"- 结论：**{total_checks - failed}/{total_checks} 项通过**"
                 + (f"，{failed} 项失败" if failed else "")
                 + f"；{len(rep.pending)} 项为“待发布步骤解决”，非缺陷。",
                 "",
                 "## A. 载荷清单（逐文件实测字节与 SHA-256）", "",
                 f"合计 **{len(files)} 个文件 / {total:,} 字节**（声明值 "
                 f"{EXPECTED_FILES} / {EXPECTED_BYTES:,}）。", "",
                 "| 文件 | 字节 | SHA-256 |", "|---|---:|---|"]
        for row in inventory:
            lines.append(f"| `{row['name']}` | {row['bytes']:,} | `{row['sha256']}` |")
        lines += ["", "## B–I. 核验项", "",
                  "| 编号 | 核验项 | 结果 | 备注 |", "|---|---|---|---|"]
        for r in rep.rows:
            mark = "待发布" if r.get("pending") else ("通过" if r["ok"] else "**失败**")
            lines.append(f"| {r['id']} | {r['label']} | {mark} | {r.get('detail','') or ''} |")
        args.md_out.parent.mkdir(parents=True, exist_ok=True)
        args.md_out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        print(f"report -> {args.md_out}")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
