#!/usr/bin/env python
"""Audit the Scientific Data Data Descriptor against the frozen deposit and the live record.

Every numeric claim and every cited filename in ``resource_paper/Data_Descriptor_v1_EN.md``
is re-derived from artifacts rather than trusted: the four TSV dimensions come from the
files themselves, the package inventory from ``deposit/manifest.json`` plus the deposit
listing, the ledger counts from the ledger itself, the DOIs and the version from
``resource_paper/dataset_doi.txt`` / ``resource_paper/code_doi.txt`` / ``deposit/manifest.json``,
and the registry facts from the DOI registry report
(``results/dataset_doi_verification_registry.json``) when present.

Two modes, one script — because the thing that has to be checked is different before and
after publication:

``--pre-publish``  the **hard gate for ``scripts/publish_dataset.py`` (stage 0)**. Nothing is
    public yet, so there is no registry report and the manuscript still cites the previous
    version's DOI. This mode therefore audits the *payload* only: inventory, checksums,
    version agreement across four files, no placeholder DOIs, the code-record link, the
    two manifest artifacts, and the machine-checked language-only regression against the
    frozen v1.0.0 payload. It answers one question: *is this byte-for-byte what we intend
    to publish?*

default (full)     everything above **plus** the manuscript cross-check (every cited
    filename exists publicly, every quoted number matches the payload, the three DOIs are
    present and distinct) **plus** the DOI registry report. Run it after publishing and
    after every manuscript edit.

Non-zero exit on any failed check.

Usage
-----
    python scripts/audit_data_descriptor.py --pre-publish
    python scripts/audit_data_descriptor.py
    python scripts/audit_data_descriptor.py --md resource_paper/Data_Descriptor_v1_EN.md
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RP = ROOT / "resource_paper"
DEPOSIT = RP / "deposit"
RELEASE = ROOT / "release"
SCRIPTS = ROOT / "scripts"
REGISTRY_REPORT = ROOT / "results" / "dataset_doi_verification_registry.json"
LANGUAGE_REPORT = ROOT / "results" / "language_normalisation.json"
SYNC_REPORT = ROOT / "results" / "documentation_sync.json"
FROZEN_DEPOSIT = RP / "frozen_deposit_v1.0.0"
REFERENCE_DEPOSIT = RP / "frozen_deposit_v1.0.1"
DATASET_DOI_FILE = RP / "dataset_doi.txt"
CODE_DOI_FILE = RP / "code_doi.txt"

sys.path.insert(0, str(SCRIPTS))

TABLES = {
    "partD_resfeas_libwide_datasets.tsv": (2238, 27),
    "partD_resfeas_libwide_collections.tsv": (397, 32),
    "partD_resfeas_geo_series.tsv": (423, 21),
    "partD_resfeas_ledger.tsv": (29, 16),
}

# The two artifacts cannot appear inside MANIFEST.sha256 (a file cannot hash itself), so the
# 14-file payload always splits 12 + 2.
N_FILES = 14
N_CHECKSUMMED = 12
PAYLOAD_VERSION = "1.0.2"
FROZEN_VERSION = "1.0.0"
PREVIOUS_VERSION = "1.0.1"

results: list[tuple[bool, str, str]] = []


def check(ok: bool, label: str, detail: str = "") -> None:
    results.append((bool(ok), label, detail))
    print(f"  [{'ok' if ok else 'FAIL'}] {label}" + (f"\n        {detail}" if detail else ""))


def note(label: str, detail: str = "") -> None:
    """A deferred check: printed, explained, and deliberately not counted as a failure."""
    print(f"  [--] {label}" + (f"\n        {detail}" if detail else ""))


def section(title: str) -> None:
    print("\n" + title)
    print("-" * len(title))


def sha256_of(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ------------------------------------------------------------------ DOI resolution

def _read_doi_file(path: Path) -> dict:
    """Parse the ``<doi>`` + ``# key: value`` convention used by the two DOI files."""
    out: dict[str, str] = {"doi": "", "concept_doi": "", "version": ""}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            m = re.match(r"#\s*([a-z ]+?)\s*:\s*(\S+)", line, flags=re.I)
            if m:
                k = m.group(1).strip().lower().replace(" ", "_")
                if k in ("concept_doi", "version", "landing_page", "record_id", "published"):
                    out[k] = m.group(2)
            continue
        out["doi"] = line
    return out


def resolve_dois(args) -> dict:
    """Read every DOI the manuscript must cite from artifacts, never from a constant here.

    Hardcoding them is how a manuscript ends up citing the previous version after a
    re-release: the constant silently wins and the audit agrees with itself. Reading the
    recorded files means a version bump cannot leave the audit checking the old numbers.

    Precedence for the dataset DOI is ``--data-doi`` > the DOI the payload README itself
    carries > ``resource_paper/dataset_doi.txt``. The payload comes before the record file
    on purpose: between reserving a DOI (stage 1) and publishing (stage 5) the payload
    already cites the reserved DOI, while ``dataset_doi.txt`` still names the previous
    version — and it is the *payload* the gate has to bless.
    """
    d = _read_doi_file(DATASET_DOI_FILE)
    c = _read_doi_file(CODE_DOI_FILE)
    payload_doi = ""
    readme = DEPOSIT / "README.md"
    if readme.exists():
        code_dois = {x for x in (c["doi"], c["concept_doi"]) if x}
        found = [x for x in re.findall(r"10\.5281/zenodo\.\d{6,}",
                                      readme.read_text(encoding="utf-8"))
                 if x not in code_dois and x != d["doi"]]
        payload_doi = found[0] if found else ""
    return {
        "dataset": args.data_doi or payload_doi or d["doi"],
        "concept": args.concept_doi or d["concept_doi"],
        "code": args.code_doi or c["doi"],
        "code_concept": c["concept_doi"],
        "code_version": c["version"],
        "recorded_dataset": d["doi"],
        "recorded_dataset_version": d["version"],
    }


# --------------------------------------------------------------------- payload facts

def table_dims(path: Path) -> tuple[int, int]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        n = sum(1 for _ in rd)
    return n, len(header)


def audit_payload(dois: dict) -> dict:
    section("1. Deposit — tables, inventory, checksums")

    dims: dict[str, tuple[int, int]] = {}
    for name, expected in TABLES.items():
        p = DEPOSIT / name
        if not p.exists():
            check(False, f"{name} present", "missing from deposit/")
            continue
        got = table_dims(p)
        dims[name] = got
        check(got == expected, f"{name} is {expected[0]} x {expected[1]}",
              f"measured {got[0]} x {got[1]}")

    files = sorted(p for p in DEPOSIT.iterdir() if p.is_file())
    total = sum(p.stat().st_size for p in files)
    check(len(files) == N_FILES, f"deposit holds {N_FILES} files", f"found {len(files)}")

    manifest_p = DEPOSIT / "manifest.json"
    manifest = json.loads(manifest_p.read_text(encoding="utf-8"))
    entries = manifest.get("files", [])
    check(len(entries) == N_CHECKSUMMED,
          f"MANIFEST.sha256 covers {N_CHECKSUMMED} payload files",
          f"manifest.json lists {len(entries)}")
    check(manifest.get("dataset_doi") == dois["dataset"],
          "manifest.json quotes the dataset version DOI",
          f"manifest {manifest.get('dataset_doi')!r} vs recorded {dois['dataset']!r}")
    check(manifest.get("code_doi") == dois["code"],
          "manifest.json quotes the code DOI",
          f"manifest {manifest.get('code_doi')!r} vs recorded {dois['code']!r}")
    check(bool(manifest.get("transformations")),
          "manifest.json records the verdict_text transformation",
          f"{len(manifest.get('transformations', []))} entr(y/ies)")

    ok = bad = 0
    for f in entries:
        name = f.get("name") or f.get("path") or ""
        p = DEPOSIT / Path(name).name
        if p.exists() and sha256_of(p) == f.get("sha256"):
            ok += 1
        else:
            bad += 1
    check(bad == 0, "every manifest entry re-hashes against the file on disk",
          f"{ok} ok / {bad} bad")

    sha_all = (DEPOSIT / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines()
    sha_entries = [l for l in sha_all if l.strip() and not l.lstrip().startswith("#")]
    check(len(sha_entries) == N_CHECKSUMMED,
          f"MANIFEST.sha256 has {N_CHECKSUMMED} checksum lines",
          f"found {len(sha_entries)} entries in {len(sha_all)} lines")

    # datapackage field counts must equal the real table widths
    dp = json.loads((DEPOSIT / "datapackage.json").read_text(encoding="utf-8"))
    declared = {r.get("path"): len(r.get("schema", {}).get("fields", []))
                for r in dp.get("resources", [])}
    check(all(declared.get(n) == d[1] for n, d in dims.items()),
          "datapackage.json field counts match the tables",
          "; ".join(f"{n}: schema {declared.get(n)} vs data {d[1]}"
                    for n, d in dims.items() if declared.get(n) != d[1]) or "all agree")

    # ledger-derived counts used in Technical Validation / Usage Notes
    with (DEPOSIT / "partD_resfeas_ledger.tsv").open(encoding="utf-8-sig", newline="") as fh:
        led = list(csv.DictReader(fh, delimiter="\t"))
    n_usable = sum(1 for r in led if str(r["USABLE"]).strip().lower() == "true")
    n_fail = len(led) - n_usable
    check((n_usable, n_fail) == (2, 27),
          "ledger splits 29 rows into 2 USABLE + 27 failing at least one gate",
          f"USABLE={n_usable} failing={n_fail}")
    gse = [r for r in led if r["resource"].upper().startswith("GSE")]
    bulk_usable = [r for r in gse if str(r["USABLE"]).strip().lower() == "true"]
    check(bulk_usable == [r for r in led if r["resource"] == "GSE133060"],
          "GSE133060 is the only GEO (bulk) candidate that is USABLE",
          "; ".join(r["resource"] for r in bulk_usable) or "none")

    with (DEPOSIT / "partD_resfeas_geo_series.tsv").open(encoding="utf-8-sig", newline="") as fh:
        geo = {r["gse"]: r for r in csv.DictReader(fh, delimiter="\t")}
    row = geo.get("GSE133060", {})
    check(row.get("n_paired_donors") == "8" and row.get("n_paired_same_site") == "8"
          and row.get("molecule") == "mRNA",
          "GSE133060 row supports 8 same-donor, same-site pairs on an mRNA assay",
          f"n_paired_donors={row.get('n_paired_donors')} "
          f"n_paired_same_site={row.get('n_paired_same_site')} molecule={row.get('molecule')}")

    with (DEPOSIT / "partD_resfeas_libwide_collections.tsv").open(encoding="utf-8-sig",
                                                                 newline="") as fh:
        coll = list(csv.DictReader(fh, delimiter="\t"))
    n_contrast = sum(1 for r in coll if str(r["has_contrast"]).strip().lower() == "true")
    n_sep = sum(1 for r in coll if str(r["has_arm_contrast"]).strip().lower() == "true")
    check((n_contrast, n_sep) == (140, 15),
          "140 collections carry a contrast, 15 have metadata-separable arms",
          f"has_contrast={n_contrast} has_arm_contrast={n_sep} "
          f"({n_sep / n_contrast * 100:.1f}%, undetermined {n_contrast - n_sep})")

    summary = json.loads(
        (DEPOSIT / "partD_resfeas_libwide_summary.json").read_text(encoding="utf-8"))
    check(summary.get("cells_total") == 292_800_864.0 and summary.get("n_datasets_total") == 2238
          and summary.get("n_collections_total") == 397,
          "summary.json headline counts are 397 / 2,238 / 292,800,864",
          f"{summary.get('n_collections_total')} / {summary.get('n_datasets_total')} / "
          f"{summary.get('cells_total'):,.0f}")

    return {"n_usable": n_usable, "n_fail": n_fail,
            "n_contrast": n_contrast, "n_sep": n_sep, "dims": dims,
            "n_files": len(files), "total_bytes": total}


# ------------------------------------------------- pre-publish internal consistency

def audit_internal(dois: dict, facts: dict) -> None:
    """Checks that only make sense *before* publication — the release is not public yet."""
    section(f"1b. Release integrity — version {PAYLOAD_VERSION}, no placeholders")

    manifest = json.loads((DEPOSIT / "manifest.json").read_text(encoding="utf-8"))
    dp = json.loads((DEPOSIT / "datapackage.json").read_text(encoding="utf-8"))

    check(manifest.get("version") == PAYLOAD_VERSION,
          f"manifest.json version is {PAYLOAD_VERSION}", str(manifest.get("version")))
    check(dp.get("version") == PAYLOAD_VERSION,
          f"datapackage.json version is {PAYLOAD_VERSION}", str(dp.get("version")))

    readme = (DEPOSIT / "README.md").read_text(encoding="utf-8")
    changelog = (DEPOSIT / "CHANGELOG.md").read_text(encoding="utf-8")
    check(f"Version {PAYLOAD_VERSION}" in readme,
          f"deposit README states Version {PAYLOAD_VERSION}")
    check(f"## [{PAYLOAD_VERSION}]" in changelog,
          f"CHANGELOG.md carries a [{PAYLOAD_VERSION}] section")
    check(f"## [{PREVIOUS_VERSION}]" in changelog,
          f"CHANGELOG.md keeps the [{PREVIOUS_VERSION}] entry (history, not overwritten)")
    check(f"## [{FROZEN_VERSION}]" in changelog,
          f"CHANGELOG.md keeps the [{FROZEN_VERSION}] entry (history, not overwritten)")
    check("synchronised with the companion Data Descriptor" in changelog
          and "gate outcome, numeric value, field name, or row membership changed" in changelog,
          f"CHANGELOG.md states the v{PAYLOAD_VERSION} scope in the required words",
          "looked for the data-dictionary-synchronisation clause and the no-change clause")
    check("English-language normalisation only" in changelog,
          f"CHANGELOG.md keeps the historical v{PREVIOUS_VERSION} scope statement")

    # --- the code record this data release was built from ---------------------------
    check(bool(dois["code"]), "a code DOI is recorded in resource_paper/code_doi.txt",
          dois["code"] or "empty")
    check(bool(dois["code_version"]),
          "resource_paper/code_doi.txt records the code version it points at",
          dois["code_version"] or "missing '# version:' line")
    code_version_ok = dois["code_version"] == f"v{PAYLOAD_VERSION}"
    check(code_version_ok,
          f"the cross-linked code record self-describes as v{PAYLOAD_VERSION}",
          f"code_doi.txt says {dois['code_version']!r}, the payload is v{PAYLOAD_VERSION}"
          + ("" if code_version_ok else
             ". Re-release the code first, then set its DOI and version here — otherwise "
             f"the published record silently cites the v{FROZEN_VERSION} code"))
    check(dois["code"] in readme,
          "deposit README cross-links the recorded code version DOI", dois["code"])
    check(dois["code"] in manifest.get("code_doi", ""),
          "deposit manifest.json cross-links the recorded code version DOI",
          str(manifest.get("code_doi")))
    if dois["code_concept"]:
        check(dois["code"] != dois["code_concept"],
              "the code version DOI differs from the code concept DOI",
              f"{dois['code']} / {dois['code_concept']}")

    # --- the data DOI that will be minted -------------------------------------------
    placeholders = []
    for p in sorted(DEPOSIT.iterdir()):
        if p.is_file() and p.suffix.lower() in {".md", ".json", ".txt", ".sha256"}:
            txt = p.read_text(encoding="utf-8", errors="replace")
            if "XXXXXXX" in txt or "YYYYYYY" in txt or "TBD" in txt:
                placeholders.append(p.name)
    check(not placeholders,
          "no DOI placeholder left anywhere in the payload (stage 2 backfilled it)",
          f"still holding XXXXXXX/YYYYYYY/TBD: {placeholders}")
    lic = (DEPOSIT / "LICENSE-DATA.md").read_text(encoding="utf-8")
    code_dois = {x for x in (dois["code"], dois["code_concept"]) if x}
    readme_dois = set(re.findall(r"10\.5281/zenodo\.\d{6,}", readme)) - code_dois
    lic_dois = set(re.findall(r"10\.5281/zenodo\.\d{6,}", lic)) - code_dois
    check(readme_dois == lic_dois == {dois["dataset"]},
          "README and LICENSE-DATA.md name exactly one dataset DOI: the one being published",
          f"README {sorted(readme_dois)}; LICENSE {sorted(lic_dois)}; "
          f"expected {{{dois['dataset']!r}}} — a superseded DOI left in either file is the "
          f"classic re-release slip")
    check(manifest.get("dataset_doi") in readme_dois,
          "manifest.json dataset_doi is the DOI the payload cites",
          f"manifest {manifest.get('dataset_doi')!r} vs payload {sorted(readme_dois)}")
    if dois["recorded_dataset"] and dois["recorded_dataset"] != dois["dataset"]:
        note(f"the reserved dataset DOI {dois['dataset']} is not yet in "
             f"resource_paper/dataset_doi.txt ({dois['recorded_dataset']})",
             "expected between stages 1 and 5; publish_dataset.py rewrites that file on "
             "publication, after which the full audit checks the two agree")

    # --- no Han characters anywhere in the payload ----------------------------------
    cjk = re.compile(r"[\u4e00-\u9fff]")
    han = []
    for p in sorted(DEPOSIT.iterdir()):
        if p.is_file() and p.suffix.lower() in {".md", ".json", ".tsv", ".txt", ".cff"}:
            if cjk.search(p.read_text(encoding="utf-8", errors="replace")):
                han.append(p.name)
    check(not han, "the payload contains no Han characters (English-only deposit)",
          f"still Chinese: {han}")

    # --- the machine-checked documentation-only claim --------------------------------
    if REFERENCE_DEPOSIT.is_dir():
        r = subprocess.run(
            [sys.executable, str(SCRIPTS / "test_documentation_sync.py"),
             "--json-out", str(SYNC_REPORT)],
            cwd=str(ROOT), capture_output=True, text=True)
        verdict = (SYNC_REPORT.read_text(encoding="utf-8")
                   if SYNC_REPORT.exists() else "{}")
        rep = json.loads(verdict) if verdict.strip().startswith("{") else {}
        check(r.returncode == 0 and rep.get("failed") == 0,
              f"documentation-only regression vs frozen v{PREVIOUS_VERSION} passes every check",
              f"exit {r.returncode}, {rep.get('total')} checks, {rep.get('failed')} failed")
        check(rep.get("total", 0) >= 40,
              "the regression actually ran (expected >= 40 checks)", str(rep.get("total")))
    else:
        check(False, f"frozen v{PREVIOUS_VERSION} reference payload present",
              f"missing {REFERENCE_DEPOSIT.relative_to(ROOT)} — the documentation-only claim "
              f"cannot be verified without it")


def audit_size_in_manuscript(md: str, facts: dict) -> None:
    """The stated package size must be the measured one, not a remembered one."""
    section("5b. Manuscript — package size re-measured")
    m = re.search(r"contains\s+(\d+)\s+files\s*\(([\d,]+)\s+bytes in total\)", md)
    if not m:
        check(False, "manuscript states the package inventory as 'N files (X bytes in total)'",
              "pattern not found — the stated size would otherwise go unchecked")
        return
    check(int(m.group(1)) == facts["n_files"],
          f"manuscript states {facts['n_files']} files", f"stated {m.group(1)}")
    stated = int(m.group(2).replace(",", ""))
    check(stated == facts["total_bytes"],
          f"manuscript states the measured package size {facts['total_bytes']:,} bytes",
          f"stated {stated:,} vs measured {facts['total_bytes']:,} "
          f"(delta {stated - facts['total_bytes']:+,})")


# ------------------------------------------------------------------- manuscript facts

def published_names() -> set[str]:
    """Every filename the manuscript is allowed to cite: the deposit + the code archive."""
    names = {p.name for p in DEPOSIT.iterdir() if p.is_file()}
    if RELEASE.is_dir():
        try:
            out = subprocess.run(["git", "-C", str(RELEASE), "ls-files"],
                                 capture_output=True, text=True, check=True).stdout
            for line in out.splitlines():
                names.add(Path(line).name)
        except Exception:
            pass
    return names


def audit_manuscript(path: Path, facts: dict, dois: dict) -> None:
    md = path.read_text(encoding="utf-8")
    body = md.split("## References", 1)[0]
    body_prose = "\n".join(l for l in body.splitlines() if not l.startswith("|"))

    section("2. Manuscript — citations")
    refs = re.findall(r"^\s*(\d+)\.\s", md.split("## References", 1)[1], flags=re.M)
    intext = sorted({int(n) for n in re.findall(r"\[(\d+)(?:\s*[-,]\s*\d+)*\]", body_prose)})
    check(len(refs) == 3, "reference list holds 3 entries", f"found {len(refs)}")
    check(intext == [1, 2, 3], "every reference is cited at least once in the text",
          f"in-text markers found: {intext or 'NONE'}")
    check("meta-analysis" not in md.lower(),
          "genre red line — the manuscript never calls itself a meta-analysis")

    section("3. Manuscript — DOIs")
    check(md.count(dois["dataset"]) >= 2, "version DOI appears repeatedly",
          f"{md.count(dois['dataset'])} occurrence(s) of {dois['dataset']}")
    check(dois["concept"] in md, "concept DOI is stated",
          f"{dois['concept']} ({md.count(dois['concept'])} occurrence(s))")
    check(dois["code"] in md, "code DOI is stated",
          f"{dois['code']} ({md.count(dois['code'])} occurrence(s))")
    trio = (dois["dataset"], dois["concept"], dois["code"])
    check(all(d in md for d in trio) and len(set(trio)) == 3,
          "the three DOIs are distinct and all present")
    check("XXXXXXX" not in md and "TBD" not in md, "no placeholder DOIs remain")
    stale = [d for d in ("10.5281/zenodo.23297703", "10.5281/zenodo.23291759")
             if d in md and d not in trio]
    check(not stale, "manuscript cites no superseded version DOI",
          f"previous-version DOI(s) still present: {stale}")
    check(f"version {PAYLOAD_VERSION}" in md,
          f"manuscript states 'version {PAYLOAD_VERSION}'")

    section("4. Manuscript — cited files are actually published")
    allowed = published_names()
    cited = sorted(set(re.findall(r"`([A-Za-z0-9_.\-]+\.(?:tsv|json|zip|md|cff|sha256))`", md)))
    unknown = [c for c in cited if c not in allowed]
    check(not unknown, f"all {len(cited)} cited filenames exist in the deposit or the code archive",
          f"not found anywhere public: {unknown}" if unknown else
          "; ".join(cited))
    for forbidden in ("PROVENANCE_MAP.tsv", "VERSIONING.md"):
        check(forbidden not in md, f"manuscript does not cite the unpublished {forbidden}")

    section("5. Manuscript — numeric claims re-derived from the payload")
    for name, (nrows, ncols) in facts["dims"].items():
        pat = re.compile(re.escape(name))
        for line in md.splitlines():
            if pat.search(line) and "|" in line:
                check(f"{nrows:,} x {ncols}".replace(",", ",") in line
                      or f"{nrows} x {ncols}" in line,
                      f"Table 1 row for {name} states {nrows} x {ncols}", line.strip()[:110])
    for token, label in [
        ("292,800,864", "cells_total"),
        ("2,238", "datasets"),
        ("397", "collections"),
        ("423", "GEO series parsed"),
        ("521", "GEO series matched by the queries"),
        ("98", "GEO series not parsed"),
        (f"{facts['n_files']} files", "package inventory"),
        (f"{facts['total_bytes']:,}", "package size in bytes"),
        (str(N_CHECKSUMMED), "checksum coverage"),
        (str(facts["n_usable"]), "USABLE candidates"),
        (str(facts["n_fail"]), "failing candidates"),
        ("140", "contrast-bearing collections"),
        ("15", "arm-separable collections"),
        (str(facts["n_contrast"] - facts["n_sep"]), "undetermined collections"),
        ("10.7%", "arm-separable share"),
    ]:
        check(token in md, f"manuscript states {label} ({token})")

    audit_size_in_manuscript(md, facts)

    section("6. Abstract budget")
    abs_txt = md.split("## Abstract", 1)[1].split("## Background", 1)[0]
    n_words = len(re.findall(r"[A-Za-z0-9][A-Za-z0-9.\-/]*", abs_txt))
    check(n_words <= 150, f"abstract is {n_words} words (Scientific Data guidance ~150)",
          f"{n_words} words")

    section("7. Gate thresholds quoted from the generator")
    gen = SCRIPTS / "partD_resfeas_ledger.py"
    if gen.exists():
        src = gen.read_text(encoding="utf-8")
        floor = {}
        for key in ("PAIR_MIN", "CTRL_MIN"):
            m = re.search(rf"^{key}\s*=\s*(\d+)", src, flags=re.M)
            floor[key] = int(m.group(1)) if m else None
        words = {6: "six", 5: "five", 4: "four", 7: "seven", 8: "eight"}
        if floor["PAIR_MIN"] and floor["CTRL_MIN"]:
            check(words[floor["PAIR_MIN"]] in md and words[floor["CTRL_MIN"]] in md
                  and "90%" in md,
                  "manuscript states the gate-3 thresholds "
                  f"({floor['PAIR_MIN']} pairs / {floor['CTRL_MIN']} controls / 90% primary)",
                  f"partD_resfeas_ledger.py: PAIR_MIN={floor['PAIR_MIN']} "
                  f"CTRL_MIN={floor['CTRL_MIN']}")
        else:
            check(False, "could parse PAIR_MIN / CTRL_MIN from partD_resfeas_ledger.py",
                  str(floor))
    else:
        check(False, "partD_resfeas_ledger.py present", "generator script not found")


def audit_registry(dois: dict) -> None:
    section("8. Published record — recorded DOI and the DataCite/doi.org channel")
    check(dois["recorded_dataset"] == dois["dataset"],
          "resource_paper/dataset_doi.txt names the DOI the payload cites",
          f"recorded {dois['recorded_dataset']!r} vs payload {dois['dataset']!r}")
    check(dois["recorded_dataset_version"] == PAYLOAD_VERSION,
          f"resource_paper/dataset_doi.txt records version {PAYLOAD_VERSION}",
          f"recorded {dois['recorded_dataset_version']!r} — if this is the previous "
          f"version, the record was not republished and the manuscript would cite bytes "
          f"that are not the ones deposited")
    if not REGISTRY_REPORT.exists():
        check(False, "registry verification report present",
              f"missing {REGISTRY_REPORT.relative_to(ROOT)} — run "
              f"scripts/verify_dataset_doi_registry.py first")
        return
    rep = json.loads(REGISTRY_REPORT.read_text(encoding="utf-8"))
    checks = rep.get("checks") or {}
    if isinstance(checks, dict):
        items = [(k, bool(v.get("pass"))) for k, v in checks.items() if isinstance(v, dict)]
    else:  # tolerate a list-shaped report from a future revision
        items = [(str(c.get("name", i)), bool(c.get("pass", c.get("ok"))))
                 for i, c in enumerate(checks)]
    failed = [name for name, ok in items if not ok]
    check(rep.get("doi") == dois["dataset"], "registry report is for the expected DOI",
          f"report {rep.get('doi')!r} vs recorded {dois['dataset']!r}")
    check(bool(items) and not failed,
          "registry report shows every check passing",
          f"{len(items)} checks ({', '.join(n for n, _ in items)}), {len(failed)} failing")
    if rep.get("failures"):
        check(False, "registry report records no failures", str(rep["failures"])[:160])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--md", default=str(RP / "Data_Descriptor_v1_EN.md"))
    ap.add_argument("--pre-publish", action="store_true",
                    help="payload-only hard gate for publish_dataset.py stage 0: skip the "
                         "manuscript and registry checks, which cannot pass before publication")
    ap.add_argument("--data-doi", help="override the dataset version DOI")
    ap.add_argument("--concept-doi", help="override the dataset concept DOI")
    ap.add_argument("--code-doi", help="override the code version DOI")
    args = ap.parse_args()

    dois = resolve_dois(args)
    mode = "pre-publish (payload only)" if args.pre_publish else "full"
    print(f"mode     : {mode}")
    print(f"auditing : {Path(args.md).relative_to(ROOT) if Path(args.md).is_absolute() else args.md}")
    print(f"deposit  : {DEPOSIT}")
    print(f"DOIs     : dataset={dois['dataset'] or '(none)'}  "
          f"concept={dois['concept'] or '(none)'}  code={dois['code'] or '(none)'}")

    facts = audit_payload(dois)
    audit_internal(dois, facts)
    if args.pre_publish:
        note("manuscript cross-check deferred to the full audit (pre-publication)")
        note("DOI registry report deferred to the full audit (nothing is public yet)")
    else:
        audit_manuscript(Path(args.md), facts, dois)
        audit_registry(dois)

    bad = [r for r in results if not r[0]]
    print("\n" + "=" * 72)
    if bad:
        print(f"RESULT: {len(bad)} of {len(results)} checks FAILED")
        for _, label, detail in bad:
            print(f"  - {label}" + (f"  ({detail})" if detail else ""))
        return 1
    print(f"RESULT: all {len(results)} checks passed"
          + (" — the payload is internally consistent and may be published"
             if args.pre_publish else " — the manuscript matches the frozen payload"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
