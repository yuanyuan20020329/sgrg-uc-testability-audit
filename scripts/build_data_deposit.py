#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build the frozen deposit of the audit *tables* (the dataset DOI record).

Current payload version: **v1.0.1**, which differs from the immutable v1.0.0 release
(`10.5281/zenodo.23297703`) only in the language of its free-text fields — see
`CHANGELOG.md` in the deposit and `scripts/test_language_normalisation.py`, which asserts
that difference file by file against `resource_paper/frozen_deposit_v1.0.0/`.

This is deliberately separate from ``scripts/build_archive.py`` (the *code* record; currently
v1.0.1, cited from ``resource_paper/code_doi.txt``). Scientific Data requires the dataset itself to
be deposited in a public repository with its own persistent identifier, so the two
records are minted and cited independently. Never merge them.

Output layout (``resource_paper/deposit/``):

    README.md                                   deposit-facing, resource-only
    LICENSE-DATA.md                             CC BY 4.0 (tables)
    CHANGELOG.md                                Added / Changed / Fixed / Source snapshot
    SCOPE_NOTE.md                               scope guard for the interpretive field
    datapackage.json                            Frictionless schema of the four tables
    MANIFEST.sha256                             sha256 over every redistributed file
    <four tables> + <two summary/verdict JSONs>
    cxg_coll_detail_snapshot_v1.0.0.zip         frozen source snapshot (397 JSON, offline re-run;
                                                name keeps the v1.0.0 tag on purpose — the bytes are
                                                unchanged between releases, so the v1.0.0 tag
                                                records the snapshot's own era, not the payload's)

Plus, next to the deposit directory:
    resource_paper/testability-audit-tables-v1.0.1.zip   single-unit download

One documented transformation is applied on the way in: the interpretive free-text field
`verdict_text` in `partD_resfeas_verdict.json` is replaced by a counts-only annotation in
the public payload. The original is retained in `results/` and in the code archive, and both
hashes are recorded in `manifest.json` under `transformations`. No other value is altered.

Usage
-----
    python scripts/build_data_deposit.py
    python scripts/build_data_deposit.py --no-snapshot     # skip the 12 MB snapshot zip
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import zipfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RP = ROOT / "resource_paper"
DEPOSIT = RP / "deposit"
RESULTS = ROOT / "results"
SNAPSHOT_SRC = ROOT / "data" / "raw" / "cxg_coll_detail"

VERSION = "1.0.1"
RELEASE_DATE = date(2026, 10, 11).isoformat()
# v1.0.0 shipped on the same calendar day; the two releases are kept distinct by version
# DOI, not by date, so both carry today's date and the changelog lists both.
FIRST_RELEASE_DATE = "2026-10-11"
SNAPSHOT_FROZEN = "2026-10-09"

# The companion **code** record. v1.0.1 re-releases the code too (the generators changed),
# so the code DOI is read from `resource_paper/code_doi.txt` rather than baked in here: the
# data record must cross-link to the *exact* code version it was built from, and the
# manuscript has to cite the same DOI. One file, one source of truth.
CODE_DOI_FILE = RP / "code_doi.txt"
CODE_DOI_V1_0_0 = "10.5281/zenodo.23291759"          # the v1.0.0 code record, still immutable
DATA_DOI_PLACEHOLDER = "10.5281/zenodo.XXXXXXX"
# A second, parallel token for the code record. `LICENSE-DATA.md` is a template that is
# *rendered* into the deposit, so both DOIs have to be substituted at build time: a real DOI
# written into that file once would silently survive every later re-release.
CODE_DOI_PLACEHOLDER = "10.5281/zenodo.YYYYYYY"

# A zip stores each member's mtime. The convenience zip is rebuilt on every run, and its
# members are the freshly written deposit files, so inheriting "now" made every rebuild
# produce different bytes. Pin it to the release date instead. This zip is **not** uploaded
# to Zenodo, so the deposit payload and its manifest are unaffected either way.
ZIP_DATE_TIME = tuple(int(x) for x in RELEASE_DATE.split("-")) + (0, 0, 0)

# --- files redistributed flat into the deposit (name in deposit -> source) -------
FLAT_FILES: list[tuple[str, Path]] = [
    ("partD_resfeas_libwide_datasets.tsv", RESULTS / "partD_resfeas_libwide_datasets.tsv"),
    ("partD_resfeas_libwide_collections.tsv", RESULTS / "partD_resfeas_libwide_collections.tsv"),
    ("partD_resfeas_geo_series.tsv", RESULTS / "partD_resfeas_geo_series.tsv"),
    ("partD_resfeas_ledger.tsv", RESULTS / "partD_resfeas_ledger.tsv"),
    ("partD_resfeas_libwide_summary.json", RESULTS / "partD_resfeas_libwide_summary.json"),
    ("partD_resfeas_verdict.json", RESULTS / "partD_resfeas_verdict.json"),
    ("datapackage.json", RP / "datapackage.json"),
]

# Rendered rather than copied, so the DOI can be substituted on the way in.
LICENSE_SRC = RP / "LICENSE-DATA.md"

# One field in one table is interpretive rather than descriptive: `verdict_text` in the
# ledger verdict file states a conclusion. A Data Descriptor published by Scientific Data
# must not carry hypothesis-level claims, and reviewers read the deposit itself — a scope
# note does not neutralise a shipped sentence. So the public payload carries a counts-only
# annotation instead, and the transformation is recorded in the manifest. The original text
# is NOT lost: it stays byte-identical in `results/` and in the companion code archive
# (the code DOI recorded in `resource_paper/code_doi.txt`), the internal provenance record.
NEUTRALISED = {"partD_resfeas_verdict.json": "verdict_text"}

NEUTRAL_VERDICT_TEXT = (
    "Ledger summary - counts only. The candidate resources screened by this audit are the "
    "rows of the ledger table; the totals above count those candidates and the subset that "
    "passes all four gates (G1_mRNA, G2_UC, G3_powered_control, G4_in_vivo). Sources: "
    "CELLxGENE Discover collections and NCBI GEO series. Audit status: each candidate's "
    "per-gate outcome is recorded in the gate columns and its failure reason in the "
    "`reason` field of the ledger table. This field carries no biological interpretation. "
    "It replaces a free-text note written during the audit; the original text is retained "
    "in the internal provenance record shipped with the analysis code (see CHANGELOG.md "
    "and SCOPE_NOTE.md)."
)

SNAPSHOT_ZIP_NAME = "cxg_coll_detail_snapshot_v1.0.0.zip"
SINGLE_UNIT_ZIP = RP / f"testability-audit-tables-v{ VERSION }.zip"

TEXT_SUFFIXES = {".md", ".json", ".tsv", ".txt", ".cff"}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_text(path: Path, text: str) -> None:
    """Write text with LF newlines and a UTF-8 BOM-free encoding."""
    path.write_text(text.replace("\r\n", "\n"), encoding="utf-8", newline="\n")


def read_code_doi_record(override: str | None = None) -> dict:
    """The companion **code** record this release cross-links to.

    Returns ``{"doi", "concept_doi", "version"}``. Priority: ``--code-doi`` >
    ``resource_paper/code_doi.txt`` > the v1.0.0 code DOI. The file outranks the constant
    on purpose: re-releasing the code becomes a one-line edit to a data file instead of a
    code change, and the deposit, the metadata block, the manuscript and
    ``audit_data_descriptor.py`` can all be pointed at one string.

    ``version`` matters as much as the DOI: linking a v1.0.1 data record to a code record
    that still self-describes as v1.0.0 is exactly the kind of silent mismatch this file
    exists to make impossible, and ``audit_data_descriptor.py`` fails the build on it.
    """
    if override:
        return {"doi": override.strip(), "concept_doi": None, "version": None}
    if CODE_DOI_FILE.exists():
        rec = {"doi": None, "concept_doi": None, "version": None}
        for line in CODE_DOI_FILE.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                m = re.match(r"#\s*(concept doi|version)\s*:\s*(\S+)", line, flags=re.I)
                if m:
                    key = "concept_doi" if m.group(1).lower().startswith("concept") else "version"
                    rec[key] = m.group(2)
                continue
            rec["doi"] = line
        if rec["doi"]:
            return rec
    return {"doi": CODE_DOI_V1_0_0, "concept_doi": "10.5281/zenodo.23291758", "version": "v1.0.0"}


def read_code_doi(override: str | None = None) -> str:
    """Just the code version DOI (see :func:`read_code_doi_record`)."""
    return read_code_doi_record(override)["doi"]


def neutralise_field(src: Path, dst: Path, field: str, code_doi: str) -> dict:
    """Copy ``src`` to ``dst`` with one interpretive field replaced by counts-only text.

    Returns the transformation record that goes into ``manifest.json`` so a reader can
    verify the change instead of taking it on trust. ``src`` is never modified.

    ``code_doi`` is the **current** code version DOI, not the one from the release that
    first recorded the original text: the record has to point a reader at the archive that
    ships alongside this payload, otherwise a re-release leaves the deposit citing a code
    version that no longer matches it.
    """
    raw = src.read_bytes()
    obj = json.loads(raw.decode("utf-8"))
    original = obj.get(field)
    if not isinstance(original, str):
        raise SystemExit(f"{src.name}: no string field '{field}' to neutralise")
    obj[field] = NEUTRAL_VERDICT_TEXT
    text = json.dumps(obj, indent=1, ensure_ascii=False) + "\n"
    write_text(dst, text)
    return {
        "file": dst.name,
        "field": field,
        "action": "replaced in the public payload",
        "reason": ("the original free-text annotation states an interpretive conclusion; the "
                   "public dataset payload carries counts, sources and audit status only"),
        "original_source": str(src.relative_to(ROOT)).replace("\\", "/"),
        "original_file_sha256": hashlib.sha256(raw).hexdigest(),
        "original_field_sha256": hashlib.sha256(original.encode("utf-8")).hexdigest(),
        "deposited_file_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "original_retained_in": ("internal provenance record: results/ and the code archive, "
                                 f"Zenodo {code_doi}"),
    }


def build_readme(n_rows: dict[str, tuple[int, int]], cells_total: int, data_doi: str,
                 code_doi: str) -> str:
    return f"""# Audit resource — public single-cell and bulk records for ulcerative colitis

**Version {VERSION}** · released {RELEASE_DATE} · data licence **CC BY 4.0**

Derived, citable annotation of public single-cell and bulk resources. For each
CELLxGENE Discover collection and dataset the tables record *where the cells came
from*, *whether control and disease arms can be separated in the metadata layer*,
*a sequencing-depth proxy*, and the outcome of a four-gate usability ledger.

This deposit contains **derived metadata records only**. It redistributes **no
primary omics data**. Every source collection and dataset stays under its own
depositor's licence; accessions and download scripts live in the companion code
repository.

## Files

| File | Rows × cols | Content |
|---|---|---|
| `partD_resfeas_libwide_datasets.tsv` | {n_rows['libwide_datasets'][0]} × {n_rows['libwide_datasets'][1]} | One row per CELLxGENE dataset (all 397 collections). Primary provenance unit. |
| `partD_resfeas_libwide_collections.tsv` | {n_rows['libwide_collections'][0]} × {n_rows['libwide_collections'][1]} | One row per collection: reuse fractions, arm separability, reuse asymmetry, depth ratio. |
| `partD_resfeas_geo_series.tsv` | {n_rows['geo_series'][0]} × {n_rows['geo_series'][1]} | One row per screened GEO series: molecular gate, in-vitro flag, same-donor pairing counts. |
| `partD_resfeas_ledger.tsv` | {n_rows['gate_ledger'][0]} × {n_rows['gate_ledger'][1]} | The four-gate ledger with a per-row human-readable `reason` string. |
| `datapackage.json` | — | Frictionless Data Package: field names, types, enum constraints, semantics for all four tables. |
| `MANIFEST.sha256` | — | SHA-256 of every redistributed file, so any published count can be re-checked without re-running the pipeline. |
| `{SNAPSHOT_ZIP_NAME}` | 397 files | Frozen CELLxGENE collection-detail snapshot ({SNAPSHOT_FROZEN}). The pipeline re-runs against this offline, so a release reproduces without depending on the live API. |

## What the coverage looks like

`cells_total` = {cells_total:,} cells across {n_rows['libwide_collections'][0]} collections and
{n_rows['libwide_datasets'][0]} datasets. These are properties of the tables, not findings:

* `arm_class` ∈ {{`ctl_only`, `dis_only`, `mixed`}}. A `mixed` dataset carries both control
  and disease labels inside **the same dataset**, so its arms cannot be split at the
  metadata layer. Counts: `ctl_only` 1,537 / `dis_only` 355 / `mixed` 346.
* 140 collections carry a case–control contrast; **15 (10.7%)** have metadata-separable
  arms. The other 125 are recorded as *undetermined*.

## The rule that must not be lost

**Not auditable is not the same as clean.** Where arms cannot be separated, the correct
record is *cannot be determined* — never *passes*. Any downstream use of these tables
must preserve that distinction.

## One field is counts-only

`partD_resfeas_verdict.json` originally carried a free-text annotation that stated a
conclusion. In this public payload that single field is replaced with a **counts-only
annotation** — counts, sources and audit status. The original text is retained in the
internal provenance record shipped with the code archive (Zenodo {code_doi}), and both
hashes are recorded under `transformations` in `manifest.json`. Every table otherwise
copies through unchanged. See `CHANGELOG.md` and `SCOPE_NOTE.md`.

## Version history

| Version | Released | What changed |
|---|---|---|
| {VERSION} | {RELEASE_DATE} | **English-language normalisation only.** The free-text fields that were written in Chinese in v1.0.0 — the per-row `reason` string of `partD_resfeas_ledger.tsv`, the `metric_definitions` glossary and the `caveat` note of `partD_resfeas_libwide_summary.json`, one example string in `partD_resfeas_verdict.json`, and two prose field descriptions in `datapackage.json` — are now English. No gate outcome, numeric value, field semantics or row membership changed: field names, table structure, column order, the four gates, `USABLE` and every number are identical to v1.0.0. |
| 1.0.0 | {FIRST_RELEASE_DATE} | First public release. Retained, unaltered, at its own version DOI. |

## Citation

Yuan, Y., Chen, H. *Testability audit of public single-cell and bulk resources for
ulcerative colitis: collection- and dataset-level provenance, arm separability and gate
ledger.* Zenodo. https://doi.org/{data_doi}

See `LICENSE-DATA.md` for attribution requirements.

## Related records

* **Code archive** (analysis code + key result tables, MIT): {code_doi}
  — a separate, independently citable record. The tables here are the dataset; the code
  record is the software. This release of the data was built from that code version.
* Source data: CELLxGENE Discover (curation API v1) and NCBI GEO. See the code
  repository's `DATA_ACCESS.md`.

_Released {RELEASE_DATE} as v{VERSION}; supersedes v1.0.0 ({FIRST_RELEASE_DATE}), which
stays available and immutable under its own version DOI._
"""


def build_changelog(code_doi: str) -> str:
    return f"""# Changelog

All notable changes to the audit resource are documented here. This file follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the resource uses semantic
versioning, and every released version keeps its own version DOI while the concept DOI
always resolves to the latest one. `README.md` carries the same history as a table.

## [{VERSION}] — {RELEASE_DATE}

### Changed
- **English-language normalisation only; no gate outcome, numeric value, field semantics,
  or row membership changed.** In v1.0.0 the free-text fields were written in Chinese.
  They are now written in English, because the deposit is cited by an English-language
  Data Descriptor and a reader must be able to see *why a given resource failed* without a
  translation layer. Specifically:
  - `partD_resfeas_ledger.tsv` — all 29 values of the per-row `reason` column.
  - `partD_resfeas_libwide_summary.json` — the `metric_definitions` field glossary
    (13 entries) and the `caveat` note.
  - `partD_resfeas_verdict.json` — one example string (`geo_max_pairs_mirna_example`).
  - `datapackage.json` — two prose field descriptions (the `geo_series` resource
    description and the `reason` field example) and the `version` field.
  - `README.md` and this file — reissued for the new version.
- What did **not** change, item by item: field names, table structure and column order,
  row count and row order, the four gates (`G1_mRNA`, `G2_UC`, `G3_powered_control`,
  `G4_in_vivo`), `USABLE`, every numeric value, every boolean, every accession, the
  `verdict_text` counts-only transformation, and every byte of the frozen CELLxGENE
  snapshot and of the three provenance-free tables. This is asserted, not merely intended:
  `scripts/test_language_normalisation.py` diffs this payload against the frozen v1.0.0
  reference file by file and allows only the fields listed above to differ.
- The companion **code archive** is re-released as v1.0.1 (`{code_doi}`) because the
  generator scripts changed; this data record cites that exact code version DOI as a
  related identifier. v1.0.0 of both records stays published and immutable.

### Added
- Nothing new is deposited. `results/partD_resfeas_ledger_reason_i18n.tsv` (in the code
  archive, not here) records the v1.0.0 → v1.0.1 reason mapping row by row.

### Fixed
- Nothing semantic.

### Source snapshot
- Unchanged: CELLxGENE Discover curation API v1, collection-detail endpoint, frozen
  **{SNAPSHOT_FROZEN}** (397 JSON files). Shipped as
  `{SNAPSHOT_ZIP_NAME}`; the pipeline re-runs against it offline with `--no-fetch`. The
  file keeps its v1.0.0 name because the snapshot itself is the v1.0.0 snapshot.

## [1.0.0] — {FIRST_RELEASE_DATE}

### Added
- First public release of the audit tables, frozen and deposited under their own
  persistent identifier (distinct from the code archive DOI).
- `partD_resfeas_libwide_datasets.tsv` — dataset-level provenance across all 397
  CELLxGENE Discover collections (2,238 datasets, 292,800,864 cells).
- `partD_resfeas_libwide_collections.tsv` — collection-level aggregates plus arm
  typology (`arm_class`), arm separability (`has_arm_contrast` / `has_mixed_contrast`),
  reuse asymmetry and the depth ratio.
- `partD_resfeas_geo_series.tsv` — 423 screened GEO series with the molecular gate
  derived from `!Series_type`, the in-vitro flag and same-donor pairing counts.
- `partD_resfeas_ledger.tsv` — the four-gate usability ledger (29 candidates) with a
  per-row `reason` string.
- `datapackage.json` — machine-readable Frictionless schema for all four tables.
- `MANIFEST.sha256` — SHA-256 of every redistributed file.
- `SCOPE_NOTE.md` — explains the one field that is counts-only in this payload, and where
  the interpretation belongs instead.

### Changed
- `partD_resfeas_verdict.json` — the free-text field `verdict_text` is **replaced with a
  counts-only annotation** in this public payload. The original text stated an interpretive
  conclusion, which a Data Descriptor must not carry; reviewers read the deposited files, so
  a scope note alone is not enough. The field now records counts, sources and audit status
  only. **The original free text is not lost:** it is retained byte-identical in the internal
  provenance record (`results/partD_resfeas_verdict.json`) and in the companion **code**
  archive (Zenodo `{code_doi}`). Both hashes — original and deposited — are
  recorded under `transformations` in `manifest.json`, so the change is verifiable rather
  than asserted. No other field or value in any table was altered; all four tables copy
  through unchanged.

### Fixed
- Nothing (first release).

### Source snapshot
- CELLxGENE Discover curation API v1, collection-detail endpoint, frozen
  **{SNAPSHOT_FROZEN}** (397 JSON files). Shipped as
  `{SNAPSHOT_ZIP_NAME}`; the pipeline re-runs against it offline with `--no-fetch`.
"""


def build_scope_note(code_doi: str) -> str:
    return f"""# Scope note — interpretive text and this resource

This resource is a **value-added aggregation of existing public metadata**. It describes
what the public single-cell and bulk records contain, how they were retrieved, how the
fields were standardised, and how the fields were technically validated. It states no
biological conclusion.

## The one field that needed handling

One free-text field, `verdict_text` in `partD_resfeas_verdict.json`, was written during the
audit and stated an **interpretive conclusion** rather than a description. Because reviewers
read the deposited files themselves, a scope note alone would not have been enough to keep
that sentence out of scope.

**It has therefore been replaced in this public payload with a counts-only annotation** —
counts, sources and audit status, with no interpretive content. The change is recorded under
`transformations` in `manifest.json` (original and deposited hashes, both), in `CHANGELOG.md`,
and it applies to that single field: every table copies through unchanged.

**The original text has not been discarded.** It remains byte-identical in the internal
provenance record shipped with the analysis code, which is a separate, independently citable
record (Zenodo `{code_doi}`).

## Where the interpretation lives

The argument built on top of these tables — what the field structure does and does not
permit, and the calibration and power analysis behind the gate thresholds — belongs to the
companion methods manuscript, which cites **this deposit's version DOI** so that its claims
stay pinned to the exact field set they were argued from.

## Reading the resource

Every remaining field is a count, a boolean, a fraction, or a verbatim metadata string.
Read them as descriptions of the resource — and remember the rule that must not be lost:
**not auditable is not the same as clean.** Where the control and disease arms cannot be
separated at the metadata layer, the correct record is *cannot be determined*, never
*passes*.
"""


def _strip_internal_notes(txt: str) -> str:
    """Drop ``<!-- internal note: not shipped -->`` … ``<!-- /internal note -->`` blocks.

    ``LICENSE-DATA.md`` doubles as a build template and as a repository document, so the
    instructions to us live between two markers instead of being written as prose a reader
    would have to skip past. Anything between the markers is stripped on the way into the
    deposit; the markers themselves never ship.
    """
    out, skipping = [], False
    for line in txt.split("\n"):
        if "internal note: not shipped" in line:
            skipping = True
            continue
        if "/internal note" in line:
            skipping = False
            continue
        if not skipping:
            out.append(line)
    return "\n".join(out)


def build_license_data(data_doi: str, code_doi: str) -> str:
    """The licence text as it ships in the deposit.

    Three deviations from the repository copy: both DOI tokens are substituted, and the
    internal build notes are dropped — they are instructions to us, not part of the licence
    a reader receives.

    The code DOI is substituted too, not just the data DOI: the companion code record is
    re-released alongside the data whenever a generator changes, and a licence that still
    points readers at the previous code version is exactly the kind of stale citation a
    reviewer notices.
    """
    txt = LICENSE_SRC.read_text(encoding="utf-8").replace("\r\n", "\n")
    txt = _strip_internal_notes(txt)
    txt = txt.replace(DATA_DOI_PLACEHOLDER, data_doi)
    txt = txt.replace(CODE_DOI_PLACEHOLDER, code_doi)
    # Belt and braces: the older instruction paragraph, if it is ever put back verbatim.
    kept = [p for p in txt.split("\n\n") if "publishing that record" not in p]
    txt = "\n\n".join(kept)
    return re.sub(r"\n{3,}", "\n\n", txt).rstrip("\n") + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-snapshot", action="store_true",
                    help="skip the frozen-snapshot zip (faster, non-depositable)")
    ap.add_argument("--data-doi", default=DATA_DOI_PLACEHOLDER,
                    help="the dataset DOI to print in README/manifest; use the reserved DOI "
                         "shown by scripts/zenodo_deposit_data.py before publishing the record")
    ap.add_argument("--code-doi", default=None,
                    help="the companion code version DOI to cross-link (default: "
                         "resource_paper/code_doi.txt, else the v1.0.0 code DOI)")
    args = ap.parse_args()
    data_doi = args.data_doi
    code_doi = read_code_doi(args.code_doi)

    # ---- sanity: everything we need exists -----------------------------------
    missing = [str(src) for _, src in FLAT_FILES if not src.exists()]
    if missing:
        print("MISSING source files:\n  " + "\n  ".join(missing), file=sys.stderr)
        return 2
    if not args.no_snapshot:
        if not SNAPSHOT_SRC.is_dir():
            print(f"MISSING snapshot dir: {SNAPSHOT_SRC}", file=sys.stderr)
            return 2
        n_json = len(list(SNAPSHOT_SRC.glob("*.json")))
        if n_json != 397:
            print(f"snapshot has {n_json} JSON files, expected 397", file=sys.stderr)
            return 2

    # ---- row/col counts for the README ---------------------------------------
    def dims(p: Path):
        with p.open(encoding="utf-8-sig", newline="") as fh:
            first = fh.readline().rstrip("\n").split("\t")
            return sum(1 for _ in fh), len(first)

    n_rows = {k: dims(RESULTS / f) for k, f in [
        ("libwide_datasets", "partD_resfeas_libwide_datasets.tsv"),
        ("libwide_collections", "partD_resfeas_libwide_collections.tsv"),
        ("geo_series", "partD_resfeas_geo_series.tsv"),
        ("gate_ledger", "partD_resfeas_ledger.tsv"),
    ]}
    summary = json.loads((RESULTS / "partD_resfeas_libwide_summary.json").read_text(encoding="utf-8"))
    cells_total = int(summary["cells_total"])

    # ---- rebuild the deposit directory (never touch the zip outside it) ------
    if DEPOSIT.exists():
        for child in DEPOSIT.iterdir():
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
    DEPOSIT.mkdir(parents=True, exist_ok=True)

    transformations: list[dict] = []
    for name, src in FLAT_FILES:
        if name in NEUTRALISED:
            transformations.append(neutralise_field(src, DEPOSIT / name,
                                                     NEUTRALISED[name], code_doi))
        else:
            shutil.copy2(src, DEPOSIT / name)

    write_text(DEPOSIT / "README.md", build_readme(n_rows, cells_total, data_doi, code_doi))
    write_text(DEPOSIT / "CHANGELOG.md", build_changelog(code_doi))
    write_text(DEPOSIT / "SCOPE_NOTE.md", build_scope_note(code_doi))
    write_text(DEPOSIT / "LICENSE-DATA.md", build_license_data(data_doi, code_doi))

    if not args.no_snapshot:
        # zf.write() inherits the source mtime. That is deliberate: the 397 JSON files under
        # data/raw/cxg_coll_detail/ are frozen (2026-10-09), so this zip — which *is* uploaded —
        # stays byte-identical across rebuilds. Do not switch it to a pinned timestamp now: its
        # hash is quoted in MANIFEST.sha256, manifest.json and the release docs.
        with zipfile.ZipFile(DEPOSIT / SNAPSHOT_ZIP_NAME, "w", zipfile.ZIP_DEFLATED) as zf:
            for j in sorted(SNAPSHOT_SRC.glob("*.json")):
                zf.write(j, arcname=f"cxg_coll_detail/{j.name}")

    # ---- manifest ------------------------------------------------------------
    entries = []
    for p in sorted(DEPOSIT.iterdir()):
        if p.is_file() and p.name != "MANIFEST.sha256":
            entries.append((p.name, p.stat().st_size, sha256_of(p)))

    lines = [f"{h}  {name}" for name, _size, h in entries]
    write_text(DEPOSIT / "MANIFEST.sha256",
               "# SHA-256 manifest — testability audit tables v" + VERSION + "\n"
               "# verify:  sha256sum -c MANIFEST.sha256\n"
               + "\n".join(lines) + "\n")

    manifest_json = {
        "version": VERSION,
        "released": RELEASE_DATE,
        "licence": "CC-BY-4.0",
        "code_doi": code_doi,
        "dataset_doi": None if data_doi == DATA_DOI_PLACEHOLDER else data_doi,
        "source_snapshot": {
            "repository": "CELLxGENE Discover curation API v1",
            "endpoint": "https://api.cellxgene.cziscience.com/curation/v1/collections",
            "frozen": SNAPSHOT_FROZEN,
            "n_files": 397 if not args.no_snapshot else 0,
        },
        "files": [{"name": n, "bytes": s, "sha256": h} for n, s, h in entries],
        "transformations": transformations,
    }
    write_text(DEPOSIT / "manifest.json", json.dumps(manifest_json, indent=2, ensure_ascii=False) + "\n")

    # ---- single-unit convenience zip -----------------------------------------
    # Written with a pinned mtime (see ZIP_DATE_TIME) so a rebuild is byte-for-byte a no-op.
    if SINGLE_UNIT_ZIP.exists():
        SINGLE_UNIT_ZIP.unlink()
    with zipfile.ZipFile(SINGLE_UNIT_ZIP, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(DEPOSIT.rglob("*")):
            if p.is_file():
                arc = f"testability-audit-tables-v{VERSION}/{p.relative_to(DEPOSIT).as_posix()}"
                info = zipfile.ZipInfo(arc, date_time=ZIP_DATE_TIME)
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o644 << 16
                zf.writestr(info, p.read_bytes())

    # ---- report --------------------------------------------------------------
    total = sum(s for _n, s, _h in entries)
    print(f"deposit dir : {DEPOSIT}")
    print(f"version     : {VERSION}  (released {RELEASE_DATE})")
    print(f"dataset doi : {data_doi}")
    print(f"code doi    : {code_doi}  (cross-linked in README and manifest.json)")
    print(f"files       : {len(entries)} ({total/1024:.1f} KB, snapshot excluded from this figure)")
    for n, s, _h in entries:
        print(f"  {s:>12,}  {n}")
    print(f"single unit : {SINGLE_UNIT_ZIP.name}  ({SINGLE_UNIT_ZIP.stat().st_size/1024/1024:.2f} MB)")
    print(f"manifest    : {len(entries)} entries -> MANIFEST.sha256 + manifest.json")
    for t in transformations:
        print(f"transformed : {t['file']}::{t['field']}  "
              f"(original {t['original_file_sha256'][:12]}… -> deposited "
              f"{t['deposited_file_sha256'][:12]}…)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
