#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Propagate the minted **v1.0.2** DOIs into the documents that cite them.

Companion to :mod:`backfill_release_dois`, which stays byte-for-byte the script that
served the v1.0.0 -> v1.0.1 transition. This one covers v1.0.1 -> v1.0.2 and shares the
same :class:`Report` machinery, so an edit is accepted in either state (pending, or
already applied) and anything else is a hard stop rather than a silent no-op.

What changed in v1.0.2: the **dataset payload's documentation** was synchronised with the
Data Descriptor that cites it (``G3_powered_control`` is a preset design-eligibility gate,
not a "power gate"; the candidate-pairing and control-vocabulary warnings; the three
per-table UC conventions). No table byte, field name, gate value or number moved. The code
archive is therefore re-released as v1.0.2 too, because the payload cross-links it.

**The two manuscripts are deliberately not touched.** The English Data Descriptor and the
methods manuscript still carry v1.0.1-era content, and stamping a v1.0.2 version string on
them without the v1.0.2 content would be worse than leaving them: they are brought up to
date, and their DOIs bumped, in a single later pass.

Usage
-----
    # rehearse: report what is still stale, write nothing
    python scripts/backfill_v1p0p2_dois.py --check

    # after the GitHub release for tag v1.0.2 has been archived by Zenodo
    python scripts/backfill_v1p0p2_dois.py --code-doi 10.5281/zenodo.<id>

    # after publish_dataset.py has written resource_paper/dataset_doi.txt
    python scripts/backfill_v1p0p2_dois.py
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from backfill_release_dois import (  # noqa: E402
    DEPOSIT, Report, _fmt, _nl_of, payload_totals, read_first_doi,
)

VERSION = "1.0.2"
RELEASE_DATE = "2026-10-16"

PREV_VERSION = "1.0.1"
PREV_DATE = "2026-10-11"

CODE_CONCEPT = "10.5281/zenodo.23291758"
DATA_CONCEPT = "10.5281/zenodo.23297702"

PREV_CODE = "10.5281/zenodo.23298914"   # the v1.0.1 code record (now the previous one)
PREV_DATA = "10.5281/zenodo.23299341"   # the v1.0.1 dataset record (now the previous one)
V100_CODE = "10.5281/zenodo.23291759"
V100_DATA = "10.5281/zenodo.23297703"

CODE_DOI_FILE = ROOT / "resource_paper" / "code_doi.txt"
DATA_DOI_FILE = ROOT / "resource_paper" / "dataset_doi.txt"


# ------------------------------------------------------------------ archive side (code)

def write_code_doi_record(new_code: str) -> None:
    text = (
        f"{new_code}\n"
        f"# concept doi: {CODE_CONCEPT}\n"
        f"# version: v{VERSION}\n"
        "#\n"
        "# Single source of truth for the COMPANION CODE record that the dataset payload, its\n"
        "# manifest.json, the deposit README and the manuscript all cross-link to.\n"
        "#\n"
        "#   line 1            the code VERSION DOI (what the data record cites)\n"
        "#   \"# concept doi:\"  the code CONCEPT DOI (always resolves to the newest code version)\n"
        "#   \"# version:\"      the code version string, e.g. v1.0.2\n"
        "#\n"
        "# This file is the ONLY place the code DOI is recorded. The payload build reads it, so a\n"
        "# re-release becomes a one-line change here rather than a hunt through the generators.\n"
        "# scripts/audit_data_descriptor.py fails when \"# version:\" disagrees with the payload\n"
        "# version, and that check is wired in as the pre-publish gate of\n"
        "# scripts/publish_dataset.py.\n"
    )
    CODE_DOI_FILE.write_text(text, encoding="utf-8", newline="\n")
    print(f"  [ok] {CODE_DOI_FILE.relative_to(ROOT)} -> {new_code} (v{VERSION})")


def edit_archive_side(rep: Report, new_code: str) -> None:
    print("\n[1/4] code archive documents")

    rep.edit("archive_src/CITATION.cff",
             "# v1.0.1. The version DOI below is the record Zenodo minted for this release; the concept\n"
             f"# DOI {CODE_CONCEPT} always resolves to the newest version of this code archive.\n",
             f"# v{VERSION}. The version DOI below is the record Zenodo minted for this release; the concept\n"
             f"# DOI {CODE_CONCEPT} always resolves to the newest version of this code archive.\n")
    rep.edit("archive_src/CITATION.cff",
             'version: "v1.0.1"', f'version: "v{VERSION}"')
    rep.edit("archive_src/CITATION.cff",
             'date-released: "2026-10-11"', f'date-released: "{RELEASE_DATE}"')
    rep.edit("archive_src/CITATION.cff",
             f"license: MIT\ndoi: {PREV_CODE}\n", f"license: MIT\ndoi: {new_code}\n")
    rep.edit(
        "archive_src/CITATION.cff",
        "identifiers:\n"
        f"  - type: doi\n    value: {PREV_CODE}\n"
        f'    description: "Zenodo archival DOI (version v{PREV_VERSION})"\n'
        f"  - type: doi\n    value: {CODE_CONCEPT}\n"
        '    description: "Zenodo concept DOI of this code archive (resolves to the newest version)"\n'
        f"  - type: doi\n    value: {V100_CODE}\n"
        '    description: "Zenodo archival DOI (version v1.0.0, published and immutable)"\n',
        "identifiers:\n"
        f"  - type: doi\n    value: {new_code}\n"
        f'    description: "Zenodo archival DOI (version v{VERSION})"\n'
        f"  - type: doi\n    value: {PREV_CODE}\n"
        f'    description: "Zenodo archival DOI (version v{PREV_VERSION}, published and immutable)"\n'
        f"  - type: doi\n    value: {CODE_CONCEPT}\n"
        '    description: "Zenodo concept DOI of this code archive (resolves to the newest version)"\n'
        f"  - type: doi\n    value: {V100_CODE}\n"
        '    description: "Zenodo archival DOI (version v1.0.0, published and immutable)"\n',
    )

    rep.edit("archive_src/zenodo/deposit_metadata.md",
             f"## Version\n\nv{PREV_VERSION}\n", f"## Version\n\nv{VERSION}\n")
    rep.edit("archive_src/zenodo/deposit_metadata.md",
             f"## Publication date\n\n{PREV_DATE}\n", f"## Publication date\n\n{RELEASE_DATE}\n")
    rep.edit(
        "archive_src/zenodo/deposit_metadata.md",
        f"* This archive (v{PREV_VERSION}): {PREV_CODE}\n"
        f"* Concept DOI (resolves to the newest version): {CODE_CONCEPT}\n"
        f"* Previous version (v1.0.0): {V100_CODE} (published and immutable)\n",
        f"* This archive (v{VERSION}): {new_code}\n"
        f"* Concept DOI (resolves to the newest version): {CODE_CONCEPT}\n"
        f"* Previous version (v{PREV_VERSION}): {PREV_CODE} (published and immutable)\n"
        f"* Previous version (v1.0.0): {V100_CODE} (published and immutable)\n",
    )
    rep.edit(
        "archive_src/zenodo/deposit_metadata.md",
        f"* \u2705 Done (v{PREV_VERSION}, {PREV_DATE}): version DOI `{PREV_CODE}` backfilled into\n"
        "  `resource_paper/code_doi.txt`, `CITATION.cff`, `README.md`, this file and the manuscript's\n"
        f"  *Data and code availability* statement. The data payload cites the **v{PREV_VERSION}** code DOI,\n"
        "  not the v1.0.0 one \u2014 `scripts/audit_data_descriptor.py --pre-publish` fails otherwise.\n",
        f"* \u2705 Done (v{PREV_VERSION}, {PREV_DATE}): version DOI `{PREV_CODE}` backfilled into\n"
        "  `resource_paper/code_doi.txt`, `CITATION.cff`, `README.md`, this file and the manuscript's\n"
        f"  *Data and code availability* statement. The data payload cites the **v{PREV_VERSION}** code DOI,\n"
        "  not the v1.0.0 one \u2014 `scripts/audit_data_descriptor.py --pre-publish` fails otherwise.\n"
        f"* \u2705 Done (v{VERSION}, {RELEASE_DATE}): version DOI `{new_code}` backfilled into\n"
        "  `resource_paper/code_doi.txt`, `CITATION.cff`, `README.md` and this file. The v1.0.2 data\n"
        f"  payload cites the **v{VERSION}** code DOI. The manuscripts are updated in a later pass \u2014\n"
        "  this release is data-dictionary synchronisation only.\n",
    )

    rep.edit(
        "archive_src/README.md",
        "## Release notes\n\n**v1.0.1** \u2014 no change to any analysis",
        f"""## Release notes

**v{VERSION}** \u2014 no change to any analysis, any number, or any conclusion. The deposited
dataset payload's **documentation** was synchronised with the Data Descriptor that cites it:

1. **`G3_powered_control` is described as a preset design-eligibility gate, not a "power
   gate".** The v1.0.1 `datapackage.json` used the word *power*, which contradicted the
   semantics the manuscript explains. The field itself is unchanged; only its description is.
2. **The warnings now ship with the data.** `n_paired_donors` / `n_paired_same_site` are
   marked *candidate* counts that need manual verification, and the known control-vocabulary
   misclassification (`normal pressure hydrocephalus`) plus the three per-table UC-selection
   conventions are written into `README.md` and the affected field descriptions. Field names,
   table structure, the four gates, `USABLE`, every numeric value and the row order are
   unchanged. `scripts/test_documentation_sync.py` asserts exactly that by diffing the deposit
   against the frozen v1.0.1 payload file by file, and
   `scripts/audit_data_descriptor.py --pre-publish` runs it as a hard pre-publication gate.
3. **This code archive was re-released** so that it is the archive that built the v1.0.2
   payload. New in v1.0.2: `scripts/test_documentation_sync.py` and the updated
   `scripts/build_data_deposit.py` templates.

**v1.0.1** \u2014 no change to any analysis""",
    )


# ----------------------------------------------------------- working docs (code + data)

def edit_working_docs(rep: Report, new_code: str, new_data: str | None) -> None:
    print("\n[2/4] resource_paper working docs")

    # ---- depends only on the code DOI ---------------------------------------
    rep.edit("resource_paper/README.md",
             f"**The frozen, upload-ready v{PREV_VERSION} payload for the dataset record**",
             f"**The frozen, upload-ready v{VERSION} payload for the dataset record**")
    rep.edit(
        "resource_paper/README.md",
        "* **The code archive** (analysis pipeline + key result tables, MIT) \u2014\n"
        f"  **{PREV_CODE}** \u2014 https://doi.org/{PREV_CODE}\n"
        f"  (v{PREV_VERSION}; the concept DOI {CODE_CONCEPT} points at all versions.)\n",
        "* **The code archive** (analysis pipeline + key result tables, MIT) \u2014\n"
        f"  **{new_code}** \u2014 https://doi.org/{new_code}\n"
        f"  (v{VERSION}; the concept DOI {CODE_CONCEPT} points at all versions.)\n",
    )
    rep.edit(
        "resource_paper/README.md",
        f"_First release 2026-10-11 (v1.0.0). **Current version v{PREV_VERSION}** \u2014 English-language "
        f"normalisation of the free-text fields only; no gate outcome, numeric value or row changed. "
        f"See `deposit/CHANGELOG.md`. v1.0.0 stays published and immutable._",
        f"_First release 2026-10-11 (v1.0.0). **Current version v{VERSION}** \u2014 the payload's "
        f"documentation is synchronised with the Data Descriptor (`G3_powered_control` is a design-"
        f"eligibility gate, not a power gate; the candidate-pairing and control-vocabulary warnings "
        f"now ship with the data). No field name, gate outcome, numeric value or row changed. "
        f"See `deposit/CHANGELOG.md`. Earlier versions stay published and immutable._",
    )
    rep.edit("resource_paper/SUBMISSION_READINESS.md",
             "A single-unit convenience zip is built at "
             "`resource_paper/testability-audit-tables-v1.0.1.zip`.",
             "A single-unit convenience zip is built at "
             f"`resource_paper/testability-audit-tables-v{VERSION}.zip`.")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             "| **Version** | 1.0.1 |", f"| **Version** | {VERSION} |")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             "| **Publication date** | 2026-10-11 |", f"| **Publication date** | {RELEASE_DATE} |")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"| `{PREV_CODE}` | is supplemented by | Software |",
             f"| `{new_code}` | is supplemented by | Software |")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"the **code** archive (`{PREV_CODE}`), and `manifest.json` records **both** hashes:",
             f"the **code** archive (`{new_code}`), and `manifest.json` records **both** hashes:")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"(4) the record's own DOI, the dataset DOI in `README.md`, and the code DOI\n"
             f"`{PREV_CODE}` are all present and **distinct**.",
             f"(4) the record's own DOI, the dataset DOI in `README.md`, and the code DOI\n"
             f"`{new_code}` are all present and **distinct**.")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             "## 9. As published \u2014 v1.0.0, 2026-10-11 (superseded by v1.0.1; kept as history)\n",
             f"## 9. As published \u2014 v1.0.0, 2026-10-11 (superseded by v{VERSION}; kept as history)\n")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             "## 10. As published \u2014 v1.0.1\n",
             f"## 10. As published \u2014 v{PREV_VERSION}, {PREV_DATE} "
             f"(superseded by v{VERSION}; kept as history)\n")

    if not new_data:
        print("  [..] the dataset-DOI edits are deferred until dataset_doi.txt is written")
        return

    # ---- depends on the dataset DOI ----------------------------------------
    rep.edit(
        "resource_paper/README.md",
        "* **The tables (this resource, dataset record)** \u2014 "
        f"**{PREV_DATA}**\n"
        f"  \u2014 https://doi.org/{PREV_DATA}\n"
        f"  (version of record, v{PREV_VERSION}; concept DOI {DATA_CONCEPT} points at all versions.)\n",
        "* **The tables (this resource, dataset record)** \u2014 "
        f"**{new_data}**\n"
        f"  \u2014 https://doi.org/{new_data}\n"
        f"  (version of record, v{VERSION}; concept DOI {DATA_CONCEPT} points at all versions.)\n",
    )
    rep.edit("resource_paper/SUBMISSION_READINESS.md",
             f"\u2705 **v{PREV_VERSION} published {PREV_DATE}** | The payload is live as Zenodo record "
             f"{PREV_DATA.rsplit('.', 1)[-1]} (v{PREV_VERSION}; see gap 3).",
             f"\u2705 **v{VERSION} published {RELEASE_DATE}** | The payload is live as Zenodo record "
             f"{new_data.rsplit('.', 1)[-1]} (v{VERSION}; see gap 3).")
    rep.edit(
        "resource_paper/SUBMISSION_READINESS.md",
        f"| 3 | Data must be **live** and directly linked to the paper | \u2705 **v{PREV_VERSION} minted "
        f"{PREV_DATE}** | **Dataset DOI `{PREV_DATA}`** (version of record) \u2014 record live at "
        f"https://zenodo.org/records/{PREV_DATA.rsplit('.', 1)[-1]} (status `published`); "
        f"v1.0.0 `{V100_DATA}` remains published and immutable. Concept DOI `{DATA_CONCEPT}`.",
        f"| 3 | Data must be **live** and directly linked to the paper | \u2705 **v{VERSION} minted "
        f"{RELEASE_DATE}** | **Dataset DOI `{new_data}`** (version of record) \u2014 record live at "
        f"https://zenodo.org/records/{new_data.rsplit('.', 1)[-1]} (status `published`); "
        f"v{PREV_VERSION} `{PREV_DATA}` and v1.0.0 `{V100_DATA}` remain published and immutable. "
        f"Concept DOI `{DATA_CONCEPT}`.",
    )
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"version DOI: both quote `{PREV_DATA}`, so the superseded `{V100_DATA}` must not be\n"
             f"   left behind).",
             f"version DOI: both quote `{new_data}`, so the superseded `{PREV_DATA}` must not be\n"
             f"   left behind).")


# ------------------------------------------------------------------- byte-dependent

def refresh_byte_counts(rep: Report) -> None:
    print("\n[3/4] byte counts (recomputed from the rebuilt payload)")
    n_files, total, mb = payload_totals()
    print(f"  payload: {n_files} files, {total:,} B ({mb:.2f} MB)")

    rel = "resource_paper/zenodo_deposit_datasets.md"
    path = ROOT / rel
    raw = path.read_text(encoding="utf-8")
    new = re.sub(r"\(14 files, [\d,]+ B = [\d.]+ MB total, of which",
                 f"({n_files} files, {_fmt(total)} B = {mb:.2f} MB total, of which", raw)
    if new != raw:
        if rep.check:
            rep.pending += 1
            print(f"  [pending] {rel}  ::  refresh the payload byte total")
        else:
            path.write_text(new, encoding="utf-8", newline=_nl_of(raw))
            rep.applied += 1
            if rel not in rep.touched:
                rep.touched.append(rel)
            print(f"  [ok] payload total refreshed -> {_fmt(total)} B")
    else:
        print("  [ok] payload byte total already current")

    row = re.compile(r"^(\| `([A-Za-z0-9_.\-]+)` \| )([\d,]+)( \|.*)$", re.M)

    def sub(m: re.Match) -> str:
        f = DEPOSIT / m.group(2)
        return f"{m.group(1)}{_fmt(f.stat().st_size)}{m.group(4)}" if f.is_file() else m.group(0)

    raw = path.read_text(encoding="utf-8")
    new, n = row.subn(sub, raw)
    if n and new != raw:
        if rep.check:
            rep.pending += 1
            print(f"  [pending] {rel}  ::  {n} byte rows to refresh")
        else:
            path.write_text(new, encoding="utf-8", newline=_nl_of(raw))
            rep.applied += 1
            if rel not in rep.touched:
                rep.touched.append(rel)
        print(f"  [ok] {n} rows in the payload table carry current sizes")
    else:
        print("  [ok] payload table already current")


def append_workbook_section(rep: Report, new_code: str, new_data: str, total: int) -> None:
    """Append the v1.0.2 'as published' record to the workbook. Idempotent."""
    rel = "resource_paper/zenodo_deposit_datasets.md"
    path = ROOT / rel
    raw = path.read_text(encoding="utf-8")
    if f"## 11. As published \u2014 v{VERSION}" in raw:
        print(f"  [ok] {rel} already carries the v{VERSION} record")
        return
    if rep.check:
        rep.pending += 1
        print(f"  [pending] {rel}  ::  append the v{VERSION} \"as published\" record")
        return
    nl = _nl_of(raw)
    section = f"""
---

## 11. As published \u2014 v{VERSION}

| | |
|---|---|
| Dataset DOI (version of record) | **{new_data}** |
| Concept DOI (all versions) | {DATA_CONCEPT} |
| Code DOI (version of record) | **{new_code}** |
| Code concept DOI (all versions) | {CODE_CONCEPT} |
| Previous dataset version | {PREV_DATA} (v{PREV_VERSION}, published and immutable) |
| Previous code version | {PREV_CODE} (v{PREV_VERSION}, published and immutable) |
| Files | 14 ({_fmt(total)} B), resource type **Dataset**, licence **CC BY 4.0**, version {VERSION} |

**What changed from v1.0.1: the payload's documentation.** The dictionary described
`G3_powered_control` as a *"power gate"*, while the Data Descriptor that cites the deposit
explains it as a **preset design-eligibility rule, not a power calculation** \u2014 a reader who
downloads the package must not be told a different semantics than the paper. The same release
marks `n_paired_donors` / `n_paired_same_site` as **candidate counts** needing manual
verification, records the one known control-vocabulary misclassification
(`normal pressure hydrocephalus`), and writes the three per-table UC-selection conventions
into `README.md` and the affected field descriptions.

**No table byte, field name, column order, gate value, boolean or number moved**, and that is
asserted rather than promised: `scripts/test_documentation_sync.py` diffs this payload against
the frozen v1.0.1 copy in `resource_paper/frozen_deposit_v1.0.1/` and only the allow-listed
descriptions and prose may differ. The test runs as a hard pre-publication gate
(`audit_data_descriptor.py --pre-publish`, wired into `publish_dataset.py`), and the payload's
`CHANGELOG.md` carries the same statement.

**Why the code archive was re-released too.** The deposit builder's templates and the release
gate changed, so the scripts in v1.0.1 no longer reproduce these bytes. Because the payload
cross-links the code DOI, a documentation-only data release with a stale code link would be a
citation error; `audit_data_descriptor.py` fails the build when `code_doi.txt`'s version
disagrees with the payload version, so it cannot ship silently.
"""
    path.write_text(raw.rstrip("\r\n") + nl + section, encoding="utf-8", newline=nl)
    rep.applied += 1
    if rel not in rep.touched:
        rep.touched.append(rel)
    print(f"  [ok] appended the v{VERSION} record to {rel}")


# ------------------------------------------------------------------------- pipeline

def run(cmd: list[str], label: str) -> int:
    print(f"\n$ {' '.join(cmd)}")
    r = subprocess.run(cmd, cwd=str(ROOT))
    print(f"  exit {r.returncode}  ({label})")
    return r.returncode


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--code-doi", help="the v1.0.2 CODE version DOI just minted by Zenodo")
    ap.add_argument("--data-doi", help="override the dataset version DOI (else dataset_doi.txt)")
    ap.add_argument("--check", action="store_true", help="report what is still stale; touch nothing")
    ap.add_argument("--no-build", action="store_true",
                    help="skip the payload rebuild and the gate (documents only)")
    args = ap.parse_args()

    recorded_code = read_first_doi(CODE_DOI_FILE)
    recorded_data = read_first_doi(DATA_DOI_FILE)
    new_code = args.code_doi or (recorded_code if recorded_code not in (None, PREV_CODE, V100_CODE)
                                 else None)
    new_data = args.data_doi or (recorded_data if recorded_data not in (None, PREV_DATA, V100_DATA)
                                 else None)

    print("=" * 74)
    print(f"backfill the v{VERSION} DOIs")
    print(f"  code  : {new_code or '(not minted yet - pass --code-doi)'}")
    print(f"  data  : {new_data or '(not minted yet - publish_dataset.py writes dataset_doi.txt)'}")
    print(f"  mode  : {'check only' if args.check else 'apply'}")
    print("=" * 74)

    if not new_code:
        print(f"\nno v{VERSION} code DOI yet.\n"
              "  The code archive is released first, because the data payload cites it:\n"
              f"    1. publish the GitHub release for tag v{VERSION}\n"
              f"    2. python scripts/backfill_v1p0p2_dois.py --code-doi 10.5281/zenodo.<id>\n"
              f"  Pass --code-doi now or let the minted DOI land in "
              f"{CODE_DOI_FILE.relative_to(ROOT)}.")
        return 2
    if not new_data:
        print(f"\n  [note] the dataset v{VERSION} DOI is not recorded yet, so only the code-side\n"
              "         documents and the payload's code cross-link are updated. Re-run this\n"
              "         after publish_dataset.py writes dataset_doi.txt.")
    print("\n  [note] the two manuscripts are NOT touched by this pass: they still carry\n"
          "         v1.0.1-era content, and are brought up to date (version, DOIs, content)\n"
          "         in a single later pass.")

    rep = Report(args.check)

    if not args.check:
        print("\n[0/4] record the code DOI")
        write_code_doi_record(new_code)

    edit_archive_side(rep, new_code)
    edit_working_docs(rep, new_code, new_data)

    if args.check:
        if new_data:
            refresh_byte_counts(rep)
            append_workbook_section(rep, new_code, new_data, payload_totals()[1])
    else:
        if not args.no_build and new_data:
            print("\n[3a/4] rebuild the payload with the new DOIs")
            rc = run([sys.executable, str(ROOT / "scripts" / "build_data_deposit.py"),
                      "--data-doi", new_data, "--code-doi", new_code], "build payload")
            if rc:
                return rc
        if new_data:
            refresh_byte_counts(rep)
            append_workbook_section(rep, new_code, new_data, payload_totals()[1])
        else:
            print("\n[3/4] byte counts: skipped - the payload has not been rebuilt with the "
                  f"v{VERSION} dataset DOI yet")

    print("\n" + "-" * 74)
    print(f"edits applied : {rep.applied}")
    print(f"already done  : {rep.uptodate}")
    if rep.check:
        print(f"still pending : {rep.pending}")
    for f in rep.touched:
        print(f"  changed: {f}")
    if rep.problems:
        print("\nPROBLEMS (nothing was written for these):")
        for p in rep.problems:
            print(f"  [FAIL] {p}")

    if args.check:
        return 1 if (rep.pending or rep.problems) else 0
    if rep.problems:
        return 1

    if not args.no_build and new_data:
        print("\n[4/4] gate + consistency audits")
        rc = run([sys.executable, str(ROOT / "scripts" / "audit_data_descriptor.py"),
                  "--pre-publish"], "pre-publish gate")
        if rc:
            print("\nthe pre-publish gate is red - fix before publishing.")
            return rc
        for script in ("audit_citations.py", "reconcile_numbers.py"):
            run([sys.executable, str(ROOT / "scripts" / script)], script)

    print("\nnext: commit `release/` (the code archive) and re-run the payload verifiers.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
