# -*- coding: utf-8 -*-
"""Stamp the two real v1.0.2 version DOIs into the Chinese manuscript.

The manuscript may not go to the journal with an unresolvable identifier in it, and it may
not carry a DOI that does not exist yet either. Until the Zenodo records are published the
manuscript therefore cites the **concept** DOIs (which always resolve) plus one interim
sentence — ``1.0.2 版已完成修订，待正式归档后以版本 DOI 引用`` — at five places. Once
``publish_dataset.py`` has minted the dataset version DOI and the companion code archive has
minted its own, this script swaps those five interim sites plus the two reference entries
over to the version DOIs, and nothing else.

It refuses to run while:

* a DOI placeholder (``zenodo.XXXXXXX`` / ``zenodo.YYYYYYY``) is still present;
* a supplied DOI is the *previous* version DOI — a re-release must mint a new one;
* an anchor is missing or matches more than once, so a silent partial edit is impossible.

    python scripts/fill_v1p0p2_manuscript_dois.py --check \
        --data-doi 10.5281/zenodo.NNNNNNNN --code-doi 10.5281/zenodo.MMMMMMMM
    python scripts/fill_v1p0p2_manuscript_dois.py \
        --data-doi 10.5281/zenodo.NNNNNNNN --code-doi 10.5281/zenodo.MMMMMMMM

The English manuscript is intentionally **not** touched: it has not yet been brought up to
date with the fourth- and fifth-round revisions, so stamping it with a v1.0.2 DOI would
misdescribe it. Extend ANCHORS when it is synchronised.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

CN = ROOT / "resource_paper" / "Data_Descriptor_v1_CN.md"

VERSION = "1.0.2"
PREVIOUS_DATA_DOI = "10.5281/zenodo.23299341"
PREVIOUS_CODE_DOI = "10.5281/zenodo.23298914"
DATA_CONCEPT = "10.5281/zenodo.23297702"
CODE_CONCEPT = "10.5281/zenodo.23291758"

DOI_RE = re.compile(r"^10\.5281/zenodo\.\d{6,}$")
PLACEHOLDER_RE = re.compile(r"zenodo\.(?:XXXXXXX|YYYYYYY)")

INTERIM = "1.0.2 版已完成修订，待正式归档后以版本 DOI 引用"

# (site id, old, new).  {DATA} / {CODE} are filled from the command line.
# Every `old` is the *whole* sentence fragment that disappears; every `new` keeps the
# concept DOI where it was, so only the version identifier is added.
ANCHORS: list[tuple[str, str, str]] = [
    (
        "abstract",
        f"（概念 DOI {DATA_CONCEPT}，指向全部版本）；**{INTERIM}。**",
        f"（{VERSION} 版版本 DOI {{DATA}}；概念 DOI {DATA_CONCEPT} 指向全部版本）。",
    ),
    (
        "software",
        f"（概念 DOI {CODE_CONCEPT}；**{INTERIM}**）",
        f"（{VERSION} 版代码版本 DOI {{CODE}}；概念 DOI {CODE_CONCEPT}）",
    ),
    (
        "data-records",
        f"其概念 DOI {DATA_CONCEPT} 指向该资源的全部版本。**{INTERIM}。**",
        f"其 {VERSION} 版版本 DOI 为 {{DATA}}；概念 DOI {DATA_CONCEPT} 指向该资源的全部版本。",
    ),
    (
        "data-availability",
        f"概念 DOI {DATA_CONCEPT} 指向该资源的全部版本。**{INTERIM}。**该记录包含四个表格",
        f"概念 DOI {DATA_CONCEPT} 指向该资源的全部版本；本文所述 {VERSION} 版的版本 DOI 为 "
        f"{{DATA}}。该记录包含四个表格",
    ),
    (
        "code-availability",
        f"概念 DOI {CODE_CONCEPT}。**{INTERIM}。**代码记录与数据集记录独立版本化并相互链接；",
        f"概念 DOI {CODE_CONCEPT}；{VERSION} 版代码版本 DOI 为 {{CODE}}。"
        f"代码记录与数据集记录独立版本化并相互链接；",
    ),
    (
        "reference-14",
        f"*Zenodo* https://doi.org/{CODE_CONCEPT} (2026).",
        f"*Zenodo* https://doi.org/{{CODE}} (2026).",
    ),
    (
        "reference-17",
        f"*Zenodo* https://doi.org/{DATA_CONCEPT} (2026).",
        f"*Zenodo* https://doi.org/{{DATA}} (2026).",
    ),
]

EXPECTED_INTERIM_SITES = 5
EXPECTED_CONCEPT_DATA_SITES = 4   # abstract, data-records, data-availability, ref 17
EXPECTED_CONCEPT_CODE_SITES = 3   # software, code-availability, ref 14


class Report:
    def __init__(self) -> None:
        self.rows: list[tuple[bool, str, str]] = []

    def check(self, ok: bool, label: str, detail: str = "") -> bool:
        self.rows.append((bool(ok), label, detail))
        print(f"  [{'ok' if ok else 'FAIL'}] {label}"
              + (f"  -- {detail}" if detail and not ok else ""))
        return bool(ok)

    def note(self, label: str, detail: str = "") -> None:
        print(f"  [--] {label}" + (f"  -- {detail}" if detail else ""))

    @property
    def failed(self) -> int:
        return sum(1 for ok, _, _ in self.rows if not ok)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-doi", required=True, help="v1.0.2 dataset version DOI")
    ap.add_argument("--code-doi", required=True, help="v1.0.2 companion code version DOI")
    ap.add_argument("--md", type=Path, default=CN)
    ap.add_argument("--check", action="store_true", help="report only, write nothing")
    args = ap.parse_args()

    data, code = args.data_doi.strip(), args.code_doi.strip()
    rep = Report()

    print("=" * 74)
    print(f"stamp v{VERSION} version DOIs into {args.md.name}"
          + ("   [CHECK ONLY]" if args.check else ""))
    print("=" * 74)

    rep.check(bool(DOI_RE.match(data)), f"--data-doi is a Zenodo DOI ({data})")
    rep.check(bool(DOI_RE.match(code)), f"--code-doi is a Zenodo DOI ({code})")
    rep.check(data != code, "the two DOIs differ from each other")
    rep.check(data != PREVIOUS_DATA_DOI,
              "--data-doi is not the previous dataset version DOI", f"{data}")
    rep.check(code != PREVIOUS_CODE_DOI,
              "--code-doi is not the previous code version DOI", f"{code}")
    if rep.failed:
        return 1

    if not args.md.exists():
        rep.check(False, f"{args.md} exists")
        return 1

    text = args.md.read_text(encoding="utf-8")
    original = text

    rep.check(not PLACEHOLDER_RE.search(text),
              "no DOI placeholder token in the manuscript",
              f"{PLACEHOLDER_RE.findall(text)[:4]}")

    applied = [site for site, _, new in ANCHORS
               if text.count(new.format(DATA=data, CODE=code)) == 1]
    if len(applied) == len(ANCHORS) and text.count(INTERIM) == 0:
        rep.note("the manuscript already carries both version DOIs — nothing to do")
        print("\nRESULT: already up to date, nothing written")
        return 0

    rep.check(text.count(INTERIM) == EXPECTED_INTERIM_SITES,
              f"interim sentence present {EXPECTED_INTERIM_SITES}x before the swap",
              f"found {text.count(INTERIM)}")
    rep.check(text.count(DATA_CONCEPT) == EXPECTED_CONCEPT_DATA_SITES,
              f"data concept DOI present {EXPECTED_CONCEPT_DATA_SITES}x",
              f"found {text.count(DATA_CONCEPT)}")
    rep.check(text.count(CODE_CONCEPT) == EXPECTED_CONCEPT_CODE_SITES,
              f"code concept DOI present {EXPECTED_CONCEPT_CODE_SITES}x",
              f"found {text.count(CODE_CONCEPT)}")
    if rep.failed:
        print("\nnothing written: the manuscript is not in the expected pre-publication "
              "state.")
        return 1

    pending = 0
    for site, old, new in ANCHORS:
        new_filled = new.format(DATA=data, CODE=code)
        n_old, n_new = text.count(old), text.count(new_filled)
        if n_new == 1 and n_old == 0:
            rep.note(f"{site}: already applied, left alone")
            continue
        if n_old != 1:
            rep.check(False, f"{site}: anchor matches exactly once", f"matched {n_old}x")
            continue
        pending += 1
        if args.check:
            rep.note(f"{site}: would be updated")
        else:
            text = text.replace(old, new_filled, 1)
            rep.note(f"{site}: updated")

    if rep.failed:
        return 1

    if args.check:
        print(f"\n{pending} of {len(ANCHORS)} anchor(s) still to apply.")
        print("RESULT: dry run only, no file written")
        return 0

    if text == original:
        print("\nRESULT: already up to date, nothing written")
        return 0

    args.md.write_text(text, encoding="utf-8", newline="\n")

    # ---- read the file back and re-assert the post-conditions
    back = args.md.read_text(encoding="utf-8")
    print("\npost-conditions")
    rep.check(back.count(INTERIM) == 0,
              "no interim sentence left in the manuscript",
              f"{back.count(INTERIM)} left")
    rep.check(back.count(data) == 4,
              "the dataset version DOI is cited 4x "
              "(abstract / data records / data availability / reference 17)",
              f"found {back.count(data)}")
    rep.check(back.count(code) == 3,
              "the code version DOI is cited 3x "
              "(software / code availability / reference 14)",
              f"found {back.count(code)}")
    for concept, n in ((DATA_CONCEPT, 3), (CODE_CONCEPT, 2)):
        rep.check(back.count(concept) == n,
                  f"concept DOI {concept} still cited {n}x", f"found {back.count(concept)}")
    rep.check(not PLACEHOLDER_RE.search(back), "no placeholder token introduced")

    print(f"\nwrote {args.md}")
    print("RESULT: " + ("version DOIs stamped" if not rep.failed else "verify the checks"))
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main())
