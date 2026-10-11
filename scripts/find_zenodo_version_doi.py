#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Find the newest Zenodo version DOI in a record family, via DataCite.

Why this exists
---------------
A new version's DOI is minted *by Zenodo* at publication and is therefore not known before
it happens — but the payload has to cite it, the manuscript has to cite it, and the code
record cross-link has to point at the right one. The remaining DOI is discovered here
rather than copied off a web page, so a transcription slip cannot enter the release.

The lookup deliberately goes through **DataCite** (`api.datacite.org`), not zenodo.org:
Zenodo mirrors every registered DOI and its `IsVersionOf` relation to DataCite, so the
family can be enumerated even when zenodo.org itself is unreachable from this machine.

Usage
-----
    # every version DOI in the code record family (concept DOI 10.5281/zenodo.23291758)
    python scripts/find_zenodo_version_doi.py 10.5281/zenodo.23291758

    # just the newest one, and only if it is newer than a known version
    python scripts/find_zenodo_version_doi.py 10.5281/zenodo.23291758 \\
            --newer-than 10.5281/zenodo.23291759

    # machine-readable, for the next step in a runbook
    python scripts/find_zenodo_version_doi.py 10.5281/zenodo.23291758 --json

Exit codes: 0 found, 1 nothing found (or nothing newer than --newer-than), 2 network error.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

DATACITE = "https://api.datacite.org/dois"
UA = "sgrg-uc-testability-audit/find_zenodo_version_doi.py"


def _get(url: str, attempts: int = 4) -> dict:
    """GET with retries, because api.datacite.org occasionally drops a TLS handshake."""
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if 400 <= e.code < 500:
                raise SystemExit(f"[FAIL] DataCite answered HTTP {e.code} {e.reason} for {url}")
            last = e
        except Exception as e:  # noqa: BLE001
            last = e
        if i < attempts - 1:
            print(f"  [retry {i + 1}/{attempts - 1}] {type(last).__name__}: {str(last)[:70]}")
            time.sleep(3)
    raise SystemExit(f"[FAIL] could not reach DataCite after {attempts} attempts: "
                     f"{type(last).__name__}: {last}")


def family(anchor: str, page_size: int = 100) -> list[dict]:
    """Every DOI whose related identifiers point at ``anchor``, plus ``anchor`` itself."""
    out: dict[str, dict] = {}

    def add(item: dict) -> None:
        a = item.get("attributes") or {}
        doi = (a.get("doi") or "").lower()
        if not doi:
            return
        rels = {r.get("relatedIdentifier", "").lower()
                for r in (a.get("relatedIdentifiers") or [])}
        out[doi] = {
            "doi": a.get("doi"),
            "version": a.get("version"),
            "published": (a.get("published") or a.get("publicationYear") or ""),
            "titles": [t.get("title") for t in (a.get("titles") or [])],
            "related": sorted(rels),
            "url": a.get("url"),
        }

    # 1. the anchor record itself
    try:
        add(_get(f"{DATACITE}/{urllib.parse.quote(anchor, safe='')}")["data"])
    except Exception:  # noqa: BLE001
        print(f"  [warn] {anchor} itself did not resolve on DataCite; continuing")

    # 2. everyone who claims a relation to it (IsVersionOf / IsNewVersionOf / IsPreviousVersionOf)
    q = urllib.parse.quote(f"relatedIdentifiers.relatedIdentifier:{anchor}")
    url = f"{DATACITE}?query={q}&page%5Bsize%5D={page_size}"
    doc = _get(url)
    for item in doc.get("data") or []:
        add(item)

    def rank(d: dict) -> tuple:
        """Order by Zenodo record id, not by date.

        DataCite reports `published` as a bare year for Zenodo records, so a date sort is a
        coin flip between two versions minted the same day. Zenodo's record ids increase
        monotonically, so the numeric suffix is the reliable ordering — with the date as a
        tie-breaker.
        """
        tail = str(d.get("doi") or "").rsplit(".", 1)[-1]
        return (int(tail) if tail.isdigit() else 0, str(d.get("published") or ""))

    return sorted(out.values(), key=rank)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("anchor", help="a concept DOI or any version DOI of the family")
    ap.add_argument("--newer-than", default=None,
                    help="a DOI already known to you; succeed only if something newer exists")
    ap.add_argument("--json", action="store_true", help="print JSON for the next tool")
    args = ap.parse_args()

    print(f"family anchor : {args.anchor}")
    members = family(args.anchor)
    print(f"found         : {len(members)} DOI(s)")

    if not members:
        print("[FAIL] nothing found. If a release was *just* created, Zenodo may not have "
              "archived it yet - a few minutes is normal; try again.")
        return 1

    for m in members:
        print(f"  {m['doi']:<28} version={str(m['version'] or '?'):<8} "
              f"published={str(m['published'])[:10]}  {(m['titles'] or [''])[0][:60]}")

    known = sum(1 for m in members if args.newer_than
                and m["doi"].lower() == args.newer_than.lower())
    newest_by_date = members[-1]
    if args.newer_than:
        idx = next((i for i, m in enumerate(members)
                    if m["doi"].lower() == args.newer_than.lower()), None)
        if idx is None:
            print(f"[warn] {args.newer_than} is not in the DataCite family listing "
                  f"({len(members)} members) - falling back to the newest by date")
        elif idx == len(members) - 1:
            print(f"[FAIL] nothing newer than {args.newer_than} yet ({known} member(s) known). "
                  f"Check that the GitHub release was actually published, and give Zenodo a "
                  f"few minutes.")
            return 1
        else:
            newest_by_date = members[-1]

    print()
    print(f"NEWEST        : {newest_by_date['doi']}  (version "
          f"{newest_by_date['version'] or '?'}, published "
          f"{str(newest_by_date['published'])[:10]})")
    if args.json:
        print(json.dumps(newest_by_date, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
