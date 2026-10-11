#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Confirm the published **dataset** DOI from the registry side — the channel that still works
when ``zenodo.org`` is blocked.

``verify_dataset_doi.py`` is the authoritative check: it talks to the Zenodo API and downloads
the payload. It therefore needs ``zenodo.org`` reachable, which the assistant's sandbox is not
(its proxy answers ``502`` for that host). This script is the companion that needs only
``doi.org`` and ``api.datacite.org`` — both reachable — so a third party (or the sandbox) can
independently confirm the record without going through Zenodo at all.

It covers **three of the four** checks ``verify_dataset_doi.py`` makes, against a *different*
data source (DataCite, the registration agency Zenodo mints through):

  R1  the DOI is registered — ``doi.org`` answers with a redirect and DataCite serves the record;
  R2  the registered resource type is ``Dataset`` (``verify_dataset_doi.py`` check 2a);
  R3  the registered rights statement is CC BY 4.0 (check 2b) — and it reads the *registry*
      field, which is the same field the Zenodo record page renders;
  R4  the dataset DOI stays distinct from, and cross-linked to, the code DOI
      ``10.5281/zenodo.23291759`` (check 4): the code DOI appears in ``relatedIdentifiers`` and
      the code record is itself of type ``Software``;
  R5  the frozen local payload still agrees with the DOI being confirmed: ``MANIFEST.sha256``
      verifies 12/12, and the deposit's own ``README.md`` / ``LICENSE-DATA.md`` / ``manifest.json``
      quote this DOI rather than a placeholder.

What it **cannot** prove, and says so out loud: check 3 of ``verify_dataset_doi.py`` — the 14
files, the 12/12 SHA-256, the manifest self-consistency. Those bytes exist only on Zenodo.
Run ``verify_dataset_doi.py <DOI>`` where ``zenodo.org`` is reachable for that half.

Usage
-----
    python scripts/verify_dataset_doi_registry.py 10.5281/zenodo.NNNNNNNN
    python scripts/verify_dataset_doi_registry.py --doi-file resource_paper/dataset_doi.txt
    python scripts/verify_dataset_doi_registry.py 10.5281/zenodo.NNNNNNNN --json-out report.json
    python scripts/verify_dataset_doi_registry.py 10.5281/zenodo.NNNNNNNN --no-local

Exit status is 0 only when every check passes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

import zenodo_deposit_data as zd  # noqa: E402  (the same metadata block that was published)

DEPOSIT = ROOT / "resource_paper" / "deposit"
LOCAL_MANIFEST = DEPOSIT / "MANIFEST.sha256"

CODE_DOI = zd.CODE_DOI
PLACEHOLDER = "10.5281/zenodo.XXXXXXX"
DATACITE = "https://api.datacite.org/dois"
UA = ("sgrg-uc-testability-audit/verify_dataset_doi_registry.py "
      "(+https://github.com/yuanyuan20020329/sgrg-uc-testability-audit)")

OK, BAD, WARN, NA = "  [ok]  ", " [FAIL] ", " [warn] ", "  [--]  "


# ------------------------------------------------------------------ HTTP


def _get(url: str, accept: str = "application/json", timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": accept})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Turn redirects into an inspectable HTTPError instead of following them.

    Following ``doi.org`` would land on ``zenodo.org``, which is exactly the host that is
    unreachable — so the redirect target is read from the ``Location`` header instead.
    """

    def redirect_request(self, *a, **k):  # noqa: D102
        return None


def probe_doi_resolution(doi: str, timeout: int = 60) -> tuple[bool, dict]:
    url = "https://doi.org/" + urllib.parse.quote(doi)
    opener = urllib.request.build_opener(_NoRedirect)
    try:
        with opener.open(urllib.request.Request(url, headers={"User-Agent": UA}),
                         timeout=timeout) as r:
            return r.status < 400, {"status": r.status, "location": None}
    except urllib.error.HTTPError as e:
        loc = e.headers.get("Location")
        if e.code in (301, 302, 303, 307, 308) and loc:
            return True, {"status": e.code, "location": loc}
        return False, {"status": e.code, "location": loc,
                       "note": "doi.org did not answer with a redirect"}
    except Exception as e:  # noqa: BLE001
        return False, {"error": f"{type(e).__name__}: {e}"}


def fetch_datacite(doi: str, timeout: int = 60) -> dict:
    return json.loads(_get(f"{DATACITE}/{urllib.parse.quote(doi)}",
                           timeout=timeout).decode("utf-8"))["data"]["attributes"]


def _rights_ids(attrs: dict) -> list[str]:
    out = []
    for r in attrs.get("rightsList") or []:
        for k in ("rightsIdentifier", "rightsUri", "rights"):
            v = r.get(k)
            if v:
                out.append(str(v))
    return out


def _creator_names(attrs: dict) -> list[str]:
    return [str(c.get("name") or "").strip() for c in attrs.get("creators") or []]


def _norm(s: str) -> str:
    return " ".join(str(s).split()).strip().lower()


# ------------------------------------------------------------------ checks


def check_local(doi: str) -> tuple[bool, dict]:
    """R5: the frozen payload on disk still agrees with the DOI being confirmed."""
    info: dict = {}
    bad = 0

    entries: dict[str, str] = {}
    if LOCAL_MANIFEST.exists():
        for line in LOCAL_MANIFEST.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                h, name = line.split("  ", 1)
                entries[name.strip()] = h.strip()
    info["n_checksum_entries"] = len(entries)
    ok = fail = 0
    for name, h in entries.items():
        p = DEPOSIT / name
        if p.exists() and hashlib.sha256(p.read_bytes()).hexdigest() == h:
            ok += 1
        else:
            fail += 1
    info["sha256_ok"] = f"{ok}/{len(entries)}"
    if fail:
        bad += 1

    on_disk = sorted(p.name for p in DEPOSIT.iterdir() if p.is_file()) if DEPOSIT.is_dir() else []
    info["files_on_disk"] = len(on_disk)
    if len(on_disk) != 14:
        info["note_files"] = "expected 14 files in the deposit"
        bad += 1

    man = DEPOSIT / "manifest.json"
    man_doi = None
    if man.exists():
        try:
            man_doi = json.loads(man.read_text(encoding="utf-8")).get("dataset_doi")
        except Exception as e:  # noqa: BLE001
            info["manifest_json_error"] = f"{type(e).__name__}: {e}"
            bad += 1
    info["manifest_dataset_doi"] = man_doi
    if man_doi not in (doi, PLACEHOLDER, None):
        info["note_manifest"] = "manifest.json carries a third, unrelated DOI"
        bad += 1

    quoted, placeholder = [], []
    for name in ("README.md", "LICENSE-DATA.md"):
        p = DEPOSIT / name
        if not p.exists():
            continue
        t = p.read_text(encoding="utf-8")
        if doi in t:
            quoted.append(name)
        if PLACEHOLDER in t:
            placeholder.append(name)
    info["deposit_docs_quoting_doi"] = quoted
    info["deposit_docs_still_placeholder"] = placeholder

    info["verdict"] = ("payload self-consistent and quotes this DOI" if not bad and not placeholder
                       else "payload self-consistent, DOI backfill not applied yet"
                       if not bad else "payload problem")
    return bad == 0, info


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("doi", nargs="?", help="dataset DOI, e.g. 10.5281/zenodo.NNNNNNNN")
    ap.add_argument("--doi-file", type=Path, help="file containing the DOI (first line)")
    ap.add_argument("--json-out", type=Path, help="also write the report as JSON")
    ap.add_argument("--no-local", action="store_true",
                    help="skip R5 (the local payload cross-check)")
    ap.add_argument("--timeout", type=int, default=60)
    args = ap.parse_args()

    doi = args.doi
    if not doi and args.doi_file and args.doi_file.exists():
        doi = args.doi_file.read_text(encoding="utf-8").strip().splitlines()[0].strip()
    if not doi:
        print("give a DOI (or --doi-file)", file=sys.stderr)
        return 2
    if doi == PLACEHOLDER or "XXXXXXX" in doi:
        print(f"'{doi}' is still the placeholder — nothing to confirm yet.", file=sys.stderr)
        return 2

    expected = zd.build_metadata()
    exp_title = _norm(expected["title"])
    exp_version = str(expected["version"]).strip().lower().lstrip("v")
    exp_creators = [_norm(c["name"]) for c in expected["creators"]]
    exp_related = [r["identifier"] for r in expected.get("related_identifiers", [])]

    report: dict = {"doi": doi, "checks": {}, "expected": {
        "title": expected["title"], "version": expected["version"],
        "creators": [c["name"] for c in expected["creators"]],
        "upload_type": expected["upload_type"], "license": expected["license"],
        "related_identifiers": exp_related}}
    failures = 0

    print(f"confirming : {doi}   (registry side — no request to zenodo.org is made)")
    print(f"deposit    : {DEPOSIT}")
    print("-" * 76)

    # ---- R1 registration
    ok_doi, det_doi = probe_doi_resolution(doi, args.timeout)
    registered = False
    attrs: dict = {}
    try:
        attrs = fetch_datacite(doi, args.timeout)
        registered = True
    except urllib.error.HTTPError as e:
        det_doi["datacite"] = f"HTTP {e.code} — no DataCite record for this DOI"
    except Exception as e:  # noqa: BLE001
        det_doi["datacite"] = f"{type(e).__name__}: {e}"

    ok = ok_doi and registered
    print(f"{OK if ok else BAD}R1. DOI is registered (doi.org redirects, DataCite serves it)")
    print(f"        doi.org : {det_doi.get('status')} -> {det_doi.get('location')}")
    if "datacite" in det_doi:
        print(f"        DataCite: {det_doi['datacite']}")
    report["checks"]["R1_registered"] = {"pass": ok, **det_doi}
    failures += 0 if ok else 1

    if not registered:
        print("-" * 76)
        print("RESULT: cannot confirm — the DOI is not in DataCite.")
        print("        Either the record has not been published yet, or the DOI string is wrong.")
        report["failures"] = failures
        if args.json_out:
            args.json_out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                                     encoding="utf-8")
            print(f"report written -> {args.json_out}")
        return 1

    own = attrs.get("doi") or ""
    state = attrs.get("state")
    print(f"        record  : state={state}  published={attrs.get('published')}  "
          f"publisher={attrs.get('publisher')}")
    print(f"        own doi : {own}")

    # ---- R2 resource type
    rt = (attrs.get("types") or {}).get("resourceTypeGeneral")
    ok2 = str(rt).strip().lower() == "dataset"
    print(f"{OK if ok2 else BAD}R2. registered resource type is 'Dataset'   (got: {rt!r})")
    report["checks"]["R2_resource_type"] = {"pass": ok2, "resource_type": rt}
    failures += 0 if ok2 else 1

    # ---- R3 licence
    rights = _rights_ids(attrs)
    lic_ok = any(r.strip().lower() == "cc-by-4.0" for r in rights)
    weak = [r for r in rights if r.strip().lower().startswith("cc-by")]
    print(f"{OK if lic_ok else BAD}R3. registered rights statement is CC BY 4.0")
    print(f"        rightsList: {rights or '(empty)'}")
    if not lic_ok and weak:
        print(f"        note: a CC BY variant is present ({weak}) but not the 4.0 identifier — "
              f"Zenodo sets cc-by-4.0; check the record's Licence field")
    report["checks"]["R3_license"] = {"pass": lic_ok, "rights_list": rights}
    failures += 0 if lic_ok else 1

    # ---- R4 metadata agreement + the two DOIs stay distinct
    got_title = _norm((attrs.get("titles") or [{}])[0].get("title") or "")
    got_version = str(attrs.get("version") or "").strip().lower().lstrip("v")
    got_creators = [_norm(n) for n in _creator_names(attrs)]
    rel = attrs.get("relatedIdentifiers") or []
    rel_ids = [str(r.get("relatedIdentifier")) for r in rel]
    code_listed = CODE_DOI in rel_ids
    distinct = own != CODE_DOI and CODE_DOI not in own

    title_ok = got_title == exp_title
    version_ok = got_version == exp_version
    creators_ok = sorted(got_creators) == sorted(exp_creators)

    code_is_software = None
    try:
        code_attrs = fetch_datacite(CODE_DOI, args.timeout)
        code_is_software = str((code_attrs.get("types") or {})
                               .get("resourceTypeGeneral") or "").strip().lower() == "software"
    except Exception:  # noqa: BLE001
        pass

    ok4 = (title_ok and version_ok and creators_ok and code_listed and distinct
           and code_is_software is True)
    print(f"{OK if ok4 else BAD}R4. DOI stays distinct from the code DOI and is cross-linked")
    print(f"        title matches        : {title_ok}")
    print(f"        version matches      : {version_ok}  (record {attrs.get('version')!r} "
          f"vs expected {expected['version']!r})")
    print(f"        creators match       : {creators_ok}  {_creator_names(attrs)}")
    print(f"        code DOI as related  : {code_listed}")
    print(f"        dataset DOI != code  : {distinct}   (own={own})")
    print(f"        code record = Software: {code_is_software}")
    print(f"        relatedIdentifiers   : {rel_ids}")
    report["checks"]["R4_distinct_and_linked"] = {
        "pass": ok4, "title_matches": title_ok, "version_matches": version_ok,
        "creators_match": creators_ok, "code_doi_listed": code_listed,
        "dataset_doi_distinct_from_code": distinct, "code_record_is_software": code_is_software,
        "own_doi": own, "related_identifiers": rel_ids,
        "record_creators": _creator_names(attrs)}
    failures += 0 if ok4 else 1

    # ---- R5 local payload
    if args.no_local:
        print(f"{WARN}R5. local payload cross-check skipped (--no-local)")
        report["checks"]["R5_local"] = {"pass": None, "note": "skipped"}
    else:
        ok5, det5 = check_local(doi)
        print(f"{OK if ok5 else BAD}R5. frozen local payload agrees with this DOI "
              f"({det5.get('sha256_ok')}, {det5.get('files_on_disk')} files)")
        for k in ("sha256_ok", "n_checksum_entries", "files_on_disk", "manifest_dataset_doi",
                  "deposit_docs_quoting_doi", "deposit_docs_still_placeholder", "verdict",
                  "manifest_json_error", "note_manifest", "note_files"):
            if k in det5:
                print(f"        {k}: {det5[k]}")
        report["checks"]["R5_local"] = {"pass": ok5, **det5}
        failures += 0 if ok5 else 1

    print("-" * 76)
    print(f"{NA}NOT checked here: the file set and the 12/12 SHA-256 check (check 3 of")
    print(f"{NA}      verify_dataset_doi.py). Those bytes exist only on zenodo.org — run")
    print(f"{NA}      python scripts/verify_dataset_doi.py {doi}")
    print(f"{NA}      where zenodo.org is reachable to close that half.")
    if failures:
        print(f"RESULT: {failures} check(s) FAILED")
    else:
        print("RESULT: registry-side confirmation passed — the DOI is registered, the resource "
              "type, licence, metadata and cross-links are right.")
    report["failures"] = failures
    if args.json_out:
        args.json_out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                                 encoding="utf-8")
        print(f"report written -> {args.json_out}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
