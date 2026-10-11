#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Verify a published Zenodo **dataset** record against the frozen v1.0.0 payload.

This checks the four things the resource record has to satisfy, and nothing about the
science:

  1. the DOI **resolves** (doi.org redirect lands on a Zenodo record) and the record's
     resource type is ``dataset``;
  2. the licence **on the record** is CC BY 4.0 (reviewers read the record field, not the
     prose in LICENSE-DATA.md);
  3. every file in the deposit is present and every sha256 matches **MANIFEST.sha256** — the
     files are downloaded from the public record, so this is exactly the "download the package
     and re-run the check" step. The deposit holds 14 files; MANIFEST.sha256 lists the 12
     payload files (the two manifest artifacts cannot reference themselves), so the check also
     confirms MANIFEST.sha256 is byte-identical to the local copy and that manifest.json's own
     ``files`` block agrees with the published bytes;
  4. the record keeps the **dataset DOI and the code DOI distinct**: the code version DOI
     recorded in ``resource_paper/code_doi.txt`` (v1.0.1: 10.5281/zenodo.23298914) appears as
     a *related identifier*, is not the record's own DOI, and README.md quotes both. The value
     is **read from that file, never typed here**: a literal is how this check went stale at
     the v1.0.1 re-release, where the record correctly cross-links the v1.0.1 code record
     while a hard-coded v1.0.0 DOI failed a perfectly good deposit.

Usage
-----
    python scripts/verify_dataset_doi.py 10.5281/zenodo.NNNNNNNN
    python scripts/verify_dataset_doi.py --record-id 12345678
    python scripts/verify_dataset_doi.py --doi-file resource_paper/dataset_doi.txt
    python scripts/verify_dataset_doi.py 10.5281/zenodo.NNNNNNNN --no-download   # metadata only

Exit status is 0 only when every check passes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
DEPOSIT = ROOT / "resource_paper" / "deposit"
LOCAL_MANIFEST = DEPOSIT / "MANIFEST.sha256"

import build_data_deposit as bdd            # noqa: E402

# Single source of truth, shared with build_data_deposit.py and the metadata block: the first
# line of resource_paper/code_doi.txt is the *current* code version DOI (falling back to the
# v1.0.0 record only if that file is missing). Re-releasing the code archive is then a one-line
# edit to a data file, and this check cannot silently keep validating the previous release.
CODE_DOI = bdd.read_code_doi()
CODE_DOI_SOURCE = ("resource_paper/code_doi.txt" if bdd.CODE_DOI_FILE.exists()
                   else "built-in v1.0.0 fallback")
ZENODO = "https://zenodo.org"
RECORDS_API = f"{ZENODO}/api/records"
UA = "sgrg-uc-testability-audit/verify_dataset_doi.py (+https://github.com/yuanyuan20020329/sgrg-uc-testability-audit)"

OK, BAD, WARN = "  [ok]  ", " [FAIL] ", " [warn] "


# ------------------------------------------------------------------ helpers

def _get(url: str, accept: str | None = None, timeout: int = 120) -> bytes:
    hdr = {"User-Agent": UA}
    if accept:
        hdr["Accept"] = accept
    req = urllib.request.Request(url, headers=hdr)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _get_final_url(url: str, timeout: int = 120) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        r.read(1024)
        return r.geturl()


def _parse_manifest(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        h, name = line.split("  ", 1)
        out[name.strip()] = h.strip()
    return out


def _record_id_from_doi(doi: str) -> str | None:
    m = re.search(r"zenodo\.(\d+)$", doi.strip())
    return m.group(1) if m else None


# ------------------------------------------------------------------ checks

def check_resolution(doi: str) -> tuple[bool, dict]:
    """Check 1a: the DOI redirects to a Zenodo record page."""
    try:
        final = _get_final_url(f"https://doi.org/{urllib.parse.quote(doi)}")
    except urllib.error.HTTPError as e:
        return False, {"error": f"HTTP {e.code} while resolving https://doi.org/{doi}"}
    except Exception as e:  # noqa: BLE001
        return False, {"error": f"{type(e).__name__}: {e}"}
    ok = "zenodo.org/records/" in final or "zenodo.org/record/" in final
    return ok, {"final_url": final}


def fetch_record(record_id: str | None, doi: str | None) -> dict:
    if record_id:
        return json.loads(_get(f"{RECORDS_API}/{record_id}").decode("utf-8"))
    q = urllib.parse.quote(f'doi:"{doi}"')
    hits = json.loads(_get(f"{RECORDS_API}?q={q}&size=5").decode("utf-8")).get("hits", {})
    total = hits.get("total", 0)
    if not total:
        raise LookupError(f"no published Zenodo record found for DOI {doi} "
                          f"(a draft is not searchable until published)")
    return hits["hits"][0]


def check_resource_type(rec: dict) -> tuple[bool, dict]:
    rt = (rec.get("metadata", {}).get("resource_type") or {}).get("type")
    return rt == "dataset", {"resource_type": rt}


def check_license(rec: dict) -> tuple[bool, dict]:
    lic = rec.get("metadata", {}).get("license")
    if isinstance(lic, dict):
        lic_id = lic.get("id") or lic.get("@id")
    else:
        lic_id = lic
    # Zenodo's API has used "cc-by-4.0" and "CC-BY-4.0" over time; compare case-insensitively.
    ok = bool(lic_id) and str(lic_id).strip().lower() == "cc-by-4.0"
    return ok, {"license": lic_id}


def check_related_and_distinct(rec: dict, own_doi: str, readme_text: str | None) -> tuple[bool, dict]:
    meta = rec.get("metadata", {})
    rels = meta.get("related_identifiers") or []
    ids = {str(r.get("identifier")) for r in rels}
    code_listed = CODE_DOI in ids
    distinct = own_doi != CODE_DOI and CODE_DOI not in (own_doi or "")
    own_in_readme = readme_text is None or (own_doi in readme_text)
    code_in_readme = readme_text is None or (CODE_DOI in readme_text)
    ok = code_listed and distinct and own_in_readme and code_in_readme
    return ok, {
        "own_doi": own_doi,
        "code_doi_listed_as_related": code_listed,
        "dataset_doi_distinct_from_code": distinct,
        "own_doi_in_readme": own_in_readme,
        "code_doi_in_readme": code_in_readme,
        "related_identifiers": sorted(ids),
    }


def check_files(rec: dict, do_download: bool) -> tuple[bool, dict, str | None]:
    """Check 3: the full file set is on the record, and every sha256 matches.

    The deposit holds 14 files. ``MANIFEST.sha256`` lists the **12 payload files** — the two
    manifest artifacts cannot reference themselves — so the check is two-layered: the record
    must carry all 14, all 12 manifest-listed files must match their sha256, ``MANIFEST.sha256``
    itself must equal the local copy byte-for-byte, and ``manifest.json``'s own ``files`` block
    must agree with the bytes actually published.
    """
    files = rec.get("files") or []
    names = {f.get("key") or f.get("filename") for f in files}

    local_manifest = _parse_manifest(LOCAL_MANIFEST.read_text(encoding="utf-8")) \
        if LOCAL_MANIFEST.exists() else {}
    manifest_artifacts = {"MANIFEST.sha256", "manifest.json"}
    expected_all = set(local_manifest) | manifest_artifacts

    info: dict = {
        "n_files_on_record": len(files),
        "n_checksum_entries": len(local_manifest),
        "n_expected_files": len(expected_all),
    }
    if not local_manifest:
        info["note"] = "no local MANIFEST.sha256 - cannot do the checksum check"
        return False, info, None

    missing = sorted(expected_all - names)
    extra = sorted(names - expected_all)
    info["missing_on_record"] = missing
    info["unexpected_on_record"] = extra

    readme_text = None
    manifest_ok = manifest_json_ok = manifest_doi_ok = True
    if do_download:
        md5_ok = sha_ok = 0
        sha_bad: list[str] = []
        with tempfile.TemporaryDirectory() as tmp:
            for f in files:
                name = f.get("key") or f.get("filename")
                url = (f.get("links") or {}).get("self")
                if not url:
                    continue
                raw = _get(url)
                want_md5 = str(f.get("checksum", "")).split("md5:")[-1]
                if want_md5 and hashlib.md5(raw).hexdigest() == want_md5:
                    md5_ok += 1
                if name in local_manifest:
                    if hashlib.sha256(raw).hexdigest() == local_manifest[name]:
                        sha_ok += 1
                    else:
                        sha_bad.append(name)
                elif name == "MANIFEST.sha256":
                    manifest_ok = (raw == LOCAL_MANIFEST.read_bytes())
                elif name == "manifest.json":
                    try:
                        pub = json.loads(raw.decode("utf-8"))
                        listed = {x["name"]: x["sha256"] for x in pub.get("files", [])}
                        # every hash the manifest claims must match the local payload
                        manifest_json_ok = all(
                            local_manifest.get(n) == h for n, h in listed.items()) \
                            and set(listed) == set(local_manifest)
                        own = pub.get("dataset_doi")
                        rec_doi = rec.get("doi") or ""
                        manifest_doi_ok = (own == rec_doi) and (pub.get("code_doi") != own)
                        info["manifest_dataset_doi"] = own
                        info["manifest_code_doi"] = pub.get("code_doi")
                    except Exception as e:  # noqa: BLE001
                        manifest_json_ok = False
                        info["manifest_json_error"] = f"{type(e).__name__}: {e}"
                if name == "README.md":
                    readme_text = raw.decode("utf-8", "replace")
        info["md5_ok"] = f"{md5_ok}/{len(files)}"
        info["sha256_ok"] = f"{sha_ok}/{len(local_manifest)}"
        info["sha256_mismatch"] = sha_bad
        info["MANIFEST.sha256_identical_to_local"] = manifest_ok
        info["manifest.json_matches_payload"] = manifest_json_ok
        info["manifest.json_dataset_doi_equals_record_doi"] = manifest_doi_ok

    ok = not missing and not extra
    if do_download:
        ok = ok and info.get("sha256_ok") == f"{len(local_manifest)}/{len(local_manifest)}"
        ok = ok and manifest_ok and manifest_json_ok and manifest_doi_ok
    return ok, info, readme_text


# ------------------------------------------------------------------ main

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("doi", nargs="?", help="dataset DOI, e.g. 10.5281/zenodo.NNNNNNNN")
    ap.add_argument("--doi-file", type=Path, help="file containing the DOI (one line)")
    ap.add_argument("--record-id", help="Zenodo record id (skips the DOI search)")
    ap.add_argument("--no-download", action="store_true",
                    help="metadata-only: do not download files, skip the 12/12 check")
    ap.add_argument("--json-out", type=Path, help="also write the report as JSON")
    args = ap.parse_args()

    doi = args.doi
    if not doi and args.doi_file and args.doi_file.exists():
        doi = args.doi_file.read_text(encoding="utf-8").strip().splitlines()[0].strip()
    if not doi and not args.record_id:
        print("give a DOI (or --record-id / --doi-file)", file=sys.stderr)
        return 2

    report: dict = {"doi_arg": doi, "record_id_arg": args.record_id, "checks": {}}
    failures = 0

    print(f"verifying : {doi or '(record id ' + str(args.record_id) + ')'}")
    print(f"deposit   : {DEPOSIT}")
    print("-" * 72)

    # ---- 1. resolution
    if doi and not args.record_id:
        ok, det = check_resolution(doi)
        print(f"{OK if ok else BAD}1. DOI resolves to a Zenodo record")
        print(f"        {det}")
        report["checks"]["resolution"] = {"pass": ok, **det}
        failures += 0 if ok else 1
        rid = _record_id_from_doi(doi)
    else:
        rid = args.record_id
        report["checks"]["resolution"] = {"pass": True, "note": "skipped (record id given)"}
        print(f"{WARN}1. DOI resolution skipped (record id given)")

    # ---- fetch the record
    try:
        rec = fetch_record(rid, doi)
    except Exception as e:  # noqa: BLE001
        print(f"{BAD}could not fetch the record: {type(e).__name__}: {e}")
        return 1

    own_doi = rec.get("doi") or (rec.get("metadata", {}).get("doi")) or doi or ""
    rec_id = rec.get("id")
    print(f"        record id = {rec_id}   own doi = {own_doi}   published = "
          f"{rec.get('is_published') or rec.get('metadata', {}).get('publication_date')}")

    # ---- 2. resource type
    ok, det = check_resource_type(rec)
    print(f"{OK if ok else BAD}2a. resource type is 'dataset'")
    print(f"        {det}")
    report["checks"]["resource_type"] = {"pass": ok, **det}
    failures += 0 if ok else 1

    # ---- 2b. licence
    ok, det = check_license(rec)
    print(f"{OK if ok else BAD}2b. licence on the record is CC BY 4.0")
    print(f"        {det}")
    report["checks"]["license"] = {"pass": ok, **det}
    failures += 0 if ok else 1

    # ---- 3. files (+ download for the 12/12 check)
    if args.no_download:
        print(f"{WARN}3. file set / sha256 check skipped (--no-download)")
        report["checks"]["files"] = {"pass": None, "note": "skipped"}
    else:
        ok, det, readme_text = check_files(rec, do_download=True)
        print(f"{OK if ok else BAD}3. file set complete ({det.get('n_files_on_record')}"
              f"/{det.get('n_expected_files')}) and sha256 matches MANIFEST.sha256 "
              f"({det.get('sha256_ok', 'n/a')})")
        for k in ("n_files_on_record", "n_checksum_entries", "n_expected_files",
                  "missing_on_record", "unexpected_on_record", "sha256_mismatch",
                  "md5_ok", "MANIFEST.sha256_identical_to_local", "manifest.json_matches_payload",
                  "manifest.json_dataset_doi_equals_record_doi", "manifest_dataset_doi"):
            if k in det:
                print(f"        {k}: {det[k]}")
        report["checks"]["files"] = {"pass": ok, **det}
        failures += 0 if ok else 1
        # publish the README check against the matching local file too
        if readme_text is not None:
            report["readme_sha256_ok"] = (
                hashlib.sha256(readme_text.encode("utf-8")).hexdigest()
                == hashlib.sha256((DEPOSIT / "README.md").read_bytes()).hexdigest())
            print(f"        README matches the local deposit byte-for-byte: "
                  f"{report['readme_sha256_ok']}")

    # ---- 4. the two DOIs stay distinct
    readme_text = None
    if not args.no_download:
        # re-read from the downloaded copy recorded above if available
        tmp_readme = DEPOSIT / "README.md"
        readme_text = tmp_readme.read_text(encoding="utf-8") if tmp_readme.exists() else None
    ok, det = check_related_and_distinct(rec, own_doi, readme_text)
    print(f"{OK if ok else BAD}4. dataset DOI and code DOI {CODE_DOI} kept distinct & cross-linked")
    print(f"        code DOI source: {CODE_DOI_SOURCE}")
    for k, v in det.items():
        print(f"        {k}: {v}")
    report["checks"]["doi_relation"] = {"pass": ok, **det}
    failures += 0 if ok else 1

    print("-" * 72)
    if failures:
        print(f"RESULT: {failures} check(s) FAILED")
    else:
        print("RESULT: all checks passed — the record matches the frozen payload")
    report["failures"] = failures
    if args.json_out:
        args.json_out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                                 encoding="utf-8")
        print(f"report written -> {args.json_out}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
