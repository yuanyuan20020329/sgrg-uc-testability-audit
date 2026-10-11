#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Mint the **dataset** DOI for the audit tables (the resource record).

The code archive already has its own Zenodo record (10.5281/zenodo.23291759, created
through Zenodo's GitHub integration on the repo's v1.0.0 release). Scientific Data
requires the *dataset* to be deposited in a public repository with its **own** persistent
identifier, so this script creates a **separate** Zenodo record by direct upload. The
two records are cross-linked but never merged.

Typical use
-----------
    # 1. see exactly what would be uploaded (no network, no token needed)
    python scripts/zenodo_deposit_data.py --dry-run

    # 2. write the paste-ready metadata block on its own
    python scripts/zenodo_deposit_data.py --emit-metadata resource_paper/zenodo_deposit_datasets.json

    # 3. rehearse on the sandbox first (sandbox token, separate account)
    set ZENODO_TOKEN=<sandbox-token>
    python scripts/zenodo_deposit_data.py --sandbox

    # 4. the real thing; leaves a DRAFT you can eyeball, then --publish
    set ZENODO_TOKEN=<real-token>
    python scripts/zenodo_deposit_data.py
    python scripts/zenodo_deposit_data.py --publish --deposition-id <id>

Zenodo does **not** assign the DOI at publication only — a draft can *pre-reserve* one
(``prereserve_doi``), it just does not resolve until the record is published. That is what
makes it possible to backfill the DOI into the archived README before publishing, so the
deposited bytes cite their own DOI. Two caveats from Zenodo's own docs:

* the pre-reserved DOI is **internal until publication and not guaranteed** — never delete
  the draft, or it is lost (Zenodo cannot reassign it);
* the *legacy* deposit API is the one that preserves the reservation across the metadata
  ``PUT``; the v1.4 draft API has a known DOI deadlock.

Token: a Zenodo personal access token with scopes ``deposit:write`` (+ ``deposit:actions``
to publish). Pass it via the ``ZENODO_TOKEN`` environment variable or ``--token-file``.
Never commit it.
"""
from __future__ import annotations

import argparse
import http.client
import json
import os
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_data_deposit as bdd          # noqa: E402  (single source for version + code DOI)

DEPOSIT = ROOT / "resource_paper" / "deposit"
META_OUT_DEFAULT = ROOT / "resource_paper" / "zenodo_deposit_datasets.json"

VERSION = bdd.VERSION
RELEASE_DATE = bdd.RELEASE_DATE
# The companion code record is read from resource_paper/code_doi.txt, so re-releasing the
# code never needs an edit here — only the data file does.
CODE_DOI = bdd.read_code_doi()
GITHUB = "https://github.com/yuanyuan20020329/sgrg-uc-testability-audit"

PROD = "https://zenodo.org/api"
SANDBOX = "https://sandbox.zenodo.org/api"

# Files inside the deposit are uploaded in a deliberate order (readme first); the
# ordering lives in upload_list() below.


def upload_list() -> list[Path]:
    """Every file in the deposit directory, README first, then tables, then snapshot."""
    if not DEPOSIT.is_dir():
        return []
    files = [p for p in DEPOSIT.iterdir() if p.is_file()]
    priority = {"README.md": 0, "SCOPE_NOTE.md": 1, "CHANGELOG.md": 2, "LICENSE-DATA.md": 3,
                "datapackage.json": 4, "MANIFEST.sha256": 5, "manifest.json": 6}
    return sorted(files, key=lambda p: (priority.get(p.name, 10), p.name))


def build_metadata() -> dict:
    """The Zenodo ``metadata`` block. HTML is allowed in ``description``."""
    description = f"""<p><strong>Derived, citable annotation of public single-cell and bulk
resources for ulcerative colitis.</strong> The tables describe, for every CELLxGENE Discover
collection and dataset, where the cells came from, whether the control and disease arms can
be separated in the metadata layer, a sequencing-depth proxy, and the outcome of a four-gate
usability ledger. They are a value-added aggregation of existing public metadata.</p>

<p><strong>Contents.</strong></p>
<ul>
<li><code>partD_resfeas_libwide_datasets.tsv</code> &mdash; 2,238 rows, one per CELLxGENE dataset
across all 397 collections (27 fields).</li>
<li><code>partD_resfeas_libwide_collections.tsv</code> &mdash; 397 rows, one per collection, with
reuse fractions, arm typology and separability, reuse asymmetry and the depth ratio (32 fields).</li>
<li><code>partD_resfeas_geo_series.tsv</code> &mdash; 423 rows, one per screened GEO series, with the
molecular gate derived from <code>!Series_type</code>, an in-vitro flag and same-donor pairing
counts (21 fields).</li>
<li><code>partD_resfeas_ledger.tsv</code> &mdash; 29 rows: the four-gate usability ledger with a
per-row human-readable <code>reason</code> string (16 fields).</li>
<li><code>datapackage.json</code> &mdash; a Frictionless Data Package describing field names,
types, enum constraints and semantics for all four tables.</li>
<li><code>cxg_coll_detail_snapshot_v1.0.0.zip</code> &mdash; the frozen CELLxGENE curation-API
snapshot (397 JSON files, 2026-10-09) against which the pipeline re-runs offline, so the release
reproduces without depending on the live API.</li>
<li><code>MANIFEST.sha256</code> &mdash; SHA-256 of every redistributed file.</li>
</ul>

<p><strong>No primary omics data are redistributed.</strong> Every source collection and dataset
remains under its own depositor's licence; accessions and retrieval scripts are in the companion
code archive. Coverage counts in the tables are properties of the resource, not findings.</p>

<p><strong>One documented transformation.</strong> A free-text field that stated an interpretive
conclusion (<code>verdict_text</code> in <code>partD_resfeas_verdict.json</code>) has been replaced
in this payload by a counts-only annotation, so that the deposited files describe the resource
rather than argue about it. The original text is retained in the internal provenance record
shipped with the analysis code; original and deposited hashes are both recorded under
<code>transformations</code> in <code>manifest.json</code>, and explained in
<code>CHANGELOG.md</code> and <code>SCOPE_NOTE.md</code>. No other value was altered.</p>

<p><strong>Usage rule.</strong> Not auditable is not the same as clean: where the control and
disease arms cannot be separated at the metadata layer, the correct record is
<em>cannot be determined</em>, never <em>passes</em>.</p>

<p><strong>This version (v{VERSION}).</strong> English-language normalisation only. The free-text
fields written in Chinese in v1.0.0 &mdash; the per-row <code>reason</code> string, the field
glossary and <code>caveat</code> note of <code>partD_resfeas_libwide_summary.json</code>, one
example string and two schema descriptions &mdash; are now English. No gate outcome, numeric
value, field semantics or row membership changed: field names, table structure, column order, the
four gates, <code>USABLE</code> and every number are identical to v1.0.0, which stays available
and unaltered at its own version DOI.</p>"""

    return {
        "title": ("Testability audit of public single-cell and bulk resources for ulcerative "
                  "colitis: collection- and dataset-level provenance, arm separability and "
                  "gate ledger"),
        "upload_type": "dataset",
        "description": description,
        "creators": [
            {
                "name": "Yuan, Yike",
                "affiliation": "Hubei University of Chinese Medicine, Wuhan, Hubei, China",
                "orcid": "0009-0009-8391-4134",
            },
            {
                "name": "Chen, Hongci",
                "affiliation": ("Hubei Provincial Hospital of Traditional Chinese Medicine, "
                                "Wuhan, Hubei, China"),
                "orcid": "0009-0000-3335-8998",
            },
        ],
        "access_right": "open",
        "license": "cc-by-4.0",
        "version": VERSION,
        "publication_date": RELEASE_DATE,
        "language": "eng",
        "keywords": [
            "data provenance", "data re-use", "metadata annotation", "arm separability",
            "single-cell genomics", "ulcerative colitis", "CELLxGENE Discover", "GEO",
            "FAIR data", "reproducibility",
        ],
        "related_identifiers": [
            {"identifier": CODE_DOI, "relation": "isSupplementedBy", "resource_type": "software",
             "scheme": "doi"},
            {"identifier": GITHUB, "relation": "isDerivedFrom", "resource_type": "software",
             "scheme": "url"},
        ],
        "notes": ("Companion code archive (separate, independently citable record): "
                  f"https://doi.org/{CODE_DOI}. A Data Descriptor describing this resource "
                  "and a methods manuscript arguing its interpretation cite this version DOI. "
                  "v1.0.1 is an English-language normalisation of the free-text fields; "
                  "v1.0.0 remains available and unaltered."),
    }


# --------------------------------------------------------------------------- HTTP

# Failures that mean "the link hiccuped", not "the request was wrong". Zenodo drops a
# connection every so often on a multi-file upload - the first real run of this deposit
# died on the 14th of 14 files with RemoteDisconnected after 13 clean uploads. Retrying is
# the correct response to these; an HTTP 4xx is not, and is raised on the spot.
TRANSIENT = (
    http.client.RemoteDisconnected,
    http.client.IncompleteRead,
    http.client.BadStatusLine,
    ConnectionResetError,
    ConnectionAbortedError,
    BrokenPipeError,
    socket.timeout,
    TimeoutError,
    urllib.error.URLError,
)


def _retry(call, what: str, attempts: int = 4, delay: float = 4.0):
    """Run ``call`` until it succeeds, retrying only *transient* failures.

    Nothing here is allowed to change the request: reuse it only for calls that can be
    repeated without side effects (a file PUT, a GET). Creating or publishing a deposition
    is never retried, because a lost response would not tell us whether it happened.
    """
    for attempt in range(1, attempts + 1):
        try:
            return call()
        except urllib.error.HTTPError as e:
            if not (500 <= e.code < 600) or attempt == attempts:
                raise
            print(f"  [retry {attempt}/{attempts - 1}] {what}: HTTP {e.code} {e.reason} "
                  f"- waiting {delay:.0f}s")
        except TRANSIENT as e:  # noqa: PERF203
            if attempt == attempts:
                raise
            print(f"  [retry {attempt}/{attempts - 1}] {what}: {type(e).__name__}: "
                  f"{str(e)[:80]} - waiting {delay:.0f}s")
        time.sleep(delay)
        delay *= 2
    raise AssertionError("unreachable")


def _req(url: str, token: str, method: str = "GET", data: bytes | None = None,
         ctype: str = "application/json") -> dict:
    headers = {"Authorization": f"Bearer {token}"}
    if data is not None:
        headers["Content-Type"] = ctype

    def call() -> dict:
        r = urllib.request.Request(url, data=data, headers=headers, method=method)
        with urllib.request.urlopen(r, timeout=120) as resp:
            body = resp.read()
        return json.loads(body) if body else {}

    # GET is safe to repeat; POST/PUT/DELETE on the deposit API are not (a lost response
    # would leave us unable to tell whether it took effect).
    return _retry(call, f"{method} {url.rsplit('/', 1)[-1]}") if method == "GET" else call()


def _upload_file(bucket: str, token: str, path: Path) -> dict:
    payload = path.read_bytes()
    url = f"{bucket.rstrip('/')}/{urllib.parse.quote(path.name)}"

    def call() -> dict:
        r = urllib.request.Request(url, data=payload, headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/octet-stream",
            "Content-Length": str(len(payload)),
        }, method="PUT")
        with urllib.request.urlopen(r, timeout=600) as resp:
            return json.loads(resp.read())

    # A PUT of the same bytes is idempotent, so this one is worth retrying.
    return _retry(call, f"upload {path.name}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org")
    ap.add_argument("--token-file", type=Path, help="file holding the Zenodo token")
    ap.add_argument("--dry-run", action="store_true", help="print the plan, do nothing")
    ap.add_argument("--emit-metadata", nargs="?", const=str(META_OUT_DEFAULT), default=None,
                    help="write the metadata JSON and exit")
    ap.add_argument("--deposition-id", help="operate on an existing deposition id")
    ap.add_argument("--publish", action="store_true", help="publish the deposition")
    ap.add_argument("--replace", action="store_true",
                    help="delete same-named files already on the deposition and re-upload "
                         "(use after backfilling the reserved DOI into README.md)")
    args = ap.parse_args()

    meta = build_metadata()

    if args.emit_metadata:
        out = Path(args.emit_metadata)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps({"metadata": meta}, indent=2, ensure_ascii=False) + "\n",
                       encoding="utf-8")
        print(f"metadata written -> {out}")
        return 0

    files = upload_list()
    if not files:
        print(f"nothing to upload: {DEPOSIT} is empty — run scripts/build_data_deposit.py first",
              file=sys.stderr)
        return 2

    token = os.environ.get("ZENODO_TOKEN", "").strip()
    if not token and args.token_file:
        token = args.token_file.read_text(encoding="utf-8").strip()

    base = SANDBOX if args.sandbox else PROD

    if args.dry_run or not token:
        print(f"target      : {base}")
        print(f"record      : {'DRAFT + PUBLISH' if args.publish else 'DRAFT only'}")
        print("prereserve  : True (the draft hands out a reserved DOI at creation; "
              "legacy deposit API preserves it across the metadata PUT)")
        print(f"upload_type : {meta['upload_type']}   licence: {meta['license']}   version: {meta['version']}")
        print(f"title       : {meta['title']}")
        print(f"related     : {[r['identifier'] for r in meta['related_identifiers']]}")
        print(f"files       : {len(files)}")
        for p in files:
            print(f"  {p.stat().st_size:>12,}  {p.name}")
        if not token:
            print("\nno token found (ZENODO_TOKEN / --token-file) -> plan only, nothing sent.")
        return 0

    try:
        if args.deposition_id:
            dep = _req(f"{base}/deposit/depositions/{args.deposition_id}", token)
        else:
            # prereserve_doi must be requested AT CREATION for the legacy API to hand one out.
            dep = _req(f"{base}/deposit/depositions", token, "POST",
                       json.dumps({"metadata": {**meta, "prereserve_doi": True}}).encode("utf-8"))
        dep_id, bucket, self_url = dep["id"], dep["links"]["bucket"], dep["links"]["self"]
        print(f"deposition id = {dep_id}")

        reserved = (dep.get("metadata", {}).get("prereserve_doi") or {}).get("doi") or dep.get("doi")
        if reserved:
            print(f"reserved doi  = {reserved}   (resolves once published; do NOT delete the draft)")
            print("        -> backfill it:  python scripts/build_data_deposit.py --data-doi "
                  f"{reserved}  then re-run with --replace")
        else:
            print("no reserved doi returned — the record will still get one at publication, "
                  "but the archived README cannot quote it in advance.")

        existing = {f["filename"]: f for f in dep.get("files", [])}
        for p in files:
            if p.name in existing:
                if not args.replace:
                    print(f"  already present, skipped: {p.name}")
                    continue
                fid = existing[p.name]["id"]
                _req(f"{self_url}/files/{fid}", token, "DELETE")
                print(f"  replaced (old copy removed): {p.name}")
            _upload_file(bucket, token, p)
            print(f"  uploaded: {p.name}")

        _req(self_url, token, "PUT",
             json.dumps({"metadata": {**meta, "prereserve_doi": True}}).encode("utf-8"))
        print("metadata set")

        if args.publish:
            pub = _req(f"{base}/deposit/depositions/{dep_id}/actions/publish", token, "POST")
            print(f"PUBLISHED  doi = {pub.get('doi')}   url = {pub.get('links', {}).get('record_html')}")
        else:
            print(f"draft ready: check it, then re-run with --publish --deposition-id {dep_id}")
        return 0
    except urllib.error.HTTPError as e:
        print(f"HTTP {e.code}: {e.read().decode('utf-8', 'replace')[:1500]}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
