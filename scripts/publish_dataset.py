#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Publish the *dataset* record (the audit tables), first release or new version.

Runs the whole documented sequence with a checkpoint between every stage, so a partial
run can always be resumed instead of restarted:

    0. preflight   - local payload integrity + connectivity probe (no side effects)
    1. draft       - create the deposition with a pre-reserved DOI
    1n. newversion - OR: derive a new-version draft from the published record (--new-version-of)
    2. backfill    - rebuild the payload so README/manifest quote the reserved DOI
    2g. gate       - HARD GATE: audit_data_descriptor.py --pre-publish must pass (14 checks
                     of inventory, checksums, version agreement, code-DOI link, no
                     placeholders, no Han characters, and the language-only regression)
    3. upload      - every file in resource_paper/deposit/
    4. metadata    - set the dataset metadata block
    5. publish     - make the record public (and write resource_paper/dataset_doi.txt)
    6. verify      - download the public files and re-run the SHA-256 check

The gate sits **between backfill and upload**, not at the very start, because that is the
first moment the payload is the one that would actually be published: the reserved DOI is
in it, so `README`, `LICENSE-DATA.md`, `manifest.json` and `MANIFEST.sha256` can all be
checked for mutual agreement. Nothing is uploaded until it passes.

The two Zenodo records stay separate: this script only ever mints the **dataset** record.
The code record is cited from it as a related identifier, and which code *version* it points
at is read from `resource_paper/code_doi.txt` (see `build_data_deposit.read_code_doi`).

Why a script and not a handful of pasted commands: the DOI has to exist *before* the bytes
that cite it are uploaded, otherwise the archived README quotes a placeholder. Stage 2 is
exactly that backfill, and it must run between the draft and the upload.

Usage
-----
    # rehearse - nothing is sent, no token needed
    python scripts/publish_dataset.py --dry-run

    # rehearse on the sandbox (separate account, throw-away records)
    python scripts/publish_dataset.py --sandbox --token-file %TEMP%\\zenodo_token

    # the real thing, first release, all stages
    python scripts/publish_dataset.py --token-file %TEMP%\\zenodo_token

    # the real thing, NEW VERSION of the published record (v1.0.0 -> v1.0.1)
    python scripts/publish_dataset.py --token-file %TEMP%\\zenodo_token \\
            --new-version-of 10.5281/zenodo.23297703

    # ...and take the token off the disk the moment the record goes public
    python scripts/publish_dataset.py --token-file %TEMP%\\zenodo_token --forget-token

    # resume a half-finished run
    python scripts/publish_dataset.py --token-file %TEMP%\\zenodo_token \\
            --deposition-id 12345678 --stop-after upload

    # resume a new-version draft: the DOI Zenodo reserved for it is read from the draft
    # itself, so it never has to be typed (and a stale one cannot be typed by mistake)
    python scripts/publish_dataset.py --token-file %TEMP%\\zenodo_token \\
            --deposition-id 12345678

The token is read from the ZENODO_TOKEN environment variable or --token-file, and is never
echoed, logged or written anywhere. Scopes needed: deposit:write, deposit:actions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

try:  # keep non-ASCII paths printable on a cp936 console
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

import zenodo_deposit_data as zd          # noqa: E402
import verify_dataset_doi as vd            # noqa: E402
import build_data_deposit as bdd           # noqa: E402

DEPOSIT = ROOT / "resource_paper" / "deposit"
SINGLE_UNIT = ROOT / "resource_paper" / f"testability-audit-tables-v{bdd.VERSION}.zip"
DOI_RECORD = ROOT / "resource_paper" / "dataset_doi.txt"
STAGES = ["preflight", "draft", "newversion", "backfill", "gate", "upload",
          "metadata", "publish", "verify"]
UA = "sgrg-uc-testability-audit/publish_dataset.py"


def rule(title: str) -> None:
    print()
    print("=" * 74)
    print(title)
    print("=" * 74)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rel(p: Path) -> str:
    """A path for printing: relative to the project when possible, absolute otherwise."""
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def _read_token(args) -> str:
    if getattr(args, "token", None):
        return args.token.strip()
    tf = getattr(args, "token_file", None)
    if tf and Path(tf).exists():
        return Path(tf).read_text(encoding="utf-8").strip()
    import os
    return os.environ.get("ZENODO_TOKEN", "").strip()


def _manifest_check(verbose: bool = True) -> tuple[int, int]:
    ok = bad = 0
    for line in (DEPOSIT / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        h, name = line.split("  ", 1)
        name = name.strip()
        p = DEPOSIT / name
        if p.exists() and _sha(p) == h:
            ok += 1
        else:
            bad += 1
            if verbose:
                print(f"  [FAIL] checksum mismatch or missing: {name}")
    return ok, bad


def _probe_api(base: str) -> tuple[bool, str]:
    """Is the Zenodo REST API reachable *and* actually answering?

    Probing ``<base>/`` is a trap: Zenodo's application server answers that bare path with a
    plain ``404 NOT FOUND``, which urllib raises as HTTPError and which is indistinguishable
    from a blocked network unless the status code is read. Ask a real endpoint instead, and
    read both the status and the shape of the body:

      * 2xx + the record list's ``hits`` key  -> the Zenodo API itself answered;
      * 4xx                                    -> something answered, but it is not the API;
      * 5xx / timeout / DNS failure             -> unreachable.

    Read the **whole** body, never a fixed-size prefix: ``/records`` carries the search
    aggregations and runs to tens of kilobytes, so reading the first few KB truncates the
    JSON, the parse fails, and a perfectly good connection gets reported as a captive portal.
    (That false alarm cost one extra round trip; hence the size and content-type in the
    message below - if the body really is not JSON, we want to see what it is instead.)
    """
    url = f"{base}/records?size=1"

    def call():
        req = urllib.request.Request(
            url, headers={"User-Agent": UA, "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=25) as r:
            return (r.status, (r.headers.get("Content-Type") or "?").split(";")[0],
                    r.read(8_000_000))

    # A single DNS/connect blip must not be reported as "the network is down". This probe
    # passed with HTTP 200, the full publish then ran over the same link for eight minutes,
    # and the next probe died on one getaddrinfo. _retry treats URLError and 5xx as
    # transient and lets a 4xx (a real rejection) through immediately.
    try:
        status, ctype, body = zd._retry(call, "reachability probe", attempts=3, delay=2.0)
    except urllib.error.HTTPError as e:
        if 400 <= e.code < 500:
            return False, (f"[warn] {url} answered HTTP {e.code} {e.reason} - the host is "
                           f"reachable but this is not Zenodo's API "
                           f"(proxy or captive portal?)")
        return False, f"[FAIL] {url} -> HTTP {e.code} {e.reason}"
    except Exception as e:  # noqa: BLE001
        text = f"{type(e).__name__}: {str(e)[:110]}"
        if any(k in text for k in ("getaddrinfo", "11004", "Name or service not known")):
            text += ("\n         the name zenodo.org did not resolve (DNS). That is normally "
                     "a transient blip - run it again.\n         If it keeps failing, it is "
                     "the DNS or the network, not the API: try another network.")
        return False, f"[FAIL] cannot reach {url} -> {text}"
    try:
        payload = json.loads(body)
    except Exception:  # noqa: BLE001
        head = body[:60].decode("utf-8", "replace").replace("\n", " ").replace("\r", " ")
        return False, (f"[warn] {url} answered HTTP {status}, {ctype}, {len(body):,} bytes, "
                       f"but the body is not JSON: {head!r} - a proxy or captive portal is "
                       f"intercepting the API")
    if isinstance(payload, dict) and "hits" in payload:
        return True, f"[ok] reachable: {base} (HTTP {status}, records endpoint confirmed)"
    return True, f"[ok] reachable: {base} (HTTP {status})"


# --------------------------------------------------------------------- stage 0

def stage_preflight(args) -> dict:
    """Local integrity + connectivity. Returns a dict of booleans."""
    rule("stage 0 - preflight")
    res = {"local": False, "remote": False, "token": False}

    if not DEPOSIT.is_dir():
        print(f"  [FAIL] deposit directory missing: {DEPOSIT}")
        print("         run: python scripts/build_data_deposit.py")
        return res
    if not (DEPOSIT / "MANIFEST.sha256").exists():
        print("  [FAIL] MANIFEST.sha256 missing")
        return res

    ok, bad = _manifest_check()
    entries = {ln.split("  ", 1)[1].strip()
               for ln in (DEPOSIT / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines()
               if ln.strip() and not ln.startswith("#")}
    print(f"  [{'ok' if not bad else 'FAIL'}] manifest: {ok}/{len(entries)} entries verify")

    on_disk = {p.name for p in DEPOSIT.iterdir() if p.is_file()}
    payload = on_disk - {"MANIFEST.sha256", "manifest.json"}
    only_disk, only_man = sorted(payload - entries), sorted(entries - payload)
    if only_disk or only_man:
        print(f"  [FAIL] manifest does not cover the payload: "
              f"unlisted={only_disk} missing={only_man}")
        bad += 1
    else:
        print(f"  [ok] manifest covers all {len(payload)} payload files "
              f"(of {len(on_disk)} uploaded; the two manifest artifacts cannot self-reference)")

    for z in sorted(DEPOSIT.glob("*.zip")) + [SINGLE_UNIT]:
        if not z.exists():
            continue
        try:
            with zipfile.ZipFile(z) as zf:
                n, badfile = len(zf.namelist()), zf.testzip()
        except zipfile.BadZipFile:
            badfile, n = "not a zip", -1
        print(f"  [{'ok' if badfile is None else 'FAIL'}] {z.name}: "
              f"{n} members, integrity {'OK' if badfile is None else badfile}")
        if badfile is not None:
            bad += 1

    has_ph = "zenodo.XXXXXXX" in (DEPOSIT / "README.md").read_text(encoding="utf-8")
    print(f"  [{'warn' if has_ph else 'ok'}] README DOI placeholder: "
          f"{'present - stage 2 replaces it' if has_ph else 'none'}")

    res["local"] = bad == 0

    base = zd.SANDBOX if args.sandbox else zd.PROD
    reachable, msg = _probe_api(base)
    print(f"  {msg}")
    if not reachable:
        print("         the Zenodo API is required; run this where it is reachable.")
    res["remote"] = reachable

    token = _read_token(args)
    res["token"] = bool(token)
    print(f"  [{'ok' if token else 'warn'}] token: "
          f"{'found (len %d, value never shown)' % len(token) if token else 'not found'}")

    files = zd.upload_list()
    total = sum(p.stat().st_size for p in files)
    print(f"  files to upload: {len(files)}  ({total/1024/1024:.2f} MB)")
    for p in files:
        print(f"    {p.stat().st_size:>12,}  {p.name}")
    return res


# --------------------------------------------------------------------- stages

def stage_draft(token: str, base: str) -> tuple[dict, str | None]:
    rule("stage 1 - create the draft (pre-reserve a DOI)")
    meta = zd.build_metadata()
    dep = zd._req(f"{base}/deposit/depositions", token, "POST",
                  json.dumps({"metadata": {**meta, "prereserve_doi": True}}).encode("utf-8"))
    reserved = (dep.get("metadata", {}).get("prereserve_doi") or {}).get("doi") or dep.get("doi")
    print(f"  deposition id = {dep['id']}")
    print(f"  reserved doi  = {reserved or '(none returned)'}")
    print("  the reservation is only guaranteed while this draft lives - never delete it.")
    return dep, reserved


def _record_id_of(spec: str) -> str:
    """Pull the record id out of a DOI, a landing-page URL or a bare id."""
    m = re.search(r"(\d{6,})", spec or "")
    if not m:
        raise SystemExit(f"cannot read a record id out of {spec!r} "
                         f"(expected e.g. 10.5281/zenodo.23297703 or the bare id 23297703)")
    return m.group(1)


def _published_parent_id(base: str, token: str, spec: str) -> str:
    """Resolve to the id of the **latest published version**, which is what Zenodo wants.

    ``actions/newversion`` refuses the concept id — and refuses it quietly enough to waste a
    round trip, because the obvious thing to pass is the DOI a reader would copy. Passing the
    concept record's id returns the published *child*, so we look up the newest version and
    use that instead.
    """
    rid = _record_id_of(spec)
    try:
        rec = zd._req(f"{base}/deposit/depositions/{rid}", token)
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
        rec = zd._req(f"{base}/records/{rid}", token)
    if rec.get("state") == "done" or rec.get("submitted") is True:
        return str(rec.get("id", rid))
    # A concept record has no files of its own; follow its newest version instead.
    versions = ((rec.get("links") or {}).get("versions")) or []
    if versions:
        newest = versions[-1]
        return _record_id_of(newest.get("href", "") if isinstance(newest, dict) else str(newest))
    latest = (rec.get("links") or {}).get("latest")
    if latest:
        return _record_id_of(latest)
    return rid


def stage_newversion(base: str, token: str, spec: str) -> tuple[dict, str | None]:
    """Derive a new-version draft from a published record (v1.0.0 -> v1.0.1).

    Three behaviours of the legacy deposit API are load-bearing here, and each of them is
    silent when you get it wrong:

    * ``actions/newversion`` takes the id of the **latest published version**, not the
      concept id;
    * the response describes the *previous* record — the new draft is only reachable through
      ``links.latest_draft``;
    * the draft inherits a **snapshot of the previous version's files**, so uploading without
      clearing them first would publish both the v1.0.0 and the v1.0.1 copies of every table.

    The new version's DOI is read from the draft. If Zenodo hands one out, the payload can be
    backfilled with it before upload (the whole point); if it does not, we stop here and say
    so rather than publish bytes that cite a placeholder.
    """
    rule("stage 1n - create a new version of the published record")
    parent = _published_parent_id(base, token, spec)
    print(f"  previous version (deposition id) = {parent}")
    resp = zd._req(f"{base}/deposit/depositions/{parent}/actions/newversion", token, "POST")
    draft_url = (resp.get("links") or {}).get("latest_draft")
    if not draft_url:
        raise SystemExit("newversion returned no links.latest_draft - nothing was created")
    draft = zd._req(draft_url, token)
    new_id = str(draft["id"])
    print(f"  new draft id = {new_id}   (the published record is untouched)")

    inherited = draft.get("files") or []
    for f in inherited:
        zd._req(f"{draft['links']['self']}/files/{f['id']}", token, "DELETE")
    print(f"  cleared {len(inherited)} inherited file(s) so the deposit holds only v"
          f"{bdd.VERSION} bytes")

    # Which DOI this draft will be published under.
    #
    # Zenodo's published docs say the response body of ``actions/newversion`` is the
    # *original* resource, so the "previous" DOI was read from it. On production
    # (InvenioRDM behind the legacy API, v14) the body is the **new draft**: reading both
    # sides out of it compares the draft with itself and rejects every honest reservation -
    # which is exactly what happened on the first real v1.0.1 run. The previous version's
    # DOI is derived from ``parent`` instead, and a reservation is trusted only when its
    # suffix is this draft's own record id: Zenodo mints ``10.5281/zenodo.<recid>`` and
    # reserves it as soon as the new version exists, so a DOI carrying any other suffix is
    # an inherited or stale value that must never be baked into the payload.
    meta = draft.get("metadata") or {}
    reserved = (meta.get("prereserve_doi") or {}).get("doi") or draft.get("doi")
    prev_doi = f"10.5281/zenodo.{parent}"
    suffix = reserved.rsplit(".", 1)[-1] if reserved else ""
    if reserved and suffix == new_id and reserved != prev_doi:
        print(f"  reserved doi (new version) = {reserved}")
    else:
        if reserved:
            print(f"  [warn] the draft reports doi={reserved!r}, which is not a new "
                  f"reservation for draft {new_id} (the previous version is {prev_doi}).")
        else:
            print("  [warn] the draft reports no DOI at all.")
        print("         Zenodo may only mint the new version DOI at publication, so "
              "publishing now would archive a README that quotes a placeholder.")
        print("         Do NOT publish yet. Open the draft, reserve the DOI, then re-run "
              "with")
        print(f"           --deposition-id {new_id} --data-doi <the reserved DOI>")
        reserved = None
    return draft, reserved


def stage_gate(doi: str | None, args) -> bool:
    """The pre-publish hard gate.

    ``audit_data_descriptor.py --pre-publish`` re-derives the payload's own facts (file
    inventory, checksums, version agreement across manifest / datapackage / README /
    CHANGELOG, the code-DOI link, no leftover placeholder, no Han characters) and runs the
    language-only regression against the frozen v1.0.0 payload. It is deliberately the last
    thing before the upload: after this point everything is a network operation, and a
    published record cannot be made to forget a mistake.
    """
    rule("stage 2g - pre-publish gate (hard)")
    cmd = [sys.executable, str(SCRIPTS / "audit_data_descriptor.py"), "--pre-publish"]
    if doi and "XXXXXXX" not in doi:
        cmd += ["--data-doi", doi]
    print("  $ " + " ".join(Path(c).name if i == 1 else c for i, c in enumerate(cmd)))
    r = subprocess.run(cmd, cwd=str(ROOT))
    if r.returncode:
        print("\n  [FAIL] the gate did not pass. Nothing was uploaded; the draft is intact.")
        print("         Fix the payload, rebuild it (scripts/build_data_deposit.py "
              f"--data-doi {doi or '<doi>'}), then resume with --deposition-id <id>.")
        return False
    print("  [ok] gate passed - the payload is the one we intend to publish")
    return True


# Only the *repository note* is rewritten in place. `LICENSE-DATA.md` is deliberately NOT
# here: it is a build-time template (`10.5281/zenodo.XXXXXXX` / `…YYYYYYY`), so writing a real
# DOI into it would make that DOI outlive the release that minted it — the exact bug this
# re-release would otherwise have shipped.
REPO_DOCS = [ROOT / "resource_paper" / "README.md"]

# The dataset-DOI placeholder and any real dataset DOI look alike, so a blanket regex would
# also rewrite the *code* DOI and the dataset *concept* DOI sitting in the same two files.
# Substitution is therefore targeted at exactly one string: the dataset version DOI this
# repository currently records, plus the placeholder.
_DATASET_DOI_RX = re.compile(r"10\.5281/zenodo\.(?:XXXXXXX|\d{6,})")


def _code_dois() -> set[str]:
    rec = bdd.read_code_doi_record()
    return {x for x in (rec.get("doi"), rec.get("concept_doi")) if x}


def _recorded_dataset_concept() -> str:
    """The dataset **concept** DOI, which must never be rewritten."""
    if not DOI_RECORD.exists():
        return ""
    for line in DOI_RECORD.read_text(encoding="utf-8").splitlines():
        m = re.match(r"#\s*concept doi:\s*(\S+)", line.strip(), flags=re.I)
        if m:
            return m.group(1)
    return ""


def _recorded_dataset_doi() -> str:
    """The dataset **version** DOI currently recorded in ``dataset_doi.txt``."""
    if not DOI_RECORD.exists():
        return ""
    for line in DOI_RECORD.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            return s
    return ""


def _payload_dataset_doi() -> str:
    """The dataset DOI the payload on disk was last built with ("" if it cannot be read)."""
    try:
        man = json.loads((DEPOSIT / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ""
    return str(man.get("dataset_doi") or "")


def stage_backfill(doi: str | None, args) -> str | None:
    """Replace the dataset DOI in the repository docs, then rebuild the payload."""
    rule("stage 2 - backfill the DOI (repo docs, then the payload)")
    doi = args.data_doi or doi
    if not doi:
        # Resuming a half-finished run: stage 2 has normally already happened, so look for
        # the backfilled DOI on disk instead of assuming the placeholder is still there.
        # Crucially, do NOT rebuild here - rebuilding would rewrite manifest.json and leave
        # the local bytes out of step with the ones already uploaded.
        repo_txt = REPO_DOCS[0].read_text(encoding="utf-8")
        found = re.search(r"10\.5281/zenodo\.\d{6,}", repo_txt)
        if found and "zenodo.XXXXXXX" not in repo_txt:
            print(f"  [ok] already backfilled with {found.group(0)} - payload left untouched")
            return found.group(0)
        print("  [warn] no DOI available and no --data-doi: README keeps the placeholder.")
        return None
    if "XXXXXXX" in doi:
        print(f"  [warn] '{doi}' is a placeholder; skipping the backfill.")
        return None

    prev = _recorded_dataset_doi()
    superseded = {prev} - {"", doi}
    if superseded:
        print(f"  superseding the recorded dataset version DOI: {', '.join(sorted(superseded))}")

    olds = sorted(superseded | {"10.5281/zenodo.XXXXXXX"})
    for p in REPO_DOCS:
        txt = p.read_text(encoding="utf-8")
        new = txt
        for old in olds:
            new = new.replace(old, doi)
        if new == txt:
            print(f"  [warn] {_rel(p)} cites no dataset DOI — left alone")
            continue
        p.write_text(new, encoding="utf-8", newline="\n")
        print(f"  [ok] {_rel(p)} now cites {doi}")

    # Anything else still sitting in these two files is either the concept DOI, a code DOI,
    # or a stale dataset DOI we were not told about — the last of which is exactly the slip
    # a re-release invites, so name it rather than pass it through.
    allowed = {doi, _recorded_dataset_concept()} | _code_dois()
    stray = sorted({d for p in REPO_DOCS if p.exists()
                    for d in _DATASET_DOI_RX.findall(p.read_text(encoding="utf-8"))
                    if d not in allowed and "XXXXXXX" not in d})
    if stray:
        print(f"  [warn] repository docs still cite DOI(s) that are neither the new version, "
              f"the concept DOI nor the code record: {stray}")
        print("         check resource_paper/README.md and LICENSE-DATA.md by hand.")

    # A resume can arrive with the DOI already baked into the payload (stage 2 ran, the
    # upload died). Rebuilding then is not a no-op: manifest.json gets rewritten, and a byte
    # that changed while keeping its length would be skipped by the upload stage on a size
    # comparison - leaving the deposit holding the earlier copy. So when the payload already
    # cites this DOI, the repository docs above are still updated but nothing is rebuilt; the
    # stage-2g gate re-checks the payload immediately after this.
    if _payload_dataset_doi() == doi:
        print(f"  [ok] the payload already cites {doi} - no rebuild; the gate re-checks it")
        return doi

    cmd = [sys.executable, str(SCRIPTS / "build_data_deposit.py"), "--data-doi", doi]
    print(f"  $ build_data_deposit.py --data-doi {doi}")
    r = subprocess.run(cmd, cwd=str(ROOT))
    if r.returncode:
        print(f"  [FAIL] rebuild exited {r.returncode}")
        return None
    man = json.loads((DEPOSIT / "manifest.json").read_text(encoding="utf-8"))
    print(f"  [ok] manifest.json dataset_doi = {man.get('dataset_doi')}")
    lic = (DEPOSIT / "LICENSE-DATA.md").read_text(encoding="utf-8")
    stale = [d for d in superseded if d in lic]
    print(f"  [{'ok' if doi in lic and 'XXXXXXX' not in lic and not stale else 'FAIL'}] "
          f"deposit LICENSE-DATA.md cites {doi} only (internal note removed)")
    if stale:
        print(f"         still carries the superseded {stale}")
    ok, bad = _manifest_check()
    print(f"  [{'ok' if not bad else 'FAIL'}] post-backfill manifest: {ok} ok / {bad} bad")
    return doi if (not bad and not stale and "XXXXXXX" not in lic) else None


def stage_upload(dep: dict, token: str, args) -> bool:
    """Upload the payload so the deposition ends up holding *exactly* these bytes.

    The danger this stage is written against: a new version of a Zenodo record inherits a
    snapshot of the previous version's files, and "the file is already there" is not evidence
    that it is *ours*. Skipping same-named files is how a v1.0.1 record quietly ships v1.0.0
    content — the deposit says v1.0.1, the bytes say v1.0.0, and nothing in the run's output
    looks wrong. So the rule here is that the local payload wins:

    * a same-named file whose **size differs** is always replaced — a size mismatch is proof
      the deposit copy is stale, no flag needed;
    * ``--replace`` (set automatically by ``--new-version-of``) forces replacement even when
      the sizes agree, for the case where only the content moved;
    * a file in the deposit that is **not in the payload** is deleted, so a file removed or
      renamed between versions cannot survive into the new one.

    The one case that is deliberately cheap is resuming: a same-named file with an identical
    size and no ``--replace`` is left alone, which is what makes a crashed upload resumable.
    """
    rule("stage 3 - upload every file")
    self_url, bucket = dep["links"]["self"], dep["links"]["bucket"]
    payload = zd.upload_list()
    existing = {f["filename"]: f for f in dep.get("files", [])}
    force = bool(getattr(args, "replace", False) or getattr(args, "new_version_of", None))

    n = skipped = 0
    for p in payload:
        cur = existing.get(p.name)
        if cur:
            size = cur.get("filesize") or cur.get("size")
            mine = p.stat().st_size
            if force or size != mine:
                zd._req(f"{self_url}/files/{cur['id']}", token, "DELETE")
                if force and size == mine:
                    print(f"  replaced on request (same size, content may differ): {p.name}")
                else:
                    print(f"  replaced stale copy ({size} -> {mine} bytes): {p.name}")
            else:
                print(f"  already present at {mine:,} bytes, skipped (resume): {p.name}")
                skipped += 1
                continue
        zd._upload_file(bucket, token, p)
        print(f"  uploaded {p.stat().st_size:>12,}  {p.name}")
        n += 1

    stray = [name for name in existing if name not in {p.name for p in payload}]
    for name in stray:
        zd._req(f"{self_url}/files/{existing[name]['id']}", token, "DELETE")
        print(f"  removed file that is not part of this payload: {name}")

    print(f"  [ok] {n} uploaded, {skipped} skipped, {len(stray)} stray removed")
    if force:
        print("       replacement was forced, so the deposit now holds this payload's bytes "
              "and nothing else")
    return True


def stage_metadata(dep: dict, token: str) -> bool:
    rule("stage 4 - set the metadata block")
    meta = zd.build_metadata()
    zd._req(dep["links"]["self"], token, "PUT",
            json.dumps({"metadata": {**meta, "prereserve_doi": True}}).encode("utf-8"))
    print(f"  [ok] upload_type={meta['upload_type']}  licence={meta['license']}  "
          f"version={meta['version']}")
    print(f"  [ok] related: {[r['identifier'] for r in meta['related_identifiers']]}")
    return True


def _record_doi(doi: str, record_url: str | None, concept_doi: str | None = None) -> None:
    """Leave a machine-readable trace of the minted DOI next to the payload.

    ``verify_dataset_doi.py --doi-file`` reads the first line, so the DOI goes first and the
    provenance follows as ``#`` comments. This is the one place the dataset DOI is written
    down after publication: the manuscript, ``audit_data_descriptor.py`` and later
    re-verifications all take it from here rather than from a console transcript.

    The concept DOI and the payload version are recorded alongside it, because both are
    checked later: the concept DOI is what a future version has to keep, and the version is
    what tells a re-release apart from a no-op republish of the same bytes.
    """
    lines = [doi]
    if concept_doi:
        lines.append(f"# concept doi: {concept_doi}")
    if record_url:
        lines.append(f"# landing page: {record_url}")
    lines.append(f"# record id: {doi.rsplit('.', 1)[-1]}")
    lines.append(f"# published: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}")
    lines.append(f"# version: {bdd.VERSION}")
    DOI_RECORD.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    try:
        shown = DOI_RECORD.relative_to(ROOT)
    except ValueError:
        shown = DOI_RECORD
    print(f"  [ok] DOI recorded -> {shown}  (version {bdd.VERSION})")
    print(f"       verify later with: verify_dataset_doi.py --doi-file {shown}")


def stage_publish(base: str, rid: str, token: str, record_doi: bool = True,
                  expect: str | None = None) -> str | None:
    rule("stage 5 - publish")
    pub = zd._req(f"{base}/deposit/depositions/{rid}/actions/publish", token, "POST")
    doi = pub.get("doi")
    concept = pub.get("conceptdoi")
    print(f"  [ok] PUBLISHED  doi = {doi}")
    if concept:
        print(f"       concept doi = {concept} (unchanged across versions)")
    print(f"       record page: {pub.get('links', {}).get('record_html') or pub.get('links', {}).get('record_url')}")
    # Last line of defence: the uploaded README/manifest quote the DOI we backfilled, and
    # those bytes are immutable from here on. If Zenodo minted a different one, say so now -
    # while the remedy is still "publish another new version", not "explain the citation".
    if doi and expect and doi != expect:
        print(f"  [FAIL] the record was published as {doi}, but the uploaded payload cites "
              f"{expect}.")
        print("         the deposited bytes cannot be edited; publish a further new version "
              "that cites the right DOI before citing this one anywhere.")
    if doi and record_doi:
        _record_doi(doi, (pub.get("links") or {}).get("record_html")
                    or (pub.get("links") or {}).get("record_url"), concept)
    elif doi:
        print("  [note] sandbox record - the DOI was not written to dataset_doi.txt")
    return doi


def stage_verify(doi: str) -> bool:
    rule("stage 6 - verify the public record")
    r = subprocess.run([sys.executable, str(SCRIPTS / "verify_dataset_doi.py"), doi],
                       cwd=str(ROOT))
    return r.returncode == 0


def _under_temp(p: Path) -> bool:
    """True only for a file that really lives in the system temp directory."""
    try:
        tmp = Path(tempfile.gettempdir()).resolve()
        rp = p.resolve()
        return rp == tmp or tmp in rp.parents
    except OSError:
        return False


def token_hygiene(args) -> None:
    """The token has done its job. Get it off local disk if the user asked for that.

    Deleting is deliberately conservative: ``--forget-token`` only ever removes a file
    that sits in the system temp directory, so it cannot nuke a token store you keep
    somewhere else by accident. ``--force-forget`` overrides that guard.
    """
    tf = getattr(args, "token_file", None)
    if not tf:
        return
    p = Path(tf)
    if not getattr(args, "forget_token", False):
        print(f"  [note] token file still on disk: {p}")
        print("         remove it now, or re-run with --forget-token.")
        return
    if not _under_temp(p) and not getattr(args, "force_forget", False):
        print(f"  [warn] --forget-token only removes files inside the temp dir "
              f"({Path(tempfile.gettempdir())}).")
        print(f"         {p} is outside it - delete it by hand, or pass --force-forget "
              f"if you are certain.")
        return
    try:
        p.unlink()
        print(f"  [ok] token file removed: {p}")
    except OSError as e:  # noqa: BLE001
        print(f"  [warn] could not remove {p}: {e}")


def _draft_own_doi(dep: dict) -> str | None:
    """The DOI *this* deposition will be published under, or None if that cannot be told.

    Prefer the deposition's own answer over any DOI typed on the command line or inferred
    from the repository. Only a DOI whose suffix is this deposit's own record id is
    returned: Zenodo mints ``10.5281/zenodo.<recid>``, so anything else is an inherited or
    stale value and baking it into the payload would archive a README citing a DOI this
    record will never have - which is unfixable once published.
    """
    doi = ((dep.get("metadata") or {}).get("prereserve_doi") or {}).get("doi") or dep.get("doi")
    rid = str(dep.get("id") or dep.get("record_id") or "")
    if doi and rid and doi.rsplit(".", 1)[-1] == rid:
        return str(doi)
    return None


# --------------------------------------------------------------------- main

def _stage_guard(stage: str, rid: str, fn):
    """Run a stage, turning a dead connection into a message you can act on.

    An unhandled ``RemoteDisconnected`` mid-upload prints a thirty-line traceback and says
    nothing about the state of the deposit. What matters at that moment is: nothing is lost,
    the draft is intact, and here is how to pick it up.
    """
    try:
        return fn()
    except urllib.error.HTTPError as e:  # a subclass of URLError - keep this branch first
        detail = " " + e.read().decode("utf-8", "replace")[:400]
        print(f"\n  [FAIL] {stage} -> HTTP {e.code}: {e.reason}{detail}")
    except (urllib.error.URLError, OSError) as e:
        print(f"\n  [FAIL] {stage} -> {type(e).__name__}: {str(e)[:160]}")
    print(f"  nothing is lost: the draft is intact. Resume with --deposition-id {rid}")
    print("       (files already uploaded are skipped, so it continues where it stopped.")
    print("        If resource_paper/dataset_doi.txt exists, it already published - stop.)")
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org")
    ap.add_argument("--token", help="token on the command line (prefer --token-file or env)")
    ap.add_argument("--token-file", type=Path, help="file holding the token")
    ap.add_argument("--dry-run", action="store_true", help="stage 0 only; nothing is sent")
    ap.add_argument("--deposition-id", help="resume an existing draft")
    ap.add_argument("--new-version-of",
                    help="publish a NEW VERSION of a published record: the version DOI, the "
                         "landing-page URL or the bare record id (e.g. 10.5281/zenodo.23297703)")
    ap.add_argument("--data-doi", help="use this DOI instead of the pre-reserved one")
    ap.add_argument("--replace", action="store_true", help="re-upload same-named files")
    ap.add_argument("--forget-token", action="store_true",
                    help="delete the --token-file once published (temp-dir files only)")
    ap.add_argument("--force-forget", action="store_true",
                    help="allow --forget-token to delete a token file outside the temp dir")
    ap.add_argument("--stop-after", choices=STAGES, help="stop once this stage has finished")
    ap.add_argument("--assume-reachable", action="store_true",
                    help="skip the reachability gate (only if the probe false-alarms: the "
                         "first real API call then fails loudly instead)")
    ap.add_argument("--skip-gate", action="store_true",
                    help="NOT recommended: bypass the stage-2g pre-publish gate")
    args = ap.parse_args()

    at = STAGES.index
    base = zd.SANDBOX if args.sandbox else zd.PROD
    print(f"target : {base}")
    print(f"deposit: {DEPOSIT}")
    print(f"version: {bdd.VERSION}")
    print(f"code   : {bdd.read_code_doi()}  (from resource_paper/code_doi.txt)")

    stop = at(args.stop_after) if args.stop_after else len(STAGES) - 1

    res = stage_preflight(args)
    if args.dry_run:
        print(f"\ndry run: local={'ok' if res['local'] else 'FAIL'} "
              f"remote={'ok' if res['remote'] else 'unreachable'} "
              f"token={'ok' if res['token'] else 'missing'} - nothing was sent.")
        return 0 if res["local"] else 1
    if stop == at("preflight"):
        print("\nstop-after=preflight: nothing further attempted.")
        return 0 if res["local"] else 1
    if not res["local"]:
        print("\npreflight failed locally - fix the payload first. Nothing was sent.")
        return 1
    if not res["remote"] and not args.assume_reachable:
        print("\nthe Zenodo API is unreachable from here; nothing was sent.")
        print("       if you know the network is fine and the probe is wrong, re-run with")
        print("       --assume-reachable: the first real API call will then fail loudly.")
        return 3

    token = _read_token(args)
    if not token:
        print("\nno token (ZENODO_TOKEN / --token-file). Nothing was sent.")
        return 2

    rid, doi, dep = args.deposition_id, args.data_doi, None

    if not rid and args.new_version_of:
        try:
            dep, reserved = stage_newversion(base, token, args.new_version_of)
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            detail = (" " + e.read().decode("utf-8", "replace")[:400]
                      if isinstance(e, urllib.error.HTTPError) else "")
            print(f"\n  [FAIL] could not create the new version -> {type(e).__name__}: "
                  f"{e}{detail}")
            print("  nothing was published; the previous version is untouched, and any "
                  "half-made draft can be discarded from the web interface.")
            return 1
        rid = str(dep["id"])
        doi = doi or reserved
        if not doi:
            print(f"\nstopped: no new-version DOI is available yet, so nothing was uploaded.")
            print(f"         the draft id {rid} is reserved for you - finish it with:")
            print(f"           python scripts/publish_dataset.py --deposition-id {rid} "
                  f"--data-doi <reserved DOI>")
            return 0
    elif stop >= at("draft") and not rid:
        try:
            dep, reserved = stage_draft(token, base)
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            detail = ""
            if isinstance(e, urllib.error.HTTPError):
                detail = " " + e.read().decode("utf-8", "replace")[:400]
            print(f"\n  [FAIL] could not create the draft -> {type(e).__name__}: {e}{detail}")
            print("  nothing was created on Zenodo; check the network and the token, then "
                  "re-run. (A 401/403 here means the token or its scopes, not the network.)")
            return 1
        rid = str(dep["id"])
        doi = doi or reserved
    elif rid:
        print(f"\nresuming existing deposition {rid}")

    if stop == at("draft"):
        print(f"\nstopped after draft. resume with --deposition-id {rid}")
        return 0
    if stop == at("newversion"):
        print(f"\nstopped after newversion. resume with --deposition-id {rid}")
        return 0
    if dep is None:
        dep = _stage_guard(f"stage 1b - fetch draft {rid}", rid,
                           lambda: zd._req(f"{base}/deposit/depositions/{rid}", token))
        if dep is None:
            return 1
    if dep.get("submitted") is True or dep.get("state") == "done":
        print(f"\n  [FAIL] deposition {rid} is already published - there is nothing to do.")
        print("         to publish a new version, use --new-version-of <record id>")
        return 1

    # Resuming: the deposition is the authority on the DOI it will be published under, and
    # that is the one number the operator would otherwise have to supply by hand. The
    # fallback in stage_backfill (reuse the DOI the repository README already cites) is only
    # a guess, and on a re-release it guesses the *previous* version's DOI.
    if args.deposition_id:
        own = _draft_own_doi(dep)
        if own and args.data_doi and args.data_doi != own:
            print(f"\n  [FAIL] --data-doi {args.data_doi}, but this draft will be published "
                  f"as {own}.")
            print("         refusing to bake a DOI that is not the one Zenodo will mint.")
            return 1
        if own and not doi:
            doi = own
            print(f"  [ok] the draft carries its own DOI: {doi} (read from Zenodo, not typed)")
        if not own and not doi:
            print("\n  [FAIL] this draft does not report a DOI of its own and no --data-doi "
                  "was given.")
            print("         refusing to fall back on the DOI the repository README happens to "
                  "cite:")
            print("         on a re-release that is the PREVIOUS version's DOI, and it would "
                  "be baked")
            print("         into bytes that cannot be corrected after publication.")
            print("         reserve the DOI on the draft page, or pass "
                  "--data-doi <the DOI it will get>.")
            return 1

    if stop >= at("backfill"):
        out = stage_backfill(doi, args)
        if out:
            doi = out
    if stop == at("backfill"):
        print(f"\nstopped after backfill. resume with --deposition-id {rid}")
        return 0

    if stop >= at("gate"):
        if args.skip_gate:
            print("\n  [warn] --skip-gate: the pre-publish gate was bypassed by request.")
        elif not stage_gate(doi, args):
            return 1
    if stop == at("gate"):
        print(f"\nstopped after gate. resume with --deposition-id {rid}")
        return 0

    if stop >= at("upload"):
        if _stage_guard("stage 3 - upload", rid,
                        lambda: stage_upload(dep, token, args)) is None:
            return 1
    if stop == at("upload"):
        print(f"\nstopped after upload. resume with --deposition-id {rid}")
        return 0

    if stop >= at("metadata"):
        if _stage_guard("stage 4 - metadata", rid,
                        lambda: stage_metadata(dep, token)) is None:
            return 1
    if stop == at("metadata"):
        print(f"\nstopped after metadata. resume with --deposition-id {rid}")
        return 0

    published = None
    if stop >= at("publish"):
        try:
            published = stage_publish(base, rid, token, record_doi=not args.sandbox,
                                      expect=doi)
        except urllib.error.HTTPError as e:
            print(f"  [FAIL] HTTP {e.code}: {e.read().decode('utf-8', 'replace')[:800]}")
            print(f"  the draft is intact - fix and re-run with --deposition-id {rid}")
            return 1
    if published:
        token_hygiene(args)
        if doi and published != doi:
            print("\n  the published DOI and the DOI quoted by the payload disagree "
                  "(see above). Not verifying; treat this release as unusable.")
            return 1
    if stop == at("publish"):
        return 0

    target = published or doi
    if not target:
        print("\nno DOI to verify against.")
        return 1
    if args.sandbox:
        print("\nsandbox record published; production verification is a separate run.")
        return 0
    return 0 if stage_verify(target) else 1


if __name__ == "__main__":
    raise SystemExit(main())
