# -*- coding: utf-8 -*-
"""Regression test for the new-version DOI gate (``publish_dataset.stage_newversion``).

Context: the first real run of ``--new-version-of 10.5281/zenodo.23297703`` (v1.0.0 ->
v1.0.1) created draft ``23299341``, cleared the inherited files, and then *stopped*,

    [warn] the draft reports doi='10.5281/zenodo.23299341', the same as the previous version

which is wrong on both counts. Zenodo's published docs say the response body of
``actions/newversion`` is the *original* resource; production Zenodo (InvenioRDM behind the
legacy API) returns the **new draft**, so the code compared the draft's DOI against itself
and rejected a perfectly good reservation every time. The consequence is not cosmetic: the
documented one-command release cannot finish, and the operator is pushed towards typing a
DOI by hand into an irreversible publication.

What is pinned here:

  * a reservation whose suffix is the draft's own record id is accepted, in both the
    production shape (body = draft) and the documented shape (body = original resource);
  * the previous version's DOI derived from the parent id is what it is compared against,
    so a draft that reports *that* DOI is still refused (the guard keeps its teeth);
  * ``_draft_own_doi`` only ever returns a DOI carrying the deposit's own record id, so an
    inherited or stale DOI cannot be backfilled into a payload in the resume path either.

Run:  python scripts/test_newversion_doi.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import publish_dataset as pd          # noqa: E402

NEW_ID = "23299341"
PARENT = "23297703"
NEW_DOI = f"10.5281/zenodo.{NEW_ID}"
OLD_DOI = f"10.5281/zenodo.{PARENT}"

results: list[tuple[bool, str, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append((ok, name, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))


def _draft(body_is_draft: bool, draft_doi: str | None) -> dict:
    """A canned transport for stage_newversion, in one of the two response shapes."""
    original = {"id": PARENT, "state": "done", "submitted": True,
                "metadata": {"doi": OLD_DOI},
                "links": {"latest_draft": f"https://zenodo.org/api/deposit/depositions/{NEW_ID}"}}
    draft = {
        "id": NEW_ID, "state": "unsubmitted", "submitted": False,
        "doi": draft_doi,
        "metadata": {"prereserve_doi": ({"doi": draft_doi} if draft_doi else {})},
        "files": [{"id": "f1", "filename": "README.md", "filesize": 4722},
                  {"id": "f2", "filename": "SCOPE_NOTE.md", "filesize": 2021}],
        "links": {"self": f"https://zenodo.org/api/deposit/depositions/{NEW_ID}",
                  "latest_draft": f"https://zenodo.org/api/deposit/depositions/{NEW_ID}"},
    }
    body = dict(draft) if body_is_draft else dict(original)
    if body_is_draft:
        body["links"] = {**draft["links"], **original["links"]}

    deleted: list[str] = []

    def fake(url, token, method="GET", data=None, ctype="application/json"):
        if method == "DELETE":
            deleted.append(url.rsplit("/", 1)[-1])
            return {}
        if url.endswith("/actions/newversion"):
            return body
        if url.endswith(f"/depositions/{PARENT}"):
            return original
        if url.endswith(f"/depositions/{NEW_ID}"):
            return draft
        raise AssertionError(f"unexpected call: {method} {url}")

    return {"req": fake, "deleted": deleted}


def _run(body_is_draft: bool, draft_doi: str | None) -> tuple[str | None, list[str]]:
    fake = _draft(body_is_draft, draft_doi)
    real = pd.zd._req
    pd.zd._req = fake["req"]
    try:
        _dep, reserved = pd.stage_newversion(pd.zd.PROD, "token", OLD_DOI)
    finally:
        pd.zd._req = real
    return reserved, fake["deleted"]


def test_production_shape() -> None:
    print("production shape (newversion returns the draft, as Zenodo v14 actually does)")
    reserved, deleted = _run(body_is_draft=True, draft_doi=NEW_DOI)
    check("the draft's own reserved DOI is accepted", reserved == NEW_DOI, f"got {reserved!r}")
    check("the inherited file snapshot is still cleared before upload",
          deleted == ["f1", "f2"], f"deleted={deleted}")


def test_documented_shape() -> None:
    print("documented shape (newversion returns the original resource)")
    reserved, _ = _run(body_is_draft=False, draft_doi=NEW_DOI)
    check("the draft's own reserved DOI is accepted", reserved == NEW_DOI, f"got {reserved!r}")


def test_guard_still_has_teeth() -> None:
    print("the guard still refuses a draft that has no new reservation")
    for label, doi in [("reports the previous version's DOI", OLD_DOI),
                       ("reports nothing", None),
                       ("reports a DOI from another record", "10.5281/zenodo.99999999")]:
        reserved, _ = _run(body_is_draft=True, draft_doi=doi)
        check(f"refused when the draft {label}", reserved is None, f"got {reserved!r}")


def test_draft_own_doi() -> None:
    print("_draft_own_doi (only a DOI carrying this deposit's own record id is trusted)")
    check("own reservation is returned",
          pd._draft_own_doi({"id": NEW_ID, "doi": NEW_DOI}) == NEW_DOI)
    check("an inherited/stale DOI is refused",
          pd._draft_own_doi({"id": NEW_ID, "doi": OLD_DOI}) is None)
    check("a missing DOI is refused", pd._draft_own_doi({"id": NEW_ID, "doi": None}) is None)
    check("prereserve_doi is preferred over the doi field",
          pd._draft_own_doi({"id": NEW_ID, "doi": OLD_DOI,
                             "metadata": {"prereserve_doi": {"doi": NEW_DOI}}}) == NEW_DOI)


def test_backfill_does_not_rebuild() -> None:
    print("stage_backfill with a known DOI the payload already cites (must not rebuild)")
    import argparse
    import json as _json
    import tempfile
    manifest = pd.DEPOSIT / "manifest.json"
    before = manifest.read_bytes()
    stated = str(_json.loads(before.decode("utf-8")).get("dataset_doi") or "")
    check("_payload_dataset_doi() mirrors deposit/manifest.json (never guesses)",
          pd._payload_dataset_doi() == stated, f"manifest {stated!r} vs {pd._payload_dataset_doi()!r}")

    # Never point this test at the repository's own README: stage_backfill rewrites REPO_DOCS
    # in place, and an earlier version of this test left a fake DOI in
    # resource_paper/README.md. The docs it writes to are redirected to a temp file.
    sentinel = "10.5281/zenodo.42424242"
    with tempfile.TemporaryDirectory() as td:
        fake_doc = Path(td) / "README.md"
        fake_doc.write_text(f"dataset record - {pd._recorded_dataset_doi()}\n", encoding="utf-8")
        real_docs, real_payload = pd.REPO_DOCS, pd._payload_dataset_doi
        pd.REPO_DOCS = [fake_doc]
        pd._payload_dataset_doi = lambda: sentinel
        try:
            out = pd.stage_backfill(sentinel, argparse.Namespace(data_doi=None))
        finally:
            pd.REPO_DOCS, pd._payload_dataset_doi = real_docs, real_payload
        check("the temp repo doc was rewritten to the new DOI",
              sentinel in fake_doc.read_text(encoding="utf-8"))
    check("the DOI comes back unchanged", out == sentinel, f"returned {out!r}")
    check("deposit/manifest.json is byte-identical (resume did not rebuild the payload)",
          manifest.read_bytes() == before)
    check("resource_paper/README.md was not touched by this test",
          "42424242" not in (pd.ROOT / "resource_paper" / "README.md").read_text(encoding="utf-8"))


def test_verifier_reads_the_recorded_code_doi() -> None:
    print("verify_dataset_doi.CODE_DOI is read from code_doi.txt, not hard-coded")
    import tempfile
    import build_data_deposit as bdd
    import verify_dataset_doi as vd
    recorded = bdd.read_code_doi()
    check("the verifier's code DOI equals the single recorded value",
          vd.CODE_DOI == recorded, f"verifier {vd.CODE_DOI!r} vs recorded {recorded!r}")
    # The v1.0.1 record cross-links the v1.0.1 code DOI. A literal pinned to the v1.0.0 code
    # record is exactly what made check 4 fail a perfectly correct deposit on 2026-10-11.
    check("the verifier no longer pins the superseded v1.0.0 code DOI",
          vd.CODE_DOI != bdd.CODE_DOI_V1_0_0, f"got {vd.CODE_DOI!r}")
    check("the verifier can name where its value came from",
          vd.CODE_DOI_SOURCE in ("resource_paper/code_doi.txt", "built-in v1.0.0 fallback"),
          f"got {vd.CODE_DOI_SOURCE!r}")
    with tempfile.TemporaryDirectory() as td:
        fake = Path(td) / "code_doi.txt"
        fake.write_text("10.5281/zenodo.11112222\n# concept doi: 10.5281/zenodo.11110000\n"
                        "# version: v9.9.9\n", encoding="utf-8")
        real = bdd.CODE_DOI_FILE
        bdd.CODE_DOI_FILE = fake
        try:
            check("re-releasing the code is a one-line edit: the reader follows the file",
                  bdd.read_code_doi() == "10.5281/zenodo.11112222",
                  f"got {bdd.read_code_doi()!r}")
        finally:
            bdd.CODE_DOI_FILE = real


def main() -> int:
    test_production_shape()
    print()
    test_documented_shape()
    print()
    test_guard_still_has_teeth()
    print()
    test_draft_own_doi()
    print()
    test_backfill_does_not_rebuild()
    print()
    test_verifier_reads_the_recorded_code_doi()
    bad = [r for r in results if not r[0]]
    print()
    print(f"SUMMARY: {len(results) - len(bad)}/{len(results)} passed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
