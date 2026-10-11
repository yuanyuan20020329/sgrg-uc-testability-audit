# -*- coding: utf-8 -*-
"""Regression test for the *resume* path and for the upload retry.

Context: the first real run of ``publish_dataset.py`` died on the 14th of 14 files with
``http.client.RemoteDisconnected`` - a hiccup of the link, not a bad request. Two things
follow from that, and both are pinned here:

  * uploading must survive a transient drop (``zenodo_deposit_data._retry``), while a
    4xx - a request that was actually rejected - must NOT be retried;
  * a failure mid-run must say "the draft is intact, resume with --deposition-id N"
    rather than dump a traceback (``publish_dataset._stage_guard``), and resuming must
    not rebuild the payload (that would put the local bytes out of step with the ones
    already uploaded).

Run:  python scripts/test_resume_upload.py
"""
from __future__ import annotations

import argparse
import http.client
import sys
import urllib.error
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import publish_dataset as pd          # noqa: E402
import zenodo_deposit_data as zd      # noqa: E402

results: list[tuple[bool, str, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append((ok, name, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))


# --------------------------------------------------------------------- _retry

def test_retry() -> None:
    print("retry (transient failures are retried, rejected requests are not)")

    calls = {"n": 0}

    def flaky():
        calls["n"] += 1
        if calls["n"] < 3:
            raise http.client.RemoteDisconnected("Remote end closed connection")
        return "ok"

    got = zd._retry(flaky, "unit-test", attempts=4, delay=0.01)
    check("flaky call succeeds on the 3rd attempt", got == "ok" and calls["n"] == 3,
          f"attempts={calls['n']}")

    def always():
        raise http.client.RemoteDisconnected("nope")

    try:
        zd._retry(always, "unit-test", attempts=3, delay=0.01)
        check("gives up after the last attempt and re-raises", False, "no exception raised")
    except http.client.RemoteDisconnected:
        check("gives up after the last attempt and re-raises", True)

    rejected = {"n": 0}

    def bad_request():
        rejected["n"] += 1
        raise urllib.error.HTTPError("http://x", 400, "BAD REQUEST", {}, None)

    try:
        zd._retry(bad_request, "unit-test", attempts=4, delay=0.01)
        check("a 400 is raised immediately, never retried", False, "no exception raised")
    except urllib.error.HTTPError:
        check("a 400 is raised immediately, never retried", rejected["n"] == 1,
              f"attempts={rejected['n']}")

    server = {"n": 0}

    def bad_gateway():
        server["n"] += 1
        if server["n"] < 2:
            raise urllib.error.HTTPError("http://x", 503, "SERVICE UNAVAILABLE", {}, None)
        return "ok"

    got = zd._retry(bad_gateway, "unit-test", attempts=4, delay=0.01)
    check("a 503 is retried", got == "ok" and server["n"] == 2, f"attempts={server['n']}")


# --------------------------------------------------------------------- _stage_guard

def test_guard() -> None:
    print("stage guard (a dead link becomes a resume instruction, not a traceback)")

    def dies():
        raise http.client.RemoteDisconnected("Remote end closed connection")

    check("RemoteDisconnected -> None (handled, not raised)",
          pd._stage_guard("stage 3 - upload", "23297703", dies) is None)

    def rejected():
        raise urllib.error.HTTPError("http://x", 401, "UNAUTHORIZED", {}, None)

    check("HTTPError -> None (handled, not raised)",
          pd._stage_guard("stage 4 - metadata", "23297703", rejected) is None)
    check("a successful stage passes its return value through",
          pd._stage_guard("stage 3 - upload", "23297703", lambda: True) is True)


# --------------------------------------------------------------------- resume skips the rebuild

def test_resume_keeps_payload() -> None:
    print("resume (stage 2 must not rebuild an already-backfilled payload)")
    manifest = pd.DEPOSIT / "manifest.json"
    before = manifest.read_bytes()
    before_n = len(list(pd.DEPOSIT.iterdir()))
    out = pd.stage_backfill(None, argparse.Namespace(data_doi=None))
    after = manifest.read_bytes()
    check("the backfilled DOI is recovered from the repo README", out == "10.5281/zenodo.23297703"
          or (out or "").startswith("10.5281/zenodo."), f"returned {out!r}")
    check("deposit/manifest.json is byte-identical afterwards (no rebuild happened)",
          before == after)
    check("no file was added to or removed from the payload",
          len(list(pd.DEPOSIT.iterdir())) == before_n)


def main() -> int:
    test_retry()
    print()
    test_guard()
    print()
    test_resume_keeps_payload()
    bad = [r for r in results if not r[0]]
    print()
    print(f"SUMMARY: {len(results) - len(bad)}/{len(results)} passed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
