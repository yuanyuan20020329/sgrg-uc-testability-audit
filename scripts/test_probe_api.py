# -*- coding: utf-8 -*-
"""Regression test for ``publish_dataset._probe_api()``.

Two false alarms came out of this probe before it was right, and both said the network was
broken when it was not:

  1. it requested ``<base>/`` (a bare path Zenodo answers with 404), so a healthy connection
     surfaced as ``HTTP Error 404: NOT FOUND``;
  2. it read only the first 4,096 bytes of ``/records``, whose body runs to ~56 KB because of
     the search aggregations - the truncated JSON failed to parse and was reported as "a proxy
     or captive portal is intercepting the API".

Each branch is exercised here against a fake ``urlopen``, so the distinction the probe has to
make - *blocked* (5xx / no route) versus *answered, just not by Zenodo* (2xx-non-JSON) versus
*fine* - stays pinned.

Run:  python scripts/test_probe_api.py
"""
from __future__ import annotations

import json
import socket
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import publish_dataset as pd  # noqa: E402

BASE = pd.zd.PROD
URL = f"{BASE}/records?size=1"
REAL_URLOPEN = urllib.request.urlopen
results: list[tuple[bool, str, str]] = []


class FakeResp:
    def __init__(self, body: bytes, status: int = 200, ctype: str = "application/json"):
        self._body, self.status = body, status
        self.headers = {"Content-Type": ctype}

    def read(self, n: int = -1) -> bytes:
        return self._body if n is None or n < 0 else self._body[:n]

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _boom(exc):
    def f(req, *a, **k):
        raise exc
    return f


def with_urlopen(fn):
    urllib.request.urlopen = fn
    try:
        return pd._probe_api(BASE)
    finally:
        urllib.request.urlopen = REAL_URLOPEN


def _ok(body: bytes, status: int = 200, ctype: str = "application/json"):
    """A urlopen stand-in: swallows the timeout kwarg the probe passes."""
    def f(req, *a, **k):
        return FakeResp(body, status, ctype)
    return f


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append((ok, name, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))


def main() -> int:
    print("probe url under test:", URL)
    print()

    # 1. the real shape: a body far larger than the 4,096-byte prefix
    big = json.dumps({
        "hits": {"total": 7_458_241, "hits": []},
        "aggregations": {"subject": {"buckets": [{"key": "x" * 40, "doc_count": i}
                                                 for i in range(800)]}},
    }).encode()
    print(f"case 1 body size: {len(big):,} bytes (the old probe read only the first 4,096)")
    truncated = False
    try:
        json.loads(big[:4096])
    except Exception:
        truncated = True
    check("case 1: a 4,096-byte prefix really does truncate this body (reproduces bug 2)",
          truncated and len(big) > 4096)
    ok, msg = with_urlopen(_ok(big))
    check("case 1: the full-body JSON is accepted as reachable",
          ok and "records endpoint confirmed" in msg, msg)

    # 2. HTTP 200 but HTML: a genuine captive portal
    ok, msg = with_urlopen(_ok(b"<html><body>Sign in to WiFi</body></html>", ctype="text/html"))
    check("case 2: HTTP 200 HTML is rejected, and the message says what it actually saw",
          (not ok) and "not JSON" in msg and "text/html" in msg and "bytes" in msg, msg)

    # 3. any 4xx: something answered, but not Zenodo
    ok, msg = with_urlopen(
        _boom(urllib.error.HTTPError(URL, 404, "NOT FOUND", {}, None)))
    check("case 3: HTTP 404 -> warn (not FAIL), and the host is called reachable",
          (not ok) and msg.startswith("[warn]") and "404" in msg, msg)

    # 4. 5xx: the proxy really is blocking (what this sandbox emits)
    ok, msg = with_urlopen(
        _boom(urllib.error.HTTPError(URL, 502, "Bad Gateway", {}, None)))
    check("case 4: HTTP 502 -> hard FAIL",
          (not ok) and msg.startswith("[FAIL]") and "502" in msg, msg)

    # 5. no route at all
    ok, msg = with_urlopen(_boom(urllib.error.URLError("Tunnel connection failed")))
    check("case 5: URLError -> hard FAIL",
          (not ok) and msg.startswith("[FAIL]") and "cannot reach" in msg, msg)

    # 6. 200 with JSON that is not the record list
    ok, msg = with_urlopen(_ok(json.dumps({"message": "ok"}).encode()))
    check("case 6: 200 JSON without 'hits' is still accepted (the host answered)",
          ok and "records endpoint confirmed" not in msg, msg)

    # 7. one DNS blip must not be reported as "the network is down"
    blip = {"n": 0}

    def dns_blip(req, *a, **k):
        blip["n"] += 1
        if blip["n"] == 1:
            raise urllib.error.URLError(socket.gaierror(11004, "getaddrinfo failed"))
        return FakeResp(json.dumps({"hits": {"total": 1, "hits": []}}).encode())

    ok, msg = with_urlopen(dns_blip)
    check("case 7: a single getaddrinfo failure is retried, not reported as unreachable",
          ok and blip["n"] == 2, f"attempts={blip['n']}")

    # 8. persistent DNS failure: fail, but say it is DNS
    def always_dns(req, *a, **k):
        raise urllib.error.URLError(socket.gaierror(11004, "getaddrinfo failed"))

    ok, msg = with_urlopen(always_dns)
    check("case 8: persistent DNS failure fails with a DNS-specific hint",
          (not ok) and "getaddrinfo" in msg and "DNS" in msg, msg)

    bad = [r for r in results if not r[0]]
    print()
    print(f"SUMMARY: {len(results) - len(bad)}/{len(results)} passed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
