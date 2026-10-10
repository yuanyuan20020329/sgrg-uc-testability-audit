# -*- coding: utf-8 -*-
"""Probe Gemma datasets: pull all samples (cursor paging) and summarise disease x tissue.

Usage: python gemma_cohort_probe.py <dataset_id> [<dataset_id> ...]
Writes data/tmp/samp/<id>.json (cached) and prints composition.
"""
import collections
import json
import os
import subprocess
import sys

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
CACHE = os.path.join(ROOT, "data", "tmp", "samp")
os.makedirs(CACHE, exist_ok=True)
BASE = "https://gemma.msl.ubc.ca/rest/v2/datasets/%s/samples"


def _curl(url, out):
    subprocess.run(["curl", "-s", "-m", "120", "--compressed", url, "-o", out], check=False)


def fetch_all(did):
    cache = os.path.join(CACHE, "%s.json" % did)
    if os.path.exists(cache):
        return json.load(open(cache, encoding="utf-8"))
    all_s, cursor = [], None
    for _ in range(60):
        url = BASE % did + "?limit=100" + ("&cursor=%s" % cursor if cursor else "")
        tmp = os.path.join(CACHE, "_p.json")
        _curl(url, tmp)
        try:
            d = json.load(open(tmp, encoding="utf-8"))
        except Exception:
            break
        data = d.get("data") or []
        if not data:
            break
        all_s += data
        cur = d.get("nextCursor")
        if not cur:
            break
        cursor = cur
    json.dump(all_s, open(cache, "w", encoding="utf-8"), ensure_ascii=False)
    return all_s


def summarise(samples):
    dis, part = collections.Counter(), collections.Counter()
    for e in samples:
        m = {}
        for c in ((e.get("sample") or {}).get("characteristics") or []):
            m.setdefault((c.get("category") or "").lower(), set()).add(c.get("value"))
        for k in ("disease", "disease state", "diagnosis", "condition", "tumor grading"):
            if k in m:
                dis[";".join(sorted(x for x in m[k] if x))] += 1
        for k in ("organism part", "tissue", "tissue type", "isolated tissue", "organ"):
            if k in m:
                part[";".join(sorted(x for x in m[k] if x))] += 1
    return dis, part


if __name__ == "__main__":
    for did in sys.argv[1:]:
        ss = fetch_all(did)
        dis, part = summarise(ss)
        print("=" * 72)
        print("Gemma %s : %d samples" % (did, len(ss)))
        print("  disease:", dict(dis.most_common(10)))
        print("  tissue :", dict(part.most_common(10)))
