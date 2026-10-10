# -*- coding: utf-8 -*-
"""A2 探测：确认 Gemma 元数据与「后果」类端点的可用字段与规模。

1) /datasets?query=... 全量分页 → 队列总数
2) /datasets/{id}/analyses/differential → 存算 DE 分析清单
3) /datasets/{id}/expressions/differential → DE 概要（DE 基因数 / 效应量）
4) /datasets/{id}/qc-metrics, /geeq/public → 质控分项
"""
import json
import os
import ssl
import time
import urllib.parse
import urllib.request

BASE = "https://gemma.msl.ubc.ca/rest/v2"
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
CACHE = "data/raw/gemma_a2"
os.makedirs(CACHE, exist_ok=True)


def get(url, tag=None):
    if tag:
        p = os.path.join(CACHE, tag + ".json")
        if os.path.exists(p):
            return json.load(open(p, encoding="utf-8"))
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    for k in range(3):
        try:
            with urllib.request.urlopen(req, context=CTX, timeout=120) as r:
                d = json.loads(r.read().decode("utf-8"))
            if tag:
                json.dump(d, open(p, "w", encoding="utf-8"), ensure_ascii=False)
            return d
        except Exception as e:
            if k == 2:
                print("  ERR", type(e).__name__, e)
                return None
            time.sleep(2)


print("=" * 70)
print("[1] 队列搜索分页（query=ulcerative colitis）")
for q in ["ulcerative colitis", "ulcerative+colitis", "colitis", "inflammatory bowel disease"]:
    d = get(BASE + "/datasets?query=" + urllib.parse.quote(q) + "&limit=1")
    tot = (d or {}).get("totalElements")
    print("  %-28s totalElements=%s" % (q, tot))

print("=" * 70)
print("[2] 存算 DE 分析清单 /datasets/5873/analyses/differential")
d = get(BASE + "/datasets/5873/analyses/differential", "probe_an_5873")
if d:
    print("  keys:", list(d.keys()))
    for a in (d.get("data") or [])[:3]:
        print(json.dumps(a, ensure_ascii=False)[:900])

print("=" * 70)
print("[3] DE 概要 /datasets/5873/expressions/differential")
d = get(BASE + "/datasets/5873/expressions/differential", "probe_de_5873")
if d is not None:
    print("  type:", type(d).__name__, "keys:", list(d.keys()) if isinstance(d, dict) else len(d))
    print(json.dumps(d, ensure_ascii=False)[:1200])

print("=" * 70)
print("[4] qc-metrics / geeq / sample-correlation / svd")
for ep in ["qc-metrics", "geeq/public", "sample-correlation", "svd", "publications"]:
    d = get("%s/datasets/5873/%s" % (BASE, ep), "probe_%s_5873" % ep.replace("/", "_"))
    if d is None:
        continue
    s = json.dumps(d, ensure_ascii=False)
    print("  [%s] len=%d  %s" % (ep, len(s), s[:400]))
