# -*- coding: utf-8 -*-
"""
Part B · Step 2 — GSE19738 rs 540049（control vs MDD，单对照）全基因组翻页。
提取：探针→基因符号、行级 p、corrected q、contrast lfc。
产出：
  results/partB_gse19738_genomewide.tsv   （全部探针行，紧凑版）
用于：基因全集（universe）、MDD-DEG 清单（q<0.05）、SGRG 面板富集检验。
"""
import os, sys, json, time, gzip, io, urllib.request, urllib.error

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT  = os.path.join(BASE, "results")
CACHE = os.path.join(BASE, "data", "partB_rs540049_cache")
os.makedirs(CACHE, exist_ok=True)

RS = 540049
PAGE = 100
TOTAL = None  # 从第一页取

def log(*a): print(*a, flush=True)

def http_get_json(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept-Encoding": "gzip", "User-Agent": "curl/8.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                raw = r.read()
                if r.headers.get("Content-Encoding") == "gzip":
                    raw = gzip.GzipFile(fileobj=io.BytesIO(raw)).read()
                return json.loads(raw)
        except Exception as e:
            last = e
            time.sleep(2 + 3 * i)
    raise RuntimeError(f"GET failed: {url} :: {last}")

def fetch_page(offset):
    cache = os.path.join(CACHE, f"page_{offset:06d}.json")
    if os.path.exists(cache):
        with open(cache, encoding="utf-8") as f:
            return json.load(f)
    url = f"https://gemma.msl.ubc.ca/rest/v2/resultSets/{RS}?limit={PAGE}&offset={offset}"
    d = http_get_json(url)
    rows = (d.get("data") or {}).get("results", [])
    total = d.get("totalElements")
    with open(cache, "w", encoding="utf-8") as f:
        json.dump({"rows": rows, "total": total}, f, ensure_ascii=False)
    return {"rows": rows, "total": total}

def main():
    global TOTAL
    first = fetch_page(0)
    TOTAL = first["total"]
    log(f"total probe rows: {TOTAL}")
    offsets = list(range(0, TOTAL, PAGE))
    out_path = os.path.join(OUT, "partB_gse19738_genomewide.tsv")
    cols = ["probe", "symbols", "p", "q", "lfc", "tstat"]
    n_written = 0
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for off in offsets:
            page = fetch_page(off)
            rows = page["rows"]
            if not rows:
                log(f"  offset {off}: EMPTY (stopping)")
                break
            for r in rows:
                syms = "|".join(g.get("officialSymbol") or "" for g in r.get("genes", []) if g.get("officialSymbol"))
                c = (r.get("contrasts") or [{}])[0]
                f.write("\t".join([
                    r.get("probeName") or "NA",
                    syms or "NA",
                    repr(r.get("pValue")) if r.get("pValue") is not None else "NA",
                    repr(r.get("correctedPvalue")) if r.get("correctedPvalue") is not None else "NA",
                    repr(c.get("logFoldChange")) if c.get("logFoldChange") is not None else "NA",
                    repr(c.get("tStat")) if c.get("tStat") is not None else "NA",
                ]) + "\n")
                n_written += 1
            if (off // PAGE) % 20 == 0:
                log(f"  {off}/{TOTAL} ({n_written} rows)")
    log(f"done: {n_written} rows -> {out_path}")

if __name__ == "__main__":
    main()
