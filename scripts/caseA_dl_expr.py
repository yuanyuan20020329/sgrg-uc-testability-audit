# -*- coding: utf-8 -*-
"""
下载 Gemma 表达矩阵。
正确契约（从 openapi.json 提取）：
  /datasets/{id}/data/processed   Accept: text/tab-separated-values; charset=UTF-8  -> TSV
  /datasets/{id}/data/raw         Accept: text/tab-separated-values; charset=UTF-8  -> TSV
  /datasets/{id}/expressions/genes/{genes}  Accept: application/json -> JSON（可指定基因，快）
"""
import urllib.request, gzip, sys, os, time

ACC_TSV = "text/tab-separated-values; charset=UTF-8"

def fetch(url, accept, timeout=1800):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": accept})
    raw = urllib.request.urlopen(req, timeout=timeout).read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw


def main():
    did = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "processed"
    outdir = r"C:\Users\H-P\WorkBuddy\应激颗粒\data\raw\caseA\gemma"
    os.makedirs(outdir, exist_ok=True)
    if mode == "processed":
        out = os.path.join(outdir, "dataset_%s_processed.tsv" % did)
        url = "https://gemma.msl.ubc.ca/rest/v2/datasets/%s/data/processed" % did
        acc = ACC_TSV
    elif mode == "raw":
        out = os.path.join(outdir, "dataset_%s_raw.tsv" % did)
        url = "https://gemma.msl.ubc.ca/rest/v2/datasets/%s/data/raw" % did
        acc = ACC_TSV
    else:
        out = os.path.join(outdir, "dataset_%s_genes_%s.json" % (did, mode.replace(",", "_")))
        url = "https://gemma.msl.ubc.ca/rest/v2/datasets/%s/expressions/genes/%s" % (did, mode)
        acc = "application/json"
    t0 = time.time()
    try:
        raw = fetch(url, acc)
    except Exception as e:
        print("ERR", repr(e)[:200], flush=True)
        sys.exit(1)
    with open(out, "wb") as f:
        f.write(raw)
    print("OK %d bytes in %.1fs -> %s" % (len(raw), time.time() - t0, out), flush=True)
    print("head:", raw[:350].decode("utf-8", "replace").replace("\n", " | "), flush=True)


if __name__ == "__main__":
    main()
