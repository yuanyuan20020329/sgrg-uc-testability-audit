# -*- coding: utf-8 -*-
"""核对 `gdsType` 输入未生效的影响：按正确键名 `gdstype` 重取并重算。

背景
----
`partD_resfeas_geo.py` 用 `d.get("gdsType")` 读 esummary，而接口返回的键名是
`gdstype`，故 `assay_class()` 的第二个输入在归档流程中**始终为空**。这是实现缺陷，
但它是否改变了任何已发布取值，须用实际取回的 `gdstype` 重算才能回答。

做法
----
1. 对归档 raw 表的 521 个 uid 重新请求 esummary（`db=gds`），只保留 `gdstype`；
2. 用**同一个** `assay_class()`（`!Series_type` ＋ `gdstype` ＋ 标题＋摘要前 400 字符）
   对 521 条重算 `molecule` / `resolution`；
3. 与发布表 423 行的取值逐行比对，把 521 行写成一张对照表。

输出列
------
`gse, gdstype, molecule_published, resolution_published, molecule_recomputed,
 resolution_recomputed, changed, in_published_table`

- `changed`：仅对发布表内的 423 行有意义；98 行未进入发布表，`molecule_published` 等留空，
  `changed` 亦留空。
- 需要可访问 NCBI E-utilities；esummary 结果缓存到 `--cache`（默认 `out/_gdstype_cache.json`）。
"""
import argparse
import csv
import json
import pathlib
import sys
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from partD_resfeas_geo import assay_class  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def read_tsv(p):
    with open(p, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="gdstype_recompute_20261017.tsv")
    ap.add_argument("--cache", default="out/_gdstype_cache.json")
    ap.add_argument("--batch", type=int, default=100)
    args = ap.parse_args()

    raw = read_tsv(ROOT / "results" / "partD_resfeas_geo_series_raw.tsv")
    pub = {r["gse"]: r for r in read_tsv(ROOT / "resource_paper" / "deposit" / "partD_resfeas_geo_series.tsv")}
    uids = [r["uid"] for r in raw if r.get("uid")]

    cpath = ROOT / args.cache
    cache = json.load(open(cpath, encoding="utf-8")) if cpath.exists() else {}
    todo = [u for u in uids if u not in cache]
    print("uids=%d cached=%d to_fetch=%d" % (len(uids), len(cache), len(todo)))
    for i in range(0, len(todo), args.batch):
        chunk = todo[i:i + args.batch]
        url = "%s/esummary.fcgi?db=gds&retmode=json&id=%s" % (EUTILS, ",".join(chunk))
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as fh:
                    js = json.load(fh)
                break
            except Exception as e:  # noqa: BLE001
                print("retry", chunk[0], e)
                time.sleep(3)
        else:
            raise SystemExit("esummary fetch failed")
        res = js.get("result", {})
        for u in chunk:
            cache[u] = res.get(u, {}).get("gdstype", "")
        cpath.parent.mkdir(parents=True, exist_ok=True)
        cpath.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")
        print("fetched", len(cache))

    rows, n_flip = [], 0
    for r in raw:
        gd = cache.get(r.get("uid", ""), "")
        mol, res = assay_class(r.get("series_type", ""), gd,
                               (r.get("title") or "") + " " + (r.get("summary") or ""))
        p = pub.get(r["gse"])
        changed = ""
        if p:
            changed = str(mol != p["molecule"] or res != p["resolution"])
            n_flip += int(changed == "True")
        rows.append({"gse": r["gse"], "gdstype": gd,
                     "molecule_published": p["molecule"] if p else "",
                     "resolution_published": p["resolution"] if p else "",
                     "molecule_recomputed": mol, "resolution_recomputed": res,
                     "changed": changed, "in_published_table": bool(p)})

    print("rows=%d published=%d changed=%d" % (len(rows), len(pub), n_flip))
    with open(args.out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
