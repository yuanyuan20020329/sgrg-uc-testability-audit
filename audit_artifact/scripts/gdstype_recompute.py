# -*- coding: utf-8 -*-
"""核对 `gdsType` 输入未生效的影响：按正确键名 `gdstype` 重取，并用原 SOFT 的
`!Series_type` 重算，只改变一个输入。

背景
----
`partD_resfeas_geo.py`（第 320 行）用 `d.get("gdsType")` 读 esummary，而接口返回的
键名是 `gdstype`，故 `assay_class()` 的第二个输入在归档流程中**始终为空**。这是实现
缺陷，但它是否改变了任何已发布取值，须用实际取回的 `gdstype` 重算才能回答。

做法（三个输入逐一对齐归档时的实际取值，只把空的 `gdstype` 换成真值）
--------------------------------------------------------------------
1. `series_type`：调 `partD_resfeas_geo.parse_soft_series()` 读系列级 SOFT 缓存里的
   `!Series_type`（与归档流程同一个函数、同一份缓存）。
2. `gdstype`：对 raw 表的 521 个 uid 重新请求 esummary（`db=gds`）。
3. `blob`：直接用 raw 表的 `blob` 列（归档时写成 `title + " " + summary`，已逐字比对相等）。
4. 用**同一个** `assay_class()` 对 521 条重算 `molecule` / `resolution`，与发布表 423 行
   逐行比对，并把 521 行写成一张对照表。

★ 为什么必须走 SOFT 而不能从 raw 表取 `series_type`
    `results/partD_resfeas_geo_series_raw.tsv` **没有** `series_type` 列（它的第 6 列叫
    `gdsType`，且在归档中恒为空）。早期版本的本脚本写成 `r.get("series_type", "")`，
    于是 `series_type` 也被传成空串——**同时改变了两个输入**，不再是"只改变一个输入"的
    对照实验。v1.0.2 起改为从 SOFT 缓存取真值，并把取到的 `series_type` 写进输出表，
    以便读者自行核对第三个输入确实非空。

输出列
------
`gse, series_type, gdstype, molecule_published, resolution_published,
 molecule_recomputed, resolution_recomputed, changed, in_published_table`

- `changed`：仅对发布表内的 423 行有意义；98 行未进入发布表，`*_published` 留空，
  `changed` 亦留空。
- 需要可访问 NCBI E-utilities 以补取 `gdstype`；已缓存则不再请求，缓存默认
  `out/_gdstype_cache.json`。

用法
----
在**代码仓库工作树根**执行（需要 `results/`、`scripts/`、`data/raw/geo_soft_series/`）：

    python audit_artifact/scripts/gdstype_recompute.py --repo . --out gdstype_recompute_20261017.tsv
"""
import argparse
import csv
import json
import os
import pathlib
import re
import sys
import time
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def read_tsv(p):
    with open(p, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".",
                    help="代码仓库工作树根（含 results/、scripts/、data/raw/geo_soft_series/）")
    ap.add_argument("--out", default="gdstype_recompute_20261017.tsv")
    ap.add_argument("--cache", default="out/_gdstype_cache.json")
    ap.add_argument("--batch", type=int, default=100)
    args = ap.parse_args()

    repo = pathlib.Path(args.repo).resolve()
    # parse_soft_series() 内部用相对路径 data/raw/geo_soft_series，须先把 cwd 切到仓库根。
    os.chdir(repo)
    sys.path.insert(0, str(repo / "scripts"))
    from partD_resfeas_geo import assay_class, parse_soft_series  # noqa: E402

    raw = read_tsv(repo / "results" / "partD_resfeas_geo_series_raw.tsv")
    pub = {r["gse"]: r for r in read_tsv(
        repo / "resource_paper" / "deposit" / "partD_resfeas_geo_series.tsv")}
    uids = [r["uid"] for r in raw if r.get("uid")]

    cpath = repo / args.cache
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
    n_empty_st = n_empty_gd = 0
    for r in raw:
        gd = cache.get(r.get("uid", ""), "")
        # 第三个输入口径 1：系列级 SOFT 的 !Series_type（与归档流程同一函数、同一缓存）
        st = parse_soft_series(r["gse"]).get("series_type", "")
        blob = r.get("blob") or ((r.get("title") or "") + " " + (r.get("summary") or ""))
        n_empty_st += int(not st.strip())
        n_empty_gd += int(not gd.strip())
        mol, res_cls = assay_class(st, gd, blob)
        p = pub.get(r["gse"])
        changed = ""
        if p:
            changed = str(mol != p["molecule"] or res_cls != p["resolution"])
            n_flip += int(changed == "True")
        rows.append({"gse": r["gse"], "series_type": st, "gdstype": gd,
                     "molecule_published": p["molecule"] if p else "",
                     "resolution_published": p["resolution"] if p else "",
                     "molecule_recomputed": mol, "resolution_recomputed": res_cls,
                     "changed": changed, "in_published_table": bool(p)})

    print("rows=%d published=%d changed=%d" % (len(rows), len(pub), n_flip))
    print("empty series_type=%d empty gdstype=%d" % (n_empty_st, n_empty_gd))
    # 逐行核对两个输入是否承载同一信息（仅分隔符 `| ` 与 `; ` 不同）。
    # 521/521 相等即解释了 `changed = 0`：不是因为 `gdstype` 无关紧要，
    # 而是 `series_type` 已经带着完全相同的内容。
    def _norm(s):
        return re.sub(r"\s*[|;]\s*", " ; ", s.strip()).lower()
    n_same = sum(1 for r in rows if _norm(r["series_type"]) == _norm(r["gdstype"]))
    print("series_type == gdstype (sep-normalised): %d/%d" % (n_same, len(rows)))
    with open(args.out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
