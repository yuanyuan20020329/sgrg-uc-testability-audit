#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
案例 A 第 1 步：获取铁死亡基因集（多来源并列，供口径对照）
==========================================================
原则（计划书 §6.1 第 2 条）：**使用原论文公布的基因集，不根据复现结果删改基因**。
因此本脚本把多个稳定来源并列下载、单独保存，不合并、不取交集——
待选定目标论文后再决定用哪一版作为主口径，其余作为敏感性分析。

来源：
  MSigDB C2:CP（含 KEGG/Reactome 等通路，其中 FERROPTOSIS 来自 KEGG hsa04216）
  MSigDB H（hallmark，作为该平台的参照全集，用于随机集抽样的表达分层）

输出：
  data/raw/caseA/genesets/*.gmt
  results/caseA_genesets.tsv   每个基因集一行：名称 / 来源 / 基因数 / 示例基因
"""
import os, sys, re, json
import urllib.request, gzip

sys.stdout.reconfigure(encoding="utf-8")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "caseA", "genesets")
OUT = os.path.join(ROOT, "results")
os.makedirs(RAW, exist_ok=True)

BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
SOURCES = [
    ("2023.2.Hs", "h.all.v2023.2.Hs.symbols.gmt"),
    ("2023.1.Hs", "c2.cp.v2023.1.Hs.symbols.gmt"),
]


def fetch(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "Mozilla/5.0", "Accept-Encoding": "gzip"})
            with urllib.request.urlopen(req, timeout=120) as r:
                raw = r.read()
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            return raw.decode("utf-8", "ignore")
        except Exception as e:
            print("  [retry %d] %s" % (i, repr(e)[:90]))
    return None


def load_gmt(text):
    sets = {}
    for line in (text or "").splitlines():
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        sets[p[0]] = [g for g in p[2:] if g.strip()]
    return sets


def main():
    fetched = {}
    for ver, fn in SOURCES:
        url = "%s/%s/%s" % (BASE, ver, fn)
        print("下载 %s ..." % url)
        t = fetch(url)
        if t is None:
            print("  FAIL")
            continue
        p = os.path.join(RAW, fn)
        with open(p, "w", encoding="utf-8", newline="") as f:
            f.write(t)
        fetched[fn] = load_gmt(t)
        print("  OK  %d 个基因集 -> %s" % (len(fetched[fn]), p))

    # ---- 定位铁死亡相关基因集
    rows = []
    for fn, sets in fetched.items():
        for name, genes in sets.items():
            if re.search(r"FERROPTOSIS", name, re.I):
                rows.append({"gmt": fn, "set_name": name, "n_genes": len(genes),
                             "example": ";".join(genes[:8])})
    print("\n★ 铁死亡相关基因集:")
    for r in rows:
        print("   %-46s n=%-4d  %s" % (r["set_name"], r["n_genes"], r["example"]))
    if not rows:
        print("   ⚠️ 未找到——需改用原论文附表")

    # ---- 焦亡 / 自噬（为案例 B 预取，顺手）
    for kw in ["PYROPTOSIS", "AUTOPHAGY", "FERRITIN"]:
        hit = [(fn, n, len(g)) for fn, sets in fetched.items()
               for n, g in sets.items() if re.search(kw, n, re.I)]
        print("   [%s] %s" % (kw, hit[:6] if hit else "无"))

    import pandas as pd
    df = pd.DataFrame(rows)
    p = os.path.join(OUT, "caseA_genesets.tsv")
    df.to_csv(p, sep="\t", index=False)
    print("\n写出 %s (%d 行)" % (p, len(df)))


if __name__ == "__main__":
    main()
