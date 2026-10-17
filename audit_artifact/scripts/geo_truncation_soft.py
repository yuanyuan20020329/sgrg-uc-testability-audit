# -*- coding: utf-8 -*-
"""对「被 retmax=600 截断丢弃」的系列，用与发布表**完全相同**的 SOFT 规则复核。

为什么要单独一步
----------------
`out/geo_retrieval_audit.py` 的 molecule 列只用 esummary 的 `gdsType` + 摘要文本，
而发布表的 `molecule`/`resolution` 用的是 `!Series_type` + 平台（SOFT）。两者不是同一条规则，
不能直接拿来与归档的 423 行表对比。本脚本对 60 条被截断系列拉系列级 SOFT，
用 `partD_resfeas_geo.assay_class` 重新判分子类别与分辨率，并抽样级 SOFT 复核候选配对计数。

输出：out/geo_retrieval_audit/truncated_soft_audit.tsv / .json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import partD_resfeas_geo as G  # noqa: E402

OUT = ROOT / "out" / "geo_retrieval_audit"
# ★ 不污染归档 SOFT 缓存：样本级 SOFT 落到本审计目录的 scratch；系列级只读既有缓存
#   （60 条被截断系列的系列级 SOFT 在 2026-10-15 全口径复核时已抓取，存在则只读）。
G.CACHE = str(OUT / "soft_cache_scratch")
(OUT / "soft_cache_scratch").mkdir(parents=True, exist_ok=True)
T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)


def pairing_of(gse: str, max_samples: int = 400) -> dict:
    """与 partD_resfeas_geo.main() 逐字一致的配对计数规则。

    ★ 口径必须与流水线完全相同：流水线用
        paired = {k: v for k, v in donors.items() if {"inflamed", "noninflamed"} <= v}
      即**必须同时出现 inflamed 与 noninflamed 两种状态**才计为一对；
      健康的「healthy」状态**不**能顶替 noninflamed。配对键是供体 token（抽不到时用占位键）。
      同部位口径 pair_site 以 (供体, 部位) 为键，同样要求两种状态齐备且部位不是 na。
    """
    try:
        samples = G.parse_soft(gse, max_samples)
    except Exception as e:  # noqa: BLE001
        return {"err": type(e).__name__}
    if samples is None:
        return {"err": "no_soft"}
    if isinstance(samples, dict) and "too_many" in samples:
        return {"err": "too_many", "n": samples["too_many"]}
    donors: dict[str, set] = {}
    pair_site: dict[tuple, set] = {}
    states, dsets, sites = [], set(), set()
    n_unknown = n_nodonor = n_nosite = 0
    for smp in samples:
        txt = " ".join([smp["title"], smp["source"], smp["chars"]])
        if G.RE_NONINFL.search(txt):
            st = "noninflamed"
        elif G.RE_INFLAM.search(txt):
            st = "inflamed"
        elif G.RE_HEALTHY.search(txt):
            st = "healthy"
        else:
            st = "unknown"
        dn = G.donor_token(smp["title"], smp["source"])
        dis = ("UC" if G.RE_UC.search(txt) else "") or ("CD" if G.RE_CD.search(txt) else "") \
              or ("IBD" if G.RE_IBD.search(txt) else "") or ("" if st == "healthy" else "")
        site = G.site_of(txt)
        if st == "unknown":
            n_unknown += 1
        if not dn:
            n_nodonor += 1
        if not site:
            n_nosite += 1
        if dis:
            dsets.add(dis)
        if site:
            sites.add(site)
        states.append(st)
        donors.setdefault(dn or ("<%s#%s>" % (st, smp["gsm"])), set()).add(st)
        pair_site.setdefault((dn or smp["gsm"], site or "na"), set()).add(st)

    wanted = {"inflamed", "noninflamed"}
    n_paired = sum(1 for v in donors.values() if wanted <= v)
    n_paired_site = sum(1 for v in pair_site.values() if wanted <= v)
    return {"n_samples": len(samples), "n_paired_donors": n_paired,
            "n_paired_same_site": n_paired_site,
            "n_inflamed": states.count("inflamed"),
            "n_noninflamed": states.count("noninflamed"),
            "n_healthy": states.count("healthy"),
            "diseases": ";".join(sorted(dsets)), "sites": ";".join(sorted(sites)),
            "n_unknown_state": n_unknown, "n_no_donor": n_nodonor, "n_no_site": n_nosite}


def main() -> None:
    src = OUT / "gse_truncated_by_retmax600.tsv"
    T = pd.read_csv(src, sep="\t")
    log("被截断系列 = %d" % len(T))

    rows = []
    for i, r in enumerate(T.itertuples(), 1):
        try:
            ser = G.parse_soft_series(r.gse)
        except Exception as e:  # noqa: BLE001
            ser = {"series_type": "ERR:%s" % type(e).__name__}
        mol, res = G.assay_class(ser.get("series_type", ""), r.gdsType, r.title + " " + (r.summary if hasattr(r, "summary") else ""))
        rec = {
            "gse": r.gse, "n_samples": r.n_samples, "pdat": r.pdat,
            "series_type": ser.get("series_type", ""),
            "molecule_soft": mol, "resolution_soft": res,
            "molecule_esummary_proxy": r.molecule,
        }
        # 与归档同规则：只对 n_samples<=60（本例 is_sc_or_spatial 全为 0）解析样本级 SOFT
        if 0 < r.n_samples <= 60:
            rec.update(pairing_of(r.gse))
        else:
            rec["pairing"] = "not_parsed (n_samples>60)"
        rows.append(rec)
        if i % 15 == 0:
            log("  SOFT %d/%d" % (i, len(T)))

    D = pd.DataFrame(rows)
    D.to_csv(OUT / "truncated_soft_audit.tsv", sep="\t", index=False)

    mol_counts = D.molecule_soft.value_counts().to_dict()
    res_counts = D.resolution_soft.value_counts().to_dict()
    sc_sp = int(D.resolution_soft.isin(["single-cell", "spatial"]).sum())
    parsed = D[D.get("n_paired_donors").notna()] if "n_paired_donors" in D else D.iloc[0:0]
    summary = {
        "n_truncated": len(D),
        "molecule_soft": mol_counts,
        "resolution_soft": res_counts,
        "n_singlecell_or_spatial": sc_sp,
        "n_parsed_for_pairing": int(len(parsed)),
        "n_with_any_pairing": int((parsed.get("n_paired_donors", pd.Series(dtype=float)) > 0).sum()) if len(parsed) else 0,
        "note": "分子类别按发布表同源的 !Series_type + 平台规则判定（不是 esummary 代理）",
    }
    (OUT / "truncated_soft_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print("\n===== 60 条被截断系列（SOFT 同源规则复核）=====")
    print("分子类别:", mol_counts)
    print("分辨率  :", res_counts)
    print("被归为单细胞/空间的系列数:", sc_sp)
    if len(parsed):
        print("样本级 SOFT 已解析系列数:", len(parsed))
        print("其中有任意配对的:", summary["n_with_any_pairing"])
    print("\n逐条（前 20）:")
    cols = ["gse", "n_samples", "pdat", "molecule_soft", "resolution_soft", "series_type"]
    print(D[cols].head(20).to_string(index=False, max_colwidth=48))
    log("完成")


if __name__ == "__main__":
    main()
