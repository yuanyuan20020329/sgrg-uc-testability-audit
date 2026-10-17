# -*- coding: utf-8 -*-
"""第十轮审稿应答用：GEO 检索截断的对比分析与可归档工件生成。

目的（对应审稿人 R2-1、R2-2、综合意见「必须修改 3」）
---------------------------------------------------------
1) 用与归档完全相同的检索式（直接 import scripts/partD_resfeas_geo.py 的 QUERIES），
   在同一时点分别以 retmax = 600 与 retmax = 5000 取回并集，因而「截断丢弃集」
   是**同一快照内**两次请求的纯差集，不含源库增长项。
2) 对被丢弃系列与已纳入系列做逐指标对比（年份、样本数、UC/IBD/CD 标签、
   组织部位、分子模态、是否单细胞/空间）。
3) 落盘可归档的机器可读产物：三组 UID 清单、双向差集、逐式命中数、截断系列清单。

输出目录：out/geo_retrieval_audit/
"""
from __future__ import annotations

import json
import re
import statistics
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import partD_resfeas_geo as G  # noqa: E402  复用同一份 QUERIES / esearch / esummary

OUT = ROOT / "out" / "geo_retrieval_audit"
OUT.mkdir(parents=True, exist_ok=True)
TODAY = time.strftime("%Y-%m-%d")
T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)


def union_uids(retmax: int) -> tuple[set[str], dict[str, int]]:
    """四条检索式的 UID 并集去重；返回 (uid 集合, 逐式 count)"""
    uids: set[str] = set()
    counts: dict[str, int] = {}
    for k, q in G.QUERIES.items():
        ids, cnt = G.esearch(q, retmax=retmax)
        counts[k] = cnt
        uids.update(ids)
        log("  retmax=%-5d %-8s 命中 %5d 取回 %5d" % (retmax, k, cnt, len(ids)))
    return uids, counts


def human_gse_table(uids: set[str]) -> pd.DataFrame:
    """esummary -> 仅保留 entrytype=GSE 且含 Homo sapiens 的系列"""
    summ = G.esummary(sorted(uids))
    recs = []
    for uid, d in summ.items():
        if uid == "uids":
            continue
        if (d.get("entrytype") or "") != "GSE":
            continue
        if "Homo sapiens" not in (d.get("taxon") or ""):
            continue
        acc = d.get("accession")
        if not acc:
            continue
        recs.append({
            "uid": uid,
            "gse": acc,
            "title": d.get("title") or "",
            "summary": d.get("summary") or "",
            "n_samples": int(d.get("n_samples") or 0),
            "gdsType": d.get("gdsType") or "",
            "pdat": d.get("pdat") or "",
            "pubmed": ";".join(str(x) for x in (d.get("pubmedids") or [])),
        })
    R = pd.DataFrame(recs).drop_duplicates("gse")
    R["blob"] = (R.title.fillna("") + " " + R.summary.fillna(""))
    return R


def annotate(R: pd.DataFrame) -> pd.DataFrame:
    R = R.copy()
    R["is_sc_or_spatial"] = R.blob.apply(lambda x: bool(G.SC_HINT.search(x)))
    R["lbl_UC"] = R.blob.apply(lambda x: bool(G.RE_UC.search(x)))
    R["lbl_CD"] = R.blob.apply(lambda x: bool(G.RE_CD.search(x)))
    R["lbl_IBD"] = R.blob.apply(lambda x: bool(G.RE_IBD.search(x)))
    R["site"] = R.blob.apply(G.site_of)
    R["in_vitro_hint"] = R.blob.apply(lambda x: bool(G.RE_INVITRO.search(x)))
    R["molecule"] = [
        G.assay_class("", gt, bl)[0] for gt, bl in zip(R.gdsType, R.blob)
    ]
    R["year"] = R.pdat.astype(str).str.slice(0, 4).replace("", None)
    return R


def summarize(R: pd.DataFrame) -> dict:
    n = len(R)
    if n == 0:
        return {"n_series": 0}

    def med(col):
        v = pd.to_numeric(R[col], errors="coerce").dropna()
        return None if v.empty else {
            "median": float(v.median()),
            "q1": float(v.quantile(0.25)),
            "q3": float(v.quantile(0.75)),
            "min": float(v.min()),
            "max": float(v.max()),
        }

    yr = pd.to_numeric(R["year"], errors="coerce").dropna()
    latest = yr.max() if not yr.empty else None
    return {
        "n_series": n,
        "n_samples": med("n_samples"),
        "n_samples_gt400": int((pd.to_numeric(R.n_samples, errors="coerce") > 400).sum()),
        "year": {
            "median": float(yr.median()) if not yr.empty else None,
            "min": float(yr.min()) if not yr.empty else None,
            "max": float(yr.max()) if not yr.empty else None,
            "n_ge_recent5": int((yr >= (latest - 4)).sum()) if latest else 0,
            "top_year": (yr.value_counts().head(5).to_dict() if not yr.empty else {}),
        },
        "label_UC": int(R.lbl_UC.sum()),
        "label_CD": int(R.lbl_CD.sum()),
        "label_IBD": int(R.lbl_IBD.sum()),
        "label_UC_or_CD_or_IBD": int((R.lbl_UC | R.lbl_CD | R.lbl_IBD).sum()),
        "site_colon": int((R.site == "colon").sum()),
        "site_ileum": int((R.site == "ileum").sum()),
        "site_rectum": int((R.site == "rectum").sum()),
        "site_none": int((R.site == "").sum()),
        "molecule": R.molecule.value_counts().to_dict(),
        "is_sc_or_spatial": int(R.is_sc_or_spatial.sum()),
        "in_vitro_hint": int(R.in_vitro_hint.sum()),
    }


def main() -> None:
    log("== 1) 双 retmax 检索（同一时点）==")
    uids600, cnt600 = union_uids(600)
    uidsfull, cntfull = union_uids(5000)
    log("并集 UID：cap600 = %d，full = %d" % (len(uids600), len(uidsfull)))

    log("== 2) esummary ==")
    R600 = human_gse_table(uids600)
    Rfull = human_gse_table(uidsfull)
    log("人类 GSE：cap600 = %d，full = %d" % (len(R600), len(Rfull)))

    R600 = annotate(R600)
    Rfull = annotate(Rfull)

    trunc = Rfull[~Rfull.gse.isin(set(R600.gse))].copy()
    extra_in_600 = R600[~R600.gse.isin(set(Rfull.gse))].copy()
    log("被截断丢弃 = %d；仅在 600 口径出现 = %d" % (len(trunc), len(extra_in_600)))

    log("== 3) 归档运行自身产出的原始系列表（2026-10-09，521 行）==")
    raw = pd.read_csv(ROOT / "results" / "partD_resfeas_geo_series_raw.tsv", sep="\t")
    arch_gse = set(raw.gse.astype(str))
    log("归档原始表系列 = %d" % len(arch_gse))

    # ---- 双向差集（归档 521 vs 本次 cap600；本次 cap600 vs 本次 full）
    diff_arch_minus_600 = sorted(arch_gse - set(R600.gse))
    diff_600_minus_arch = sorted(set(R600.gse) - arch_gse)
    diff_600_minus_full = sorted(set(R600.gse) - set(Rfull.gse))
    diff_full_minus_600 = sorted(set(Rfull.gse) - set(R600.gse))
    log("归档521 − 本次cap600 = %d；本次cap600 − 归档521 = %d"
        % (len(diff_arch_minus_600), len(diff_600_minus_arch)))

    log("== 4) 对比汇总 ==")
    summary = {
        "run_date": TODAY,
        "note": "同一时点两次请求的差集即纯截断效应；归档 521 行来自 2026-10-09 运行自身产出",
        "query_counts": {"cap600": cnt600, "full": cntfull},
        "union_uid": {"cap600": len(uids600), "full": len(uidsfull)},
        "human_gse": {"cap600": len(R600), "full": len(Rfull)},
        "archived_run_20261009": {"n_series_raw_table": len(arch_gse)},
        "diff": {
            "archived_minus_cap600": diff_arch_minus_600,
            "cap600_minus_archived": diff_600_minus_arch,
            "cap600_minus_full": diff_600_minus_full,
            "full_minus_cap600": diff_full_minus_600,
        },
        "included_cap600": summarize(R600),
        "truncated_full_minus_cap600": summarize(trunc),
        "pairing_note": "候选配对计数需要样本级 SOFT，本次对比不评估；给出 in_vitro 与单细胞/空间线索代理",
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---- 落盘工件
    (OUT / "uids_cap600_union.txt").write_text(
        "\n".join(sorted(uids600)) + "\n", encoding="utf-8")
    (OUT / "uids_full_union.txt").write_text(
        "\n".join(sorted(uidsfull)) + "\n", encoding="utf-8")
    (OUT / "gse_cap600_human.tsv").write_text(
        R600[["uid", "gse", "title", "n_samples", "gdsType", "pdat",
              "molecule", "is_sc_or_spatial", "lbl_UC", "lbl_CD", "lbl_IBD", "site"]]
        .to_csv(sep="\t", index=False), encoding="utf-8")
    (OUT / "gse_full_human.tsv").write_text(
        Rfull[["uid", "gse", "title", "n_samples", "gdsType", "pdat",
               "molecule", "is_sc_or_spatial", "lbl_UC", "lbl_CD", "lbl_IBD", "site"]]
        .to_csv(sep="\t", index=False), encoding="utf-8")
    (OUT / "gse_truncated_by_retmax600.tsv").write_text(
        trunc[["uid", "gse", "title", "n_samples", "gdsType", "pdat",
               "molecule", "is_sc_or_spatial", "lbl_UC", "lbl_CD", "lbl_IBD",
               "site", "in_vitro_hint"]].to_csv(sep="\t", index=False), encoding="utf-8")
    (OUT / "diff_tables").mkdir(exist_ok=True)
    for name, lst in [("archived_minus_cap600", diff_arch_minus_600),
                      ("cap600_minus_archived", diff_600_minus_arch),
                      ("cap600_minus_full", diff_600_minus_full),
                      ("full_minus_cap600", diff_full_minus_600)]:
        (OUT / "diff_tables" / f"{name}.txt").write_text(
            "\n".join(lst) + ("\n" if lst else ""), encoding="utf-8")

    print("\n===== 对比表（本次 %s 复核）=====" % TODAY)
    inc, tr = summary["included_cap600"], summary["truncated_full_minus_cap600"]
    rows = [
        ("系列数", inc["n_series"], tr["n_series"]),
        ("样本数中位数", inc["n_samples"]["median"], tr["n_samples"]["median"]),
        ("样本数 IQR", f"{inc['n_samples']['q1']}–{inc['n_samples']['q3']}",
         f"{tr['n_samples']['q1']}–{tr['n_samples']['q3']}"),
        ("样本数 >400 的系列", inc["n_samples_gt400"], tr["n_samples_gt400"]),
        ("年份中位数", inc["year"]["median"], tr["year"]["median"]),
        ("年份范围", f"{inc['year']['min']}–{inc['year']['max']}",
         f"{tr['year']['min']}–{tr['year']['max']}"),
        ("近 5 年（≥max−4）系列", inc["year"]["n_ge_recent5"], tr["year"]["n_ge_recent5"]),
        ("标题/摘要含 UC", inc["label_UC"], tr["label_UC"]),
        ("标题/摘要含 CD", inc["label_CD"], tr["label_CD"]),
        ("标题/摘要含 IBD", inc["label_IBD"], tr["label_IBD"]),
        ("含任一 IBD 标签", inc["label_UC_or_CD_or_IBD"], tr["label_UC_or_CD_or_IBD"]),
        ("组织=结肠", inc["site_colon"], tr["site_colon"]),
        ("组织=回肠", inc["site_ileum"], tr["site_ileum"]),
        ("组织=直肠", inc["site_rectum"], tr["site_rectum"]),
        ("组织未识别", inc["site_none"], tr["site_none"]),
        ("模态 mRNA", inc["molecule"].get("mRNA", 0), tr["molecule"].get("mRNA", 0)),
        ("模态 chromatin", inc["molecule"].get("chromatin", 0), tr["molecule"].get("chromatin", 0)),
        ("模态 ncRNA", inc["molecule"].get("ncRNA", 0), tr["molecule"].get("ncRNA", 0)),
        ("模态 other", inc["molecule"].get("other", 0), tr["molecule"].get("other", 0)),
        ("像单细胞/空间", inc["is_sc_or_spatial"], tr["is_sc_or_spatial"]),
        ("含离体模型线索", inc["in_vitro_hint"], tr["in_vitro_hint"]),
    ]
    w = max(len(r[0]) for r in rows) + 2
    print(f"{'指标':<{w}}{'已纳入(cap600)':>16}{'被截断':>12}")
    for a, b, c in rows:
        print(f"{a:<{w}}{str(b):>16}{str(c):>12}")
    log("工件已写入 out/geo_retrieval_audit/")


if __name__ == "__main__":
    main()
