# -*- coding: utf-8 -*-
"""用**发布表同一** `assay_class()` 规则复核 98 个未进入 GEO 表的系列。

背景
----
`partD_resfeas_geo.py` 对 521 个归档人类系列做筛选，其中 423 个进入发布表，
其余 98 个因样本量或分辨率规则跳过。对早期一次复核（`out/verify_bulk_claim.py`）
只把 `!Series_type` 与标题拼接后判定，**不使用** `gdsType` 与摘要，与发布表
所用的规则不同源。本脚本改用发布表的 `assay_class()`（三项输入：
`!Series_type`、`gdsType`、标题＋摘要前 400 字符）重新判定这 98 个系列，
并同时给出旧规则的取值与两者的差异标记，供逐条核对。

输入
----
- `results/partD_resfeas_geo_series_raw.tsv`：521 条归档记录（含 `gdsType` 列）
- `data/raw/geo_soft_series/<GSE>.txt`：系列级 SOFT 缓存（**不属本包**，见 README §3）
- `resource_paper/deposit/partD_resfeas_geo_series.tsv`：已发布 423 行（用于取补集）
- `scripts/partD_resfeas_geo.py`：解析与分类函数（属代码归档，不属本包）

输出
----
- 默认写出 `skipped98_modality_recheck.tsv`（98 行 ＋ 表头）

注意
----
归档的 raw 表中 `gdsType` 一列**全部为空**：生成脚本读的是 esummary 的
`gdsType` 键，而该接口实际返回的键名是 `gdstype`，故 `assay_class()` 的
第二个输入在归档流程中从未携带信息。本脚本因此用 `gdstype` 为空值复现
发布表当时的实际输入；其影响另行评估（见 README §3 第 6 条）。
"""
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from partD_resfeas_geo import assay_class  # noqa: E402

SER = ROOT / "data" / "raw" / "geo_soft_series"
RAW = ROOT / "results" / "partD_resfeas_geo_series_raw.tsv"
PARSED = ROOT / "resource_paper" / "deposit" / "partD_resfeas_geo_series.tsv"

OLD_MOL = re.compile(r"non[\s\-]?coding rna|microrna|mirna|ncrna|lncrna", re.I)
GENE = re.compile(r"genome binding|occupancy|methylation|atac|chip[\s\-]?seq|chromatin|histone", re.I)
EXPR = re.compile(r"expression profiling|transcriptom|rna[\s\-]?seq|microarray|third[\s\-]?party reanalysis", re.I)
SC_OLD = re.compile(r"spatial transcriptom|xenium|cosmx|visium|merscope|in situ|slide[\s\-]?seq|"
                    r"stereo[\s\-]?seq|seqfish|merfish|single[\s\-]?cell|single[\s\-]?nucle|"
                    r"\bsc[\s\-]?rna|\bsn[\s\-]?rna|10x", re.I)


def read_tsv(p):
    with open(p, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def old_rule(series_type, title):
    """早期复核所用规则：只有 !Series_type ＋ 标题，且分辨率不区分阵列/测序。"""
    st = (str(series_type or "") + " " + str(title or "")).lower()
    if OLD_MOL.search(st):
        return "ncRNA", "ncRNA"
    if GENE.search(st):
        return "chromatin", "chromatin"
    if not EXPR.search(st):
        return "other", "other"
    return "mRNA", ("single-cell/spatial" if SC_OLD.search(str(series_type or "") + " " + str(title or "")) else "bulk")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="skipped98_modality_recheck.tsv")
    args = ap.parse_args()

    raw = read_tsv(RAW)
    parsed = {r["gse"] for r in read_tsv(PARSED)}
    skipped = [r for r in raw if r["gse"] not in parsed]
    print("raw=%d parsed=%d skipped=%d" % (len(raw), len(parsed), len(skipped)))

    rows, missing = [], []
    for r in skipped:
        fp = SER / (r["gse"] + ".txt")
        if fp.exists():
            st = " | ".join(m.strip() for m in re.findall(
                r"^!Series_type = (.*)$", fp.read_text(encoding="utf-8", errors="replace"), re.M))
        else:
            missing.append(r["gse"])
            st = ""
        mol, res = assay_class(st, r["gdsType"], (r.get("title") or "") + " " + (r.get("summary") or ""))
        om, ores = old_rule(st, r["title"])
        rows.append({"gse": r["gse"], "n_samples": r["n_samples"], "title": r["title"],
                     "gdsType": r["gdsType"], "series_type": st,
                     "molecule_true_rule": mol, "resolution_true_rule": res,
                     "molecule_old_rule": om, "resolution_old_rule": ores,
                     "disagree": (mol != om) or (res != ores)})

    print("missing SOFT cache:", len(missing), missing[:10])
    print("molecule:", Counter(r["molecule_true_rule"] for r in rows))
    print("resolution:", Counter(r["resolution_true_rule"] for r in rows))
    print("rows where the two rules differ:", sum(1 for r in rows if r["disagree"]))

    with open(args.out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
