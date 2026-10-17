# -*- coding: utf-8 -*-
"""检索式敏感性分析（对应审稿人 R2-3）。

做法
----
基线 = 归档四条检索式在 `retmax = 5000` 下的 UID 并集
      （读 out/geo_retrieval_audit/uids_full_union.txt，不重复请求）。
对每一个「候选补充词」构造变体检索式，取回其 UID 集合，报告：

  * 变体命中数；
  * **增量** = 变体 UID − 基线 UID，即「现有四条检索式未覆盖、该词能补到的系列」；
  * 重叠 = 变体 ∩ 基线。

对增量规模较小的**疾病/组织特异性术语**再做一步分诊：取 esummary，
统计其中「人类 GSE」条数，以及其中标题/摘要含 IBD 相关词
（`ulcerative colitis` / `crohn` / `IBD` / `UC` / `CD`）的条数——
后者才是「可能相关但被漏检」的候选。泛词（如单独 `UC`、`colorectal`）
增量过大且以无关命中为主，只报计数并显式标注未逐条分诊。

说明：GEO 的 esearch 在 db=gds 上匹配标题/摘要/描述等系列级字段，
**不检索 GSM 级样本元数据**；因此样本级表达（如 “non-inflamed mucosa”）
无法由系列级检索式覆盖，这一点在正文中单独声明。

输出：out/query_sensitivity/query_sensitivity.json 与 _increment_detail.tsv
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

OUT = ROOT / "out" / "query_sensitivity"
OUT.mkdir(parents=True, exist_ok=True)
BASE_UIDS = ROOT / "out" / "geo_retrieval_audit" / "uids_full_union.txt"
TRIAGE_MAX = 500          # 增量不超过该值才做逐条分诊
T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)

GUT = G.GUT
IBD_IN_TITLE = G.RE_UC

# 分组：A 疾病词表变体 / B 肠道部位词表变体 / C 状态与对照表达（基线完全不含）
VARIANTS = [
    ("A", "缩写 UC（无全称）", "(UC OR \"U.C.\")"),
    ("A", "临床术语 proctitis", "(proctitis OR proctocolitis)"),
    ("A", "临床术语 pancolitis", "(pancolitis OR \"pan-colitis\")"),
    ("A", "全称截断式 ulcerative colit", "(\"ulcerative colit\" OR \"ulcerative-colitis\")"),
    ("A", "泛称 colitis", "(colitis)"),
    ("B", "部位 ileocolonic", "(ileocolonic OR \"ileo-colonic\" OR ileocaecal OR \"ileocecal\")"),
    ("B", "部位 colorectal", "(colorectal)"),
    ("B", "部位 caecum/cecum/sigmoid/升结肠/降结肠",
     "(caecum OR cecum OR sigmoid OR \"ascending colon\" OR \"descending colon\")"),
    ("C", "对照侧 healthy donor", "(\"healthy donor\" OR \"healthy donors\" OR \"healthy control\")"),
    ("C", "对照侧 non-inflamed mucosa", "(noninflamed OR \"non-inflamed\" OR uninflamed OR \"unaffected mucosa\")"),
    ("C", "状态 quiescent/remission", "(quiescent OR remission OR \"inactive disease\")"),
]


def triage(new_uids: list[str]) -> dict:
    """对增量 UID 取 esummary，统计人类 GSE 与其中含 IBD 相关词者。"""
    summ = G.esummary(new_uids)
    human, ibd, examples = 0, 0, []
    for uid, d in summ.items():
        if uid == "uids":
            continue
        if (d.get("entrytype") or "") != "GSE":
            continue
        if "Homo sapiens" not in (d.get("taxon") or ""):
            continue
        human += 1
        blob = "%s %s" % (d.get("title") or "", d.get("summary") or "")
        if IBD_IN_TITLE.search(blob) or G.RE_CD.search(blob) or G.RE_IBD.search(blob):
            ibd += 1
            if len(examples) < 5:
                examples.append({"gse": d.get("accession"), "title": d.get("title"),
                                 "n_samples": d.get("n_samples"), "pdat": d.get("pdat")})
    return {"n_human_gse": human, "n_human_gse_with_ibd_term": ibd, "examples": examples}


def main() -> None:
    base = set(BASE_UIDS.read_text(encoding="utf-8").split())
    log("基线 UID（四条检索式 retmax=5000 并集）= %d" % len(base))

    rows = []
    for grp, name, term in VARIANTS:
        q = "(%s AND %s)" % (term, GUT)
        ids, cnt = G.esearch(q, retmax=5000)
        s = set(ids)
        new = sorted(s - base)
        row = {"group": grp, "variant": name, "term": term, "n_hits": cnt,
               "n_retrieved": len(s), "n_overlap": len(s & base), "n_increment": len(new),
               "increment_uids": new, "triaged": False}
        if new and len(new) <= TRIAGE_MAX:
            log("  %-30s 增量 %3d → 取 esummary 分诊" % (name, len(new)))
            row.update(triage(new))
            row["triaged"] = True
        else:
            log("  %-30s 命中 %5d | 重叠 %5d | 增量 %4d（未分诊）"
                % (name, cnt, len(s & base), len(new)))
        rows.append(row)

    slim_rows, uid_rows = [], []
    for r in rows:
        slim_rows.append({
            "group": r["group"], "variant": r["variant"], "term": r["term"],
            "n_hits": r["n_hits"], "n_overlap": r["n_overlap"],
            "n_increment": r["n_increment"], "triaged": r["triaged"],
            "n_increment_human_gse": r.get("n_human_gse") if r["triaged"] else None,
            "n_increment_human_gse_with_ibd_term":
                r.get("n_human_gse_with_ibd_term") if r["triaged"] else None,
            "human_gse_with_ibd_term":
                [e["gse"] for e in r.get("examples", [])] if r["triaged"] else [],
        })
        if r["triaged"]:
            for u in r["increment_uids"]:
                uid_rows.append({"group": r["group"], "variant": r["variant"], "uid": u})

    slim = {
        "run_date": time.strftime("%Y-%m-%d"),
        "baseline": {"source": str(BASE_UIDS.relative_to(ROOT)).replace("\\", "/"),
                     "n_uid": len(base),
                     "note": "四条归档检索式在 retmax=5000 下的 UID 并集，构成敏感性的基线"},
        "triage_max_increment": TRIAGE_MAX,
        "note": ("增量 = 变体命中 − 基线命中。泛词（单独 UC／colitis／colorectal 等）增量数千条、"
                 "以无关命中为主，未逐条分诊，其增量不得读作「漏掉的 IBD 系列数」。"
                 "仅对增量 ≤ 阈值者取 esummary 分诊，统计其中人类 GSE 及标题/摘要含 IBD 相关词者。"),
        "variants": slim_rows,
    }
    (OUT / "query_sensitivity_summary.json").write_text(
        json.dumps(slim, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(uid_rows).to_csv(OUT / "query_sensitivity_triaged_uids.tsv",
                                  sep="\t", index=False)

    print("\n===== 检索式敏感性（基线 %d 个 UID）=====" % len(base))
    w = max(len(r["variant"]) for r in rows) + 2
    print(f"{'组':<3}{'候选补充词':<{w}}{'命中':>7}{'重叠':>7}{'增量':>7}{'增量中人类GSE':>13}{'其中含IBD词':>12}")
    for r in rows:
        print(f"{r['group']:<3}{r['variant']:<{w}}{r['n_hits']:>7}{r['n_overlap']:>7}"
              f"{r['n_increment']:>7}"
              f"{(r.get('n_human_gse') if r['triaged'] else '-'):>13}"
              f"{(r.get('n_human_gse_with_ibd_term') if r['triaged'] else '-'):>12}")
    log("工件已写入 out/query_sensitivity/")


if __name__ == "__main__":
    main()
