#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
案例 A 第 0 步：队列枚举与批次标注地景（铁死亡 × 肝癌/结直肠癌）
================================================================
目的：为独立案例 A 建立可审计的候选队列池，并量化「批次-疾病混杂」在
      该疾病领域的基线分布。**零表达矩阵下载** —— 全部来自 Gemma 数据集元数据。

与 Part A2 的关系：本脚本输出结构与 partA2 的队列表同构，便于直接复用
      partA2_analyse.py 的 L2 指标（Cramér's V / batch-only AUC）。

输出：
  results/caseA_gemma_cohorts.tsv       每个队列一行，含 batchEffect 枚举与统计文本
  results/caseA_landscape.json          汇总数字（可写入论文）

用法：
  python scripts/caseA_cohort_enum.py [--query hepatocellular] [--limit 400]
"""
import os, sys, json, gzip, time, argparse, re
import urllib.request, urllib.parse

sys.stdout.reconfigure(encoding="utf-8")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results")
RAW = os.path.join(ROOT, "data", "raw", "caseA")
os.makedirs(RAW, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

GEMMA = "https://gemma.msl.ubc.ca/rest/v2"

# 与 Part A2 同一套粗四类口径，便于跨疾病比较
ANNOT_MAP = {
    "BATCH_EFFECT_FAILURE":              ("batch_artifact", "批次伪影"),
    "PROBLEMATIC_BATCH_INFO_FAILURE":    ("batch_artifact", "批次伪影"),
    "BATCH_EFFECT_UNDETERMINED_FAILURE": ("batch_artifact", "批次伪影"),
    "SINGLETON_BATCHES_FAILURE":         ("batch_artifact", "批次伪影"),
    "UNINFORMATIVE_HEADERS_FAILURE":     ("insufficient",   "信息不足"),
    "NO_BATCH_INFO":                     ("insufficient",   "信息不足"),
    "BATCH_CORRECTED_SUCCESS":           ("corrected",      "已校正"),
    "NO_BATCH_EFFECT_SUCCESS":           ("clean",          "无问题"),
    "SINGLE_BATCH_SUCCESS":              ("single_batch",   "单批次"),
}


def get_json(url, tries=4, timeout=60):
    """Gemma 返回 gzip；必须解压后再 json.loads（否则 UnicodeDecodeError）。"""
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "Mozilla/5.0", "Accept-Encoding": "gzip"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                raw = r.read()
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            return json.loads(raw.decode("utf-8"))
        except Exception as e:
            last = e
            time.sleep(1.5 + i * 2.5)
    print("  [FAIL] %s -> %r" % (url[:110], last))
    return None


def fetch_datasets(query, limit=400):
    """分页拉取 query 命中的全部数据集。列表条目自带 batchEffect 等策展字段。"""
    rows, offset, page = [], 0, 100
    while offset < limit:
        url = "%s/datasets?query=%s&limit=%d&offset=%d" % (
            GEMMA, urllib.parse.quote(query), page, offset)
        d = get_json(url)
        if d is None:
            break
        batch = d.get("data", d) if isinstance(d, dict) else d
        if not isinstance(batch, list) or not batch:
            break
        rows.extend(batch)
        total = d.get("totalElements") if isinstance(d, dict) else None
        offset += page
        if total is not None and offset >= total:
            break
    # 去重
    seen, uniq = set(), []
    for r in rows:
        k = r.get("shortName") or r.get("id")
        if k in seen:
            continue
        seen.add(k)
        uniq.append(r)
    return uniq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", default="hepatocellular carcinoma")
    ap.add_argument("--limit", type=int, default=400)
    a = ap.parse_args()

    print("=" * 74)
    print("案例 A · 队列枚举  query = %r" % a.query)
    print("=" * 74)

    rows = fetch_datasets(a.query, a.limit)
    print("去重后数据集数 = %d" % len(rows))

    recs = []
    for r in rows:
        tax = r.get("taxon") or {}
        if isinstance(tax, dict):
            taxon = tax.get("commonName") or tax.get("scientificName") or ""
        else:
            taxon = str(tax)
        be = r.get("batchEffect")
        if isinstance(be, dict):
            be = be.get("name") or be.get("value") or str(be)
        cls, cls_cn = ANNOT_MAP.get(be, ("other", "其他/无标注"))

        # ★ 修正 1：batchConfound 不可截断 —— 关键信息在句末（"Factor(s) confounded were: ..."）
        bc = (r.get("batchConfound") or "").strip()
        # ★ 修正 2：geeq 是嵌套字典，不是扁平键 "geeq.publicQualityScore"
        gq = r.get("geeq")
        geeq_score = None
        if isinstance(gq, dict):
            geeq_score = gq.get("publicQualityScore")
        else:
            geeq_score = gq

        # 从完整文本里抽混杂因素名单
        conf_factors = ""
        m = re.search(r"[Ff]actor\(?s?\)?\s*(?:that\s+were\s+)?confounded[^:]*:\s*(.+)$", bc)
        if m:
            conf_factors = m.group(1).strip()

        recs.append({
            "shortName": r.get("shortName"),
            "gemma_id": r.get("id"),
            "taxon": taxon,
            "platform": ";".join(p.get("shortName", "") for p in (r.get("platforms") or [])
                                 if isinstance(p, dict)),
            "n_assay": r.get("bioAssayCount"),
            "is_single_cell": r.get("isSingleCell"),
            "external_db": r.get("externalDatabase"),
            "batchEffect": be,
            "annot_class": cls,
            "annot_class_cn": cls_cn,
            "batchConfound_full": bc,          # 完整文本，不截断
            "confound_factors": conf_factors,  # 抽出的混杂因素名单
            "batchEffectStatistics": (r.get("batchEffectStatistics") or "")[:200],
            "geeq": geeq_score,
            "name": (r.get("name") or "")[:200],
        })

    # 只保留人类（案例 A 明确要求人类疾病队列）
    human = [x for x in recs if "human" in (x["taxon"] or "").lower()]
    print("其中人类 = %d" % len(human))

    # ---- 分布统计（人类，且有样本数）
    from collections import Counter
    dist = Counter(x["annot_class_cn"] for x in human)
    n = len(human)
    print("\n批次策展标注分布（人类 n=%d）:" % n)
    for k, v in dist.most_common():
        print("   %-12s %4d  (%.1f%%)" % (k, v, 100.0 * v / n))

    # 与 Part A2（UC）同口径的「有问题或无法判定」比例
    bad = sum(v for k, v in dist.items() if k in ("批次伪影", "信息不足"))
    print("\n★ 有问题或无法判定 = %d/%d (%.1f%%)" % (bad, n, 100.0 * bad / n))

    # 涉及疾病混杂的队列 —— 用抽出后的 confound_factors（完整文本尾段）
    conf_dis = [x for x in human
                if re.search(r"disease|tumou?r|cancer|carcinoma|condition|group|diagnos|stage",
                             x.get("confound_factors") or "", re.I)]
    print("★ batchConfound 明确提及疾病/分组 = %d 个" % len(conf_dis))
    for x in conf_dis[:12]:
        print("     %-12s %-30s %s" % (x["shortName"], (x["batchEffect"] or "")[:28],
                                        (x.get("confound_factors") or "")[:70]))

    # 混杂因素词频（给 Discussion 用）
    from collections import Counter as C2
    words = []
    for x in human:
        cf = (x.get("confound_factors") or "")
        for w in re.split(r"[,;/]| and ", cf):
            w = w.strip()
            if 2 < len(w) < 40:
                words.append(w.lower())
    wfreq = C2(words)
    print("\n混杂因素词频 top10:")
    for w, c in wfreq.most_common(10):
        print("   %-40s %d" % (w, c))

    # ---- 落盘
    import pandas as pd
    df = pd.DataFrame(human)
    p1 = os.path.join(OUT, "caseA_gemma_cohorts.tsv")
    df.to_csv(p1, sep="\t", index=False)
    print("\n写出 %s  (%d 行)" % (p1, len(df)))

    summ = {
        "query": a.query,
        "n_datasets_total": len(rows),
        "n_human": n,
        "annot_dist": dict(dist),
        "n_problem_or_unjudgeable": bad,
        "share_problem_or_unjudgeable": round(bad / n, 4) if n else None,
        "n_batchconfound_disease": len(conf_dis),
        "confound_word_freq_top15": dict(wfreq.most_common(15)),
        "n_with_confound_factors_parsed": sum(1 for x in human if x.get("confound_factors")),
        "annot_map": ANNOT_MAP,
        "source": "Gemma REST /datasets 列表条目（list endpoint 自带策展字段）",
        "fixes_applied": [
            "batchConfound 不截断（关键信息在句末 'Factor(s) confounded were: ...'）",
            "geeq 取嵌套字典的 publicQualityScore（非扁平键）",
        ],
    }
    p2 = os.path.join(OUT, "caseA_landscape.json")
    with open(p2, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=2)
    print("写出 %s" % p2)


if __name__ == "__main__":
    main()
