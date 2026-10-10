# -*- coding: utf-8 -*-
"""A2 第二步：把抓到的原始 JSON 整理成一张可用于统计的队列宽表。

输入  data/raw/gemma_a2/{cohort_index.json, an_*.json, svd_*.json, samp_*_*.json}
输出  results/partA2_cohorts.tsv        队列级宽表（标注侧 + 后果侧）
      results/partA2_build_report.json  解析诊断

标签抽取要点：Gemma 样本对象里的 factorValues[] 携带结构化的实验因子值
（characteristics[].category == 'disease' / 'block' ...），比解析自由文本特征可靠得多。
"""
import glob
import json
import os
import re

import numpy as np
import pandas as pd

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
CACHE = os.path.join(ROOT, "data", "raw", "gemma_a2")
OUT = os.path.join(ROOT, "results")

DISEASE_CATS = {"disease", "disease state", "disease status", "diagnosis",
                "clinical diagnosis", "condition", "disease outcome",
                "disease subtype", "disease activity", "disease staging",
                "disease model", "subject status"}
TISSUE_CATS = {"organism part", "tissue", "tissue type", "isolated tissue", "organ",
               "biosource", "bio source", "specimen", "anatomical site"}
# 技术性因子（批次/流程）——与生物学对比相对
BATCH_CATS = {"block", "batch", "scan date", "hybridization date", "label",
              "protocol", "time point", "timepoint", "extraction batch"}
# 主因子优先级：疾病 > 疾病模型/分期 > 治疗 > 基因型 > 组织 > 细胞类型 > 时间 > 表型
PRIMARY_PRIORITY = ["disease", "disease model", "disease staging", "disease subtype",
                    "disease activity", "treatment", "genotype", "strain",
                    "organism part", "cell type", "timepoint", "phenotype",
                    "biological sex", "age", "age group", "collection of material"]

UC_RE = re.compile(r"ulcerative colitis|\buc\b", re.I)
CASE_EXTRA = re.compile(r"reference subject role|(?<!non-)(?<!in)active disease|"
                        r"colitis|inflammatory bowel|ibd", re.I)
CTRL_RE = re.compile(r"reference subject role|\bcontrol\b|\bnormal\b|healthy|"
                     r"non-?inflam|unaffected|quiescent|remission|inactive", re.I)
EXCL_RE = re.compile(r"non-?involved|non-?neoplastic|neoplas|dysplas|cancer|carcinoma|"
                     r"adenoma|tumou?r|polyp", re.I)

diag = {"n_index": 0, "n_an": 0, "n_svd": 0, "n_samp": 0, "no_de_analysis": 0,
        "no_disease_factor": 0, "svd_auc_ok": 0, "svd_batch_ok": 0,
        "n_svd_matched": 0, "errors": []}


# ------------------------------------------------------------------ 工具
def load(p):
    try:
        d = json.load(open(p, encoding="utf-8"))
        if isinstance(d, dict) and d.get("error"):
            return None
        return d
    except Exception:
        return None


def auc(score, y):
    """秩和法 AUC（y 取 0/1）"""
    y = np.asarray(y, dtype=float)
    s = np.asarray(score, dtype=float)
    ok = np.isfinite(s)
    y, s = y[ok], s[ok]
    n1, n0 = y.sum(), len(y) - y.sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank().values
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def dir_auc(score, y):
    """方向无关可分性：max(AUC, 1-AUC)"""
    a = auc(score, y)
    return np.nan if not np.isfinite(a) else float(max(a, 1 - a))


def sample_factor_table(did):
    """bioAssayId -> {(experimentalFactorId, category): value}

    读样本对象 factorValues[] 的结构化实验因子值（比自由文本特征可靠）。
    """
    per = {}
    n_files = 0
    for p in sorted(glob.glob(os.path.join(CACHE, "samp_%s_*.json" % did))):
        d = load(p)
        if not d or not d.get("data"):
            continue
        n_files += 1
        for e in d["data"]:
            s = e.get("sample") or e
            bids = s.get("bioAssayIds") or ([e.get("id")] if e.get("id") else [])
            rec = {}
            for fv in (s.get("factorValues") or s.get("factorValueObjects") or []):
                fid = fv.get("experimentalFactorId")
                for c in (fv.get("characteristics") or []):
                    cat = (c.get("category") or "").strip().lower()
                    if cat and c.get("value") is not None:
                        rec[(fid, cat)] = str(c["value"])
                if not (fv.get("characteristics")) and fv.get("value"):
                    rec[(fid, "__raw__")] = str(fv["value"])
            for bid in bids:
                tgt = per.setdefault(bid, {})
                for k, v in rec.items():
                    tgt[k] = v
    return per, n_files


def factor_auc(V, vals, min_arm=4, min_keep=12):
    """对某因子的最常见两个水平做方向无关 AUC（PC1..3）。

    返回 (auc_pc1, auc_best3, n0, n1, lvl0, lvl1) 或 None。
    """
    vals = np.array(vals, dtype=object)
    present = [v for v in vals if v is not None]
    if len(present) < min_keep:
        return None
    mc = pd.Series(present).value_counts()
    if len(mc) < 2:
        return None
    lv0, c0 = mc.index[0], int(mc.iloc[0])
    lv1, c1 = mc.index[1], int(mc.iloc[1])
    if c0 < min_arm or c1 < min_arm:
        return None
    y = np.array([1 if x == lv0 else (0 if x == lv1 else np.nan) for x in vals], dtype=float)
    keep = np.isfinite(y)
    if keep.sum() < min_keep:
        return None
    yy = y[keep].astype(int)
    npc = min(3, V.shape[1])
    aucs = [dir_auc(V[keep, i], yy) for i in range(npc)]
    if not any(np.isfinite(a) for a in aucs):
        return None
    return (aucs[0], float(np.nanmax(aucs)), c0, c1, str(lv0)[:40], str(lv1)[:40])


def svd_metrics(did, primary_key=None):
    """用 SVD 载荷 + 结构化因子标签算主轴可分性（不下载表达矩阵）。

    核心对比：同一批主成分上，「生物学主因子」（Gemma 存算 DE 指定的因子）
    与「批次因子」谁更能分开样本 —— 即批次是否占据了主导轴。
    """
    d = load(os.path.join(CACHE, "svd_%s.json" % did))
    if not d or not d.get("data"):
        return {}
    o = d["data"]
    var = np.array(o.get("variances") or [], dtype=float)
    V = np.array(o.get("vmatrix") or [], dtype=float)
    bids = o.get("bioAssayIds") or []
    if V.ndim != 2 or len(bids) != V.shape[0] or len(var) == 0:
        return {}
    diag["n_svd"] += 1
    out = {"pc1_var_frac": float(var[0] / var.sum()) if var.sum() > 0 else np.nan,
           "pc1_2_var_frac": float(var[:2].sum() / var.sum()) if var.sum() > 0 else np.nan,
           "n_svd_samples": int(V.shape[0]), "n_pcs": int(V.shape[1])}

    per, n_files = sample_factor_table(did)
    if n_files:
        diag["n_samp"] += 1
    out["n_labels_matched"] = int(sum(1 for b in bids if b in per))

    # ---- 所有可用因子的 PC 可分性
    keys = {}
    for b in bids:
        for k in per.get(b, {}):
            keys[k] = keys.get(k, 0) + 1
    n = len(bids)
    feats = {}
    for (fid, cat), cov in keys.items():
        if cat in ("__raw__",) or cov < 0.9 * n:
            continue
        vals = [per.get(b, {}).get((fid, cat)) for b in bids]
        r = factor_auc(V, vals)
        if r:
            feats[(fid, cat)] = r

    # ---- 生物学主因子：优先用 Gemma 存算 DE 指定的因子；否则按类别优先级
    bio = None
    if primary_key and primary_key in feats:
        bio = (primary_key[1], feats[primary_key])
    if bio is None:
        for pref in PRIMARY_PRIORITY:
            hit = [(c, v) for (f, c), v in feats.items() if c == pref]
            if hit:
                hit.sort(key=lambda x: -(x[1][2] + x[1][3]))
                bio = hit[0]
                break
    if bio is None:
        cand = [(c, v) for (f, c), v in feats.items() if c not in BATCH_CATS]
        if cand:
            cand.sort(key=lambda x: -(x[1][2] + x[1][3]))
            bio = cand[0]
    if bio:
        cat, (a1, ab, c0, c1, l0, l1) = bio
        out.update({"primary_cat": cat, "auc_pc1_primary": a1, "auc_best_primary": ab,
                    "primary_n0": c0, "primary_n1": c1,
                    "primary_lvl0": l0, "primary_lvl1": l1,
                    "primary_is_disease": cat in DISEASE_CATS})
        if cat in DISEASE_CATS:
            out["auc_pc1_disease"] = a1
            out["auc_disease_best_pc"] = ab

    # ---- 批次因子（技术性）：所有 block/batch 类取覆盖最好者
    bat = [(c, v) for (f, c), v in feats.items() if c in BATCH_CATS]
    if bat:
        bat.sort(key=lambda x: -(x[1][2] + x[1][3]))
        cat, (a1, ab, c0, c1, l0, l1) = bat[0]
        out.update({"batch_cat": cat, "auc_pc1_batch": a1, "auc_best_batch": ab,
                    "batch_n0": c0, "batch_n1": c1, "batch_lvl0": l0,
                    "n_batch_labels": int(len({per.get(b, {}).get(k) for b in bids
                                               for k in [kk for kk in per.get(b, {}) if kk[1] in BATCH_CATS]}))})
        if np.isfinite(a1) and out.get("auc_pc1_primary") is not None and np.isfinite(out["auc_pc1_primary"]):
            out["pc1_batch_minus_primary"] = float(a1 - out["auc_pc1_primary"])
            out["pc1_batch_wins"] = bool(a1 > out["auc_pc1_primary"])
            diag["svd_batch_ok"] += 1
    if out.get("auc_pc1_primary") is not None and np.isfinite(out["auc_pc1_primary"]):
        diag["svd_auc_ok"] += 1

    # ---- UC vs 对照的可用臂样本数（用于评估「队列池里真正能回答 UC 问题的比例」）
    uc = ct = 0
    for b in bids:
        vals = [v for (f, c), v in per.get(b, {}).items() if c in DISEASE_CATS]
        if any(re.search(r"ulcerative colitis|\buc\b", v, re.I) and not EXCL_RE.search(v) for v in vals):
            uc += 1
        elif any(CTRL_RE.search(v) for v in vals):
            ct += 1
    out["n_uc"] = int(uc)
    out["n_ctrl"] = int(ct)
    return out


def parse_analysis(did):
    """从存算 DE 分析里取后果指标"""
    d = load(os.path.join(CACHE, "an_%s.json" % did))
    if not d or not d.get("data"):
        diag["no_de_analysis"] += 1
        return {}
    diag["n_an"] += 1
    cands = []
    n_batch_models = 0
    for a in d["data"]:
        fvu = a.get("factorValuesUsed") or {}
        kinds = {}
        for fid, fvs in fvu.items():
            cats = set()
            descs = set()
            for fv in fvs:
                descs.add(str(fv.get("description") or fv.get("summary") or "")[:40])
                if fv.get("category"):
                    cats.add(str(fv["category"]).strip().lower())
                for c in (fv.get("characteristics") or []):
                    cats.add((c.get("category") or "").strip().lower())
            kinds[fid] = {"cats": cats, "n_levels": len(fvs), "descs": sorted(descs)[:3],
                          "is_disease": bool(cats & DISEASE_CATS),
                          "is_batch": bool(cats & {"block", "batch", "scan date", "label"})}
        dis_f = [f for f, v in kinds.items() if v["is_disease"]]
        bat_f = [f for f, v in kinds.items() if v["is_batch"]]
        if bat_f:
            n_batch_models += 1
        rs = (a.get("resultSets") or [])
        if not dis_f or not rs:
            continue
        rs0 = rs[0]
        f_id = dis_f[0]
        f_cat = sorted(kinds[f_id]["cats"])[0] if kinds[f_id]["cats"] else None
        try:
            f_id_int = int(f_id)
        except Exception:
            f_id_int = None
        rec = {
            "de_primary_fid": f_id_int,
            "de_primary_cat": f_cat,
            "de_analysis_id": a.get("id"),
            "de_analysis_type": a.get("name"),
            "de_is_subset": bool(a.get("isSubset")),
            "de_factor_id": dis_f[0],
            "de_factor_levels": kinds[dis_f[0]]["n_levels"],
            "de_factor_values": " | ".join(kinds[dis_f[0]]["descs"]),
            "de_has_batch_block": bool(bat_f),
            "de_n_batch_factors": len(bat_f),
            "de_n_model_factors": len(fvu),
            "de_threshold": rs0.get("threshold"),
            "n_genes_analyzed": rs0.get("numberOfGenesAnalyzed"),
            "n_probes_analyzed": rs0.get("numberOfProbesAnalyzed"),
            "n_de_probes": rs0.get("numberOfDiffExpressedProbes"),
            "n_up_probes": rs0.get("numberOfUpregulatedProbes"),
            "n_down_probes": rs0.get("numberOfDownregulatedProbes"),
        }
        tot = rec["n_probes_analyzed"] or rec["n_genes_analyzed"]
        ud = (rec["n_up_probes"] or 0) + (rec["n_down_probes"] or 0)
        try:
            rec["de_rate"] = rec["n_de_probes"] / tot
        except Exception:
            rec["de_rate"] = np.nan
        try:
            rec["de_dir_rate"] = ud / tot
        except Exception:
            rec["de_dir_rate"] = np.nan
        cands.append(rec)
    if not cands:
        diag["no_disease_factor"] += 1
        return {"n_de_analyses_total": len(d["data"]), "n_de_models_with_batch": n_batch_models}
    cands.sort(key=lambda r: (not r["de_is_subset"], r.get("n_probes_analyzed") or 0), reverse=True)
    best = cands[0]
    best["n_de_analyses_total"] = len(d["data"])
    best["n_de_models_with_batch"] = n_batch_models
    best["n_de_candidate_analyses"] = len(cands)
    return best


# ------------------------------------------------------------------ 主流程
def main():
    idx = load(os.path.join(CACHE, "cohort_index.json")) or {}
    diag["n_index"] = len(idx)
    rows = []
    for k, e in idx.items():
        did = int(k) if str(k).isdigit() else k
        chars = e.get("characteristics") or []
        bycat = {}
        for c in chars:
            cat = (c.get("category") or "").strip().lower()
            if c.get("value"):
                bycat.setdefault(cat, set()).add(str(c.get("value")))
        geeq = e.get("geeq") or {}
        bes = e.get("batchEffectStatistics") or ""
        m_p = re.search(r"p\s*=?\s*([0-9]+\.?[0-9]*(?:[eE][-+]?\d+)?)", bes)
        m_pc = re.search(r"PC\s*(\d+)", bes)
        conf = e.get("batchConfound") or ""
        conf_names = [x.strip() for x in re.split(r"[,;]", conf.split(":")[-1]) if x.strip()] if conf else []

        r = {
            "gid": did,
            "gse": e.get("shortName"),
            "name": (e.get("name") or "")[:120],
            "n_samples": e.get("numberOfBioAssays"),
            "tech": e.get("technologyType"),
            "library_strategy": ",".join(sorted({x.get("value") or "?" for x in (e.get("libraryStrategies") or [])})),
            "n_platforms": len(e.get("platforms") or []),
            "platform": ",".join(sorted({p.get("shortName") or "" for p in (e.get("platforms") or [])})),
            "is_single_cell": bool(e.get("isSingleCell")),
            "date_created": (e.get("dateCreated") or "")[:10],
            "last_updated": (e.get("lastUpdated") or "")[:10],
            "troubled": bool(e.get("troubled")),
            "needs_attention": bool(e.get("needsAttention")),
            # ---- 标注侧
            "batch_effect": e.get("batchEffect"),
            "be_stat_text": bes,
            "be_p": float(m_p.group(1)) if m_p else np.nan,
            "be_pc": int(m_pc.group(1)) if m_pc else np.nan,
            "batch_confound": conf,
            "confound_names": ";".join(conf_names),
            "confound_disease": bool(re.search(r"disease|diagnosis|condition", conf, re.I)),
            "n_confound_factors": len(conf_names),
            "geeq_public": geeq.get("publicQualityScore"),
            "geeq_batch_info": geeq.get("qScoreBatchInfo"),
            "geeq_pub_batch_effect": geeq.get("qScorePublicBatchEffect"),
            "geeq_pub_batch_confound": geeq.get("qScorePublicBatchConfound"),
            "geeq_replicates": geeq.get("qScoreReplicates"),
            "geeq_platforms": geeq.get("qScorePlatformsTech"),
            "geeq_outliers": geeq.get("qScoreOutliers"),
            "geeq_mean_corr": geeq.get("qScoreSampleMeanCorrelation"),
            "geeq_corr_var": geeq.get("qScoreSampleCorrelationVariance"),
            "batch_corrected_flag": bool(geeq.get("batchCorrected")),
            "disease_values": " | ".join(sorted(bycat.get("disease", set()))[:6]),
            "tissue_values": " | ".join(sorted(bycat.get("organism part", set()) | bycat.get("tissue", set()))[:4]),
            "is_uc_query": bool(UC_RE.search(" ".join(sorted(bycat.get("disease", set())))) or UC_RE.search(e.get("name") or "")),
        }
        # ---- 后果侧
        de = parse_analysis(did)
        r.update(de)
        pk = None
        if de.get("de_primary_fid") is not None and de.get("de_primary_cat"):
            pk = (de["de_primary_fid"], de["de_primary_cat"])
        r.update(svd_metrics(did, primary_key=pk))
        rows.append(r)
        if len(rows) % 20 == 0:
            print("  ... %d/%d" % (len(rows), len(idx)), flush=True)

    df = pd.DataFrame(rows)
    cls_map = {
        "BATCH_EFFECT_FAILURE": "失败",
        "PROBLEMATIC_BATCH_INFO_FAILURE": "失败",
        "SINGLETON_BATCHES_FAILURE": "失败",
        "UNINFORMATIVE_HEADERS_FAILURE": "失败",
        "NO_BATCH_EFFECT_SUCCESS": "无效应",
        "SINGLE_BATCH_SUCCESS": "单批次",
        "BATCH_CORRECTED_SUCCESS": "已校正",
        "NO_BATCH_INFO": "信息缺失",
        "BATCH_EFFECT_UNDETERMINED_FAILURE": "无法判定",
        "NO_BATCH_EFFECT": "无效应",
    }
    df["batch_class"] = df["batch_effect"].map(lambda x: cls_map.get(x, "其他" if x else "未标注"))
    # 粗分类：区分「查出批次伪影」「批次信息不可用」「无批次问题」「已校正」
    coarse = {
        "BATCH_EFFECT_FAILURE": "批次伪影",
        "NO_BATCH_INFO": "信息不足",
        "UNINFORMATIVE_HEADERS_FAILURE": "信息不足",
        "PROBLEMATIC_BATCH_INFO_FAILURE": "信息不足",
        "SINGLETON_BATCHES_FAILURE": "信息不足",
        "BATCH_EFFECT_UNDETERMINED_FAILURE": "信息不足",
        "NO_BATCH_EFFECT_SUCCESS": "无问题",
        "SINGLE_BATCH_SUCCESS": "无问题",
        "BATCH_CORRECTED_SUCCESS": "已校正",
    }
    df["annot_status"] = df["batch_effect"].map(lambda x: coarse.get(x, "其他" if x else "未标注"))
    df.to_csv(os.path.join(OUT, "partA2_cohorts.tsv"), sep="\t", index=False)
    json.dump(diag, open(os.path.join(OUT, "partA2_build_report.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("=" * 70)
    print("写出 %s（%d 行 × %d 列）" % (os.path.join(OUT, "partA2_cohorts.tsv"), len(df), df.shape[1]))
    print("诊断：", diag)
    print("\n批次标注分布：")
    print(df["batch_class"].value_counts().to_string())
    print("\n有主轴可分性（主因子）的队列：%d / %d" % (df["auc_pc1_primary"].notna().sum() if "auc_pc1_primary" in df else 0, len(df)))
    print("有批次因子可分性的队列：%d / %d" % (df["auc_pc1_batch"].notna().sum() if "auc_pc1_batch" in df else 0, len(df)))
    if "auc_pc1_primary" in df:
        print(df.groupby("batch_class")[["auc_pc1_primary", "auc_pc1_batch"]].agg(["count", "median"]).to_string())
    if "pc1_batch_wins" in df:
        sub = df[df["pc1_batch_wins"].notna()]
        print("\nPC1 上「批次可分性 > 生物学可分性」的比例：")
        print(sub.groupby("batch_class")["pc1_batch_wins"].agg(["count", "mean"]).to_string())
    print("\n有存算 DE 后果数据的队列：%d / %d" % (df["de_rate"].notna().sum(), len(df)))


if __name__ == "__main__":
    main()
