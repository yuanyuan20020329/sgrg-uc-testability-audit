# -*- coding: utf-8 -*-
"""A2 第三步：把「批次质控标注」与「下游数据侧后果」连起来。

可用后果（均为不下载表达矩阵即可得到的量）：
  ① 主轴可分性（SVD 载荷 + 结构化样本因子值）
     auc_pc1_primary  队列生物学主因子在第 1 主成分上的可分性
     auc_pc1_batch    批次因子在同一主成分上的可分性
     pc1_batch_wins   批次比生物学更能分开样本（主轴上批次占了上风）
  ② 存算 DE 信号规模 de_rate（仅 44 个队列有疾病对比）

主要问题：
  Q1 标注被判「批次失败」的队列，其 PC1 是否真的由批次主导？
  Q2 这类队列的生物学信号在主轴上是否更难恢复？
  Q3 标注与数据侧后果是否**不一致**（标注漏检）？规模多大？
  Q4 阳性/阴性对照：已校正组应回落；标注不应预测队列规模。

产出 results/partA2_*.tsv + partA2_analysis.json
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT = os.path.join(ROOT, "results")
NBOOT = 4000
NPERM = 20000

LEVELS = ["无效应", "失败", "无法判定", "单批次", "已校正", "信息缺失", "未标注", "其他"]
PROBLEM = ["失败", "无法判定", "信息缺失"]
CLEAN = ["无效应", "已校正", "单批次"]
LEVELS2 = ["批次伪影", "信息不足", "无问题", "已校正"]
PROBLEM2 = ["批次伪影", "信息不足"]


# ------------------------------------------------------------------ 基础工具
def ols(X, y):
    XtXi = np.linalg.pinv(X.T @ X)
    return XtXi @ (X.T @ y), XtXi


def design(df, terms):
    """terms: ('num',name) 或 ('cat',name,ref,levels)"""
    cols = [np.ones(len(df))]
    names = ["intercept"]
    for t in terms:
        if t[0] == "num":
            v = pd.to_numeric(df[t[1]], errors="coerce").values.astype(float)
            med = np.nanmedian(v) if np.isfinite(np.nanmedian(v)) else 0.0
            cols.append(np.nan_to_num(v, nan=med))
            names.append(t[1])
        else:
            s = df[t[1]].astype(str)
            for x in [z for z in t[3] if z != t[2] and z in set(s.unique())]:
                cols.append((s == x).values.astype(float))
                names.append("%s[%s]" % (t[1], x))
    return np.column_stack(cols), names


def ols_report(df, ycol, terms, ref_levels=None, nboot=NBOOT, seed=11, nperm=3000):
    """OLS + bootstrap CI + 对第一个分类项做置换 F 检验"""
    d = df.copy()
    y = pd.to_numeric(d[ycol], errors="coerce").values.astype(float)
    d = d[np.isfinite(y)]
    y = y[np.isfinite(y)]
    if len(d) < 12:
        return None
    X, names = design(d, terms)
    beta, _ = ols(X, y)
    rss_f = float(np.sum((y - X @ beta) ** 2))
    r2 = 1 - rss_f / np.sum((y - y.mean()) ** 2) if np.sum((y - y.mean()) ** 2) > 0 else np.nan

    rng = np.random.default_rng(seed)
    B = []
    for _ in range(nboot):
        idx = rng.integers(0, len(d), len(d))
        Xi, nm = design(d.iloc[idx], terms)
        yi = pd.to_numeric(d.iloc[idx][ycol], errors="coerce").values.astype(float)
        if not np.all(np.isfinite(yi)):
            continue
        try:
            b, _ = ols(Xi, yi)
            if len(b) == len(names):
                B.append(b)
        except Exception:
            pass
    B = np.array(B)
    if len(B) >= 100 and B.shape[1] == len(names):
        lo, hi = np.percentile(B, 2.5, axis=0), np.percentile(B, 97.5, axis=0)
    else:
        lo = hi = np.full(len(names), np.nan)
    tab = pd.DataFrame({"term": names, "beta": beta, "ci_lo": lo, "ci_hi": hi})

    # 置换检验：第一个分类项（batch_class / confound 标志）的联合效应
    p_perm = np.nan
    ftest = None
    cat_terms = [t for t in terms if t[0] == "cat"]
    if cat_terms:
        t0 = cat_terms[0]
        red = [t for t in terms if t is not t0]
        Xr, _ = design(d, red)
        br, _ = ols(Xr, y)
        rss_r = float(np.sum((y - Xr @ br) ** 2))
        df1 = X.shape[1] - Xr.shape[1]
        df2 = len(d) - X.shape[1]
        if df1 > 0 and df2 > 0:
            Fobs = ((rss_r - rss_f) / df1) / (rss_f / df2)
            cnt = 0
            tot = 0
            for _ in range(nperm):
                dp = d.copy()
                dp[t0[1]] = rng.permutation(dp[t0[1]].values)
                Xp, _ = design(dp, terms)
                try:
                    bp, _ = ols(Xp, y)
                except Exception:
                    continue
                rss1 = float(np.sum((y - Xp @ bp) ** 2))
                Xrp, _ = design(dp, red)
                brp, _ = ols(Xrp, y)
                rss0 = float(np.sum((y - Xrp @ brp) ** 2))
                fp = ((rss0 - rss1) / df1) / (rss1 / df2) if rss1 > 0 else 0
                tot += 1
                cnt += fp >= Fobs
            p_perm = (cnt + 1) / (tot + 1) if tot else np.nan
            ftest = {"F": float(Fobs), "df1": int(df1), "df2": int(df2), "p_perm": float(p_perm),
                     "term": t0[1]}
    return {"outcome": ycol, "n": int(len(d)), "R2": float(r2) if np.isfinite(r2) else None,
            "terms": tab.round(4).to_dict("records"), "F": ftest}


def group_stat(df, col, group_col="annot_status", min_n=3):
    rows = []
    for g, sub in df.groupby(group_col):
        v = pd.to_numeric(sub[col], errors="coerce").dropna()
        if len(v) < min_n:
            continue
        rows.append({"group": g, "n": int(len(v)), "median": float(v.median()),
                     "q1": float(v.quantile(.25)), "q3": float(v.quantile(.75)),
                     "mean": float(v.mean())})
    return pd.DataFrame(rows).sort_values("median", ascending=False)


def two_group(df, col, a_mask, b_mask, label=""):
    a = pd.to_numeric(df.loc[a_mask, col], errors="coerce").dropna()
    b = pd.to_numeric(df.loc[b_mask, col], errors="coerce").dropna()
    if len(a) < 3 or len(b) < 3:
        return None
    u, p = stats.mannwhitneyu(a, b)
    rng = np.random.default_rng(5)
    bd = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(4000)]
    return {"label": label, "col": col, "n_a": int(len(a)), "n_b": int(len(b)),
            "med_a": float(np.median(a)), "med_b": float(np.median(b)),
            "diff": float(np.median(a) - np.median(b)),
            "ci": [float(np.percentile(bd, 2.5)), float(np.percentile(bd, 97.5))],
            "p_mwu": float(p), "auc": float(u / (len(a) * len(b)))}


def fisher(a, b, c, d):
    odds, p = stats.fisher_exact([[a, b], [c, d]])
    return {"odds_ratio": float(odds), "p": float(p)}


# ------------------------------------------------------------------ 主分析
def main():
    df = pd.read_csv(os.path.join(OUT, "partA2_cohorts.tsv"), sep="\t")
    res = {}

    # ---------- 派生变量
    df["log2n"] = np.log2(pd.to_numeric(df["n_samples"], errors="coerce"))
    df["is_gut"] = df["tissue_values"].fillna("").str.contains(
        r"colon|rectum|ileum|intestin|biopsy|cecum|sigmoid|gut", case=False).astype(int)
    df["is_rnaseq"] = df["library_strategy"].fillna("").str.lower().str.contains("rna_seq").astype(int)
    df["is_twocolor"] = (df["tech"].astype(str) == "TWOCOLOR").astype(int)
    df["log_ndeg"] = np.log2(pd.to_numeric(df["n_de_probes"], errors="coerce").clip(lower=1))
    df["year"] = pd.to_numeric(df["date_created"].astype(str).str[:4], errors="coerce")
    df["is_problem"] = df["annot_status"].isin(PROBLEM2).astype(int)
    df["is_clean"] = df["annot_status"].isin(["无问题", "已校正"]).astype(int)
    df["normalized_batch_auc"] = df["auc_pc1_batch"] - df["auc_pc1_primary"]

    # 数据侧「批次主导」定义。
    # ★ 口径修正（2026-10-09）：Gemma 官方文档明确 batch factor「almost always correlates
    #   with one of the first three principal components」，故主体口径必须用 PC1-3 的最佳值
    #   （auc_best_batch），只取 PC1 会系统性低估其标注灵敏度、人为制造"过报"假象。
    ok = df["auc_best_batch"].notna() & df["auc_best_primary"].notna()
    df["svd_batch_dominant"] = np.where(
        ok, (df["auc_best_batch"] >= 0.75) & (df["auc_best_batch"] > df["auc_best_primary"] + 0.05),
        np.nan)
    # 旧 PC1-only 口径保留为敏感性对照
    ok1 = df["auc_pc1_batch"].notna() & df["auc_pc1_primary"].notna()
    df["svd_batch_dominant_pc1only"] = np.where(
        ok1, (df["auc_pc1_batch"] >= 0.75) & (df["auc_pc1_batch"] > df["auc_pc1_primary"] + 0.05),
        np.nan)

    res["n_total"] = int(len(df))

    # ---------- 1. 队列池质控地景（描述）
    land = df.groupby("annot_status").agg(
        n=("gid", "size"),
        med_samples=("n_samples", "median"),
        n_batch_confound=("batch_confound", lambda s: int(s.fillna("").str.len().gt(0).sum())),
        n_outlier_flag=("troubled", "sum"),
        med_geeq_public=("geeq_public", "median"),
        med_geeq_batch_eff=("geeq_pub_batch_effect", "median"),
        med_geeq_batch_conf=("geeq_pub_batch_confound", "median"),
        n_svd_ok=("auc_pc1_primary", lambda s: int(s.notna().sum())),
        med_pc1_var=("pc1_var_frac", "median"),
    ).reset_index()
    land["pct_confound"] = (land["n_batch_confound"] / land["n"]).round(3)
    land["pct_of_pool"] = (land["n"] / len(df)).round(3)
    land.to_csv(os.path.join(OUT, "partA2_landscape.tsv"), sep="\t", index=False)
    res["landscape"] = land.to_dict("records")
    res["pool_summary"] = {
        "n": int(len(df)),
        "n_problem": int(df["is_problem"].sum()),
        "n_with_batch_confound_text": int(df["batch_confound"].fillna("").str.len().gt(0).sum()),
        "n_confound_disease": int(df["confound_disease"].sum()),
        "n_uc_vs_ctrl_usable": int(((df["n_uc"] >= 4) & (df["n_ctrl"] >= 4)).sum()),
        "n_with_de_disease": int(df["de_rate"].notna().sum()),
        "n_with_pc1_batch": int(df["auc_pc1_batch"].notna().sum()),
    }

    # ---------- 2. 主轴后果：按标注类别的分布
    g_batch_auc = group_stat(df, "auc_pc1_batch")
    g_bio_auc = group_stat(df, "auc_pc1_primary")
    g_norm = group_stat(df, "normalized_batch_auc")
    for g, nm in ((g_batch_auc, "auc_pc1_batch"), (g_bio_auc, "auc_pc1_primary"),
                  (g_norm, "normalized_batch_auc")):
        g.to_csv(os.path.join(OUT, "partA2_group_%s.tsv" % nm), sep="\t", index=False)
    res["group_auc_pc1_batch"] = g_batch_auc.to_dict("records")
    res["group_auc_pc1_primary"] = g_bio_auc.to_dict("records")
    res["group_normalized"] = g_norm.to_dict("records")

    # ---------- 3. 两两对比：问题组 vs 干净组
    sv = df.dropna(subset=["auc_pc1_batch", "auc_pc1_primary"]).copy()
    res["contrasts"] = []
    for col in ["auc_pc1_batch", "auc_pc1_primary", "normalized_batch_auc", "pc1_var_frac"]:
        r = two_group(sv, col, sv["is_problem"] == 1, sv["is_clean"] == 1, "问题组 vs 干净组")
        if r:
            res["contrasts"].append(r)
    for col in ["auc_pc1_batch", "auc_pc1_primary", "normalized_batch_auc"]:
        r = two_group(sv, col, sv["annot_status"] == "批次伪影", sv["annot_status"] == "已校正",
                      "批次伪影 vs 已校正")
        if r:
            res["contrasts"].append(r)
        r = two_group(sv, col, sv["confound_disease"] == True, sv["confound_disease"] == False,  # noqa: E712
                      "批次-疾病混杂 vs 无混杂")
        if r:
            res["contrasts"].append(r)

    # ---------- 4. 调整模型
    COV_SVD = [("num", "log2n"), ("num", "is_rnaseq"), ("num", "is_twocolor"),
               ("num", "is_gut"), ("num", "pc1_var_frac"), ("num", "n_svd_samples")]
    m_svd_batch = ols_report(sv, "auc_pc1_batch",
                             [("cat", "batch_class", "无效应", LEVELS)] + COV_SVD)
    m_svd_bio = ols_report(sv, "auc_pc1_primary",
                           [("cat", "batch_class", "无效应", LEVELS)] + COV_SVD)
    m_norm = ols_report(sv, "normalized_batch_auc",
                        [("cat", "batch_class", "无效应", LEVELS)] + COV_SVD)
    res["models"] = {"auc_pc1_batch": m_svd_batch, "auc_pc1_primary": m_svd_bio,
                     "normalized_batch_auc": m_norm}

    # ---------- 5. 不一致性：标注 vs 数据侧后果（2×2，McNemar）
    disc = df.dropna(subset=["svd_batch_dominant"]).copy()
    disc["svd_batch_dominant"] = disc["svd_batch_dominant"].astype(bool)
    disc["annot_problem"] = disc["is_problem"] == 1
    a = int(((disc["annot_problem"]) & (disc["svd_batch_dominant"])).sum())          # 一致：都阳性
    b = int(((disc["annot_problem"]) & (~disc["svd_batch_dominant"])).sum())         # 标注阳性/数据阴性
    c = int(((~disc["annot_problem"]) & (disc["svd_batch_dominant"])).sum())         # 标注阴性/数据阳性（漏检）
    d0 = int(((~disc["annot_problem"]) & (~disc["svd_batch_dominant"])).sum())       # 一致：都阴性
    nmis = len(disc)
    res["discordance"] = {
        "n": int(nmis), "table": {"both_pos": a, "annot_only": b, "svd_only": c, "both_neg": d0},
        "sens": float(a / (a + b)) if (a + b) else None,
        "spec": float(d0 / (c + d0)) if (c + d0) else None,
        "ppv": float(a / (a + c)) if (a + c) else None,
        "missed_n": c, "missed_frac_of_negatives": float(c / (c + d0)) if (c + d0) else None,
        "fisher": fisher(a, b, c, d0),
        "criterion": "PC1-3 最佳批次可分性（对齐 Gemma 官方口径）",
    }
    # 敏感性对照：仅用 PC1（旧口径，会低估灵敏度）
    d1 = df.dropna(subset=["svd_batch_dominant_pc1only"]).copy()
    d1["dom"] = d1["svd_batch_dominant_pc1only"].astype(bool)
    d1["pos"] = d1["is_problem"] == 1
    a1 = int((d1["pos"] & d1["dom"]).sum()); b1 = int((d1["pos"] & ~d1["dom"]).sum())
    c1 = int((~d1["pos"] & d1["dom"]).sum()); e1 = int((~d1["pos"] & ~d1["dom"]).sum())
    res["discordance_pc1only"] = {
        "n": int(len(d1)), "table": {"both_pos": a1, "annot_only": b1, "svd_only": c1, "both_neg": e1},
        "sens": float(a1 / (a1 + b1)) if (a1 + b1) else None,
        "spec": float(e1 / (c1 + e1)) if (c1 + e1) else None,
        "ppv": float(a1 / (a1 + c1)) if (a1 + c1) else None,
        "criterion": "仅 PC1（旧口径；因 Gemma 用 PC1-3 判定，此口径系统性低估灵敏度）",
    }
    # 漏检组特征
    miss = disc[(~disc["annot_problem"]) & (disc["svd_batch_dominant"])]
    res["discordance"]["missed_profile"] = {
        "n": int(len(miss)),
        "median_n": float(miss["n_samples"].median()) if len(miss) else None,
        "by_class": miss["annot_status"].value_counts().to_dict() if len(miss) else {},
        "median_batch_auc": float(miss["auc_pc1_batch"].median()) if len(miss) else None,
        "median_bio_auc": float(miss["auc_pc1_primary"].median()) if len(miss) else None,
        "pct_rnaseq": float(miss["is_rnaseq"].mean()) if len(miss) else None,
    }
    miss.to_csv(os.path.join(OUT, "partA2_discordant_missed.tsv"), sep="\t", index=False)

    # ---------- 6. pc1_batch_wins 比例（Fisher vs 干净组）
    win = df.dropna(subset=["pc1_batch_wins"]).copy()
    pa = int(((win["is_problem"] == 1) & (win["pc1_batch_wins"] == True)).sum())   # noqa: E712
    pb = int(((win["is_problem"] == 1) & (win["pc1_batch_wins"] == False)).sum())  # noqa: E712
    ca = int(((win["is_clean"] == 1) & (win["pc1_batch_wins"] == True)).sum())     # noqa: E712
    cb = int(((win["is_clean"] == 1) & (win["pc1_batch_wins"] == False)).sum())    # noqa: E712
    res["pc1_batch_wins"] = {
        "n": int(len(win)), "problem": [pa, pb], "clean": [ca, cb],
        "frac_problem": float(pa / (pa + pb)) if (pa + pb) else None,
        "frac_clean": float(ca / (ca + cb)) if (ca + cb) else None,
        "fisher": fisher(pa, pb, ca, cb),
    }

    # ---------- 7. （数据侧连续打分的交叉验证见第 15 节）

    # ---------- 8. 存算 DE 后果（次要，n 小）
    dec = df[(df["de_rate"].notna()) & (df["is_single_cell"] != True)].copy()  # noqa: E712
    dec["y_de"] = dec["de_rate"]
    dec.to_csv(os.path.join(OUT, "partA2_de_subset.tsv"), sep="\t", index=False)
    if len(dec) >= 15:
        m_de = ols_report(dec, "y_de",
                          [("cat", "annot_status", "无问题", LEVELS2),
                           ("num", "log2n"), ("num", "is_rnaseq"), ("num", "is_twocolor"),
                           ("num", "is_gut")])
        res["de_model"] = m_de
        res["de_contrasts"] = [x for x in [
            two_group(dec, "de_rate", dec["is_problem"] == 1, dec["is_clean"] == 1, "DE率 问题 vs 干净"),
            two_group(dec, "de_rate", dec["annot_status"] == "批次伪影",
                      dec["annot_status"] == "无问题", "DE率 批次伪影 vs 无问题"),
        ] if x]

    # ---------- 9. 阴性对照
    nc = {}
    nc["n_samples_vs_class"] = two_group(df, "n_samples", df["is_problem"] == 1, df["is_clean"] == 1,
                                         "样本量 问题 vs 干净（应无差）")
    nc["year_vs_class"] = two_group(df.dropna(subset=["year"]), "year",
                                    df.dropna(subset=["year"])["is_problem"] == 1,
                                    df.dropna(subset=["year"])["is_clean"] == 1,
                                    "年份 问题 vs 干净")
    nc["ncbi_pool_frac_1kg"] = None
    res["negcontrol"] = nc

    # ---------- 10. 混杂结构：标注 与 队列属性 的关联
    conf_tab = []
    for col in ["is_rnaseq", "is_twocolor", "is_single_cell", "confound_disease", "troubled"]:
        ct = pd.crosstab(df["is_problem"], df[col].astype(bool))
        if ct.shape == (2, 2):
            r = fisher(int(ct.iloc[0, 0]), int(ct.iloc[0, 1]), int(ct.iloc[1, 0]), int(ct.iloc[1, 1]))
            conf_tab.append({"var": col, "problem_true": [int(ct.iloc[0, 1]), int(ct.iloc[0, 0])],
                             "clean_true": [int(ct.iloc[1, 1]), int(ct.iloc[1, 0])], **r})
    res["confounding"] = conf_tab

    # ---------- 11. 稳健性：粗分类 + 二分类模型
    df["status_problem"] = df["annot_status"].isin(PROBLEM2).astype(int)
    sv2 = df.dropna(subset=["auc_pc1_batch", "auc_pc1_primary"]).copy()
    COVB = [("num", "log2n"), ("num", "is_rnaseq"), ("num", "is_twocolor"),
            ("num", "is_gut"), ("num", "pc1_var_frac"), ("num", "n_svd_samples")]
    res["models_binary"] = {
        "auc_pc1_batch": ols_report(sv2, "auc_pc1_batch",
                                    [("cat", "annot_status", "无问题",
                                      ["批次伪影", "信息不足", "无问题", "已校正"])] + COVB),
        "auc_pc1_primary": ols_report(sv2, "auc_pc1_primary",
                                      [("cat", "annot_status", "无问题",
                                        ["批次伪影", "信息不足", "无问题", "已校正"])] + COVB),
        "normalized_batch_auc": ols_report(sv2, "normalized_batch_auc",
                                           [("cat", "annot_status", "无问题",
                                             ["批次伪影", "信息不足", "无问题", "已校正"])] + COVB),
    }
    res["group_by_status"] = {c: group_stat(df, c, "annot_status").to_dict("records")
                              for c in ["auc_pc1_batch", "auc_pc1_primary", "normalized_batch_auc"]}
    res["status_contrasts"] = []
    for col in ["auc_pc1_batch", "auc_pc1_primary", "normalized_batch_auc"]:
        for a_name, b_name in [("批次伪影", "无问题"), ("批次伪影", "已校正"), ("信息不足", "无问题")]:
            r = two_group(sv2, col, sv2["annot_status"] == a_name, sv2["annot_status"] == b_name,
                          "%s vs %s" % (a_name, b_name))
            if r:
                res["status_contrasts"].append(r)

    # ---------- 12. 阈值敏感性（不一致表随口径变化）
    sens_rows = []
    base = df.dropna(subset=["svd_batch_dominant"]).copy()
    base["svd_batch_dominant"] = base["svd_batch_dominant"].astype(bool)
    base["annot_pos"] = base["annot_status"].isin(PROBLEM2)
    for thr in [0.65, 0.70, 0.75, 0.80, 0.85, 0.90]:
        dom = (base["auc_best_batch"] >= thr) & (base["auc_best_batch"] > base["auc_best_primary"] + 0.05)
        a = int((base["annot_pos"] & dom).sum())
        b = int((base["annot_pos"] & ~dom).sum())
        c = int((~base["annot_pos"] & dom).sum())
        d0 = int((~base["annot_pos"] & ~dom).sum())
        sens_rows.append({"thr": thr, "both_pos": a, "annot_only": b, "svd_only": c, "both_neg": d0,
                          "sens": a / (a + b) if a + b else None,
                          "spec": d0 / (c + d0) if c + d0 else None,
                          "ppv": a / (a + c) if a + c else None,
                          "missed": c, "fisher_p": fisher(a, b, c, d0)["p"]})
    res["threshold_sensitivity"] = sens_rows
    pd.DataFrame(sens_rows).to_csv(os.path.join(OUT, "partA2_threshold_sensitivity.tsv"),
                                   sep="\t", index=False)

    # ---------- 13. 匹配对照（按 log2n + pc1_var_frac 最近邻配对）
    def match_pairs(dfx, a_mask, b_mask, out_col, keys=("log2n", "pc1_var_frac")):
        A, B = dfx[a_mask].copy(), dfx[b_mask].copy()
        used, rows = set(), []
        for i, ra in A.iterrows():
            best, bd = None, np.inf
            for j, rb in B.iterrows():
                if j in used:
                    continue
                dd = sum(abs(float(ra[k]) - float(rb[k])) for k in keys)
                if dd < bd:
                    bd, best = dd, j
            if best is not None:
                used.add(best)
                rows.append({"a": ra[out_col], "b": dfx.loc[best, out_col],
                             "dist": bd, "gid_a": ra["gid"], "gid_b": dfx.loc[best, "gid"]})
        if len(rows) < 6:
            return None
        P = pd.DataFrame(rows)
        try:
            w = stats.wilcoxon(P["a"], P["b"])
            wp = float(w.pvalue)
        except Exception:
            wp = np.nan
        return {"n_pairs": int(len(P)), "med_a": float(P["a"].median()), "med_b": float(P["b"].median()),
                "median_diff": float((P["a"] - P["b"]).median()), "wilcoxon_p": wp,
                "n_win": int((P["a"] > P["b"]).sum())}
    res["matched"] = {
        "batch_auc_伪影vs无问题": match_pairs(sv2, sv2["annot_status"] == "批次伪影",
                                             sv2["annot_status"] == "无问题", "auc_pc1_batch"),
        "batch_auc_伪影vs已校正": match_pairs(sv2, sv2["annot_status"] == "批次伪影",
                                             sv2["annot_status"] == "已校正", "auc_pc1_batch"),
        "norm_伪影vs无问题": match_pairs(sv2, sv2["annot_status"] == "批次伪影",
                                        sv2["annot_status"] == "无问题", "normalized_batch_auc"),
    }

    # ---------- 14. UC 专属子集（UC vs 对照，每臂>=4）
    uc = df[(pd.to_numeric(df["n_uc"], errors="coerce") >= 4) &
            (pd.to_numeric(df["n_ctrl"], errors="coerce") >= 4)].copy()
    uc["auc_pc1_disease_use"] = uc["auc_pc1_disease"].fillna(uc["auc_pc1_primary"])
    res["uc_subset"] = {
        "n": int(len(uc)),
        "by_status": uc["annot_status"].value_counts().to_dict(),
        "disease_auc_by_status": group_stat(uc.dropna(subset=["auc_pc1_disease_use"]),
                                            "auc_pc1_disease_use", "annot_status", min_n=2).to_dict("records"),
        "list": uc[["gid", "gse", "n_samples", "annot_status", "n_uc", "n_ctrl",
                    "auc_pc1_primary", "auc_pc1_batch", "auc_pc1_disease_use"]].to_dict("records"),
    }
    uc.to_csv(os.path.join(OUT, "partA2_uc_subset.tsv"), sep="\t", index=False)

    # ---------- 15. GEEQ 综合质量分与主轴后果（数据侧交叉验证）
    for gcol in ["geeq_public", "geeq_pub_batch_confound", "geeq_mean_corr", "geeq_outliers"]:
        sub = df.dropna(subset=[gcol, "auc_pc1_batch"])
        if len(sub) >= 15 and sub[gcol].nunique() > 3:
            rho, p = stats.spearmanr(sub[gcol], sub["auc_pc1_batch"])
            res.setdefault("geeq_vs_svd", []).append(
                {"score": gcol, "n": int(len(sub)), "rho": float(rho), "p": float(p),
                 "n_unique": int(sub[gcol].nunique())})

    # ---------- 16. 配对子集（两个主轴指标都可用）—— 头对头比较的正式口径
    paird = df.dropna(subset=["auc_pc1_batch", "auc_pc1_primary"]).copy()
    paird["normalized_batch_auc"] = paird["auc_pc1_batch"] - paird["auc_pc1_primary"]
    paird["pc1_batch_wins"] = paird["auc_pc1_batch"] > paird["auc_pc1_primary"]
    res["paired_subset"] = {
        "n": int(len(paird)),
        "by_status": paird["annot_status"].value_counts().to_dict(),
        "medians": {},
    }
    for g in ["批次伪影", "信息不足", "无问题", "已校正"]:
        sub = paird[paird["annot_status"] == g]
        if len(sub) < 3:
            continue
        res["paired_subset"]["medians"][g] = {
            "n": int(len(sub)),
            "auc_batch": float(sub["auc_pc1_batch"].median()),
            "auc_bio": float(sub["auc_pc1_primary"].median()),
            "auc_diff": float(sub["normalized_batch_auc"].median()),
            "frac_batch_wins": float(sub["pc1_batch_wins"].mean()),
        }
    paird.to_csv(os.path.join(OUT, "partA2_paired_subset.tsv"), sep="\t", index=False)

    json.dump(res, open(os.path.join(OUT, "partA2_analysis.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1, default=float)

    # ---------------------------------------------------------------- 打印
    L = 78
    print("=" * L)
    print("A2：公共 UC 队列池的批次质控标注 × 下游数据侧后果")
    print("=" * L)
    print("队列池：%d 个" % res["n_total"])
    print("  标注为问题类（批次伪影/信息不足）：%d" % res["pool_summary"]["n_problem"])
    print("  有批次混杂文本：%d（其中含疾病混杂 %d）" %
          (res["pool_summary"]["n_with_batch_confound_text"], res["pool_summary"]["n_confound_disease"]))
    print("  能构成 UC vs 对照（每臂>=4）的：%d" % res["pool_summary"]["n_uc_vs_ctrl_usable"])
    print("  有批次主轴可分性的：%d" % res["pool_summary"]["n_with_pc1_batch"])
    print("\n[地景] 按标注类别：")
    print(land[["annot_status", "n", "pct_of_pool", "med_samples", "n_batch_confound",
                "med_geeq_public", "med_geeq_batch_eff", "n_svd_ok"]].to_string(index=False))

    print("\n[Q1] PC1 上批次可分性 auc_pc1_batch（按标注类别）：")
    print(g_batch_auc.to_string(index=False))
    print("\n[Q2] PC1 上生物学可分性 auc_pc1_primary：")
    print(g_bio_auc.to_string(index=False))
    print("\n[Q1/Q2] 归一化差值 = auc_pc1_batch - auc_pc1_primary：")
    print(g_norm.to_string(index=False))

    print("\n[两两对比]")
    for r in res["contrasts"]:
        print("  %-22s %-22s 中位差 %+0.3f [%+0.3f,%+0.3f] p=%.3g n=%d/%d" % (
            r["col"], r["label"], r["diff"], r["ci"][0], r["ci"][1], r["p_mwu"], r["n_a"], r["n_b"]))

    print("\n[调整模型]")
    for k, m in res["models"].items():
        if not m:
            continue
        print("  %-22s n=%-3d R2=%.3f  F=%.2f p_perm=%.4f" % (
            k, m["n"], m["R2"] or float("nan"), m["F"]["F"], m["F"]["p_perm"]))
        for t in m["terms"]:
            if t["term"].startswith("batch_class"):
                print("      %-30s %+0.4f [%+0.4f,%+0.4f]" % (t["term"], t["beta"], t["ci_lo"], t["ci_hi"]))

    dd = res["discordance"]
    print("\n[Q3] 标注 vs 数据侧主轴后果（批次主导 = PC1-3 最佳批次可分性>=0.75 且 > 生物学+0.05）")
    print("  n=%d  一致阳性=%d  仅标注阳性=%d  仅数据阳性(漏检)=%d  一致阴性=%d" % (
        dd["n"], dd["table"]["both_pos"], dd["table"]["annot_only"], dd["table"]["svd_only"],
        dd["table"]["both_neg"]))
    print("  灵敏 %.3f  特异 %.3f  PPV %.3f  漏检占阴性比例 %.3f  Fisher p=%.4g" % (
        dd["sens"] or float("nan"), dd["spec"] or float("nan"), dd["ppv"] or float("nan"),
        dd["missed_frac_of_negatives"] or float("nan"), dd["fisher"]["p"]))
    print("  漏检组：", json.dumps(dd["missed_profile"], ensure_ascii=False)[:300])

    w = res["pc1_batch_wins"]
    print("\n[Q1b] PC1 上批次胜过生物学 的比例")
    print("  问题组 %.3f (%d/%d) vs 干净组 %.3f (%d/%d)  Fisher p=%.4g" % (
        w["frac_problem"] or float("nan"), w["problem"][0], sum(w["problem"]),
        w["frac_clean"] or float("nan"), w["clean"][0], sum(w["clean"]), w["fisher"]["p"]))

    if res.get("de_contrasts"):
        print("\n[次要] 存算 DE 率（n=%d）：" % len(dec))
        for r in res["de_contrasts"]:
            print("  %-22s 中位差 %+0.3f [%+0.3f,%+0.3f] p=%.3g n=%d/%d" % (
                r["label"], r["diff"], r["ci"][0], r["ci"][1], r["p_mwu"], r["n_a"], r["n_b"]))

    print("\n[阴性对照]")
    for k, v in nc.items():
        if v:
            print("  %-20s 中位差 %+0.3f p=%.3g" % (k, v["diff"], v["p_mwu"]))
    print("\n写出：partA2_landscape.tsv, partA2_group_*.tsv, partA2_discordant_missed.tsv, partA2_analysis.json")

    # ---------------- 稳健性打印
    print("\n" + "=" * L)
    print("[稳健性] 粗分类（批次伪影 / 信息不足 / 无问题 / 已校正）")
    st = pd.DataFrame(res["group_by_status"]["auc_pc1_batch"])
    if len(st):
        print("  auc_pc1_batch：")
        print(st.to_string(index=False))
    print("\n  两两对比：")
    for r in res["status_contrasts"]:
        print("    %-22s %-16s 中位差 %+0.3f [%+0.3f,%+0.3f] p=%.3g n=%d/%d" % (
            r["col"], r["label"], r["diff"], r["ci"][0], r["ci"][1], r["p_mwu"], r["n_a"], r["n_b"]))
    print("\n  二分类调整模型（以 无问题 为参照）：")
    for k, m in res["models_binary"].items():
        if not m:
            continue
        print("    %-22s n=%-3d R2=%.3f  F=%.2f p_perm=%.4f" % (
            k, m["n"], m["R2"] or float("nan"), m["F"]["F"], m["F"]["p_perm"]))
        for t in m["terms"]:
            if t["term"].startswith("annot_status"):
                print("        %-32s %+0.4f [%+0.4f,%+0.4f]" % (t["term"], t["beta"], t["ci_lo"], t["ci_hi"]))

    print("\n[阈值敏感性] 不一致表随「批次主导」阈值变化：")
    for r in res["threshold_sensitivity"]:
        print("    thr=%.2f  一致阳性=%2d 仅标注=%2d 仅数据=%2d 一致阴性=%2d | 灵敏 %.2f 特异 %.2f PPV %.2f 漏检 %d" % (
            r["thr"], r["both_pos"], r["annot_only"], r["svd_only"], r["both_neg"],
            r["sens"] or float("nan"), r["spec"] or float("nan"), r["ppv"] or float("nan"), r["missed"]))

    print("\n[匹配对照] log2n+pc1_var_frac 最近邻配对：")
    for k, v in res["matched"].items():
        if v:
            print("    %-24s n_pairs=%d 中位差 %+0.3f Wilcoxon p=%.3g 上升侧 %d/%d" % (
                k, v["n_pairs"], v["median_diff"], v["wilcoxon_p"], v["n_win"], v["n_pairs"]))

    print("\n[UC 专属子集] 每臂>=4 的队列 n=%d  分布 %s" % (
        res["uc_subset"]["n"], json.dumps(res["uc_subset"]["by_status"], ensure_ascii=False)))
    print("    疾病主轴可分性按标注：", json.dumps(res["uc_subset"]["disease_auc_by_status"], ensure_ascii=False)[:300])

    if isinstance(res.get("geeq_vs_svd"), list):
        print("\n[GEEQ 交叉验证] 各质量分与 PC1 批次可分性的 Spearman：")
        for r in res["geeq_vs_svd"]:
            print("    %-28s rho=%+.3f p=%.3g n=%d (唯一值 %d)" % (
                r["score"], r["rho"], r["p"], r["n"], r["n_unique"]))
    print("\n[混杂结构] 标注问题 vs 队列属性（Fisher）：")
    for r in res["confounding"]:
        print("    %-20s OR=%.2f p=%.3g" % (r["var"], r["odds_ratio"], r["p"]))

    ps = res["paired_subset"]
    print("\n[配对子集] 两个主轴指标都可用 n=%d  分布 %s" % (ps["n"], json.dumps(ps["by_status"], ensure_ascii=False)))
    print("    %-10s %-4s %-12s %-12s %-12s %s" % ("类别", "n", "AUC_batch", "AUC_bio", "差值", "批次胜出比例"))
    for g, v in ps["medians"].items():
        print("    %-10s %-4d %-12.3f %-12.3f %+-12.3f %.3f" % (
            g, v["n"], v["auc_batch"], v["auc_bio"], v["auc_diff"], v["frac_batch_wins"]))
    print("    （图件与报告统一引用此口径）")


if __name__ == "__main__":
    main()
