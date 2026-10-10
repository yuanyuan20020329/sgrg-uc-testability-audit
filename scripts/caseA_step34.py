# -*- coding: utf-8 -*-
"""
案例 A · Step 3 + Step 4
 Step 3 忠实复现原文流水线：t 检验 p<0.001（无 FDR） -> 差异 FRG -> LASSO + SVM-RFE -> logistic ROC
 Step 4 审计阶梯：
   L2 批次-疾病嵌套：需要批次标签（Gemma batch info）+ 技术轴（SVD PC1 与分组相关）
   L3 技术负载：面板评分 vs 技术协变量
   L4 随机集校准：同规模随机基因集跑同一流水线，看 AUC 与「入选数」的零分布
   L5 检验力/设计：非配对两组，报效应量与最小可达 p

输出 results/caseA_step34_*.tsv / .json / .md
"""
import os, json, sys, math, itertools, collections
import numpy as np
import pandas as pd

ROOT = r"C:\Users\H-P\WorkBuddy\应激颗粒"
EXPR = os.path.join(ROOT, "results", "caseA_expr_7673.tsv")
RES = os.path.join(ROOT, "results")
RNG = np.random.default_rng(20261009)

TUMOR = "liver tumor"
NORMAL = "adjacent non-tumor liver"


def load():
    df = pd.read_csv(EXPR, sep="\t")
    genes = [c for c in df.columns if c not in ("sample", "label")]
    return df, genes


def group_masks(df):
    y = (df["label"] == TUMOR).astype(int).values
    return y


def ttest_p(X, y):
    """按基因做 Welch/Student 两样本 t 检验，返回 p 值数组（等价 R t.test 默认 var.equal=FALSE? 原文未写；两口径都算）。"""
    from scipy import stats
    a = X[y == 1]
    b = X[y == 0]
    t, p = stats.ttest_ind(a, b, axis=0, equal_var=True)
    t2, p2 = stats.ttest_ind(a, b, axis=0, equal_var=False)
    return p, p2


def cohens_d(X, y):
    a = X[y == 1]
    b = X[y == 0]
    na, nb = len(a), len(b)
    va, vb = a.var(axis=0, ddof=1), b.var(axis=0, ddof=1)
    sp = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    sp[sp == 0] = np.nan
    return (a.mean(axis=0) - b.mean(axis=0)) / sp


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    ranked = p[o] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(ranked, 0, 1)
    return out


def roc_auc(y, score):
    from sklearn.metrics import roc_auc_score
    return roc_auc_score(y, score)


def pipeline_auc(X, y, genes, p_cut=0.001, use_fdr=False, seed=0):
    """
    走原文流水线的一个简化但忠实的替身：
      1) 按 p_cut 选差异基因（可选用 BH-FDR）
      2) 若差异基因 >= 3：LASSO(logistic, 10-fold CV) 选非零系数
      3) SVM-RFE（线性核）逐步淘汰到 ~7 个
      4) 交集基因做 logistic 回归，在**同一数据**上算 AUC（与原文一致，无嵌套 CV）
    返回 dict。
    """
    from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
    from sklearn.svm import LinearSVC
    from sklearn.preprocessing import StandardScaler
    p, _ = ttest_p(X, y)
    pv = bh_fdr(p) if use_fdr else p
    sel = np.where(pv < p_cut)[0]
    out = {"n_de": int(len(sel)), "auc": np.nan, "n_lasso": 0, "n_svm": 0, "n_final": 0,
           "final_genes": []}
    if len(sel) < 3:
        return out
    Xs = StandardScaler().fit_transform(X[:, sel])
    yy = y
    # LASSO logistic
    try:
        lasso = LogisticRegressionCV(Cs=20, cv=10, penalty="l1", solver="liblinear",
                                     scoring="roc_auc", max_iter=5000, random_state=seed)
        lasso.fit(Xs, yy)
        coef = lasso.coef_.ravel()
        lsel = np.where(np.abs(coef) > 1e-8)[0]
    except Exception:
        lsel = np.arange(Xs.shape[1])
    out["n_lasso"] = int(len(lsel))
    # SVM-RFE
    try:
        keep = list(range(Xs.shape[1]))
        target = max(7, int(round(0.4 * len(keep))))
        while len(keep) > target and len(keep) > 3:
            clf = LinearSVC(C=1.0, max_iter=20000, random_state=seed)
            clf.fit(Xs[:, keep], yy)
            w = np.abs(clf.coef_.ravel())
            drop = int(np.argmin(w))
            keep.pop(drop)
        ssel = np.array(keep)
    except Exception:
        ssel = np.arange(Xs.shape[1])
    out["n_svm"] = int(len(ssel))
    inter = np.intersect1d(lsel, ssel)
    if len(inter) == 0:
        inter = lsel if len(lsel) else ssel
    out["n_final"] = int(len(inter))
    out["final_genes"] = [genes[sel[i]] for i in inter]
    # logistic on same data (原文口径：无嵌套 CV)
    try:
        lr = LogisticRegression(max_iter=5000, C=1.0)
        lr.fit(Xs[:, inter], yy)
        prob = lr.predict_proba(Xs[:, inter])[:, 1]
        out["auc"] = float(roc_auc(yy, prob))
    except Exception as e:
        out["err"] = repr(e)[:80]
    return out


def make_random_sets(genes, X, n_sets=200, stratify_bins=5, seed=1):
    """
    表达分层匹配随机集：按每个基因的平均表达分箱，从同箱基因池抽同规模集合。
    这里基因宇宙 = 提供给 X 的基因（62 个）——注意样本量小，故改为
    「从全基因组抽」不可行（未下载全矩阵）；本版做的是**在 62 基因内置换**，
    用于检验「面板是否比同规模任意子集更好」。★ 这是保守版，须在报告中注明。
    """
    rng = np.random.default_rng(seed)
    mean_expr = X.mean(axis=0)
    order = np.argsort(mean_expr)
    n = len(genes)
    # 分层：按表达排序后切成 k 层，随机集从各层按比例取
    k = stratify_bins
    bins = np.array_split(order, k)
    sets = []
    kk = int(sys.argv[1]) if len(sys.argv) > 1 else n  # 面板规模
    for _ in range(n_sets):
        pick = []
        per = max(1, kk // k)
        for b in bins:
            take = min(len(b), per)
            pick.extend(rng.choice(b, size=take, replace=False))
        while len(pick) < kk:
            pick.append(rng.integers(0, n))
        pick = list(dict.fromkeys(pick))[:kk]
        sets.append(np.array(pick))
    return sets


def main():
    df, genes = load()
    y = group_masks(df)
    X = df[genes].values.astype(float)
    print("n samples=%d  tumor=%d  normal=%d  genes=%d" % (len(df), y.sum(), (1 - y).sum(), len(genes)))

    # ---------- Step 3 忠实复现 ----------
    p, p_welch = ttest_p(X, y)
    d = cohens_d(X, y)
    fdr = bh_fdr(p)
    base = pd.DataFrame({"gene": genes, "p_student": p, "p_welch": p_welch,
                         "p_bh_fdr": fdr, "cohens_d": d})
    base = base.sort_values("p_student")
    base.to_csv(os.path.join(RES, "caseA_step3_de.tsv"), sep="\t", index=False)
    n_sig = int((p < 0.001).sum())
    n_sig_fdr = int((fdr < 0.001).sum())
    n_fdr05 = int((fdr < 0.05).sum())
    print("DE at p<0.001: %d ; at BH-FDR<0.001: %d ; BH-FDR<0.05: %d" % (n_sig, n_sig_fdr, n_fdr05))

    rep = pipeline_auc(X, y, genes, p_cut=0.001, use_fdr=False)
    print("忠实复现：DE=%d LASSO=%d SVM=%d final=%d AUC=%.4f" % (
        rep["n_de"], rep["n_lasso"], rep["n_svm"], rep["n_final"], rep["auc"]))
    print("  final genes:", rep["final_genes"])

    rep_fdr = pipeline_auc(X, y, genes, p_cut=0.05, use_fdr=True)
    print("FDR 口径：DE=%d final=%d AUC=%.4f" % (
        rep_fdr["n_de"], rep_fdr["n_final"], rep_fdr["auc"]))

    # ---------- L4 随机集校准 ----------
    kk = rep["n_de"] if rep["n_de"] >= 3 else max(3, n_sig)
    print("\nL4: 随机基因集校准（面板规模 k=%d，200 集）..." % kk)
    rng = np.random.default_rng(20261009)
    mean_expr = X.mean(axis=0)
    order = np.argsort(mean_expr)
    bins = np.array_split(order, 5)
    aurocs, ndes = [], []
    for i in range(200):
        pick = []
        per = max(1, kk // 5)
        for b in bins:
            take = min(len(b), per)
            pick.extend(rng.choice(b, size=take, replace=False))
        pick = np.array(sorted(set(pick)))[:kk]
        r = pipeline_auc(X[:, pick], y, [genes[j] for j in pick], p_cut=0.001, seed=i)
        if not math.isnan(r["auc"]):
            aurocs.append(r["auc"])
        ndes.append(r["n_de"])
    aurocs = np.array(aurocs)
    ndes = np.array(ndes)
    panel_auc = rep["auc"]
    emp_p = float((1 + (aurocs >= panel_auc).sum()) / (len(aurocs) + 1)) if len(aurocs) else np.nan
    print("  random AUC: median=%.3f  q95=%.3f  max=%.3f" % (
        np.median(aurocs), np.quantile(aurocs, 0.95), aurocs.max()))
    print("  panel AUC=%.3f -> empirical p=%.4f  (rank %d/%d)" % (
        panel_auc, emp_p, int((aurocs >= panel_auc).sum()) + 1, len(aurocs) + 1))
    print("  random DE counts: median=%.0f  (panel DE=%d)" % (np.median(ndes), rep["n_de"]))

    summary = {
        "dataset": "GSE36376 (Gemma id 7673)",
        "n_samples": int(len(df)), "n_tumor": int(y.sum()), "n_normal": int((1 - y).sum()),
        "n_genes": len(genes),
        "step3": {"n_de_p001": n_sig, "n_de_fdr001": n_sig_fdr, "n_de_fdr05": n_fdr05,
                  "replicate": rep, "replicate_fdr": rep_fdr},
        "l4": {"k": int(kk), "n_random": int(len(aurocs)),
               "random_auc_median": float(np.median(aurocs)),
               "random_auc_q95": float(np.quantile(aurocs, 0.95)),
               "random_auc_max": float(aurocs.max()),
               "panel_auc": panel_auc, "emp_p_one_sided": emp_p,
               "random_de_median": float(np.median(ndes)), "panel_de": rep["n_de"]},
    }
    json.dump(summary, open(os.path.join(RES, "caseA_step34_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("\nwrote", os.path.join(RES, "caseA_step34_summary.json"))


if __name__ == "__main__":
    main()
