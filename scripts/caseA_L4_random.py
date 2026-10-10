# -*- coding: utf-8 -*-
"""
案例 A · L4 随机基因集校准（快速版）。

科学问题：原文 AUC≈1.000 的「面板」是否**特异于铁死亡基因**，还是任何同等规模的基因子集
在「选特征 + 同集评估」的流水线下都能拿到 AUC≈1？

方法（避免 200×10 折 CV 太慢）：
  对每个随机集，执行与原文同构的简化流水线：
    (a) 组间 t 检验，取 p<0.001 的基因（如 <3 个则整集入模）
    (b) 用全部候选基因直接做 logistic 回归（等价于「不限稀疏」的上界）
    (c) 用**同一数据**算 AUC（原文口径，无嵌套 CV）
  再平行做「7 折交叉验证 AUC」作为诚实的泛化下界。
  最终把真实面板的 AUC 放进随机集 AUC 的零分布里定位。

★ 关键设计：随机集在**表达分层**内抽取（按平均表达分 5 层），否则会系统抽到高表达基因。
"""
import os, json, sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

ROOT = r"C:\Users\H-P\WorkBuddy\应激颗粒"
EXPR = os.path.join(ROOT, "results", "caseA_expr_7673.tsv")
RES = os.path.join(ROOT, "results")
N_RANDOM = 500
TUMOR = "liver tumor"


def load():
    df = pd.read_csv(EXPR, sep="\t")
    genes = [c for c in df.columns if c not in ("sample", "label")]
    y = (df["label"] == TUMOR).astype(int).values
    return df, genes, y


def de_select(X, y, cut=0.001):
    p = stats.ttest_ind(X[y == 1], X[y == 0], axis=0, equal_var=True)[1]
    return np.where(p < cut)[0], p


def auc_insample(X, y, idx):
    if len(idx) < 1:
        return np.nan, np.nan
    Xs = StandardScaler().fit_transform(X[:, idx])
    lr = LogisticRegression(max_iter=5000, C=1.0)
    lr.fit(Xs, y)
    pr = lr.predict_proba(Xs)[:, 1]
    a_in = roc_auc_score(y, pr)
    # honest 7-fold CV AUC
    try:
        cv = StratifiedKFold(n_splits=7, shuffle=True, random_state=0)
        pipe_auc = cross_val_score(LogisticRegression(max_iter=5000, C=1.0),
                                   Xs, y, cv=cv, scoring="roc_auc").mean()
    except Exception:
        pipe_auc = np.nan
    return float(a_in), float(pipe_auc)


def main():
    df, genes, y = load()
    X = df[genes].values.astype(float)
    n_genes = len(genes)
    print("n=%d tumor=%d genes=%d" % (len(df), y.sum(), n_genes), flush=True)

    # --- 真实铁死亡面板 ---
    idx_panel, p_panel = de_select(X, y)
    n_de = len(idx_panel)
    a_panel_in, a_panel_cv = auc_insample(X, y, idx_panel)
    # 全 62 基因
    a_all_in, a_all_cv = auc_insample(X, y, np.arange(n_genes))
    print("PANEL: DE=%d  in-sample AUC=%.4f  CV AUC=%.4f" % (n_de, a_panel_in, a_panel_cv), flush=True)
    print("ALL62: in-sample AUC=%.4f  CV AUC=%.4f" % (a_all_in, a_all_cv), flush=True)

    # --- 分层随机集 ---
    rng = np.random.default_rng(20261009)
    order = np.argsort(X.mean(axis=0))
    bins = np.array_split(order, 5)
    k = max(3, n_de)  # 与面板同样的规模
    print("\nL4: %d 个分层随机集 (k=%d) ..." % (N_RANDOM, k), flush=True)

    in_aucs, cv_aucs, des = [], [], []
    for i in range(N_RANDOM):
        pick = []
        per = max(1, k // 5)
        for b in bins:
            take = min(len(b), per)
            pick.extend(rng.choice(b, size=take, replace=False))
        pick = np.array(sorted(set(pick)))[:k]
        idxr, _ = de_select(X[:, pick], y)
        des.append(len(idxr))
        ai, ac = auc_insample(X[:, pick], y, idxr if len(idxr) >= 1 else np.arange(len(pick)))
        if not np.isnan(ai):
            in_aucs.append(ai)
        if not np.isnan(ac):
            cv_aucs.append(ac)
        if (i + 1) % 100 == 0:
            print("   %d/%d  (median in-sample AUC so far %.3f)" % (
                i + 1, N_RANDOM, np.median(in_aucs)), flush=True)

    in_aucs = np.array(in_aucs)
    cv_aucs = np.array(cv_aucs)
    des = np.array(des)

    def rankp(val, dist):
        return float((1 + (dist >= val).sum()) / (len(dist) + 1))

    out = {
        "dataset": "GSE36376 (Gemma 7673)", "n": int(len(df)),
        "n_tumor": int(y.sum()), "n_normal": int((1 - y).sum()), "n_genes": n_genes,
        "panel": {"n_de": int(n_de), "auc_insample": a_panel_in, "auc_cv7": a_panel_cv},
        "all_genes": {"auc_insample": a_all_in, "auc_cv7": a_all_cv},
        "random": {
            "n_sets": int(len(in_aucs)), "k": int(k),
            "in_auc_median": float(np.median(in_aucs)),
            "in_auc_q95": float(np.quantile(in_aucs, 0.95)),
            "in_auc_max": float(in_aucs.max()),
            "cv_auc_median": float(np.median(cv_aucs)),
            "cv_auc_q95": float(np.quantile(cv_aucs, 0.95)),
            "cv_auc_max": float(cv_aucs.max()),
            "de_median": float(np.median(des)),
            "frac_random_in_auc_ge_panel": float((in_aucs >= a_panel_in).mean()),
            "frac_random_in_auc_ge_099": float((in_aucs >= 0.99).mean()),
            "frac_random_in_auc_ge_095": float((in_aucs >= 0.95).mean()),
            "emp_p_panel_insample": rankp(a_panel_in, in_aucs),
            "emp_p_panel_cv": rankp(a_panel_cv, cv_aucs),
        },
    }
    json.dump(out, open(os.path.join(RES, "caseA_L4_random.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("\n=== L4 RESULT ===")
    print("panel in-sample AUC=%.4f (emp p=%.4f)" % (a_panel_in, out["random"]["emp_p_panel_insample"]))
    print("random in-sample AUC: median=%.3f q95=%.3f max=%.3f" % (
        out["random"]["in_auc_median"], out["random"]["in_auc_q95"], out["random"]["in_auc_max"]))
    print("random CV AUC:        median=%.3f q95=%.3f max=%.3f" % (
        out["random"]["cv_auc_median"], out["random"]["cv_auc_q95"], out["random"]["cv_auc_max"]))
    print("frac random in-sample AUC>=0.95 : %.1f%%" % (100 * out["random"]["frac_random_in_auc_ge_095"]))
    print("frac random in-sample AUC>=0.99 : %.1f%%" % (100 * out["random"]["frac_random_in_auc_ge_099"]))
    print("wrote results/caseA_L4_random.json")


if __name__ == "__main__":
    main()
