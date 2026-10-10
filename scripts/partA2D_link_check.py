# -*- coding: utf-8 -*-
"""A2（bulk 队列级）× Part D（单细胞 donor 级）能否接上的可行性检验

核心问题：两侧是不是共享同一条【定量】规律，而不只是【定性类比】？

  假设 H：技术污染强度 → 名义结果的虚高幅度（剂量-反应）

  Part D 侧：per cell_type，depth-disease 共线强度 R²(log_nUMI ~ disease, donor 级)
             与「原始 Δ 到抽稀 Δ 的收缩幅度」是否单调相关
  A2 侧   ：(AUC_batch − AUC_bio) 与队列的技术污染强度分是否单调相关

若两侧都成立 → 可以写成「同一规律在两个分辨率上的复现」（强）
若只有一侧 → 只能写成类比（弱），必须如实标注。

只读已有结果表，不重跑 h5ad。
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT = os.path.join(ROOT, "results")


def fit_r2(y, x):
    """单变量 OLS 的 R²（x 为 0/1 分组变量）。"""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 6 or len(np.unique(x)) < 2:
        return np.nan, np.nan, len(x)
    X = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - (r ** 2).sum() / ss_tot if ss_tot > 0 else np.nan
    # 斜率 t 检验
    n, p = X.shape
    s2 = (r ** 2).sum() / (n - p)
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(cov[1, 1])
    t = beta[1] / se if se > 0 else np.nan
    pval = 2 * stats.t.sf(abs(t), n - p) if np.isfinite(t) else np.nan
    return float(r2), float(pval), len(x)


# ============================================================
# Part D 侧
# ============================================================
print("=" * 78)
print("【Part D 侧】ds4：depth-disease 共线强度  vs  效应量收缩幅度")
print("=" * 78)

pb = pd.read_csv(os.path.join(OUT, "partD_ds4_tearle_pseudobulk_score.tsv"), sep="\t")
t_orig = pd.read_csv(os.path.join(OUT, "partD_ds4_tearle_uc_vs_hc_tests.tsv"), sep="\t")
t_ds1k = pd.read_csv(os.path.join(OUT, "partD_ds4_tearle_ds1k_uc_vs_hc_tests.tsv"), sep="\t")

pb["is_uc"] = (pb["disease"] == "ulcerative colitis").astype(int)
pb = pb[pb["disease"].isin(["ulcerative colitis", "normal"])]

rows = []
for ct, g in pb.groupby("cell_type"):
    n_uc = int(g["is_uc"].sum())
    n_hc = int((1 - g["is_uc"]).sum())
    if n_uc < 3 or n_hc < 3:
        continue
    r2, pv, n = fit_r2(g["log_nUMI"].values, g["is_uc"].values)
    rows.append({"cell_type": ct, "n_donor": n, "r2_depth_disease": r2,
                 "p_depth": pv, "med_nUMI_uc": g.loc[g.is_uc == 1, "nUMI"].median(),
                 "med_nUMI_hc": g.loc[g.is_uc == 0, "nUMI"].median()})
dep = pd.DataFrame(rows)

o = t_orig[t_orig["score"] == "SGRG"][["cell_type", "diff", "p_perm", "cohen_d"]].rename(
    columns={"diff": "diff_orig", "p_perm": "p_orig", "cohen_d": "d_orig"})
k = t_ds1k[t_ds1k["score"] == "SGRG"][["cell_type", "diff", "p_perm", "cohen_d"]].rename(
    columns={"diff": "diff_ds1k", "p_perm": "p_ds1k", "cohen_d": "d_ds1k"})

m = dep.merge(o, on="cell_type").merge(k, on="cell_type")
m["shrink"] = m["diff_orig"] - m["diff_ds1k"]          # 收缩幅度（正=原始更抬高）
m["shrink_abs"] = m["diff_orig"].abs() - m["diff_ds1k"].abs()
m["log2_depth_ratio"] = np.log2(m["med_nUMI_uc"] / m["med_nUMI_hc"])
m = m.dropna(subset=["r2_depth_disease", "shrink"])
m = m.sort_values("r2_depth_disease", ascending=False)

print("匹配到的细胞类型：%d" % len(m))
print("\n按 depth-disease 共线强度 R² 排序（前 15）：")
print(m.head(15)[["cell_type", "n_donor", "r2_depth_disease", "med_nUMI_uc", "med_nUMI_hc",
                  "diff_orig", "diff_ds1k", "shrink", "p_orig", "p_ds1k"]]
      .to_string(index=False, float_format=lambda v: "%.4g" % v))

for ycol, lab in [("shrink", "Δ原始−Δ抽稀"), ("shrink_abs", "|Δ原始|−|Δ抽稀|")]:
    rho, p = stats.spearmanr(m["r2_depth_disease"], m[ycol])
    r, pr = stats.pearsonr(m["r2_depth_disease"], m[ycol])
    print("\n[R²(depth~disease)] vs [%s]:  Spearman rho=%+.3f p=%.3g | Pearson r=%+.3f p=%.3g (n=%d)"
          % (lab, rho, p, r, pr, len(m)))

rho, p = stats.spearmanr(m["log2_depth_ratio"], m["shrink"])
print("[log2(UC/对照 深度比)] vs [收缩幅度]:  rho=%+.3f p=%.3g" % (rho, p))

# 分组对比：共线强 vs 弱
med = m["r2_depth_disease"].median()
hi = m[m["r2_depth_disease"] >= med]
lo = m[m["r2_depth_disease"] < med]
print("\n共线强度中位 = %.3f" % med)
print("  高共线组 (n=%d): 原始显著数=%d/%d, 抽稀后显著数=%d/%d, 中位收缩=%.4f"
      % (len(hi), (hi.p_orig < 0.05).sum(), len(hi), (hi.p_ds1k < 0.05).sum(), len(hi), hi.shrink.median()))
print("  低共线组 (n=%d): 原始显著数=%d/%d, 抽稀后显著数=%d/%d, 中位收缩=%.4f"
      % (len(lo), (lo.p_orig < 0.05).sum(), len(lo), (lo.p_ds1k < 0.05).sum(), len(lo), lo.shrink.median()))
u, pu = stats.mannwhitneyu(hi["shrink"], lo["shrink"])
print("  两组收缩幅度 Mann-Whitney p=%.4g" % pu)

m.to_csv(os.path.join(OUT, "partD_ds4_depth_bio_collinearity.tsv"), sep="\t", index=False)
print("\n写出 → results/partD_ds4_depth_bio_collinearity.tsv")


# ============================================================
# A2 侧
# ============================================================
print()
print("=" * 78)
print("【A2 侧】队列级：(批次占据 − 生物占据)  vs  技术污染强度分")
print("=" * 78)

a2 = pd.read_csv(os.path.join(OUT, "partA2_cohorts.tsv"), sep="\t")
a2["delta"] = a2["auc_pc1_batch"] - a2["auc_pc1_primary"]      # 用 PC1 口径
a2["delta_best"] = a2["auc_best_batch"] - a2["auc_best_primary"]  # PC1-3 最佳口径

cands = ["geeq_public", "geeq_pub_batch_effect", "geeq_pub_batch_confound",
         "geeq_batch_info", "geeq_replicates", "geeq_outliers", "geeq_platforms",
         "geeq_mean_corr", "geeq_corr_var", "n_confound_factors", "n_samples",
         "n_platforms", "pc1_var_frac", "pc1_2_var_frac"]

print("\n%-26s %6s %14s %14s" % ("污染强度分", "非空", "rho(PC1口径)", "rho(PC1-3口径)"))
print("-" * 66)
for c in cands:
    if c not in a2.columns:
        continue
    v = pd.to_numeric(a2[c], errors="coerce")
    for dcol, store in [("delta", []), ("delta_best", [])]:
        pass
    results = []
    for dcol in ["delta", "delta_best"]:
        sub = pd.DataFrame({"x": v, "y": a2[dcol]}).dropna()
        if len(sub) >= 20 and sub["x"].nunique() > 2:
            rr, pp = stats.spearmanr(sub["x"], sub["y"])
            results.append("%+.3f(p=%.1e,n=%d)" % (rr, pp, len(sub)))
        else:
            results.append("—")
    print("%-26s %6d %14s %14s" % (c, v.notna().sum(), results[0], results[1]))

# 分类口径：Gemma 明说 disease 被 batch 混杂的队列
print("\n--- 分类：confound_disease（Gemma 明标 disease 与批次混杂）---")
for dcol, lab in [("delta", "PC1"), ("delta_best", "PC1-3最佳")]:
    sub = a2.dropna(subset=[dcol, "auc_pc1_batch"])
    a = sub[sub["confound_disease"] == True][dcol]
    b = sub[sub["confound_disease"] == False][dcol]
    if len(a) >= 3 and len(b) >= 3:
        u, pp = stats.mannwhitneyu(a, b)
        print("  %-10s confounded n=%2d 中位=%+.3f | 非混杂 n=%3d 中位=%+.3f | MWU p=%.3g"
              % (lab, len(a), a.median(), len(b), b.median(), pp))

print("\n--- annot_status 分组下的 delta ---")
print(a2.groupby("annot_status")[["delta", "delta_best", "auc_pc1_batch", "auc_best_batch",
                                  "auc_pc1_primary", "auc_best_primary"]]
      .agg(["count", "median"]).round(3).to_string())
